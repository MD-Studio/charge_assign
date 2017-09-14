import hashlib
import json
import subprocess
import tempfile
import warnings
import zipfile
from collections import defaultdict, deque
from zipfile import ZipFile

import jsonschema
import msgpack
import re

import os
from rdkit import Chem
from itertools import groupby


from charge_assign.settings import NAUTY_EXC, IACM_MAP, DATA_SCHEMA


class LoadError(Exception):
    pass


class AssignmentError(Warning):
    pass


class Charger:

    def __init__(self, repository: str='atb.zip', nauty: str=NAUTY_EXC):
        if not zipfile.is_zipfile(repository):
            raise LoadError('%s is not a valid zip file.' % repository)

        if not os.path.isfile(nauty) or not os.access(nauty, os.X_OK):
            raise LoadError('Could not find dreadnaut executable at: "%s". Did you install nauty (http://users.cecs.'
                            'anu.edu.au/~bdm/nauty/)?' % nauty)

        self.__charges_iacm = defaultdict(dict)
        self.__charges_elem = defaultdict(dict)

        with ZipFile(repository, 'r') as zf:
            for name in zf.namelist():
                m = re.match('(\d+)\.json$', name)
                if not m:
                    continue
                shell = int(m.group(1))

                data = json.loads(zf.read(name).decode())
                if not data or len(data) == 0:
                    raise LoadError('%s is invalid.' % name)

                jsonschema.validate(data, DATA_SCHEMA)

                for _, obj in data.items():
                    iacm, elem = self.__canonize(self.__json_to_rdmol(obj))
                    self.__charges_iacm[shell][iacm] = (obj['charges'][0], obj['uncertainties'][0])
                    self.__charges_elem[shell][elem] = (obj['charges'][0], obj['uncertainties'][0])


    def __json_to_rdmol(self, obj: dict) -> Chem.Mol:
        rdmol = Chem.RWMol()
        mapping = {}
        for idx, iacm in obj['atoms'].items():
            idx = int(idx)
            if not iacm in IACM_MAP:
                raise LoadError('%s not known.' % iacm)

            element = IACM_MAP[iacm]
            atom = Chem.Atom(element)
            atom.SetProp('iacm', iacm)
            atom.SetBoolProp('core', idx in obj['core_ids'])
            mapping[idx] = rdmol.AddAtom(atom)
        for bond in obj['bonds']:
            rdmol.AddBond(mapping[bond[0]], mapping[bond[1]], Chem.BondType.SINGLE)
        return rdmol


    def __nauty_input(self, rdmol: Chem.Mol) -> str:
        input_str = ' n={num_atoms} g {edges}. f=[{node_partition}] cxb'.format(
            num_atoms=rdmol.GetNumAtoms(),
            edges=self.__nauty_edges(rdmol),
            node_partition=self.__nauty_node_partition(rdmol),
        )

        return input_str


    def __nauty_edges(self, rdmol: Chem.Mol) -> str:
        bonds = sorted(set(map(lambda bond: (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) \
            if bond.GetBeginAtomIdx() < bond.GetEndAtomIdx() else (bond.GetEndAtomIdx(), bond.GetBeginAtomIdx()),
                        rdmol.GetBonds())))

        return ';'.join(map(lambda bond: '{0}:{1}'.format(bond[0], bond[1]), bonds))


    def __nauty_node_partition(self, rdmol) -> str:
        atoms = map(lambda atom: (atom.GetIdx(), atom.GetSymbol()), rdmol.GetAtoms())

        return '|'.join(map(lambda group: ','.join(map(lambda atom: str(atom[0]), group[1])),
                            groupby(atoms, key=lambda atom: atom[1])))


    def __nauty_output(self, nautstr: str, rdmol: Chem.Mol) -> tuple:
        lines = [line.strip() for line in nautstr.split('seconds')[-1].strip().split('\n')]

        adj = [[int(val) for val in line[:-1].split(':')[-1].split()] for line in lines[1:]]

        iacm, elements, core = zip(*map(
            lambda idx: (rdmol.GetAtomWithIdx(idx).GetProp('iacm'),
                         rdmol.GetAtomWithIdx(idx).GetSymbol(),
                         rdmol.GetAtomWithIdx(idx).GetBoolProp('core')),
            map(int, lines[0].split())))
        return hashlib.md5(msgpack.packb([adj, iacm, core])).hexdigest(), \
               hashlib.md5(msgpack.packb([adj, elements, core])).hexdigest()


    def __canonize(self, rdmol: Chem.Mol) -> tuple:

        with tempfile.TemporaryFile(buffering=0) as tmp:
            tmp.write(self.__nauty_input(rdmol).encode())
            tmp.seek(0)

            p = subprocess.Popen(
                [NAUTY_EXC],
                shell=True,
                stdin=tmp,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )

            out, err = p.communicate()

            if len(err) > 0:
                raise Exception(err.decode())
            if len(out) == 0:
                raise Exception()

        return self.__nauty_output(out.strip().decode(), rdmol)


    def __set_canons(self, rdmol: Chem.Mol) -> Chem.Mol:

        for atom in rdmol.GetAtoms():
            for shell in sorted(self.__charges_iacm.keys(), reverse=True):
                fragment = Chem.RWMol()
                mapping = {}

                visited = set()
                depth = defaultdict(int)
                queue = deque()

                idx = atom.GetIdx()
                queue.append(idx)
                visited.add(idx)

                while len(queue) > 0:
                    current_idx = queue.popleft()
                    current_atom = rdmol.GetAtomWithIdx(current_idx)

                    copy = Chem.Atom(current_atom.GetSymbol())
                    copy.SetBoolProp('core', current_idx == idx)
                    copy.SetProp('iacm', current_atom.GetProp('iacm'))
                    mapping[fragment.AddAtom(copy)] = current_idx

                    if depth[current_idx] < shell:
                        for bond in current_atom.GetBonds():
                            other_idx = bond.GetEndAtomIdx() \
                                if bond.GetBeginAtomIdx() == current_idx else bond.GetBeginAtomIdx()
                            if not other_idx in visited:
                                visited.add(other_idx)
                                depth[other_idx] = depth[current_idx] + 1
                                queue.append(other_idx)

                for u in range(fragment.GetNumAtoms()):
                    for v in range(u+1, fragment.GetNumAtoms()):
                        if rdmol.GetBondBetweenAtoms(mapping[u], mapping[v]):
                            fragment.AddBond(u, v, Chem.BondType.SINGLE)

                iacm_key, elem_key = self.__canonize(fragment)
                atom.SetProp('fragment_%d' % shell, elem_key)
                atom.SetProp('fragment_iacm_%d' % shell, iacm_key)

        return rdmol


    def __set_partial_charges(self, rdmol: Chem.Mol) -> Chem.Mol:
        copy = Chem.Mol(rdmol)
        for atom in copy.GetAtoms():
            atom.SetProp('iacm', self.__iacm(atom))
        self.__set_canons(copy)
        for idx in range(rdmol.GetNumAtoms()):
            atom = rdmol.GetAtomWithIdx(idx)
            for shell in sorted(self.__charges_iacm.keys(), reverse=True):
                key_iacm = copy.GetAtomWithIdx(idx).GetProp('fragment_iacm_%d' % shell)
                key_elem = copy.GetAtomWithIdx(idx).GetProp('fragment_%d' % shell)
                if key_iacm in self.__charges_iacm[shell]:
                    charge, uncertainty = self.__charges_iacm[shell][key_iacm]
                    atom.SetDoubleProp('partial_charge', charge)
                    atom.SetDoubleProp('uncertainty', uncertainty)
                    break
                elif key_elem in self.__charges_elem[shell]:
                    charge, uncertainty = self.__charges_elem[shell][key_elem]
                    atom.SetDoubleProp('partial_charge', charge)
                    atom.SetDoubleProp('uncertainty', uncertainty)
                    break
            else:
                warnings.warn(AssignmentError('Could not assign charge to atom #%d.' % idx))

        return rdmol


    def __iacm(self, atom):
        def get_bonded_atoms(atom, bond_type=None):
            if not bond_type:
                return [bond.GetEndAtom() if atom.GetIdx() == bond.GetBeginAtomIdx() else bond.GetBeginAtom()
                        for bond in atom.GetBonds()]
            else:
                return [bond.GetEndAtom() if atom.GetIdx() == bond.GetBeginAtomIdx() else bond.GetBeginAtom()
                        for bond in atom.GetBonds() if bond.GetBondType() == bond_type]

        element = atom.GetSymbol()

        bas = get_bonded_atoms(atom)
        if element == 'C':
            bhs = list(filter(lambda a: a.GetSymbol() == 'H', bas))
            if len(bas) == 4 and len(bhs) == 0:
                return 'CH0'
            else:
                return 'C'
        elif element == 'H':
            if bas and bas[0].GetSymbol() == 'C':
                return 'HC'
            else:
                return 'H'
        elif element == 'O':
            if len(list(filter(lambda a: a.GetSymbol() == 'C', bas))) == len(bas) and len(bas) > 1:
                return 'OE'
            elif len(bas) > 1:
                return 'OA'
            elif bas and len(list(filter(lambda a: a.GetSymbol() == 'O' and \
                            len(get_bonded_atoms(a)) == 1, get_bonded_atoms(bas[0])))) > 1 and \
                            bas != get_bonded_atoms(atom, Chem.BondType.AROMATIC):
                return 'OM'
            else:
                return 'O'
        elif element == 'N':
            if len(bas) > 3:
                return 'NL'
            elif len(bas) == 1:
                return 'NR'
            elif len(get_bonded_atoms(atom, Chem.BondType.AROMATIC)) > 1:
                return 'NR'
            elif len(list(filter(lambda a: a.GetSymbol() == 'H', bas))) < 2:
                return 'N'
            else:
                return 'NT'
        elif element == 'S':
            if len(bas) > 2:
                return 'SDmso'
            else:
                return 'S'
        elif element == 'P':
            return 'P,SI'
        elif element == 'Si':
            return 'AR'
        elif element == 'F':
            return 'F'
        elif element == 'Cl':
            return 'CL'
        elif element == 'Br':
            return 'BR'
        else:
            return element


    def charge_inchi(self, inchi: str) -> Chem.Mol:
        rdmol = Chem.MolFromInchi(inchi, removeHs=False, sanitize=False)
        rdmol = Chem.AddHs(rdmol)
        return self.__set_partial_charges(rdmol)


    def charge_smiles(self, smiles: str) -> Chem.Mol:
        rdmol = Chem.MolFromSmiles(smiles)
        rdmol = Chem.AddHs(rdmol)
        return self.__set_partial_charges(rdmol)


    def charge_pdb(self, pdb: str) -> Chem.Mol:
        rdmol = Chem.MolFromPDBBlock(pdb, removeHs=False, sanitize=False)
        rdmol = Chem.AddHs(rdmol)
        return self.__set_partial_charges(rdmol)


    def charge_rdmol(self, rdmol: Chem.Mol) -> Chem.Mol:
        rdmol = Chem.AddHs(rdmol)
        return self.__set_partial_charges(rdmol)

