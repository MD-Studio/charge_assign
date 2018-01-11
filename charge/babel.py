from collections import defaultdict
from enum import Enum
from io import StringIO
from typing import Any

import networkx as nx

from charge.settings import IACM_MAP, IACM_ELEMENTS


class IOType(Enum):
    LGF = 1
    GML = 2
    RDKIT = 3
    OPENBABEL = 4
    PYBEL = 5


class BondType(Enum):
    SINGLE = 'SINGLE'
    DOUBLE = 'DOUBLE'
    TRIPLE = 'TRIPLE'
    AROMATIC = 'AROMATIC'
    UNKNOWN = 'UNKNOWN'


def __lgf_to_nx(obj: str) -> nx.Graph:
    graph = nx.Graph()

    nodes, edges, header = False, False, False
    idx = 0
    for line in map(lambda line: line.strip(), obj.splitlines()):
        if not line or len(line) == 0:
            continue
        if '@nodes' in line:
            nodes, edges, header = True, False, True
            continue
        if '@edges' in line:
            nodes, edges, header = False, True, True
            continue
        if nodes:
            if header:
                keys = line.split()
                header = False
                if not 'label' in keys:
                    raise ValueError('Missing attribute "label".')
                if not 'label2' in keys:
                    raise ValueError('Missing attribute "label2".')
                if not 'atomType' in keys:
                    raise ValueError('Missing attribute "atomType".')
            else:
                values = line.split()
                attr = {'atom_type': IACM_ELEMENTS[int(values[keys.index('atomType')])-1],
                        'label': values[keys.index('label2')],
                        'idx': idx}
                idx += 1
                for key in keys:
                    if key == 'label' or key == 'atomType' or key == 'label2':
                        continue
                    attr[key] = values[keys.index(key)]
                graph.add_node(int(values[keys.index('label')]), **attr)

        if edges:
            if header:
                keys = line.split()
                header = False
                if not 'label' in keys:
                    raise ValueError('Missing attribute "label".')
            else:
                values = line.split()
                attr = {}
                for key in keys:
                    if key == 'label':
                        continue
                    attr[key] = values[keys.index(key)+2]
                graph.add_edge(int(values[0]), int(values[1]), **attr)
    return graph


def __nx_to_lgf(graph: nx.Graph) -> str:
    out = StringIO()
    out.write('@nodes\n')

    idxmap = {}
    keys = set.union(*map(lambda n: set(n[1].keys()), graph.nodes(data=True)))
    keys.remove('atom_type')
    keys.remove('idx')
    if 'label' in keys:
        keys.remove('label')
    out.write('label\tlabel2\tatomType\t')
    for k in keys:
        out.write('%s\t' % k)
    out.write('\n')

    el_count = defaultdict(lambda: 1)
    for i, (v, data) in enumerate(graph.nodes(data=True)):
        idx = i + 1
        idxmap[v] = idx
        if not data['atom_type'] in IACM_MAP:
            raise ValueError('Unknown atom type: %s' % data['atom_type'])

        element = IACM_MAP[data['atom_type']]
        if 'label' in data:
            label = data['label']
        else:
            label = '%s%d' % (element, el_count[element])
        el_count[element] += 1
        iacm_num = IACM_ELEMENTS.index(data['atom_type'])+1

        out.write('%d\t%s\t%d\t' % (idx, label, iacm_num))

        for k in keys:
            if k in data:
                out.write('%s\t' % str(data[k]))
            else:
                out.write('%s\t')
        out.write('\n')

    keys = set.union(*map(lambda n: set(n[2].keys()), graph.edges(data=True)))
    if 'orig_bond_type' in keys:
        keys.remove('orig_bond_type')

    out.write('@edges\n\t\tlabel\t')
    for k in keys:
        out.write('%s\t' % k)
    out.write('\n')

    for i, (u, v, data) in enumerate(graph.edges(data=True)):
        out.write('%d\t%d\t%d\t' % (idxmap[u], idxmap[v], i))
        for k in keys:
            if k in data:
                out.write('%s\t' % str(data[k]))
            else:
                out.write('\t')
        out.write('\n')

    content = out.getvalue()

    out.close()
    return content


def __gml_to_nx(obj: str) -> nx.Graph:
    graph = nx.parse_gml(obj)
    for idx, v in enumerate(graph.nodes_iter()):
        attr = graph.node[v]
        attr['idx'] = idx
        if not 'atom_type' in attr:
            raise ValueError('Missing attribute "atom_type" for atom %s' % v)
    return graph


def __nx_to_gml(graph: nx.Graph) -> str:
    cp = graph.copy()
    for v, data in cp.nodes(data=True):
        if 'idx' in data:
            del data['idx']

    for e, data in cp.edges(data=True):
        if 'orig_bond_type' in data:
            del data['orig_bond_type']

    return nx.generate_gml(graph)


def __rdmol_to_nx(obj: Any) -> nx.Graph:
    from rdkit import Chem
    if not isinstance(obj, Chem.Mol):
        raise ValueError('Invalid input. Chem.Mol required.')
    graph = nx.Graph()
    
    for idx, atom in enumerate(obj.GetAtoms()):
        props = atom.GetPropsAsDict()
        props['idx'] = idx
        props['atom_type'] = atom.GetSymbol()
        graph.add_node(atom.GetIdx(), **props)
    
    for bond in obj.GetBonds():
        props = bond.GetPropsAsDict()
        props['orig_bond_type'] = bond.GetBondType()
        if props['orig_bond_type'] == Chem.BondType.SINGLE:
            props['bond_type'] = BondType.SINGLE
        elif props['orig_bond_type'] == Chem.BondType.DOUBLE:
            props['bond_type'] = BondType.DOUBLE
        elif props['orig_bond_type'] == Chem.BondType.TRIPLE:
            props['bond_type'] = BondType.TRIPLE
        elif props['orig_bond_type'] == Chem.BondType.AROMATIC:
            props['bond_type'] = BondType.AROMATIC
        else:
            props['bond_type'] = BondType.UNKNOWN
        graph.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), **props)
    
    return graph


def __nx_to_rdmol(graph: nx.Graph) -> Any:
    from rdkit import Chem
    mol = Chem.RWMol()
    idxmap = {}

    for v, data in graph.nodes(data=True):
        idxmap[v] = mol.AddAtom(data['atom_type'])
        atom = mol.GetAtomWithIdx(idxmap[v])
        for k, v in data:
            if k == 'atom_type':
                continue
            if isinstance(v, bool):
                atom.SetBoolProp(k, v)
            elif isinstance(v, float):
                atom.SetDoubleProp(k, v)
            elif isinstance(v, int):
                atom.SetIntProp(k, v)
            elif isinstance(v, str):
                atom.SetProp(k, v)
            else:
                raise ValueError('Invalid property type for key %s' % k)

    for (u, v), data in graph.edges(data=True):
        if 'orig_bond_type' in data and isinstance(data['orig_bond_type'], Chem.BondType):
            btype = data['orig_bond_type']
        elif data['bond_type'] == BondType.SINGLE:
            btype = Chem.BondType.SINGLE
        elif data['bond_type'] == BondType.DOUBLE:
            btype = Chem.BondType.DOUBLE
        elif data['bond_type'] == BondType.TRIPLE:
            btype = Chem.BondType.TRIPLE
        elif data['bond_type'] == BondType.AROMATIC:
            btype = Chem.BondType.AROMATIC
        else:
            btype = Chem.BondType.UNSPECIFIED

        idx = mol.AddBond(idxmap[u], idxmap[v], btype)
        bond = mol.GetBondWithIdx(idx)
        for k, v in data:
            if k == 'bond_type' or k == 'orig_bond_type':
                continue
            if isinstance(v, bool):
                bond.SetBoolProp(k, v)
            elif isinstance(v, float):
                bond.SetDoubleProp(k, v)
            elif isinstance(v, int):
                bond.SetIntProp(k, v)
            elif isinstance(v, str):
                bond.SetProp(k, v)
            else:
                raise ValueError('Invalid property type for key %s' % k)

    return Chem.Mol(mol)


def __obmol_to_nx(obj: Any) -> nx.Graph:
    # TODO openbabel
    pass


def __nx_to_obmol(graph: nx.Graph) -> Any:
    # TODO openbabel
    pass


def __pymol_to_nx(obj: Any) -> nx.Graph:
    # TODO pybel
    pass


def __nx_to_pymol(graph: nx.Graph) -> Any:
    # TODO pybel
    pass


def convert_from(molecule: Any, type: IOType) -> nx.Graph:
    if type == IOType.LGF:
        return __lgf_to_nx(molecule)
    elif type == IOType.GML:
        return __gml_to_nx(molecule)
    elif type == IOType.RDKIT:
        return __rdmol_to_nx(molecule)
    elif type == IOType.OPENBABEL:
        return __obmol_to_nx(molecule)
    elif type == IOType.PYBEL:
        return __pymol_to_nx(molecule)
    else:
        raise ValueError('Unknown data format: %s' % type)


def convert_to(graph: nx.Graph, type: IOType) -> Any:
    if type == IOType.LGF:
        return __nx_to_lgf(graph)
    elif type == IOType.GML:
        return __nx_to_gml(graph)
    elif type == IOType.RDKIT:
        return __nx_to_rdmol(graph)
    elif type == IOType.OPENBABEL:
        return __nx_to_obmol(graph)
    elif type == IOType.PYBEL:
        return __nx_to_pymol(graph)
    else:
        raise ValueError('Unknown data format: %s' % type)
