"""
Conversion of a `networkx <http://networkx.github.io/>`_ graph from/to different chemical formats.

.. moduleauthor:: Martin S. Engler <martin.engler@cwi.nl>

"""

import re
from collections import defaultdict
from enum import Enum
from io import StringIO
from typing import Any

import networkx as nx

from charge.bond_type import BondType
from charge.settings import IACM_MAP, IACM_ELEMENTS


class IOType(Enum):
    """Supported I/O formats."""
    LGF = 1
    GML = 2
    RDKIT = 3
    OPENBABEL = 4
    PYBEL = 5
    ITP = 6

    def get_extension(self) -> str:
        """Returns the corresponding filename extension, including
        the period.

        Raises ValueError for RDKIT, OPENBABEL and PYBEL, as these
        support multiple formats."""
        ext = {
                IOType.LGF: ".lgf",
                IOType.GML: ".gml",
                IOType.ITP: ".itp",
                }
        try:
            return ext[self]
        except KeyError:
            raise ValueError('Unsupported file type: {}'.format(self.name))


def __lgf_to_nx(obj: str) -> nx.Graph:
    graph = nx.Graph()
    node_keys = None
    edge_keys = None
    el_count = defaultdict(int)
    charge_groups = set()

    nodes, edges, header = False, False, False
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
                # create a dict mapping column names to column indices
                node_keys = dict((key, i) for i, key in enumerate(line.split()))
                header = False
                # label is the numerical id of the atom
                if not 'label' in node_keys:
                    raise ValueError('Missing attribute "label".')
                # atomType is the IACM atom type number
                if not 'atomType' in node_keys:
                    raise ValueError('Missing attribute "atomType".')
            elif node_keys:
                values = line.split()
                type_idx = int(values[node_keys['atomType']])-1
                if type_idx < 0 or type_idx >= len(IACM_ELEMENTS):
                    raise ValueError('Unknown atom type: %d' % type_idx)
                iacm_atom_type = IACM_ELEMENTS[type_idx]
                plain_atom_type = IACM_MAP[iacm_atom_type]
                # label2 is the name, e.g. 'H4'
                if 'label2' in node_keys:
                    label = values[node_keys['label2']]
                else:
                    # make a name if there isn't one
                    el = IACM_MAP[iacm_atom_type]
                    el_count[el] += 1
                    label = '{}{}'.format(el, el_count[el])
                if 'initColor' in node_keys:
                    charge_group = int(values[node_keys['initColor']])
                    charge_groups.add(charge_group)
                else:
                    charge_group = 0
                # graph node is an int (from 'label')
                # iacm: str is IACM symbolic atom type (from 'atomType')
                # atom_type: str is the plain atom type (derived from IACM)
                # label: str is the name of the atom (from 'label2' or generated)
                # charge_group: int is the charge group the atom belongs to (from 'initColor')
                attr = {
                        'iacm': iacm_atom_type, 'atom_type': plain_atom_type,
                        'label': label, 'charge_group': charge_group}
                for key in node_keys:
                    if key == 'label' or key == 'atomType' or key == 'label2' or key == 'initColor':
                        continue
                    attr[key] = values[node_keys[key]]
                graph.add_node(int(values[node_keys['label']]), **attr)

        if edges:
            if header:
                edge_keys = dict((key, i) for i, key in enumerate(line.split()))
                if 'label' in edge_keys:
                    del edge_keys['label']

                header = False
            else:
                values = line.split()
                attr = {}
                for key in edge_keys:
                    attr[key] = values[edge_keys[key]+2]
                if 'bond_type' in attr:
                    if int(attr['bond_type']) == 1:
                        attr['bond_type'] = BondType.SINGLE
                    elif int(attr['bond_type']) == 2:
                        attr['bond_type'] = BondType.DOUBLE
                    elif int(attr['bond_type']) == 3:
                        attr['bond_type'] = BondType.TRIPLE
                    elif int(attr['bond_type']) == 4:
                        attr['bond_type'] = BondType.AROMATIC
                else:
                    attr['bond_type'] = BondType.UNKNOWN
                # 'label' is removed, rest becomes attributes
                graph.add_edge(int(values[0]), int(values[1]), **attr)

    # Assume all charge groups are perfect
    # ATB data seem to have one-atom charge groups if they have them at all
    # Anyway, it's never actually used
    if charge_groups:
        graph.graph['group_charges'] = dict((k, 0.0) for k in charge_groups)
    else:
        graph.graph['group_charges'] = {0: 0.0}

    return graph


def __nx_to_lgf(graph: nx.Graph) -> str:
    out = StringIO()
    out.write('@nodes\n')

    idxmap = {}
    # make a set of all the unique keys in all the node attributes
    keys = set.union(*map(lambda n: set(n[1].keys()), graph.nodes(data=True)))

    if not 'atom_type' in keys and not 'iacm' in keys:
        raise ValueError('Need at least one of "atom_type" and "iacm" attributes.')

    if 'label' in keys:
        keys.remove('label')

    out.write('label\tlabel2\tatomType\t')
    for k in keys:
        if k == 'iacm' or k == 'atom_type':
            continue
        elif k == 'charge_group':
            out.write('initColor\t')
        else:
            out.write('%s\t' % k)
    out.write('\n')

    el_count = defaultdict(int)
    for i, (v, data) in enumerate(graph.nodes(data=True)):
        idx = i + 1
        idxmap[v] = idx
        if 'iacm' in keys:
            if not data['iacm'] in IACM_MAP:
                raise ValueError('Unknown atom type: %s' % data['iacm'])
        elif not data['atom_type'] in IACM_MAP:
            raise ValueError('Unknown atom type: %s' % data['atom_type'])

        if 'label' in data:
            label = data['label']
        else:
            element = IACM_MAP[data['atom_type']]
            el_count[element] += 1
            label = '%s%d' % (element, el_count[element])

        if 'iacm' in keys:
            iacm_num = IACM_ELEMENTS.index(data['iacm']) + 1
        else:
            iacm_num = IACM_ELEMENTS.index(data['atom_type'])+1

        out.write('%d\t%s\t%d\t' % (idx, label, iacm_num))

        for k in keys:
            if k == 'iacm' or k == 'atom_type':
                continue
            if k in data:
                out.write('%s\t' % str(data[k]))
            else:
                raise ValueError('Missing attribute {} for node {}.'.format(k, v))
        out.write('\n')

    keys = set.union(*map(lambda n: set(n[2].keys()), graph.edges(data=True)))
    if 'rkit_bond_type' in keys:
        keys.remove('rdkit_bond_type')

    if 'bond_type' in keys:
        bond_types = set(n[2]['bond_type'] for n in graph.edges(data=True)
                         if 'bond_type' in n[2] and isinstance(n[2]['bond_type'], BondType))
        if not bond_types or bond_types == {BondType.UNKNOWN}:
            keys.remove('bond_type')

    out.write('@edges\n\t\tlabel\t')
    for k in keys:
        if k == 'bond_type':
            out.write('bondType\t')
        else:
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
    el_count = defaultdict(int)
    group_charges = dict()

    graph = nx.parse_gml(obj)
    cp = nx.convert_node_labels_to_integers(graph, first_label=1, label_attribute='label')
    if 'name' in graph.graph:
        cp.graph['name'] = graph.graph['name']
    else:
        del cp.graph['name']

    for v, data in cp.nodes.data():
        if not isinstance(v, int):
            raise ValueError('id {0} is not int.'.format(v))
        if 'partialcharge' in data:
            data['partial_charge'] = data.pop('partialcharge')
        if not 'atomtype' in data:
            raise ValueError('Missing attribute "atomtype" for atom {0}'.format(v))
        data['iacm'] = data.pop('atomtype')
        if not data['iacm'] in IACM_ELEMENTS:
            raise ValueError('Unknown "atom_type" for atom {0}: {1}'.format(v, data['atom_type']))
        element = IACM_MAP[data['iacm']]
        data['atom_type'] = element
        el_count[element] += 1
        if not 'label' in data or not isinstance('label', str):
            data['label'] = '%s%d' % (element, el_count[element])
        if not 'chargegroup' in data:
            data['charge_group'] = 0
        else:
            data['charge_group'] = data.pop('chargegroup')
        group_charges[data['charge_group']] = 0.0

    for _, _, data in cp.edges.data():
        if not 'bondtype' in data:
            data['bond_type'] = BondType.UNKNOWN
        else:
            try:
                data['bond_type'] = BondType(data.pop('bondtype'))
            except:
                del data['bondtype']
                data['bond_type'] = BondType.UNKNOWN

    for k in list(cp.graph.keys()):
        m = re.fullmatch(r'groupcharge(\d+)', k)
        if m:
            group_idx = int(m.group(1))
            if group_idx in group_charges:
                group_charges[group_idx] = float(cp.graph[k])
            del cp.graph[k]

    for _, data in cp.nodes.data():
        if not data['charge_group'] in group_charges:
            group_charges[data['charge_group']] = 0.0
    cp.graph['group_charges'] = group_charges

    return cp


def __nx_to_gml(graph: nx.Graph) -> str:
    cp = graph.copy()
    el_count = defaultdict(int)

    for _, data in cp.nodes(data=True):
        element = IACM_MAP[data['iacm']]
        el_count[element] += 1
        if not 'label' in data or not isinstance('label', str):
            data['label'] = '%s%d' % (element, el_count[element])

    nx.relabel_nodes(cp, mapping=dict((v, data['label']) for v, data in cp.nodes.data()),  copy=False)

    for u, v, data in cp.edges(data=True):
        if 'rdkit_bond_type' in data:
            del data['rdkit_bond_type']
        if 'bond_type' in data and isinstance(data['bond_type'], BondType):
            data['bondtype'] = data.pop('bond_type').value

    for v, data in cp.nodes(data=True):
        if 'iacm' in data:
            data['atomtype'] = data.pop('iacm')
            del data['atom_type']
        else:
            data['atomtype'] = data.pop('atom_type')
        if 'charge_group' in data:
            data['chargegroup'] = data.pop('charge_group')
        if 'partial_charge' in data:
            data['partialcharge'] = data.pop('partial_charge')

    if 'group_charges' in cp.graph and isinstance(cp.graph['group_charges'], dict):
        for k, v in cp.graph['group_charges'].items():
            cp.graph['groupcharge%d' % k] = v
        del cp.graph['group_charges']

    return "\n".join(nx.generate_gml(cp))


def __rdmol_to_nx(obj: Any) -> nx.Graph:
    from rdkit import Chem
    if not isinstance(obj, Chem.Mol):
        raise ValueError('Invalid input. Chem.Mol required.')
    graph = nx.Graph()
    el_count = defaultdict(int)
    group_charges = {}

    for atom in obj.GetAtoms():
        props = atom.GetPropsAsDict()

        element = atom.GetSymbol()
        el_count[element] += 1
        if 'atom_type' in props:
            props['iacm'] = props['atom_type']
            if props['iacm'] not in IACM_MAP:
                raise ValueError('Unknown atom type: %s' % props['iacm'])
            element = IACM_MAP[props['iacm']]
            props['atom_type'] = element
        if not 'charge_group' in props or not isinstance(props['charge_group'], int):
            props['charge_group'] = 0
        group_charges[props['charge_group']] = 0.0
        if not 'label' in props or not isinstance(props['label'], str):
            props['label'] = '%s%d' % (element, el_count[element])
        graph.add_node(atom.GetIdx(), **props)

    for bond in obj.GetBonds():
        props = bond.GetPropsAsDict()
        props['rdkit_bond_type'] = bond.GetBondType()
        if props['rdkit_bond_type'] == Chem.BondType.SINGLE:
            props['bond_type'] = BondType.SINGLE
        elif props['rdkit_bond_type'] == Chem.BondType.DOUBLE:
            props['bond_type'] = BondType.DOUBLE
        elif props['rdkit_bond_type'] == Chem.BondType.TRIPLE:
            props['bond_type'] = BondType.TRIPLE
        elif props['rdkit_bond_type'] == Chem.BondType.AROMATIC:
            props['bond_type'] = BondType.AROMATIC
        else:
            props['bond_type'] = BondType.UNKNOWN
        graph.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), **props)

    for k, v in obj.GetPropsAsDict().items():
        m = re.fullmatch(r'group_charge_(\d+)', k)
        if m:
            group_idx = int(m.group(1))
            if group_idx in group_charges:
                group_charges[group_idx] = float(v)

    if len(group_charges) == 0:
        for _, data in graph.nodes.data():
            if not data['charge_group'] in group_charges:
                group_charges[data['charge_group']] = 0.0
    graph.graph['group_charges'] = group_charges

    return graph


def __nx_to_rdmol(graph: nx.Graph) -> Any:
    from rdkit import Chem
    mol = Chem.RWMol()
    idxmap = {}
    el_count = defaultdict(int)

    for v, data in graph.nodes(data=True):
        if not 'atom_type' in data:
            raise ValueError('Missing atom type for atom {}'.format(v))

        element = IACM_MAP[data['atom_type']]
        atom = Chem.Atom(element)
        el_count[element] += 1

        if not 'label' in data:
            atom.SetProp('label', '%s%d' % (element, el_count[element]))

        if 'iacm' in data and data['iacm'] != element:
            atom.SetProp('atom_type', data['iacm'])
        elif data['atom_type'] != element:
            atom.SetProp('atom_type', data['atom_type'])

        for k, val in data.items():
            if k == 'atom_type' or k == 'iacm':
                continue
            if isinstance(val, bool):
                atom.SetBoolProp(k, val)
            elif isinstance(val, float):
                atom.SetDoubleProp(k, val)
            elif isinstance(val, int):
                atom.SetIntProp(k, val)
            elif isinstance(val, str):
                atom.SetProp(k, val)
            else:
                raise ValueError('Invalid property type for key %s' % k)
        idxmap[v] = mol.AddAtom(atom)

    for u, v, data in graph.edges(data=True):
        if 'rdkit_bond_type' in data and isinstance(data['rdkit_bond_type'], Chem.BondType):
            btype = data['rdkit_bond_type']
        elif 'bond_type' in data:
            if data['bond_type'] == BondType.SINGLE:
                btype = Chem.BondType.SINGLE
            elif data['bond_type'] == BondType.DOUBLE:
                btype = Chem.BondType.DOUBLE
            elif data['bond_type'] == BondType.TRIPLE:
                btype = Chem.BondType.TRIPLE
            elif data['bond_type'] == BondType.AROMATIC:
                btype = Chem.BondType.AROMATIC
            else:
                btype = Chem.BondType.UNSPECIFIED
        else:
            btype = Chem.BondType.UNSPECIFIED

        idx = mol.AddBond(idxmap[u], idxmap[v], btype)
        bond = mol.GetBondWithIdx(idx-1)
        for k, v in data.items():
            if k == 'bond_type' or k == 'rdkit_bond_type':
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

    if 'group_charges' in graph.graph and isinstance(graph.graph['group_charges'], dict):
        for k, v in graph.graph['group_charges'].items():
            mol.SetDoubleProp('group_charge_%d' % k, float(v))

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

def __itp_to_nx(obj: Any) -> nx.Graph:
    graph = nx.Graph()
    node_keys = None
    edge_keys = None
    group_charges = dict()
    charge_groups = 0

    nodes, edges, header = False, False, False
    for line in map(lambda line: line.strip(), obj.splitlines()):
        if not line or len(line) == 0:
            continue
        if re.search(r'\[\s*atoms\s*\]', line):
            nodes, edges, header = True, False, True
            continue
        if re.search(r'\[\s*pairs\s*\]', line):
            nodes, edges, header = False, True, True
            continue
        if not header and re.match(r'\s*;', line):
            continue
        if nodes:
            if header:
                node_keys = re.search(r';([^;]+)', line)
                header = False
                if node_keys:
                    node_keys = node_keys.group(1).split()
                    if not 'nr' in node_keys:
                        raise ValueError('Missing attribute "nr".')
                    if not 'atom' in node_keys:
                        raise ValueError('Missing attribute "atom".')
                    if not 'type' in node_keys:
                        raise ValueError('Missing attribute "type".')
            elif node_keys:
                m = re.search(r'([^;]+)(;([^;]+))*', line)
                if m:
                    values = m.group(1).split()
                    atom_type = values[node_keys.index('type')]
                    if not atom_type in IACM_MAP:
                        raise ValueError('Unknown atom type: {}'.format(atom_type))
                    element = IACM_MAP[atom_type]
                    attr = {'iacm': atom_type,
                            'atom_type': element,
                            'label': values[node_keys.index('atom')]}
                    if 'charge' in node_keys:
                        attr['partial_charge'] = values[node_keys.index('charge')]
                    for key in node_keys:
                        if key == 'nr' or key == 'atom' or key == 'type' or key == 'charge':
                            continue
                        idx = node_keys.index(key)
                        if idx < len(values):
                            attr[key] = values[idx]
                    attr['charge_group'] = charge_groups
                    graph.add_node(int(values[node_keys.index('nr')]), **attr)

                    if m.group(3):
                        try:
                            group_charges[charge_groups] = float(m.group(3).strip())
                            charge_groups += 1
                        except:
                            pass

        if edges:
            if header:
                edge_keys = re.search(r';([^;]+)', line)
                header = False
                if edge_keys:
                    edge_keys = edge_keys.group(1).split()
                    if not 'ai' in edge_keys:
                        raise ValueError('Missing attribute "ai".')
                    if not 'aj' in edge_keys:
                        raise ValueError('Missing attribute "aj".')
            elif edge_keys:
                values = re.search(r'([^;]+)', line)
                if values:
                    values = values.group(1).split()
                    attr = {}
                    for key in edge_keys:
                        if key == 'ai' or key == 'aj':
                            continue
                        attr[key] = values[edge_keys.index(key)]
                    attr['bond_type'] = BondType.UNKNOWN
                    graph.add_edge(int(values[edge_keys.index('ai')]), int(values[edge_keys.index('aj')]), **attr)

    if len(group_charges) == 0:
        group_charges = {0: 0.0}
    graph.graph['group_charges'] = group_charges

    return graph

def convert_from(molecule: Any, type: IOType) -> nx.Graph:
    """Convert a molecule to a `networkx <http://networkx.github.io/>`_ graph.

    See :doc:`conversion` for more information.

    :param molecule: lgf, gml or itp string or rdkit, openbabel or pybel molecule object
    :type molecule: Any
    :param type: type to convert from
    :type type: charge.babel.IOType
    :returns: (networkx.Graph) -- the converted graph
    """
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
    elif type == IOType.ITP:
        return __itp_to_nx(molecule)
    else:
        raise ValueError('Unknown data format: %s' % type)


def convert_to(graph: nx.Graph, type: IOType) -> Any:
    """Convert a `networkx <http://networkx.github.io/>`_ graph to a molecule.

    See :doc:`conversion` for more information.

    :param graph: the graph to convert
    :type graph: networkx.Graph
    :param type: type to convert to
    :type type: charge.babel.IOType
    :returns: (lgf, gml or itp strings or rdkit, openbabel or pybel molecule object) -- the converted molecule
    """
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
    elif type == IOType.ITP:
        raise ValueError('Exporting to ITF not supported.')
    else:
        raise ValueError('Unknown data format: %s' % type)
