"""
Conversion of a `networkx <http://networkx.github.io/>`_ graph from/to different chemical formats.

.. moduleauthor:: Martin S. Engler <martin.engler@cwi.nl>

"""

from collections import defaultdict
from enum import Enum
from io import StringIO
from typing import Any

import networkx as nx
import re

from charge.settings import IACM_MAP, IACM_ELEMENTS


class IOType(Enum):
    """Supported I/O formats."""
    LGF = 1
    GML = 2
    RDKIT = 3
    OPENBABEL = 4
    PYBEL = 5
    ITP = 6


class BondType(Enum):
    """Supported chemical bond types."""

    SINGLE = 'SINGLE'
    DOUBLE = 'DOUBLE'
    TRIPLE = 'TRIPLE'
    AROMATIC = 'AROMATIC'
    UNKNOWN = 'UNKNOWN'


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
                node_keys = dict((key, i) for i, key in enumerate(line.split()))
                header = False
                if not 'label' in node_keys:
                    raise ValueError('Missing attribute "label".')
                if not 'atomType' in node_keys:
                    raise ValueError('Missing attribute "atomType".')
            elif node_keys:
                values = line.split()
                type_idx = int(values[node_keys['atomType']])-1
                if type_idx < 0 or type_idx >= len(IACM_ELEMENTS):
                    raise ValueError('Unknown atom type: %d' % type_idx)
                atom_type = IACM_ELEMENTS[type_idx]
                if 'label2' in node_keys:
                    label = values[node_keys['label2']]
                else:
                    el = IACM_MAP[atom_type]
                    el_count[el] += 1
                    label = '{}{}'.format(el, el_count[el])
                if 'initColor' in node_keys:
                    charge_group = int(values[node_keys['initColor']])
                    charge_groups.add(charge_group)
                else:
                    charge_group = 0
                attr = {'atom_type': atom_type, 'label': label, 'charge_group': charge_group}
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
                graph.add_edge(int(values[0]), int(values[1]), **attr)

    if charge_groups:
        graph.graph['charge_groups'] = dict((k, 0.0) for k in charge_groups)
    else:
        graph.graph['charge_groups'] = {0: 0.0}

    return graph


def __nx_to_lgf(graph: nx.Graph) -> str:
    out = StringIO()
    out.write('@nodes\n')

    idxmap = {}
    keys = set.union(*map(lambda n: set(n[1].keys()), graph.nodes(data=True)))

    if not 'atom_type' in keys:
        raise ValueError('Missing attribute "atom_type".')

    keys.remove('atom_type')
    if 'label' in keys:
        keys.remove('label')

    out.write('label\tlabel2\tatomType\t')
    for k in keys:
        if k == 'iacm':
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
            if k == 'iacm':
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
    # TODO check all keys and values
    # TODO group charges
    graph = nx.parse_gml(obj)
    for v in enumerate(graph.nodes_iter()):
        if not 'atom_type' in graph.node[v]:
            raise ValueError('Missing attribute "atom_type" for atom {0}'.format(v))
        if not graph.node[v]['atom_type'] in IACM_ELEMENTS:
            raise ValueError('Unknown "atom_type" for atom {0}: {1}'.format(v, graph.node[v]['atom_type']))

    return graph


def __nx_to_gml(graph: nx.Graph) -> str:
    # TODO check all keys and values
    # TODO group charges
    cp = graph.copy()

    for e, data in cp.edges(data=True):
        if 'rdkit_bond_type' in data:
            del data['rdkit_bond_type']

    return nx.generate_gml(cp)


def __rdmol_to_nx(obj: Any) -> nx.Graph:
    # TODO check all keys and values
    from rdkit import Chem
    if not isinstance(obj, Chem.Mol):
        raise ValueError('Invalid input. Chem.Mol required.')
    graph = nx.Graph()
    
    for atom in obj.GetAtoms():
        props = atom.GetPropsAsDict()
        props['atom_type'] = atom.GetSymbol()
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
    
    return graph


def __nx_to_rdmol(graph: nx.Graph) -> Any:
    # TODO check all keys and values
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
        bond = mol.GetBondWithIdx(idx)
        for k, v in data:
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
                node_keys = re.search(r';(.+);*.*', line)
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
                values = re.search(r'(.+);*.*', line)
                if values:
                    values = values.group(1).split()
                    attr = {'atom_type': IACM_ELEMENTS[int(values[node_keys.index('type')]) - 1],
                            'label': values[node_keys.index('atom')]}
                    for key in node_keys:
                        if key == 'nr' or key == 'atom' or key == 'type':
                            continue
                        idx = node_keys.index(key)
                        if idx < len(values):
                            attr[key] = values[idx]
                    graph.add_node(int(values[node_keys.index('nr')]), **attr)

        if edges:
            if header:
                edge_keys = re.search(r';(.+);*.*', line)
                header = False
                if edge_keys:
                    edge_keys = edge_keys.group(1).split()
                    if not 'ai' in edge_keys:
                        raise ValueError('Missing attribute "ai".')
                    if not 'aj' in edge_keys:
                        raise ValueError('Missing attribute "aj".')
            elif edge_keys:
                values = re.search(r'(.+);*.*', line)
                if values:
                    values = values.group(1).split()
                    attr = {}
                    for key in node_keys:
                        attr[key] = values[node_keys.index(key) + 2]
                    graph.add_edge(int(values[0]), int(values[1]), **attr)
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
