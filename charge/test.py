import unittest

from charge.babel import convert_from, IOType, convert_to, BondType


GRAPH_ATTRIBUTES = {'group_charges': {0: 0.0}}
GRAPH_NODES = [(1, {'atom_type': 'C', 'label': 'C1', 'charge_group': 0}),
               (2, {'atom_type': 'HC', 'label': 'H1', 'charge_group': 0}),
               (3, {'atom_type': 'HC', 'label': 'H2', 'charge_group': 0}),
               (4, {'atom_type': 'HC', 'label': 'H3', 'charge_group': 0}),
               (5, {'atom_type': 'HC', 'label': 'H4', 'charge_group': 0})]
GRAPH_EDGES = [(1, 2, {'bond_type': BondType.UNKNOWN}),
               (1, 3, {'bond_type': BondType.UNKNOWN}),
               (1, 4, {'bond_type': BondType.UNKNOWN}),
               (1, 5, {'bond_type': BondType.UNKNOWN})]


class TestCharger(unittest.TestCase):

    def test_babel_lgf(self):
        lgf = """@nodes
        label   label2  atomType    initColor
        1       C1      12          0
        2       H1      20          0
        3       H2      20          0
        4       H3      20          0
        5       H4      20          0
        @edges
                        label
        1       2       0
        1       3       1
        1       4       2
        1       5       3"""

        graph = convert_from(lgf, IOType.LGF)
        self.assertDictEqual(graph.graph, GRAPH_ATTRIBUTES)
        self.assertListEqual(list(graph.nodes.data()), GRAPH_NODES)
        self.assertListEqual(list(graph.edges.data()), GRAPH_EDGES)

        graph2 = convert_from(convert_to(graph, IOType.LGF), IOType.LGF)
        self.assertListEqual(list(graph.nodes.data()), list(graph2.nodes.data()))
        self.assertListEqual(list(graph.edges.data()), list(graph2.edges.data()))

    def test_babel_gml(self):
        gml = \
"""graph [
    groupcharge0 0.0
    node [
        id 0
        label "C1"
        atomtype "C"
        chargegroup 0
    ]
    node [
        id 1
        label "H1"
        atomtype "HC"
        chargegroup 0
    ]
    node [
        id 2
        label "H2"
        atomtype "HC"
        chargegroup 0
    ]
    node [
        id 3
        label "H3"
        atomtype "HC"
        chargegroup 0
    ]
    node [
        id 4
        label "H4"
        atomtype "HC"
    ]
    edge [
        source 0
        target 1
        bondtype "UNKNOWN"
    ]
    edge [
        source 0
        target 2
        bondtype "UNKNOWN"
    ]
    edge [
        source 0
        target 3
        bondtype "UNKNOWN"
    ]
    edge [
        source 0
        target 4
        bondtype "UNKNOWN"
    ]
]"""

        graph = convert_from(gml, IOType.GML)
        self.assertDictEqual(graph.graph, GRAPH_ATTRIBUTES)
        self.assertListEqual(list(graph.nodes.data()), GRAPH_NODES)
        self.assertListEqual(list(graph.edges.data()), GRAPH_EDGES)

        graph2 = convert_from(convert_to(graph, IOType.LGF), IOType.LGF)
        self.assertListEqual(list(graph.nodes.data()), list(graph2.nodes.data()))
        self.assertListEqual(list(graph.edges.data()), list(graph2.edges.data()))

    def test_babel_itp(self):
        itp = \
"""[ atoms ]
;  nr  type  atom total_charge
    1     C    C1
    2    HC    H1
    3    HC    H2
    4    HC    H3
    5    HC    H4  ;  0.000
[ pairs ]
;  ai   aj
    1    2
    1    3
    1    4
    1    5"""

        graph = convert_from(itp, IOType.ITP)
        self.assertDictEqual(graph.graph, GRAPH_ATTRIBUTES)
        self.assertListEqual(list(graph.nodes.data()), GRAPH_NODES)
        self.assertListEqual(list(graph.edges.data()), GRAPH_EDGES)

    def test_babel_rdkit(self):
        from rdkit import Chem
        rdmol = Chem.RWMol()

        rdmol.SetDoubleProp('group_charge_0', 0.0)

        c = Chem.Atom('C')
        h1 = Chem.Atom('H')
        h2 = Chem.Atom('H')
        h3 = Chem.Atom('H')
        h4 = Chem.Atom('H')

        h1.SetProp('label', 'H1')
        h2.SetProp('label', 'H2')
        h3.SetProp('label', 'H3')
        h4.SetProp('label', 'H4')

        h1.SetProp('atom_type', 'HC')
        h2.SetProp('atom_type', 'HC')
        h3.SetProp('atom_type', 'HC')
        h4.SetProp('atom_type', 'HC')

        idx = dict()
        idx[c] = rdmol.AddAtom(c)
        idx[h1] = rdmol.AddAtom(h1)
        idx[h2] = rdmol.AddAtom(h2)
        idx[h3] = rdmol.AddAtom(h3)
        idx[h4] = rdmol.AddAtom(h4)

        rdmol.AddBond(idx[c], idx[h1])
        rdmol.AddBond(idx[c], idx[h2])
        rdmol.AddBond(idx[c], idx[h3])
        rdmol.AddBond(idx[c], idx[h4])

        rdkit_graph_nodes = [(v-1, data) for v, data in GRAPH_NODES]
        rdkit_graph_edges = \
            [(u - 1, v - 1, {**{'rdkit_bond_type': Chem.BondType.UNSPECIFIED}, **data}) for u, v, data in GRAPH_EDGES]

        graph = convert_from(rdmol, IOType.RDKIT)

        self.assertDictEqual(graph.graph, GRAPH_ATTRIBUTES)
        self.assertListEqual(list(graph.nodes.data()), rdkit_graph_nodes)
        self.assertListEqual(list(graph.edges.data()), rdkit_graph_edges)

        graph2 = convert_from(convert_to(graph, IOType.RDKIT), IOType.RDKIT)
        self.assertListEqual(list(graph.nodes.data()), list(graph2.nodes.data()))
        self.assertListEqual(list(graph.edges.data()), list(graph2.edges.data()))

    # TODO add babel tests
    # TODO add charging tests

if __name__ == '__main__':
    unittest.main()