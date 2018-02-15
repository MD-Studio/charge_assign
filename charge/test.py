import unittest

from charge.babel import convert_from, IOType, convert_to, BondType


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
        self.assertDictEqual(graph.graph, {'charge_groups': {0: 0.0}})
        self.assertListEqual(list(graph.nodes.data()),
            [(1, {'atom_type': 'C', 'label': 'C1', 'charge_group': 0}),
             (2, {'atom_type': 'HC', 'label': 'H1', 'charge_group': 0}),
             (3, {'atom_type': 'HC', 'label': 'H2', 'charge_group': 0}),
             (4, {'atom_type': 'HC', 'label': 'H3', 'charge_group': 0}),
             (5, {'atom_type': 'HC', 'label': 'H4', 'charge_group': 0})])
        self.assertListEqual(list(graph.edges.data()),
            [(1, 2, {'bond_type': BondType.UNKNOWN}),
             (1, 3, {'bond_type': BondType.UNKNOWN}),
             (1, 4, {'bond_type': BondType.UNKNOWN}),
             (1, 5, {'bond_type': BondType.UNKNOWN})])

        graph2 = convert_from(convert_to(graph, IOType.LGF), IOType.LGF)
        self.assertListEqual(list(graph.nodes.data()), list(graph2.nodes.data()))
        self.assertListEqual(list(graph.edges.data()), list(graph2.edges.data()))


    # TODO add babel tests
    # TODO add charging tests

if __name__ == '__main__':
    unittest.main()