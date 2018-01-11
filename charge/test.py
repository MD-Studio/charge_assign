import unittest

from charge.babel import convert_from, IOType, convert_to


class TestCharger(unittest.TestCase):

    def test_babel_lgf(self):
        lgf = """@nodes
        label   label2  atomType    prop    
        1   H1  20   1   
        2   H2  20   2   
        3   O   0   3   
        @edges
                label   prop    
        1   3   0   0
        2   3   1   1"""
        graph = convert_from(lgf, IOType.LGF)
        self.assertListEqual([(1, {'prop': '1', 'atom_type': 'H', 'label': 'H1', 'idx': 0}),
                              (2, {'prop': '2', 'atom_type': 'H', 'label': 'H2', 'idx': 1}),
                              (3, {'prop': '3', 'atom_type': 'O', 'label': 'O', 'idx': 2})],
                      graph.nodes(data=True))
        self.assertListEqual([(1,3,{'prop': '0'}),
                              (2,3,{'prop': '1'})], graph.edges(data=True))
        graph2 = convert_from(convert_to(graph, IOType.LGF), IOType.LGF)
        self.assertListEqual(graph.nodes(data=True), graph2.nodes(data=True))
        self.assertListEqual(graph.edges(data=True), graph2.edges(data=True))


    # TODO add babel tests
    # TODO add charging tests

if __name__ == '__main__':
    unittest.main()