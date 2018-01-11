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

    #TODO GML test
    #TODO rdkit test


#     def test_inchi(self):
#         charger = Charger()
#         # molid 8
#         rdmol = charger.charge_inchi('InChI=1S/C4H10/c1-4(2)3/h4H,1-3H3')
#
#         for atom in rdmol.GetAtoms():
#             names = atom.GetPropNames()
#             if len(names) > 0:
#                 self.assertIn('partial_charge', names)
#                 self.assertIn('uncertainty', names)
#
#
#     def test_rdmol(self):
#         charger = Charger()
#         # molid 8
#         rdmol = charger.charge_rdmol(Chem.MolFromInchi('InChI=1S/C4H10/c1-4(2)3/h4H,1-3H3'))
#
#         for atom in rdmol.GetAtoms():
#             names = atom.GetPropNames()
#             if len(names) > 0:
#                 self.assertIn('partial_charge', names)
#                 self.assertIn('uncertainty', names)
#
#
#     def test_smiles(self):
#         charger = Charger()
#         # molid 8
#         rdmol = charger.charge_smiles('CC(C)C')
#
#         for atom in rdmol.GetAtoms():
#             names = atom.GetPropNames()
#             if len(names) > 0:
#                 self.assertIn('partial_charge', names)
#                 self.assertIn('uncertainty', names)
#
#
#     def test_pdb(self):
#         charger = Charger()
#         # molid 8
#         rdmol = charger.charge_pdb("""
# HETATM    1   H7 _I08    0      -0.988  -1.950   0.257  1.00  0.00           H
# HETATM    2   C3 _I08    0      -0.075  -1.460  -0.105  1.00  0.00           C
# HETATM    3   H5 _I08    0      -0.078  -1.517  -1.202  1.00  0.00           H
# HETATM    4   H6 _I08    0       0.784  -2.040   0.256  1.00  0.00           H
# HETATM    5   C2 _I08    0       0.000   0.000   0.365  1.00  0.00           C
# HETATM    6   H4 _I08    0      -0.000   0.000   1.466  1.00  0.00           H
# HETATM    7   C1 _I08    0      -1.228   0.795  -0.105  1.00  0.00           C
# HETATM    8   H1 _I08    0      -2.159   0.343   0.260  1.00  0.00           H
# HETATM    9   H2 _I08    0      -1.193   1.831   0.254  1.00  0.00           H
# HETATM   10   H3 _I08    0      -1.277   0.823  -1.202  1.00  0.00           H
# HETATM   11   C4 _I08    0       1.302   0.665  -0.105  1.00  0.00           C
# HETATM   12   H8 _I08    0       1.354   0.689  -1.202  1.00  0.00           H
# HETATM   13   H9 _I08    0       1.375   1.699   0.255  1.00  0.00           H
# HETATM   14  H10 _I08    0       2.182   0.120   0.258  1.00  0.00           H
# CONECT    1    2
# CONECT    2    1    3    4    5
# CONECT    3    2
# CONECT    4    2
# CONECT    5    2    6    7   11
# CONECT    6    5
# CONECT    7    5    8    9   10
# CONECT    8    7
# CONECT    9    7
# CONECT   10    7
# CONECT   11    5   12   13   14
# CONECT   12   11
# CONECT   13   11
# CONECT   14   11
# END""")
#         for atom in rdmol.GetAtoms():
#             names = atom.GetPropNames()
#             if len(names) > 0:
#                 self.assertIn('partial_charge', names)
#                 self.assertIn('uncertainty', names)


if __name__ == '__main__':
    unittest.main()