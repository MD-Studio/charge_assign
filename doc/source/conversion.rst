
Data Conversion
===============

A molecule in charge_assign is represented by a `networkx <http://networkx.github.io/>`_ graph.

Introduction
------------

:func:`convert_from() <charge.babel.convert_from>` converts from the supported chemical formats to a networkx graph and :func:`convert_to() <charge.babel.convert_to>` from a networkx graph to the supported chemical formats. The chemical formats can be either a str or an object. In the following, we will use this small example networkx graph::

        >> graph.graph
        {'group_charges': {0: 0.0}}
        >> graph.nodes.data()
        NodeDataView({1: {'atom_type': 'C', 'label': 'C1', 'charge_group': 0},\
            2: {'atom_type': 'HC', 'label': 'H1', 'charge_group': 0},\
            3: {'atom_type': 'HC', 'label': 'H2', 'charge_group': 0},\
            4: {'atom_type': 'HC', 'label': 'H3', 'charge_group': 0},\
            5: {'atom_type': 'HC', 'label': 'H4', 'charge_group': 0}})
        >> graph.edges.data()
        EdgeDataView([(1, 2, {'bond_type': <BondType.UNKNOWN: 'UNKNOWN'>}),\
            (1, 3, {'bond_type': <BondType.UNKNOWN: 'UNKNOWN'>}),\
            (1, 4, {'bond_type': <BondType.UNKNOWN: 'UNKNOWN'>}),\
            (1, 5, {'bond_type': <BondType.UNKNOWN: 'UNKNOWN'>})])


Supported Formats
-----------------

For each format, we list those node and edge attributes that are required and those that will be changed during conversion, any other attributes are preserved.

A word on atom types
^^^^^^^^^^^^^^^^^^^^

charge_assign has the option to convert basic elemental types to :data:`IACM types <charge.settings.IACM_ELEMENTS>` and saves it as the iacm node attribute. When converting from networkx to the chemical string formats (lgf, gml, itp), the IACM type takes precedence over the elemental type. When converting from networkx to chemical object formats (rdkit, openbabel, pybel) the atom_type is always converted to an elemental type.

lgf
^^^

`Lemon graph format <http://lemon.cs.elte.hu/pub/tutorial/a00018.html>`_ string.

+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| **required node attributes**                                                                                                                              |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| **lgf**                                                                     | **networkx**                                                                |
+=============================================================================+=============================================================================+
| label (int)                                                                 | node id (int)                                                               |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| atomType (int)                                                              | atom_type/iacm (str)                                                        |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+

+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| **optional node attributes**                                                                                                                              |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| **lgf**                                                                     | **networkx**                                                                |
+=============================================================================+=============================================================================+
| label2 (str, auto-generated if missing)                                     | label (str, auto-generated if missing)                                      |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| initColor (int)                                                             | charge_group (int, all atoms in group 0 if missing)                         |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+

+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| **optional edge attributes**                                                                                                                              |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| **lgf**                                                                     | **networkx**                                                                |
+=============================================================================+=============================================================================+
| bondType (int)                                                              | bond_type (:class:`BondType <charge.babel.BondType>`, UNKNOWN if missing)   |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| label (int, auto-generated)                                                 | None                                                                        |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+

The atom type is mapped from ints to strings according to the index (starting at 1) in :data:`IACM_ELEMENTS <charge.settings.IACM_ELEMENTS>`. group_charges are set to 0.0.

Bond types are mapped as:
    * 1, 2, 3: SINGLE, DOUBLE, TRIPLE
    * 4: AROMATIC
    * else: UNKNOWN

Bonds in the lgf file are indexed by the label attribute, which is ignored when converting from lgf and automatically added when converting to lgf.

Example::

    @nodes
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
    1       5       3

gml
^^^

`Graph modelling language <https://en.wikipedia.org/wiki/Graph_Modelling_Language>`_ string.

+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| **optional graph attributes**                                                                                                                             |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| **gml**                                                                     | **networkx**                                                                |
+=============================================================================+=============================================================================+
| groupchargei (float)                                                        | group_charges (map[int, float], 0.0 if missing)                             |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+

+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| **required node attributes**                                                                                                                              |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| **gml**                                                                     | **networkx**                                                                |
+=============================================================================+=============================================================================+
| id (int)                                                                    | node id (int)                                                               |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| atomtype (str)                                                              | atom_type/iacm (str)                                                        |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+

+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| **optional node attributes**                                                                                                                              |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| **gml**                                                                     | **networkx**                                                                |
+=============================================================================+=============================================================================+
| label (str, auto-generated if missing)                                      | label (str, auto-generated if missing)                                      |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| chargegroup (int)                                                           | charge_group (int, all atoms in group 0 if missing)                         |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| partialcharge (float)                                                       | partial_charge                                                              |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+

+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| **optional edge attributes**                                                                                                                              |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| **gml**                                                                     | **networkx**                                                                |
+=============================================================================+=============================================================================+
| bondtype (str)                                                              | bond_type (:class:`BondType <charge.babel.BondType>`, UNKNOWN if missing)   |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+

Allowed strings for the atom_type are the :data:`IACM_ELEMENTS <charge.settings.IACM_ELEMENTS>`. The bond_type string is mapped to its corresponding :class:`BondType <charge.babel.BondType>`. Missing label strings are auto-generated, missing bond_types get mapped to :class:`BondType <charge.babel.BondType>` UNKNOWN.

Example::

    graph [
        groupcharge0 0.0
        node [
            id 0
            label "C1"
            atomtype "C"
        ]
        node [
            id 1
            label "H1"
            atomtype "HC"
        ]
        node [
            id 2
            label "H2"
            atomtype "HC"
        ]
        node [
            id 3
            label "H3"
            atomtype "HC"
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
    ]

itp
^^^

`GROMOS include topology <http://manual.gromacs.org/current/online/itp.html>`_ string. **Only converting from itp to networkx is supported, not converting from networkx to itp.**

+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| **required node attributes**                                                                                                                              |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| **itp**                                                                     | **networkx**                                                                |
+=============================================================================+=============================================================================+
| nr (int)                                                                    | node id (int)                                                               |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| type (str)                                                                  | atom_type/iacm (str)                                                        |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+

+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| **optional node attributes**                                                                                                                              |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| **itp**                                                                     | **networkx**                                                                |
+=============================================================================+=============================================================================+
| atom (str)                                                                  | label (str, auto-generated if missing)                                      |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| charge (float)                                                              | partial_charge (float)                                                      |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| total_charge (float)                                                        | charge_group (int, all atoms in group 0 if missing)                         |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+

Allowed strings for the atom_type are the :data:`IACM_ELEMENTS <charge.settings.IACM_ELEMENTS>`. There is no bond type in the ITP file format, so all bonds will be of the :class:`BondType <charge.babel.BondType>` 'UNKNOWN'. The values of the total_charge attribute are used to determine the group_charges values.

Example::

    [ atoms ]
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
        1    5

rdkit
^^^^^

`rdkit <http://www.rdkit.org/>`_'s `rdchem.Mol <http://www.rdkit.org/Python_Docs/rdkit.Chem.rdchem.Mol-class.html>`_ object.

+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| **optional graph attributes**                                                                                                                             |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| **rdkit**                                                                   | **networkx**                                                                |
+=============================================================================+=============================================================================+
| `group_charge_i <http://www.rdkit.org/Python_Docs/rdkit.Chem.rdchem.Mol-    | group_charges (map[int, float], 0.0 if missing)                             |
| class.html#GetDoubleProp>`_ (float)                                         |                                                                             |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+

+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| **required node attributes**                                                                                                                              |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| **rdkit**                                                                   | **networkx**                                                                |
+=============================================================================+=============================================================================+
| `node idx <http://www.rdkit.org/Python_Docs/rdkit.Chem.rdchem.Atom-class    | node id (int)                                                               |
| .html#GetIdx>`_ (int)                                                       |                                                                             |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| `symbol <http://www.rdkit.org/Python_Docs/rdkit.Chem.rdchem.Atom-class.html | atom_type/iacm (str)                                                        |
| #GetSymbol>`_ (str) / `atom_type <http://www.rdkit.org/Python_Docs/         |                                                                             |
| rdkit.Chem.rdchem.Atom-class.html#GetProp>`_ (str)                          |                                                                             |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+

+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| **optional node attributes**                                                                                                                              |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| **rdkit**                                                                   | **networkx**                                                                |
+=============================================================================+=============================================================================+
| `label <http://www.rdkit.org/Python_Docs/rdkit.Chem.rdchem.Atom-class.html  | label (str, auto-generated if missing)                                      |
| #GetProp>`_ (str)                                                           |                                                                             |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| `charge_group <http://www.rdkit.org/Python_Docs/rdkit.Chem.rdchem.Atom-     | charge_group (int, all atoms in group 0 if missing)                         |
| class.html#GetIntProp>`_ (int)                                              |                                                                             |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+

+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| **required edge attributes**                                                                                                                              |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+
| **rdkit**                                                                   | **networkx**                                                                |
+=============================================================================+=============================================================================+
| `bond type <http://www.rdkit.org/Python_Docs/rdkit.Chem.rdchem.Bond-class   | bond_type (:class:`BondType <charge.babel.BondType>`)                       |
| .html#GetBondType>`_ (`rdchem.BondType <http://www.rdkit.org/Python_Docs/   |                                                                             |
| rdkit.Chem.rdchem.BondType-class.html>`_)                                   |                                                                             |
+-----------------------------------------------------------------------------+-----------------------------------------------------------------------------+

Missing label strings are auto-generated. SINGLE, DOUBLE, TRIPLE and AROMATIC bond types are mapped, all else are mapped to UNKNOWN. However, when converting from an rdkit molecule to a networkx graph, the original bond type is preserved in the rdkit_bond_type attribute, which is only exported when converting back to an rdkit molecule. RDKIT does not support :data:`IACM types <charge.settings.IACM_ELEMENTS>`. Therefore, when converting to rdkit, the IACM type is saved in the atom_type property, which takes precedence over the rdkit atom symbol when converting back from rdkit.



