import pytest

import networkx as nx
import networkx.algorithms.isomorphism as nxiso

from charge.babel import convert_from, IOType, convert_to

def equal_attributes(att1, att2):
     return att1 == att2

def assert_same_graph(graph, ref_graph):
    assert graph.graph == ref_graph.graph
    assert graph.nodes == ref_graph.nodes
    assert graph.edges == ref_graph.edges
    assert nxiso.is_isomorphic(graph, ref_graph,
            node_match=equal_attributes, edge_match=equal_attributes)


def test_from_lgf(ref_graph_lgf, ref_graph):
    graph = convert_from(ref_graph_lgf, IOType.LGF)
    assert_same_graph(graph, ref_graph)


def test_to_lgf(ref_graph, ref_graph_lgf):
    graph_lgf = convert_to(ref_graph, IOType.LGF)
    assert graph_lgf == ref_graph_lgf


def test_from_gml(ref_graph_gml, ref_graph):
    graph = convert_from(ref_graph_gml, IOType.GML)
    assert_same_graph(graph, ref_graph)


def test_to_gml(ref_graph, ref_graph_gml):
    graph_gml = convert_to(ref_graph, IOType.GML)
    assert graph_gml == ref_graph_gml


def test_from_itp(ref_graph_itp, ref_graph):
    graph = convert_from(ref_graph_itp, IOType.ITP)
    assert_same_graph(graph, ref_graph)


def test_to_itp(ref_graph, ref_graph_itp):
    # Not yet supported, check that it throws exception
    with pytest.raises(ValueError):
        graph_itp = convert_to(ref_graph, IOType.ITP)


def test_from_rdkit(ref_graph_rdkit, ref_graph_shifted):
    graph = convert_from(ref_graph_rdkit, IOType.RDKIT)
    assert_same_graph(graph, ref_graph_shifted)


def test_to_rdkit(ref_graph, ref_graph_rdkit):
    def assert_rdkit_seq_equal(seq, seq_ref):
        # works for both atoms and bonds
        # docs don't say what Match means though...
        for i, obj in enumerate(seq):
            assert obj.Match(seq_ref[i])

    graph_rdkit = convert_to(ref_graph, IOType.RDKIT)
    assert_rdkit_seq_equal(graph_rdkit.GetAtoms(), ref_graph_rdkit.GetAtoms())
    assert_rdkit_seq_equal(graph_rdkit.GetBonds(), ref_graph_rdkit.GetBonds())
    assert graph_rdkit.GetPropsAsDict() == ref_graph_rdkit.GetPropsAsDict()
