from pathlib import Path
import pytest
import stat

from charge.nauty import Nauty

def test_create():
    nauty = Nauty()
    assert nauty.exe != ''
    nauty_path = Path(nauty.exe)
    assert nauty_path.exists()
    assert nauty_path.is_file()
    assert nauty_path.stat().st_mode & stat.S_IXUSR
    assert not nauty._Nauty__process.poll()


def test_canonize_neighborhood_graph_1(nauty, ref_graph):
    print(list(ref_graph.nodes(data=True)))
    key = nauty.canonize_neighborhood(ref_graph, 1, 0)
    ref_key = nauty._Nauty__make_hash(
            [(True, 'C')],
            [])
    assert key == ref_key

def test_canonize_neighborhood_graph_2(nauty, ref_graph):
    key = nauty.canonize_neighborhood(ref_graph, 1, 1)
    ref_key = nauty._Nauty__make_hash(
            [(False, 'C'), (False, 'HC'), (False, 'HC'), (False, 'HC'), (True, 'HC')],
            [(0, 1), (0, 2), (0, 3), (0, 4), (1, 0), (1, 2), (1, 3), (1, 4)])

def test_canonize_neighborhood_same_graph(nauty, ref_graph, ref_graph2):
    """Tests different representations of the same graph."""
    key = nauty.canonize_neighborhood(ref_graph, 2, 1)
    key2 = nauty.canonize_neighborhood(ref_graph2, 3, 1)
    assert key == key2

def test_canonize_neighborhood_different_core(nauty, ref_graph):
    # C atom in methane
    key = nauty.canonize_neighborhood(ref_graph, 1, 3)
    # H atom in methane
    key2 = nauty.canonize_neighborhood(ref_graph, 2, 3)
    assert key != key2

    # different H atom in methane, same due to symmetry
    key3 = nauty.canonize_neighborhood(ref_graph, 3, 3)
    assert key2 == key3

def test_make_nauty_input(nauty, ref_graph):
    colors = map(lambda node: node[1]['atom_type'], ref_graph.nodes(data=True))
    nauty_input = nauty._Nauty__make_nauty_input(ref_graph, colors)
    assert 'n=5' in nauty_input
    assert 'g 0:1;0:2;0:3;0:4' in nauty_input
    assert 'f=[0|1,2,3,4]' in nauty_input
    assert 'cxb' in nauty_input

    colors = ['x', 'y', 'y', 'y', 'y']
    nauty_input = nauty._Nauty__make_nauty_input(ref_graph, colors)
    assert 'n=5' in nauty_input
    assert 'g 0:1;0:2;0:3;0:4' in nauty_input
    assert 'f=[0|1,2,3,4]' in nauty_input
    assert 'cxb' in nauty_input

    colors = ['x', 'x', 'y', 'y', 'y']
    nauty_input = nauty._Nauty__make_nauty_input(ref_graph, colors)
    assert 'n=5' in nauty_input
    assert 'g 0:1;0:2;0:3;0:4' in nauty_input
    assert 'f=[0,1|2,3,4]' in nauty_input
    assert 'cxb' in nauty_input

def test_make_nauty_edges(nauty, ref_graph):
    to_nauty_id = { 1: 3, 2: 2, 3: 4, 4: 1, 5: 5 }
    nauty_edges = nauty._Nauty__make_nauty_edges(
            ref_graph.edges(), to_nauty_id)
    assert nauty_edges == [(3, 1), (3, 2), (3, 4), (3, 5)]
