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


def test_canonize_neighborhood_1(nauty, ref_graph):
    print(list(ref_graph.nodes(data=True)))
    key = nauty.canonize_neighborhood(ref_graph, 1, 0)
    ref_key = nauty._Nauty__make_hash(
            [(True, 'C')],
            [])
    assert key == ref_key

def test_canonize_neighborhood_2(nauty, ref_graph):
    key = nauty.canonize_neighborhood(ref_graph, 1, 1)
    ref_key = nauty._Nauty__make_hash(
            [(False, 'C'), (False, 'HC'), (False, 'HC'), (False, 'HC'), (True, 'HC')],
            [(0, 1), (0, 2), (0, 3), (0, 4), (1, 0), (1, 2), (1, 3), (1, 4)])

def test_canonize_neighborhood_3(nauty, ref_graph, ref_graph2):
    key = nauty.canonize_neighborhood(ref_graph, 2, 1)
    key2 = nauty.canonize_neighborhood(ref_graph2, 3, 1)
    assert key == key2
