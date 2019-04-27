import pytest

from charge.solvers import CDPSolver, DPSolver, ILPSolver, SimpleSolver, SymmetricILPSolver, SymmetricDPSolver


def test_simple_solver(ref_graph):
    solver = SimpleSolver(2)
    charge_dists = {
            1: ([-0.516], [1.0]),
            2: ([0.129], [1.0]),
            3: ([0.129], [1.0]),
            4: ([0.129], [1.0]),
            5: ([0.129], [1.0])
            }
    solver.solve_partial_charges(
            ref_graph,
            charge_dists,
            0)

    assert ref_graph.node[1]['partial_charge'] == -0.516
    assert ref_graph.node[1]['score'] == 1.0
    assert ref_graph.node[2]['partial_charge'] == 0.129
    assert ref_graph.node[2]['score'] == 1.0
    assert ref_graph.node[3]['partial_charge'] == 0.129
    assert ref_graph.node[3]['score'] == 1.0
    assert ref_graph.node[4]['partial_charge'] == 0.129
    assert ref_graph.node[4]['score'] == 1.0
    assert ref_graph.node[5]['partial_charge'] == 0.129
    assert ref_graph.node[5]['score'] == 1.0

    assert ref_graph.graph['time'] < 0.1
    assert ref_graph.graph['total_charge'] == pytest.approx(0.0)
    assert ref_graph.graph['score'] == pytest.approx(5.0)


def test_ilp_solver(ref_graph):
    solver = ILPSolver(2)
    charge_dists = {
            1: ([-0.62, -0.52, -0.42], [0.25, 0.5, 0.25]),
            2: ([0.11, 0.13, 0.15], [0.1, 0.8, 0.1]),
            3: ([0.11, 0.13, 0.15], [0.1, 0.8, 0.1]),
            4: ([0.11, 0.13, 0.15], [0.1, 0.8, 0.1]),
            5: ([0.11, 0.13, 0.15], [0.1, 0.8, 0.1])
            }
    solver.solve_partial_charges(
            ref_graph,
            charge_dists,
            0)

    assert ref_graph.node[1]['partial_charge'] == -0.52
    assert ref_graph.node[1]['score'] == 0.5
    assert ref_graph.node[2]['partial_charge'] == 0.13
    assert ref_graph.node[2]['score'] == 0.8
    assert ref_graph.node[3]['partial_charge'] == 0.13
    assert ref_graph.node[3]['score'] == 0.8
    assert ref_graph.node[4]['partial_charge'] == 0.13
    assert ref_graph.node[4]['score'] == 0.8
    assert ref_graph.node[5]['partial_charge'] == 0.13
    assert ref_graph.node[5]['score'] == 0.8

    assert ref_graph.graph['time'] < 0.1
    assert ref_graph.graph['total_charge'] == pytest.approx(0.0)
    assert ref_graph.graph['score'] == pytest.approx(3.7)

def test_symmetricilp_solver(ref_graph):
    solver = SymmetricILPSolver(2)
    charge_dists = {
            1: ([-0.62, -0.52, -0.42], [0.25, 0.5, 0.25]),
            2: ([0.12, 0.13, 0.14, 0.15], [0.25, 0.25, 0.25, 0.25]),
            3: ([0.12, 0.13, 0.14, 0.15], [0.25, 0.25, 0.25, 0.25]),
            4: ([0.12, 0.13, 0.14, 0.15], [0.25, 0.25, 0.25, 0.25]),
            5: ([0.12, 0.13, 0.14, 0.15], [0.25, 0.25, 0.25, 0.25])
            }
    solver.solve_partial_charges(
            ref_graph,
            charge_dists,
            0,
            shells=[2, 1, 0])

    assert ref_graph.node[1]['partial_charge'] == -0.52
    assert ref_graph.node[1]['score'] == 0.5
    assert ref_graph.node[2]['partial_charge'] == 0.13
    assert ref_graph.node[2]['score'] == 0.25
    assert ref_graph.node[3]['partial_charge'] == 0.13
    assert ref_graph.node[3]['score'] == 0.25
    assert ref_graph.node[4]['partial_charge'] == 0.13
    assert ref_graph.node[4]['score'] == 0.25
    assert ref_graph.node[5]['partial_charge'] == 0.13
    assert ref_graph.node[5]['score'] == 0.25

    assert ref_graph.graph['time'] < 0.1
    assert ref_graph.graph['total_charge'] == pytest.approx(0.0)
    assert ref_graph.graph['score'] == pytest.approx(1.5)


def test_dp_solver(ref_graph):
    solver = DPSolver(2)
    charge_dists = {
            1: ([-0.62, -0.52, -0.42], [0.25, 0.5, 0.25]),
            2: ([0.11, 0.13, 0.15], [0.1, 0.8, 0.1]),
            3: ([0.11, 0.13, 0.15], [0.1, 0.8, 0.1]),
            4: ([0.11, 0.13, 0.15], [0.1, 0.8, 0.1]),
            5: ([0.11, 0.13, 0.15], [0.1, 0.8, 0.1])
            }
    solver.solve_partial_charges(
            ref_graph,
            charge_dists,
            0)

    assert ref_graph.node[1]['partial_charge'] == -0.52
    assert ref_graph.node[1]['score'] == 0.5
    assert ref_graph.node[2]['partial_charge'] == 0.13
    assert ref_graph.node[2]['score'] == 0.8
    assert ref_graph.node[3]['partial_charge'] == 0.13
    assert ref_graph.node[3]['score'] == 0.8
    assert ref_graph.node[4]['partial_charge'] == 0.13
    assert ref_graph.node[4]['score'] == 0.8
    assert ref_graph.node[5]['partial_charge'] == 0.13
    assert ref_graph.node[5]['score'] == 0.8

    assert ref_graph.graph['time'] < 0.1
    assert ref_graph.graph['total_charge'] == pytest.approx(0.0)
    assert ref_graph.graph['score'] == pytest.approx(3.7)

def test_symmetric_dp_solver(ref_graph):
    solver = SymmetricDPSolver(2)
    charge_dists = {
        1: ([-0.62, -0.52, -0.42], [0.25, 0.5, 0.25]),
        2: ([0.12, 0.13, 0.14, 0.15], [0.25, 0.25, 0.25, 0.25]),
        3: ([0.12, 0.13, 0.14, 0.15], [0.25, 0.25, 0.25, 0.25]),
        4: ([0.12, 0.13, 0.14, 0.15], [0.25, 0.25, 0.25, 0.25]),
        5: ([0.12, 0.13, 0.14, 0.15], [0.25, 0.25, 0.25, 0.25])
    }
    solver.solve_partial_charges(
            ref_graph,
            charge_dists,
            0)

    assert ref_graph.node[1]['partial_charge'] == -0.52
    assert ref_graph.node[1]['score'] == 0.5
    assert ref_graph.node[2]['partial_charge'] == 0.13
    assert ref_graph.node[2]['score'] == 0.25
    assert ref_graph.node[3]['partial_charge'] == 0.13
    assert ref_graph.node[3]['score'] == 0.25
    assert ref_graph.node[4]['partial_charge'] == 0.13
    assert ref_graph.node[4]['score'] == 0.25
    assert ref_graph.node[5]['partial_charge'] == 0.13
    assert ref_graph.node[5]['score'] == 0.25

    assert ref_graph.graph['time'] < 0.1
    assert ref_graph.graph['total_charge'] == pytest.approx(0.0)
    assert ref_graph.graph['score'] == pytest.approx(1.5)


def test_cdp_solver(ref_graph):
    solver = CDPSolver(2)
    charge_dists = {
            1: ([-0.62, -0.52, -0.42], [0.25, 0.5, 0.25]),
            2: ([0.11, 0.13, 0.15], [0.1, 0.8, 0.1]),
            3: ([0.11, 0.13, 0.15], [0.1, 0.8, 0.1]),
            4: ([0.11, 0.13, 0.15], [0.1, 0.8, 0.1]),
            5: ([0.11, 0.13, 0.15], [0.1, 0.8, 0.1])
            }
    solver.solve_partial_charges(
            ref_graph,
            charge_dists,
            0)

    assert ref_graph.node[1]['partial_charge'] == -0.52
    assert ref_graph.node[1]['score'] == 0.5
    assert ref_graph.node[2]['partial_charge'] == 0.13
    assert ref_graph.node[2]['score'] == 0.8
    assert ref_graph.node[3]['partial_charge'] == 0.13
    assert ref_graph.node[3]['score'] == 0.8
    assert ref_graph.node[4]['partial_charge'] == 0.13
    assert ref_graph.node[4]['score'] == 0.8
    assert ref_graph.node[5]['partial_charge'] == 0.13
    assert ref_graph.node[5]['score'] == 0.8

    assert ref_graph.graph['time'] < 0.1
    assert ref_graph.graph['total_charge'] == pytest.approx(0.0)
    assert ref_graph.graph['score'] == pytest.approx(3.7)

def test_atom_neighborhood_class():
    solver = SimpleSolver(2)
    atom_idx = {0: 1, 1: 2, 2: 3, 3: 4, 4: 5, 5:6, 6:7}
    keydict = {1: 'A', 2: 'B', 3: 'B', 4: 'B', 5: 'B', 6: 'C', 7: 'C'}

    classes = solver.compute_atom_neighborhood_classes(atom_idx, keydict)

    assert classes[0] == [0]
    assert classes[1] == [1, 2, 3, 4]
    assert classes[2] == [5, 6]
