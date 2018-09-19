from charge.solvers import CDPSolver, DPSolver, ILPSolver, SimpleSolver

import networkx as nx
import pytest


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
