from charge.chargers import CDPCharger, DPCharger, ILPCharger, MeanCharger, MedianCharger, ModeCharger

from math import log
import networkx as nx
import pytest


def test_mean_charger(mock_repository, ref_graph):
    charger = MeanCharger(mock_repository, 2)
    charger.charge(ref_graph, 0)

    assert ref_graph.node[1]['partial_charge'] == 0.34
    assert ref_graph.node[1]['score'] == 1.0
    assert ref_graph.node[2]['partial_charge'] == 0.34
    assert ref_graph.node[2]['score'] == 1.0
    assert ref_graph.node[3]['partial_charge'] == 0.34
    assert ref_graph.node[3]['score'] == 1.0
    assert ref_graph.node[4]['partial_charge'] == 0.34
    assert ref_graph.node[4]['score'] == 1.0
    assert ref_graph.node[5]['partial_charge'] == 0.34
    assert ref_graph.node[5]['score'] == 1.0

    assert ref_graph.graph['time'] < 0.1
    assert ref_graph.graph['total_charge'] == pytest.approx(1.7)
    assert ref_graph.graph['total_charge_redist'] == pytest.approx(0.0)
    assert ref_graph.graph['score'] == pytest.approx(5.0)


def test_iacmize(mock_repository, plain_ref_graph):
    charger = MeanCharger(mock_repository, 2)
    charger.charge(plain_ref_graph, 0, True)

    assert plain_ref_graph.node[1]['partial_charge'] == 0.34
    assert plain_ref_graph.node[1]['score'] == 1.0
    assert plain_ref_graph.node[2]['partial_charge'] == 0.34
    assert plain_ref_graph.node[2]['score'] == 1.0
    assert plain_ref_graph.node[3]['partial_charge'] == 0.34
    assert plain_ref_graph.node[3]['score'] == 1.0
    assert plain_ref_graph.node[4]['partial_charge'] == 0.34
    assert plain_ref_graph.node[4]['score'] == 1.0
    assert plain_ref_graph.node[5]['partial_charge'] == 0.34
    assert plain_ref_graph.node[5]['score'] == 1.0

    assert plain_ref_graph.graph['time'] < 0.1
    assert plain_ref_graph.graph['total_charge'] == pytest.approx(1.7)
    assert plain_ref_graph.graph['total_charge_redist'] == pytest.approx(0.0)
    assert plain_ref_graph.graph['score'] == pytest.approx(5.0)


def test_median_charger(mock_repository, ref_graph):
    charger = MedianCharger(mock_repository, 2)
    charger.charge(ref_graph, 0, True)

    assert ref_graph.node[1]['partial_charge'] == pytest.approx(0.31)
    assert ref_graph.node[1]['score'] == pytest.approx(1.0)
    assert ref_graph.node[2]['partial_charge'] == pytest.approx(0.31)
    assert ref_graph.node[2]['score'] == pytest.approx(1.0)
    assert ref_graph.node[3]['partial_charge'] == pytest.approx(0.31)
    assert ref_graph.node[3]['score'] == pytest.approx(1.0)
    assert ref_graph.node[4]['partial_charge'] == pytest.approx(0.31)
    assert ref_graph.node[4]['score'] == pytest.approx(1.0)
    assert ref_graph.node[5]['partial_charge'] == pytest.approx(0.31)
    assert ref_graph.node[5]['score'] == pytest.approx(1.0)

    assert ref_graph.graph['time'] < 0.1
    assert ref_graph.graph['total_charge'] == pytest.approx(1.55)
    assert ref_graph.graph['total_charge_redist'] == pytest.approx(0.0)
    assert ref_graph.graph['score'] == pytest.approx(5.0)


def test_mode_charger(mock_repository, ref_graph):
    charger = ModeCharger(mock_repository, 2)
    charger.charge(ref_graph, 0, True)

    assert ref_graph.node[1]['partial_charge'] == pytest.approx(0.31)
    assert ref_graph.node[1]['score'] == pytest.approx(2.0)
    assert ref_graph.node[2]['partial_charge'] == pytest.approx(0.31)
    assert ref_graph.node[2]['score'] == pytest.approx(2.0)
    assert ref_graph.node[3]['partial_charge'] == pytest.approx(0.31)
    assert ref_graph.node[3]['score'] == pytest.approx(2.0)
    assert ref_graph.node[4]['partial_charge'] == pytest.approx(0.31)
    assert ref_graph.node[4]['score'] == pytest.approx(2.0)
    assert ref_graph.node[5]['partial_charge'] == pytest.approx(0.31)
    assert ref_graph.node[5]['score'] == pytest.approx(2.0)

    assert ref_graph.graph['time'] < 0.1
    assert ref_graph.graph['total_charge'] == pytest.approx(1.55)
    assert ref_graph.graph['total_charge_redist'] == pytest.approx(0.0)
    assert ref_graph.graph['score'] == pytest.approx(10.0)


def test_ilp_charger(mock_methane_repository, ref_graph):
    charger = ILPCharger(mock_methane_repository, 2, 10)
    charger.charge(ref_graph, 0)

    assert ref_graph.node[1]['partial_charge'] == pytest.approx(-0.52)
    assert ref_graph.node[1]['score'] == pytest.approx(log(4))
    assert ref_graph.node[2]['partial_charge'] == pytest.approx(0.13)
    assert ref_graph.node[2]['score'] == pytest.approx(log(3))
    assert ref_graph.node[3]['partial_charge'] == pytest.approx(0.13)
    assert ref_graph.node[3]['score'] == pytest.approx(log(3))
    assert ref_graph.node[4]['partial_charge'] == pytest.approx(0.13)
    assert ref_graph.node[4]['score'] == pytest.approx(log(3))
    assert ref_graph.node[5]['partial_charge'] == pytest.approx(0.13)
    assert ref_graph.node[5]['score'] == pytest.approx(log(3))

    assert ref_graph.graph['time'] < 0.1
    assert ref_graph.graph['total_charge'] == pytest.approx(0.0)
    assert ref_graph.graph['total_charge_redist'] == pytest.approx(0.0)
    assert ref_graph.graph['score'] == pytest.approx(4 * log(3) + log(4))


def test_dp_charger(mock_methane_repository, ref_graph):
    charger = DPCharger(mock_methane_repository, 2)
    charger.charge(ref_graph, 0)

    assert ref_graph.node[1]['partial_charge'] == pytest.approx(-0.52)
    assert ref_graph.node[1]['score'] == pytest.approx(log(4))
    assert ref_graph.node[2]['partial_charge'] == pytest.approx(0.13)
    assert ref_graph.node[2]['score'] == pytest.approx(log(3))
    assert ref_graph.node[3]['partial_charge'] == pytest.approx(0.13)
    assert ref_graph.node[3]['score'] == pytest.approx(log(3))
    assert ref_graph.node[4]['partial_charge'] == pytest.approx(0.13)
    assert ref_graph.node[4]['score'] == pytest.approx(log(3))
    assert ref_graph.node[5]['partial_charge'] == pytest.approx(0.13)
    assert ref_graph.node[5]['score'] == pytest.approx(log(3))

    assert ref_graph.graph['time'] < 0.1
    assert ref_graph.graph['total_charge'] == pytest.approx(0.0)
    assert ref_graph.graph['total_charge_redist'] == pytest.approx(0.0)
    assert ref_graph.graph['score'] == pytest.approx(4 * log(3) + log(4))


def test_cdp_charger(mock_methane_repository, ref_graph):
    charger = CDPCharger(mock_methane_repository, 2)
    charger.charge(ref_graph, 0)

    assert ref_graph.node[1]['partial_charge'] == pytest.approx(-0.52)
    assert ref_graph.node[1]['score'] == pytest.approx(log(4))
    assert ref_graph.node[2]['partial_charge'] == pytest.approx(0.13)
    assert ref_graph.node[2]['score'] == pytest.approx(log(3))
    assert ref_graph.node[3]['partial_charge'] == pytest.approx(0.13)
    assert ref_graph.node[3]['score'] == pytest.approx(log(3))
    assert ref_graph.node[4]['partial_charge'] == pytest.approx(0.13)
    assert ref_graph.node[4]['score'] == pytest.approx(log(3))
    assert ref_graph.node[5]['partial_charge'] == pytest.approx(0.13)
    assert ref_graph.node[5]['score'] == pytest.approx(log(3))

    assert ref_graph.graph['time'] < 0.1
    assert ref_graph.graph['total_charge'] == pytest.approx(0.0)
    assert ref_graph.graph['total_charge_redist'] == pytest.approx(0.0)
    assert ref_graph.graph['score'] == pytest.approx(4 * log(3) + log(4))
