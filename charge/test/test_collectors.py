from math import log

import pytest

from charge.collectors import AssignmentError, HistogramCollector, MeanCollector, CachingCollector


def test_mean_collector(ref_graph, mock_repository):
    collector = MeanCollector(mock_repository, 2)

    means = collector.collect_values(ref_graph, False, [3, 2, 1, 0])

    assert means[1][:2] == ([0.34], [1.0])
    assert means[2][:2] == ([0.34], [1.0])
    assert means[3][:2] == ([0.34], [1.0])
    assert means[4][:2] == ([0.34], [1.0])
    assert means[5][:2] == ([0.34], [1.0])


def test_mean_collector_elem_fallback(ref_graph, mock_elem_repository):
    collector = MeanCollector(mock_elem_repository, 2)

    means = collector.collect_values(ref_graph, False, [3, 2, 1, 0])

    assert means[1][:2] == ([0.2], [1.0])
    assert means[2][:2] == ([0.34], [1.0])
    assert means[3][:2] == ([0.34], [1.0])
    assert means[4][:2] == ([0.34], [1.0])
    assert means[5][:2] == ([0.34], [1.0])


def test_iacm_data_only(ref_graph, mock_elem_repository):
    collector = MeanCollector(mock_elem_repository, 2)

    with pytest.raises(AssignmentError):
        collector.collect_values(ref_graph, True, [3, 2, 1, 0])


def test_shell_size_order(ref_graph, mock_order_repository):
    collector = MeanCollector(mock_order_repository, 2)

    means = collector.collect_values(ref_graph, False, [3, 2, 1, 0])
    assert means[1][:2] == ([0.2], [1.0])

    means = collector.collect_values(ref_graph, False, [0, 1, 2, 3])
    assert means[1][:2] == ([0.34], [1.0])


def test_histogram_collector(ref_graph, mock_repository):
    collector = HistogramCollector(mock_repository, 2)

    means = collector.collect_values(ref_graph, False, [3, 2, 1, 0])

    assert means[1][0] == pytest.approx([0.31, 0.46])
    assert means[1][1] == pytest.approx([log(2.0), log(1.0)])
    assert means[2][0] == pytest.approx([0.31, 0.46])
    assert means[2][1] == pytest.approx([log(2.0), log(1.0)])
    assert means[3][0] == pytest.approx([0.31, 0.46])
    assert means[3][1] == pytest.approx([log(2.0), log(1.0)])
    assert means[4][0] == pytest.approx([0.31, 0.46])
    assert means[4][1] == pytest.approx([log(2.0), log(1.0)])
    assert means[5][0] == pytest.approx([0.31, 0.46])
    assert means[5][1] == pytest.approx([log(2.0), log(1.0)])


def test_histogram_collector_max_bins(ref_graph, mock_repository):
    collector = HistogramCollector(mock_repository, 2, max_bins=1)

    means = collector.collect_values(ref_graph, False, [3, 2, 1, 0])

    assert means[1][0] == pytest.approx([0.31])
    assert means[1][1] == pytest.approx([log(3.0)])
    assert means[2][0] == pytest.approx([0.31])
    assert means[2][1] == pytest.approx([log(3.0)])
    assert means[3][0] == pytest.approx([0.31])
    assert means[3][1] == pytest.approx([log(3.0)])
    assert means[4][0] == pytest.approx([0.31])
    assert means[4][1] == pytest.approx([log(3.0)])
    assert means[5][0] == pytest.approx([0.31])
    assert means[5][1] == pytest.approx([log(3.0)])


def test_caching_collector0(ref_graph, mock_repository):

    collector = CachingCollector(MeanCollector(mock_repository, 2))

    means0 = collector.collect_values(ref_graph, False, [3, 2, 1, 0])

    for shell, chdct in mock_repository.charges_iacm.items():
        for key in chdct.keys():
            chdct[key].append(1)

    means1 = collector.collect_values(ref_graph, False, [3, 2, 1, 0])

    assert means0 == means1


def test_caching_collector1(ref_graph, mock_versioning_repository):

    collector = CachingCollector(MeanCollector(mock_versioning_repository, 2))

    means0 = collector.collect_values(ref_graph, False, [3, 2, 1, 0])

    for shell, chdct in mock_versioning_repository.charges_iacm.items():
        for key in chdct.keys():
            chdct[key].append(1)

    means1 = collector.collect_values(ref_graph, False, [3, 2, 1, 0])

    assert means0 != means1
