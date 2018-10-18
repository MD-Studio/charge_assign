import os
import warnings

from charge.nauty import Nauty
from charge.validation import (cross_validate_molecules, _FilteredCharges,
                               _FilteredRepository, AtomReport, MoleculeReport, strip_molecule,
                               cross_validate_molecule)


def test_filtered_charges_1(mock_traceable_charges) -> None:
    fc = _FilteredCharges(mock_traceable_charges, [])
    assert fc['key1'] == [0.1, 0.1,0.2]
    assert fc['key2'] == [0.3, 0.4, 0.4, 0.5]


def test_filtered_charges_2(mock_traceable_charges) -> None:
    fc = _FilteredCharges(mock_traceable_charges, [1])
    assert fc['key1'] == [0.2]
    assert fc['key2'] == [0.4, 0.4, 0.5]


def test_filtered_charges_3(mock_traceable_charges) -> None:
    fc = _FilteredCharges(mock_traceable_charges, [2])
    assert fc['key1'] == [0.1, 0.1]
    assert fc['key2'] == [0.3, 0.5]


def test_filtered_charges_4(mock_traceable_charges) -> None:
    fc = _FilteredCharges(mock_traceable_charges, [1, 2])
    assert 'key1' not in fc
    assert fc['key2'] == [0.5]


def test_filtered_repo_1(mock_traceable_repository) -> None:
    fr = _FilteredRepository(mock_traceable_repository, 1)
    assert fr.charges_iacm[1]['key'] == [0.129, 0.130, 0.329]
    assert fr.charges_iacm[1]['c18208da9e290c6faf8a0c58017d24d9'] == [0.129, 0.130, 0.329]


def test_filtered_repo_2(mock_traceable_repository) -> None:
    fr = _FilteredRepository(mock_traceable_repository, 2)
    assert fr.charges_iacm[1]['key'] == [-0.516, 0.129]
    assert fr.charges_iacm[1]['c18208da9e290c6faf8a0c58017d24d9'] == [-0.516, -0.321, 0.129]


def test_filtered_repo_3(mock_traceable_repository) -> None:
    fr = _FilteredRepository(mock_traceable_repository, 3)
    assert fr.charges_iacm[1]['key'] == [0.129, 0.130, 0.329]
    assert fr.charges_iacm[1]['c18208da9e290c6faf8a0c58017d24d9'] == [0.129, 0.130, 0.329]


def test_atom_report_add_atom_error() -> None:
    report = AtomReport()

    report.add_atom_error(0.25)
    assert report.total_atoms == 1
    assert report.sum_abs_atom_err == 0.25
    assert report.sum_sq_atom_err == 0.0625
    assert report.atom_errors == [0.25]

    report.add_atom_error(0.75)
    assert report.total_atoms == 2
    assert report.sum_abs_atom_err == 1.0
    assert report.sum_sq_atom_err == 0.625
    assert report.atom_errors == [0.25, 0.75]


def test_molecule_report_add_total_error() -> None:
    report = MoleculeReport()

    report.add_total_charge_error(0, 0.25)
    assert report.total_mols == 1
    assert report.sum_abs_total_err == 0.25
    assert report.sum_sq_total_err == 0.0625
    assert report.total_charge_errors == [(0, 0.25)]

    report.add_total_charge_error(1, 0.75)
    assert report.total_mols == 2
    assert report.sum_abs_total_err == 1.0
    assert report.sum_sq_total_err == 0.625
    assert report.total_charge_errors == [(0, 0.25), (1, 0.75)]


def test_atom_report_calculations() -> None:
    report = AtomReport()
    report.total_atoms = 16

    report.sum_abs_atom_err = 0.5
    assert report.mean_abs_atom_err() == 1.0 / 32.0

    report.sum_sq_atom_err = 1.0
    assert report.mean_sq_atom_err() == 0.0625
    assert report.rms_atom_err() == 0.25


def test_molecule_report_calculations() -> None:
    report = MoleculeReport()
    report.total_mols = 4

    report.sum_abs_total_err = 2.0
    assert report.mean_abs_total_err() == 0.5

    report.sum_sq_total_err = 1.0
    assert report.mean_sq_total_err() == 0.25
    assert report.rms_total_err() == 0.5

    report.solver_stats = [(0.0, 0.0, 0.125, 0.0, 0.0),
            (0.0, 0.0, 0.250, 0.0, 0.0), (0.0, 0.0, 0.375, 0.0, 0.0),
            (0.0, 0.0, 0.250, 0.0, 0.0)]
    assert report.mean_time() == 0.25


def test_atom_report_aggregation() -> None:
    report1 = AtomReport()
    report2 = AtomReport()

    report1.add_atom_error(0.25)
    report2.add_atom_error(0.75)

    report3 = report1 + report2
    assert report3.total_atoms == 2
    assert report3.sum_abs_atom_err == 1.0
    assert report3.sum_sq_atom_err == 0.625
    assert report3.atom_errors == [0.25, 0.75]

    report1 += report2
    assert report1.total_atoms == 2
    assert report1.sum_abs_atom_err == 1.0
    assert report1.sum_sq_atom_err == 0.625
    assert report1.atom_errors == [0.25, 0.75]


def test_molecule_report_aggregation() -> None:
    report1 = MoleculeReport()
    report2 = MoleculeReport()

    report1.add_total_charge_error(0, 0.25)
    report2.add_total_charge_error(1, 0.75)

    report3 = report1 + report2
    assert report3.total_mols == 2
    assert report3.sum_abs_total_err == 1.0
    assert report3.sum_sq_total_err == 0.625
    assert report3.total_charge_errors == [(0, 0.25), (1, 0.75)]

    report1 += report2
    assert report1.total_mols == 2
    assert report1.sum_abs_total_err == 1.0
    assert report1.sum_sq_total_err == 0.625
    assert report1.total_charge_errors == [(0, 0.25), (1, 0.75)]


def test_strip_molecule(ref_graph_charged) -> None:
    for atom, data in ref_graph_charged.nodes(data=True):
        assert 'partial_charge' in data
        assert 'iacm' in data
    test_graph = strip_molecule(ref_graph_charged, True)
    for atom, data in test_graph.nodes(data=True):
        assert 'partial_charge' not in data
        assert 'iacm' in data


def test_strip_molecule_iacm(ref_graph_charged) -> None:
    for atom, data in ref_graph_charged.nodes(data=True):
        assert 'partial_charge' in data
        assert 'iacm' in data
    test_graph = strip_molecule(ref_graph_charged, False)
    for atom, data in test_graph.nodes(data=True):
        assert 'partial_charge' not in data
        assert 'iacm' not in data


def test_cross_validate_molecule(mock_traceable_methane_repository, ref_graph_charged, charger_iacm):
    charger, iacm = charger_iacm
    nauty = Nauty()

    report = cross_validate_molecule(
            mock_traceable_methane_repository,
            1, ref_graph_charged, charger, [1], iacm, nauty)

    assert report.category('C').total_atoms == 1
    assert report.category('H').total_atoms == 4
    assert report.category('O').total_atoms == 0
    assert report.molecule.total_mols == 1
    assert report.category('C').sum_abs_atom_err == 0.0
    assert report.category('H').sum_abs_atom_err == 0.0
    assert report.molecule.sum_abs_total_err == 0.0


def test_cross_validate_molecules(lgf_data_dir):
    num_molecules = len([fn
                         for fn in os.listdir(str(lgf_data_dir))
                         if fn.endswith('.lgf')])

    with warnings.catch_warnings(record=True) as w:
        report = cross_validate_molecules(
                'MeanCharger', True, str(lgf_data_dir))
        all_warnings = w

    assert report.molecule.total_mols + len(all_warnings) == num_molecules
    assert report.category('C').mean_abs_atom_err() < 0.4
    assert report.category('H').mean_abs_atom_err() < 0.2
