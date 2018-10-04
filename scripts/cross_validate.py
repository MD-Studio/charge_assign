import os
import warnings

from charge.validation import cross_validate_molecules


def cross_validate(charger, iacm, shell, test_data_dir) -> None:
    num_warnings = 0
    with warnings.catch_warnings(record=True) as w:
        report = cross_validate_molecules(
                charger, iacm, test_data_dir, shell=shell)
        num_warnings = len(w)

    print('{}, IACM: {}, shell: {}'.format(charger, iacm, shell))
    print('warnings: {}'.format(num_warnings))
    print('mols: {}'.format(report.molecule.total_mols))
    print('mean time: {}'.format(report.molecule.mean_time()))
    print('total mae: {}'.format(report.molecule.mean_abs_total_err()))
    print('total mse: {}, rmse: {}'.format(report.molecule.mean_sq_total_err(), report.molecule.rms_total_err()))

    for cat in ['C', 'H', 'N', 'O', 'P', 'S', 'Other']:
        atom_rep = report.category(cat)
        print('Cat: {}, atoms: {}'.format(cat, atom_rep.total_atoms))
        print('mae: {}'.format(atom_rep.mean_abs_atom_err()))
        print('mse: {}, rmse: {}'.format(atom_rep.mean_sq_atom_err(), atom_rep.rms_atom_err()))



if __name__ == '__main__':
    test_data_dir = os.path.realpath(
            os.path.join(__file__, '..', 'cross_validation_data'))

    #cross_validate('SimpleCharger', False, 1, test_data_dir)
    #cross_validate('ILPCharger', False, 1, test_data_dir)
    #cross_validate('ILPCharger', True, 1, test_data_dir)
    cross_validate('ILPCharger', False, 2, test_data_dir)
    cross_validate('ILPCharger', True, 2, test_data_dir)
    #cross_validate('CDPCharger', True, 1, test_data_dir)
