import os
import sys
import warnings

from charge.repository import Repository
from charge.validation import cross_validate_molecules


def print_report(charger, iacm, shell, report, num_warnings):
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


def cross_validate(charger, iacm, shell, test_data_dir, repo, bucket, num_buckets) -> None:
    num_warnings = 0
    with warnings.catch_warnings(record=True) as w:
        report = cross_validate_molecules(
                charger, iacm, test_data_dir, shell=shell,
                repo=repo, bucket=bucket, num_buckets=num_buckets)
        num_warnings = len(w)

    outfile = 'cross_validation_report_{}_{}_{}.json'.format(
            charger.lower(), bucket, num_buckets)
    with open(outfile, 'w') as f:
        f.write(report.as_json())

    print_report(charger, iacm, shell, report, num_warnings)


if __name__ == '__main__':
    test_data_dir = os.path.realpath(
            os.path.join(__file__, '..', 'cross_validation_data'))

    if len(sys.argv) != 3:
        print('Usage: {} bucket num_buckets'.format(sys.argv[0]))
        print('Got {} arguments.'.format(len(sys.argv) - 1))
        quit()

    bucket = int(sys.argv[1])
    num_buckets = int(sys.argv[2])

    repo_file = 'cross_validation_repository.zip'
    repo = Repository.read(repo_file)

    cross_validate('SimpleCharger', False, 1, test_data_dir, repo, bucket, num_buckets)
    #cross_validate('ILPCharger', False, 1, test_data_dir)
    #cross_validate('ILPCharger', True, 1, test_data_dir)
    #cross_validate('ILPCharger', False, 2, test_data_dir, repo, bucket, num_buckets)
    #cross_validate('ILPCharger', False, 2, test_data_dir, repo, 1, 2)
    #cross_validate('ILPCharger', True, 2, test_data_dir)
    #cross_validate('CDPCharger', True, 1, test_data_dir)
