import os

from charge.repository import Repository

if __name__ == '__main__':

    test_data_dir = os.path.realpath(
            os.path.join(__file__, '..', 'cross_validation_data'))

    out_file = 'cross_validation_repository.zip'

    repo = Repository.create_from(test_data_dir,
            min_shell=0, max_shell=3, traceable=True)

    repo.write(out_file)
