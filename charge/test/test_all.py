import msgpack
import unittest
import zipfile

from charge.repository import Repository


class TestAll(unittest.TestCase):
    '''Create a repository and store it, as a crude regression test
    while refactoring and adding unit tests.'''

    def test_all(self):
        # sample = 'sample'
        sample = 'larger_sample'
        # sample = 'qm_1_and_2_21-07-17'

        repository = Repository.create_from(data_location='../../data/{}'.format(sample))
        repository.write('test_all_repository.zip')

        def unpack_data(zipfile, name):
            return msgpack.unpackb(zipfile.read(name), raw=False)

        def compare(result_zip, reference_zip, name):
            result_data = unpack_data(result_zip, name)
            reference_data = unpack_data(reference_zip, name)
            if type(result_data) == list:
                result_data.sort()
                reference_data.sort()
            assert result_data == reference_data

        with zipfile.ZipFile('test_all_repository.zip', 'r') as result:
            with zipfile.ZipFile('test_all_{}_repo.zip.ref'.format(sample), 'r') as reference:
                compare(result, reference, 'meta')
                compare(result, reference, 'charges_iacm')
                compare(result, reference, 'charges_elem')
                compare(result, reference, 'iso_iacm')
                compare(result, reference, 'iso_elem')


if __name__ == '__main__':
    unittest.main()
