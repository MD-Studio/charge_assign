from io import BytesIO

import pytest

from charge.repository import Repository, _VersioningList


def test_create_empty():
    repo = Repository()
    assert repo._Repository__min_shell == 1
    assert repo._Repository__max_shell == 7
    assert repo._Repository__traceable == False
    assert repo._Repository__versioning == False
    assert len(repo.charges_iacm) == 0
    assert len(repo.charges_elem) == 0
    assert not hasattr(repo, 'iso_iacm')
    assert not hasattr(repo, 'iso_elem')


def test_create_empty_traceable():
    repo = Repository(traceable=True)
    assert repo._Repository__min_shell == 1
    assert repo._Repository__max_shell == 7
    assert repo._Repository__traceable == True
    assert repo._Repository__versioning == False
    assert len(repo.charges_iacm) == 0
    assert len(repo.charges_elem) == 0
    assert len(repo.iso_iacm) == 0
    assert len(repo.iso_elem) == 0


def test_create_empty_versioning():
    repo = Repository(versioning=True)
    assert repo._Repository__min_shell == 1
    assert repo._Repository__max_shell == 7
    assert repo._Repository__traceable == False
    assert repo._Repository__versioning == True
    assert len(repo.charges_iacm) == 0
    assert len(repo.charges_elem) == 0
    assert not hasattr(repo, 'iso_iacm')
    assert not hasattr(repo, 'iso_elem')


def test_set_shell_sizes():
    repo = Repository(2, 4)
    assert repo._Repository__min_shell == 2
    assert repo._Repository__max_shell == 4


def test_create_from_dir(lgf_data_dir):
    repo = Repository.create_from(str(lgf_data_dir))
    assert repo._Repository__min_shell == 1
    assert repo._Repository__max_shell == 7
    assert repo._Repository__traceable == False
    assert repo._Repository__versioning == False

    assert len(repo.charges_iacm) == 7
    assert len(repo.charges_iacm[1]) == 10
    assert len(repo.charges_iacm[2]) == 14
    assert len(repo.charges_iacm[3]) == 15

    assert repo.charges_iacm[1]['c18208da9e290c6faf8a0c58017d24d9'] == [-0.516]
    assert repo.charges_iacm[3]['76198d87470cc1b2f871da60449146bc'] == [0.077, 0.077, 0.077]

    assert len(repo.charges_elem) == 7
    assert len(repo.charges_elem[1]) == 9
    assert len(repo.charges_elem[2]) == 14
    assert len(repo.charges_elem[3]) == 15

    assert repo.charges_elem[1]['92ed00c54b2190be94748bee34b22847'] == [-0.516]
    assert repo.charges_elem[3]['17ac3199bf634022485c145821f358d5'] == [0.077, 0.077, 0.077]

    assert not hasattr(repo, 'iso_iacm')
    assert not hasattr(repo, 'iso_elem')


def test_create_traceable_from_dir(lgf_data_dir):
    repo = Repository.create_from(str(lgf_data_dir), traceable=True)
    assert repo._Repository__min_shell == 1
    assert repo._Repository__max_shell == 7
    assert repo._Repository__traceable == True
    assert repo._Repository__versioning == False

    assert len(repo.charges_iacm) == 7
    assert len(repo.charges_iacm[1]) == 10
    assert len(repo.charges_iacm[2]) == 14
    assert len(repo.charges_iacm[3]) == 15

    assert repo.charges_iacm[1]['c18208da9e290c6faf8a0c58017d24d9'] == [(-0.516, 15610, 2)]
    assert repo.charges_iacm[3]['76198d87470cc1b2f871da60449146bc'] == [
            (0.077, 1204, 1), (0.077, 1204, 4), (0.077, 1204, 5)]

    assert len(repo.charges_elem) == 7
    assert len(repo.charges_elem[1]) == 9
    assert len(repo.charges_elem[2]) == 14
    assert len(repo.charges_elem[3]) == 15

    assert repo.charges_elem[1]['92ed00c54b2190be94748bee34b22847'] == [(-0.516, 15610, 2)]
    assert repo.charges_elem[3]['17ac3199bf634022485c145821f358d5'] == [
            (0.077, 1204, 1), (0.077, 1204, 4), (0.077, 1204, 5)]

    assert hasattr(repo, 'iso_iacm')
    assert hasattr(repo, 'iso_elem')
    assert len(repo.iso_iacm) == 2
    assert len(repo.iso_elem) == 2


def test_create_versioning_from_dir(lgf_data_dir):
    repo = Repository.create_from(str(lgf_data_dir), versioning=True)
    assert repo._Repository__min_shell == 1
    assert repo._Repository__max_shell == 7
    assert repo._Repository__traceable == False
    assert repo._Repository__versioning == True

    assert len(repo.charges_iacm) == 7
    assert len(repo.charges_iacm[1]) == 10
    assert len(repo.charges_iacm[2]) == 14
    assert len(repo.charges_iacm[3]) == 15

    assert isinstance(repo.charges_iacm[1]['c18208da9e290c6faf8a0c58017d24d9'], _VersioningList)
    assert isinstance(repo.charges_iacm[3]['76198d87470cc1b2f871da60449146bc'], _VersioningList)

    assert repo.charges_iacm[1]['c18208da9e290c6faf8a0c58017d24d9'] == [-0.516]
    assert repo.charges_iacm[3]['76198d87470cc1b2f871da60449146bc'] == [0.077, 0.077, 0.077]

    assert len(repo.charges_elem) == 7
    assert len(repo.charges_elem[1]) == 9
    assert len(repo.charges_elem[2]) == 14
    assert len(repo.charges_elem[3]) == 15

    assert isinstance(repo.charges_elem[1]['92ed00c54b2190be94748bee34b22847'], _VersioningList)
    assert isinstance(repo.charges_elem[3]['17ac3199bf634022485c145821f358d5'], _VersioningList)

    assert repo.charges_elem[1]['92ed00c54b2190be94748bee34b22847'] == [-0.516]
    assert repo.charges_elem[3]['17ac3199bf634022485c145821f358d5'] == [0.077, 0.077, 0.077]

    assert not hasattr(repo, 'iso_iacm')
    assert not hasattr(repo, 'iso_elem')

    ov0 = repo.charges_iacm[1]['c18208da9e290c6faf8a0c58017d24d9'].version
    repo.charges_iacm[1]['c18208da9e290c6faf8a0c58017d24d9'].append(1)
    ov1 = repo.charges_iacm[1]['c18208da9e290c6faf8a0c58017d24d9'].version
    repo.charges_iacm[1]['c18208da9e290c6faf8a0c58017d24d9'] += [2]
    ov2 = repo.charges_iacm[1]['c18208da9e290c6faf8a0c58017d24d9'].version
    del repo.charges_iacm[1]['c18208da9e290c6faf8a0c58017d24d9'][-1]
    ov3 = repo.charges_iacm[1]['c18208da9e290c6faf8a0c58017d24d9'].version
    repo.charges_iacm[1]['c18208da9e290c6faf8a0c58017d24d9'].remove(1)
    ov4 = repo.charges_iacm[1]['c18208da9e290c6faf8a0c58017d24d9'].version

    assert len({ov0, ov1, ov2, ov3, ov4}) == 5


def test_read_write(lgf_data_dir):
    repo0 = Repository.create_from(str(lgf_data_dir))

    tmp = BytesIO()
    repo0.write(tmp)

    repo1 = Repository.read(tmp)

    assert repo0._Repository__min_shell == repo1._Repository__min_shell
    assert repo0._Repository__max_shell == repo1._Repository__max_shell
    assert repo0._Repository__traceable == repo1._Repository__traceable
    assert repo0._Repository__versioning == repo1._Repository__versioning

    assert repo0.charges_iacm == repo1.charges_iacm
    assert repo0.charges_elem == repo1.charges_elem
    assert not hasattr(repo0, 'iso_iacm')
    assert not hasattr(repo0, 'iso_elem')
    assert not hasattr(repo1, 'iso_iacm')
    assert not hasattr(repo1, 'iso_elem')


def test_read_write_traceable(lgf_data_dir):
    repo0 = Repository.create_from(str(lgf_data_dir), traceable=True)

    tmp = BytesIO()
    repo0.write(tmp)

    repo1 = Repository.read(tmp)

    assert repo0._Repository__min_shell == repo1._Repository__min_shell
    assert repo0._Repository__max_shell == repo1._Repository__max_shell
    assert repo0._Repository__traceable == repo1._Repository__traceable
    assert repo0._Repository__versioning == repo1._Repository__versioning


    # msgpack writes tuples as lists. this is ok for our purposes, but we need to take care of that in testing!
    assert repo0.charges_iacm.keys() == repo1.charges_iacm.keys()
    for shell_size in repo0.charges_iacm.keys():
        assert repo0.charges_iacm[shell_size].keys() == repo1.charges_iacm[shell_size].keys()
        for key in repo0.charges_iacm[shell_size].keys():
            assert repo0.charges_iacm[shell_size][key] == [tuple(x) for x in repo0.charges_iacm[shell_size][key]]

    assert repo0.charges_elem.keys() == repo1.charges_elem.keys()
    for shell_size in repo0.charges_elem.keys():
        assert repo0.charges_elem[shell_size].keys() == repo1.charges_elem[shell_size].keys()
        for key in repo0.charges_elem[shell_size].keys():
            assert repo0.charges_elem[shell_size][key] == [tuple(x) for x in repo0.charges_elem[shell_size][key]]

    assert repo0.iso_iacm == repo1.iso_iacm
    assert repo0.iso_elem == repo1.iso_elem


def test_read_write_versioning(lgf_data_dir):
    repo0 = Repository.create_from(str(lgf_data_dir), versioning=True)

    tmp = BytesIO()
    repo0.write(tmp)

    repo1 = Repository.read(tmp, versioning=True)

    assert repo0._Repository__min_shell == repo1._Repository__min_shell
    assert repo0._Repository__max_shell == repo1._Repository__max_shell
    assert repo0._Repository__traceable == repo1._Repository__traceable
    assert repo0._Repository__versioning == repo1._Repository__versioning

    assert repo0.charges_iacm == repo1.charges_iacm
    assert repo0.charges_elem == repo1.charges_elem
    assert not hasattr(repo0, 'iso_iacm')
    assert not hasattr(repo0, 'iso_elem')
    assert not hasattr(repo1, 'iso_iacm')
    assert not hasattr(repo1, 'iso_elem')

    assert isinstance(repo1.charges_iacm[1]['c18208da9e290c6faf8a0c58017d24d9'], _VersioningList)
    assert isinstance(repo1.charges_iacm[3]['76198d87470cc1b2f871da60449146bc'], _VersioningList)
    assert isinstance(repo1.charges_elem[1]['92ed00c54b2190be94748bee34b22847'], _VersioningList)
    assert isinstance(repo1.charges_elem[3]['17ac3199bf634022485c145821f358d5'], _VersioningList)


def test_add_from(lgf_data_dir):
    repo = Repository.create_from(str(lgf_data_dir))
    repo.add_from(lgf_data_dir)

    assert repo._Repository__min_shell == 1
    assert repo._Repository__max_shell == 7
    assert repo._Repository__traceable == False
    assert repo._Repository__versioning == False

    assert len(repo.charges_iacm) == 7
    assert len(repo.charges_iacm[1]) == 10
    assert len(repo.charges_iacm[2]) == 14
    assert len(repo.charges_iacm[3]) == 15

    assert repo.charges_iacm[1]['c18208da9e290c6faf8a0c58017d24d9'] == [-0.516, -0.516]
    assert repo.charges_iacm[3]['76198d87470cc1b2f871da60449146bc'] == [0.077, 0.077, 0.077, 0.077, 0.077, 0.077]

    assert len(repo.charges_elem) == 7
    assert len(repo.charges_elem[1]) == 9
    assert len(repo.charges_elem[2]) == 14
    assert len(repo.charges_elem[3]) == 15

    assert repo.charges_elem[1]['92ed00c54b2190be94748bee34b22847'] == [-0.516, -0.516]
    assert repo.charges_elem[3]['17ac3199bf634022485c145821f358d5'] == [0.077, 0.077, 0.077, 0.077, 0.077, 0.077]

    assert not hasattr(repo, 'iso_iacm')
    assert not hasattr(repo, 'iso_elem')


def test_add_from_traceable(lgf_data_dir):
    repo = Repository.create_from(str(lgf_data_dir), traceable=True)

    with pytest.raises(ValueError):
        repo.add_from(lgf_data_dir)


def test_add_from_versioning(lgf_data_dir):
    repo = Repository.create_from(str(lgf_data_dir), versioning=True)
    repo.add_from(lgf_data_dir)

    assert repo._Repository__min_shell == 1
    assert repo._Repository__max_shell == 7
    assert repo._Repository__traceable == False
    assert repo._Repository__versioning == True

    assert len(repo.charges_iacm) == 7
    assert len(repo.charges_iacm[1]) == 10
    assert len(repo.charges_iacm[2]) == 14
    assert len(repo.charges_iacm[3]) == 15

    assert repo.charges_iacm[1]['c18208da9e290c6faf8a0c58017d24d9'] == [-0.516, -0.516]
    assert repo.charges_iacm[3]['76198d87470cc1b2f871da60449146bc'] == [0.077, 0.077, 0.077, 0.077, 0.077, 0.077]

    assert len(repo.charges_elem) == 7
    assert len(repo.charges_elem[1]) == 9
    assert len(repo.charges_elem[2]) == 14
    assert len(repo.charges_elem[3]) == 15

    assert repo.charges_elem[1]['92ed00c54b2190be94748bee34b22847'] == [-0.516, -0.516]
    assert repo.charges_elem[3]['17ac3199bf634022485c145821f358d5'] == [0.077, 0.077, 0.077, 0.077, 0.077, 0.077]

    assert not hasattr(repo, 'iso_iacm')
    assert not hasattr(repo, 'iso_elem')

    assert isinstance(repo.charges_iacm[1]['c18208da9e290c6faf8a0c58017d24d9'], _VersioningList)
    assert isinstance(repo.charges_iacm[3]['76198d87470cc1b2f871da60449146bc'], _VersioningList)
    assert isinstance(repo.charges_elem[1]['92ed00c54b2190be94748bee34b22847'], _VersioningList)
    assert isinstance(repo.charges_elem[3]['17ac3199bf634022485c145821f358d5'], _VersioningList)


def test_remove_from(lgf_data_dir):
    repo = Repository.create_from(str(lgf_data_dir))
    repo.remove_from(lgf_data_dir)

    assert repo._Repository__min_shell == 1
    assert repo._Repository__max_shell == 7
    assert repo._Repository__traceable == False
    assert repo._Repository__versioning == False

    assert len(repo.charges_iacm) == 0
    assert len(repo.charges_elem) == 0

    assert not hasattr(repo, 'iso_iacm')
    assert not hasattr(repo, 'iso_elem')


def test_remove_from_traceable(lgf_data_dir):
    repo = Repository.create_from(str(lgf_data_dir), traceable=True)

    with pytest.raises(ValueError):
        repo.remove_from(lgf_data_dir)
