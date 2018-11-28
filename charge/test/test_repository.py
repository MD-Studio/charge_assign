from io import StringIO, BytesIO

from charge.repository import Repository, _VersioningList


def test_create_empty():
    repo = Repository()
    assert repo._Repository__min_shell == 1
    assert repo._Repository__max_shell == 7
    assert len(repo.charges_iacm) == 0
    assert len(repo.charges_elem) == 0
    assert len(repo.iso_iacm) == 0
    assert len(repo.iso_elem) == 0


def test_set_shell_sizes():
    repo = Repository(2, 4)
    assert repo._Repository__min_shell == 2
    assert repo._Repository__max_shell == 4


def test_create_from_dir(lgf_data_dir):
    repo = Repository.create_from(str(lgf_data_dir))
    assert repo._Repository__min_shell == 1
    assert repo._Repository__max_shell == 7

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

    assert len(repo.iso_iacm) == 2
    assert len(repo.iso_elem) == 2


def test_create_traceable_from_dir(lgf_data_dir):
    repo = Repository.create_from(str(lgf_data_dir), traceable=True)
    assert repo._Repository__min_shell == 1
    assert repo._Repository__max_shell == 7

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

    assert len(repo.iso_iacm) == 2
    assert len(repo.iso_elem) == 2


def test_create_versioning_from_dir(lgf_data_dir):
    repo = Repository.create_from(str(lgf_data_dir), versioning=True)
    assert repo._Repository__min_shell == 1
    assert repo._Repository__max_shell == 7

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

    assert len(repo.iso_iacm) == 2
    assert len(repo.iso_elem) == 2

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

    assert repo0.charges_iacm == repo1.charges_iacm
    assert repo0.charges_elem == repo1.charges_elem


def test_read_write_versioning(lgf_data_dir):
    repo0 = Repository.create_from(str(lgf_data_dir), versioning=True)

    tmp = BytesIO()
    repo0.write(tmp)

    repo1 = Repository.read(tmp, versioning=True)

    assert repo0.charges_iacm == repo1.charges_iacm
    assert repo0.charges_elem == repo1.charges_elem

    assert isinstance(repo1.charges_iacm[1]['c18208da9e290c6faf8a0c58017d24d9'], _VersioningList)
    assert isinstance(repo1.charges_iacm[3]['76198d87470cc1b2f871da60449146bc'], _VersioningList)
    assert isinstance(repo1.charges_elem[1]['92ed00c54b2190be94748bee34b22847'], _VersioningList)
    assert isinstance(repo1.charges_elem[3]['17ac3199bf634022485c145821f358d5'], _VersioningList)
