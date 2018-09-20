import pytest

from charge.repository import Repository

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
    repo = Repository.create_from(str(lgf_data_dir));
    assert repo._Repository__min_shell == 1
    assert repo._Repository__max_shell == 7

    assert len(repo.charges_iacm) == 7
    assert len(repo.charges_iacm[1]) == 5
    assert len(repo.charges_iacm[2]) == 7
    assert len(repo.charges_iacm[3]) == 7

    assert repo.charges_iacm[1]['c18208da9e290c6faf8a0c58017d24d9'] == [-0.516]
    assert repo.charges_iacm[3]['76198d87470cc1b2f871da60449146bc'] == [0.077, 0.077, 0.077]

    assert len(repo.charges_elem) == 7
    assert len(repo.charges_elem[1]) == 5
    assert len(repo.charges_elem[2]) == 7
    assert len(repo.charges_elem[3]) == 7

    assert repo.charges_elem[1]['92ed00c54b2190be94748bee34b22847'] == [-0.516]
    assert repo.charges_elem[3]['17ac3199bf634022485c145821f358d5'] == [0.077, 0.077, 0.077]

    assert len(repo.iso_iacm) == 2
    assert len(repo.iso_elem) == 2


def test_create_traceable_from_dir(lgf_data_dir):
    repo = Repository.create_from(str(lgf_data_dir), traceable=True)
    assert repo._Repository__min_shell == 1
    assert repo._Repository__max_shell == 7

    assert len(repo.charges_iacm) == 7
    assert len(repo.charges_iacm[1]) == 5
    assert len(repo.charges_iacm[2]) == 7
    assert len(repo.charges_iacm[3]) == 7

    assert repo.charges_iacm[1]['c18208da9e290c6faf8a0c58017d24d9'] == [(-0.516, 15610, 2)]
    assert repo.charges_iacm[3]['76198d87470cc1b2f871da60449146bc'] == [
            (0.077, 1204, 1), (0.077, 1204, 4), (0.077, 1204, 5)]

    assert len(repo.charges_elem) == 7
    assert len(repo.charges_elem[1]) == 5
    assert len(repo.charges_elem[2]) == 7
    assert len(repo.charges_elem[3]) == 7

    assert repo.charges_elem[1]['92ed00c54b2190be94748bee34b22847'] == [(-0.516, 15610, 2)]
    assert repo.charges_elem[3]['17ac3199bf634022485c145821f358d5'] == [
            (0.077, 1204, 1), (0.077, 1204, 4), (0.077, 1204, 5)]

    assert len(repo.iso_iacm) == 2
    assert len(repo.iso_elem) == 2


# TODO: test read and write
