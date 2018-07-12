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

    assert len(repo.charges_elem) == 7
    assert len(repo.charges_elem[1]) == 5
    assert len(repo.charges_elem[2]) == 7
    assert len(repo.charges_elem[3]) == 7

    assert len(repo.iso_iacm) == 2
    assert len(repo.iso_elem) == 2
