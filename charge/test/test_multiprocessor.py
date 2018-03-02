import pytest

from charge.multiprocessor import MultiProcessor


class ProcClass:
    def __init__(self):
        pass


class ProcClassInitTester:
    def __init__(self, arg1, arg2='3'):
        assert arg1 == 1
        assert arg2 == '2'


class ProcClassAdder:
    def __init__(self, addend):
        self.__addend = addend

    def process(self, number):
        return number + self.__addend


def test_create1():
    mp = MultiProcessor(ProcClass, num_processes=1)
    mp.shutdown()


def test_create2():
    mp = MultiProcessor(ProcClass, num_processes=2)
    mp.shutdown()


def test_create_initargs1():
    mp = MultiProcessor(ProcClassInitTester, (1,), 1)
    mp.shutdown()


def test_create_initargs2():
    mp = MultiProcessor(ProcClassInitTester, 1, 1)
    mp.shutdown()


def test_create_initargs3():
    mp = MultiProcessor(ProcClassInitTester, (1, '2'), 1)
    mp.shutdown()


def test_scope_guard():
    with MultiProcessor(ProcClass, num_processes=2) as mp:
        assert mp is not None


def test_processing_serial():
    mp = MultiProcessor(ProcClassAdder, 5, num_processes=1)
    result = list(mp.processed(range(100)))
    mp.shutdown()
    assert result == list(range(5, 105))


def test_processing_parallel():
    mp = MultiProcessor(ProcClassAdder, 5, num_processes=2)
    result = list(mp.processed(range(100)))
    mp.shutdown()
    assert set(result) == set(range(5, 105))
