import inspect
import sys
from abc import ABC, abstractmethod
from types import MethodType
from typing import Iterable, Optional, Union

import networkx as nx

from charge import util
from charge.collectors import HistogramCollector, MeanCollector, ModeCollector, MedianCollector
from charge.nauty import Nauty
from charge.repository import Repository
from charge.settings import ROUNDING_DIGITS, DEFAULT_TOTAL_CHARGE, MAX_ROUNDING_DIGITS, MAX_BINS
from charge.solvers import CDPSolver, DPSolver, ILPSolver, SimpleSolver


class Charger(ABC):
    """Assigns charges to a molecule using data from a repository.

    This is a base class, derive from it and set the _collector and \
    _solver member variables to the desired objects in your derived \
    class.
    """
    @abstractmethod
    def __init__(self,
                 repository: Repository,
                 rounding_digits: int = ROUNDING_DIGITS,
                 nauty: Optional[Nauty]=None,
                 **kwargs) -> None:
        """Create a ChargerBase.

        Args:
            repository: The repository to obtain charges from
            rounding_digits: Number of significant digits to round charges to.
            nauty: An external Nauty instance to use for canonization
        """
        # These are all protected, not private
        self._nauty = nauty if nauty is not None else Nauty()
        self._repo = repository
        self._rounding_digits = min(max(rounding_digits, 0), MAX_ROUNDING_DIGITS)

        # To be assigned in derived classes
        self._collector = None      # type: Collector
        self._solver = None         # type: Solver

    def charge(
            self,
            graph: nx.Graph,
            total_charge: int=DEFAULT_TOTAL_CHARGE,
            iacmize: bool=False,
            iacm_data_only: bool=False,
            shell: Union[None, int, Iterable[int]]=None,
            **kwargs
            ) -> None:
        """Add charges to a molecule.

        The iacmize option will try to assign IACM atom types to a \
        plain element molecule, and then match to the IACM data only.

        Args:
            graph: The molecule grapth to add charges to
            iacmize: Whether to convert to IACM types and use IACM data
            iacm_data_only: Use only IACM data in the repository
            shell: Shell sizes to use.
            total_charge: Total charge of the molecule
        """
        if not shell:
            shells = sorted(self._repo.charges_iacm.keys(), reverse=True)
        elif isinstance(shell, int):
            shells = [shell]
        elif isinstance(shell, Iterable):
            shells = shell
        else:
            raise TypeError('shell must be int or List[int]')

        if iacmize:
            graph = util.iacmize(graph)

        values = self._collector.collect_values(graph, iacm_data_only or iacmize, shells, **kwargs)
        self._solver.solve_partial_charges(graph, values, total_charge, **kwargs)
        self.__add_redistributed_charges(graph, total_charge)

    def __add_redistributed_charges(
            self,
            graph: nx.Graph,
            total_charge: int
            ) -> None:
        """Add adjusted charges that fit the exact desired total charge.

        This function adds a 'partial_charge_redist' field to the graph \
        which contains the adjusted charge. Summing these fields over \
        the entire molecule will give the desired total charge.

        This uses the scores associated with each charge to make \
        larger adjustments to charges we're less certain about.

        Args:
            graph: The molecule graph to adjust
            total_charge: The total charge to match.
        """
        total_score = graph.graph['score']
        total_error = total_charge - graph.graph['total_charge']

        total_charge_redist = 0
        min_score = float('inf')
        min_atom_data = None
        proportions = {}
        total_prop = 0
        for v, data in graph.nodes().data():
            proportions[v] = total_score / data['score'] if data['score'] > 0 else 0
            total_prop += proportions[v]

        for v, data in graph.nodes().data():
            delta = round((proportions[v] * total_error) / total_prop, self._rounding_digits)
            data['partial_charge_redist'] = round(data['partial_charge'] + delta, self._rounding_digits)
            total_charge_redist += data['partial_charge_redist']
            if data['score'] < min_score:
                min_score = data['score']
                min_atom_data = data

        if min_atom_data:
            min_atom_data['partial_charge_redist'] = round(min_atom_data['partial_charge_redist'] +
                    total_charge - total_charge_redist, self._rounding_digits)
            total_charge_redist += total_charge - total_charge_redist
            graph.graph['total_charge_redist'] = round(total_charge_redist, self._rounding_digits)


class MeanCharger(Charger):
    """A simple reliable charger.

    This charger assigns the mean of all matching charges in the \
    repository. As long as there is at least one matching charge for \
    each atom, it will give an answer.
    """
    def __init__(
            self,
            repository: Repository,
            rounding_digits: int,
            nauty: Optional[Nauty]=None
            ) -> None:
        """Create a MeanCharger.

        Nauty instances manage an external process, so they're \
        somewhat expensive to create. If you have multiple chargers, \
        you could consider sharing one between them.

        Args:
            repository: The repository to get charges from
            rounding_digits: Number of digits to round charges to
            nauty: An external Nauty instance to use
        """
        super().__init__(repository, rounding_digits, nauty)
        self._collector = MeanCollector(repository, rounding_digits, self._nauty)
        self._solver = SimpleSolver(rounding_digits)


class MedianCharger(Charger):
    """A simple reliable charger.

    This charger assigns the median of all matching charges in the \
    repository. As long as there is at least one matching charge for \
    each atom, it will give an answer.
    """
    def __init__(
            self,
            repository: Repository,
            rounding_digits: int,
            nauty: Optional[Nauty]=None
            ) -> None:
        """Create a MeanCharger.

        Nauty instances manage an external process, so they're \
        somewhat expensive to create. If you have multiple chargers, \
        you could consider sharing one between them.

        Args:
            repository: The repository to get charges from
            rounding_digits: Number of digits to round charges to
            nauty: An external Nauty instance to use
        """
        super().__init__(repository, rounding_digits, nauty)
        self._collector = MedianCollector(repository, rounding_digits, self._nauty)
        self._solver = SimpleSolver(rounding_digits)


class ModeCharger(Charger):
    """A simple reliable charger.

    This charger assigns the mode of the histogram of all matching \
    charges in the repository. If the histogram is multimodal, it \
    chooses the mode closest to the median. As long as there is at \
    least one matching charge for each atom, it will give an answer.
    """
    def __init__(
            self,
            repository: Repository,
            rounding_digits: int,
            nauty: Optional[Nauty]=None
            ) -> None:
        """Create a ModeCharger.

        Nauty instances manage an external process, so they're \
        somewhat expensive to create. If you have multiple chargers, \
        you could consider sharing one between them.

        Args:
            repository: The repository to get charges from
            rounding_digits: Number of digits to round charges to
            nauty: An external Nauty instance to use
        """
        super().__init__(repository, rounding_digits, nauty)
        self._collector = ModeCollector(repository, rounding_digits, self._nauty)
        self._solver = SimpleSolver(rounding_digits)


class ILPCharger(Charger):
    """A charger that uses Integer Linear Programming.

    This charger calculates an optimal charge distribution given the \
    charges found in a repository, using Integer Linear Programming.
    """
    def __init__(
            self,
            repository: Repository,
            rounding_digits: int,
            max_seconds: int,
            nauty: Optional[Nauty]=None,
            scoring: Optional[MethodType]=None,
            max_bins: Optional[int] = MAX_BINS
            ) -> None:
        """Create an ILPCharger.

        Nauty instances manage an external process, so they're \
        somewhat expensive to create. If you have multiple chargers, \
        you could consider sharing one between them.

        Args:
            repository: The repository to get charges from
            rounding_digits: Number of digits to round charges to
            max_seconds: A limit for the run time. If no solution is \
                    found within this limit, an exception will be \
                    raised.
            nauty: An external Nauty instance to use
            scoring: A scoring function for the histogram. See \
             :func:`~charge.collectors.HistogramCollector.score_histogram_count`, \
             :func:`~charge.collectors.HistogramCollector.score_histogram_log`, and \
             :func:`~charge.collectors.HistogramCollector.score_histogram_martin`.
        """
        super().__init__(repository, rounding_digits, nauty)
        self._collector = HistogramCollector(repository, rounding_digits, self._nauty, scoring, max_bins)
        self._solver = ILPSolver(rounding_digits, max_seconds)


class DPCharger(Charger):
    """A charger that uses Dynamic Programming, Python version.

    This charger calculates an optimal charge distribution given the \
    charges found in a repository, using Dynamic Programming. The C \
    implementation in CDPCharger does the same, but runs faster.
    """
    def __init__(
            self,
            repository: Repository,
            rounding_digits: int,
            nauty: Optional[Nauty]=None,
            scoring: Optional[MethodType]=None,
            max_bins: Optional[int] = MAX_BINS
            ) -> None:
        """Create an DPCharger.

        Nauty instances manage an external process, so they're \
        somewhat expensive to create. If you have multiple chargers, \
        you could consider sharing one between them.

        Args:
            repository: The repository to get charges from
            rounding_digits: Number of digits to round charges to
            nauty: An external Nauty instance to use
            scoring: A scoring function for the histogram. See \
             :func:`~charge.collectors.HistogramCollector.score_histogram_count`, \
             :func:`~charge.collectors.HistogramCollector.score_histogram_log`, and \
             :func:`~charge.collectors.HistogramCollector.score_histogram_martin`.
        """
        super().__init__(repository, rounding_digits, nauty)
        self._collector = HistogramCollector(repository, rounding_digits, self._nauty, scoring, max_bins)
        self._solver = DPSolver(rounding_digits)


class CDPCharger(Charger):
    """A charger that uses Dynamic Programming, C version.

    This charger calculates an optimal charge distribution given the \
    charges found in a repository, using Dynamic Programming. This is \
    a faster C implementation than the Python one in DPCharger.
    """
    def __init__(
            self,
            repository: Repository,
            rounding_digits: int,
            nauty: Optional[Nauty]=None,
            scoring: Optional[MethodType]=None,
            max_bins: Optional[int] = MAX_BINS
            ) -> None:
        """Create a CDPCharger.

        Nauty instances manage an external process, so they're \
        somewhat expensive to create. If you have multiple chargers, \
        you could consider sharing one between them.

        Args:
            repository: The repository to get charges from
            rounding_digits: Number of digits to round charges to
            nauty: An external Nauty instance to use
            scoring: A scoring function for the histogram collector. See \
             :func:`~charge.collectors.HistogramCollector.score_histogram_count`, \
             :func:`~charge.collectors.HistogramCollector.score_histogram_log`, and \
             :func:`~charge.collectors.HistogramCollector.score_histogram_martin`.
        """
        super().__init__(repository, rounding_digits, nauty)
        self._collector = HistogramCollector(repository, rounding_digits, self._nauty, scoring, max_bins)
        self._solver = CDPSolver(rounding_digits)


def make_charger(
        name: str,
        repository: Repository,
        rounding_digits: int,
        max_seconds: int,
        nauty: Optional[Nauty] = None,
        scoring: Optional[MethodType] = None,
        max_bins: Optional[int] = MAX_BINS) -> Charger:
    # get function parameter names and values:
    # since we did not declare any new variables yet, this equals the local variables
    local_vars = locals().copy()
    # get all classes in this module
    clsmembers = inspect.getmembers(sys.modules[__name__], inspect.isclass)

    for cls_name, cls in clsmembers:
        # find class with the matching name
        if cls_name == name:
            parameters = dict()
            # match cls.__init__() parameters to make_charger() parameters
            for param_name in inspect.signature(cls).parameters:
                if param_name in local_vars:
                    parameters[param_name] = local_vars[param_name]
            # return instance of cls
            return cls(**parameters)

    # class not found
    raise ValueError('Invalid charger name {}'.format(name))
