from abc import ABC, abstractmethod
from math import log
from types import MethodType
from typing import Any, Dict, List, Optional, Tuple

import networkx as nx
import numpy as np

from charge.charge_types import ChargeList, WeightList
from charge.nauty import Nauty
from charge.repository import Repository, EitherChargeSet, _VersioningList
from charge.settings import MAX_BINS
from charge.util import AssignmentError, third_quartile, round_to, median, first_quartile

Atom = Any  # TODO: define this for the whole library


class Collector(ABC):
    """Base class for collectors.

    Collectors query the repository for possible charges for each atom \
    in a molecule graph, and return a histogram describing the \
    distribution of the obtained charges for each atom.
    Args:
        repo: The Repository to collect charges from.
        rounding_digits: The number of decimals to round charges to.
    """
    @abstractmethod
    def __init__(
            self,
            repository: Repository,
            rounding_digits: int,
            nauty: Optional[Nauty]=None
            ) -> None:
        """Create a SimpleCollector.

        Args:
            repository: The repository to collect charges from
            rounding_digits: Number of digits to round charges to
            nauty: An external Nauty instance to use
        """
        self._repository = repository
        self._rounding_digits = rounding_digits
        self._nauty = nauty if nauty is not None else Nauty()

    def collect_values(
            self,
            graph: nx.Graph,
            iacm_data_only: bool,
            shells: List[int]
            ) -> Dict[Atom, Tuple[ChargeList, WeightList]]:
        """Collect charges for a graph's atoms.

        For each atom in the graph, return a list of possible \
        charges.

        Args:
            graph: The graph to collect charges for
            iacm_data_only: If true, do not fall back to plain elements
            shells: A list of shell sizes to try, in order, until a \
                    match is found.

        Raises:
            AssignmentError: If no charges could be found for at least \
                    one of the atoms in the molecule.

        Returns:
            A dictionary mapping atoms (nodes) in graph to a tuple of \
                    lists, the first with charges, the second with \
                    weights.
        """
        charges = dict()
        no_vals = list()
        keys = dict()

        for atom in graph.nodes():
            for shell_size in shells:
                atom_has_iacm = 'iacm' in graph.node[atom]

                if atom_has_iacm:
                    if shell_size in self._repository.charges_iacm:
                        key = self._nauty.canonize_neighborhood(graph, atom, shell_size, 'iacm')
                        if key in self._repository.charges_iacm[shell_size]:
                            charges[atom] = self._collect(self._repository.charges_iacm, shell_size, key)
                            keys[atom] = key

                if not atom_has_iacm or (not atom in charges and not iacm_data_only):
                    if shell_size in self._repository.charges_elem:
                        key = self._nauty.canonize_neighborhood(graph, atom, shell_size, 'atom_type')
                        if key in self._repository.charges_elem[shell_size]:
                            charges[atom] = self._collect(self._repository.charges_elem, shell_size, key)
                            keys[atom] = key

                if atom in charges:
                    break
            else:
                no_vals.append(atom)

        self._handle_error(no_vals, shells)
        return charges, keys

    @abstractmethod
    def _collect(
            self,
            chargeset: EitherChargeSet,
            shell_size: int,
            key: str
            ) -> Tuple[ChargeList, WeightList]:
        """Collect charges for a single atom.

        Queries the given chargeset with the current shell_size (k) and key \
        (hash of the atom's k-neighborhood graph) and returns \
        processed charge values.

        Args:
            chargeset: A dictonary index by shell_size and key that holds the charge values
            shell_size: Shell size k
            key: Hash of the atom's k-neighborhood

        Returns:
            A tuple of the processed charges and associated weights.
        """
        pass

    def _handle_error(
            self,
            no_vals: List[Atom],
            shells: List[int]
            ) -> None:
        """Raises an AssignmentError if no_vals is not empty.

        Args:
            no_vals: List of atoms for which no charges could be found.
            shells: List of tried shell sizes.

        Raises:
            AssignmentError: If no_vals is not empty.
        """
        if len(no_vals) > 0:
            err = 'Could not find charges for atoms {0}.'.format(', '.join(map(str, no_vals)))
            if not 0 in shells:
                err += ' Please retry with a smaller "shell" parameter.'
            raise AssignmentError(err)


class MeanCollector(Collector):
    """A collector that returns the mean of all charges found.

    For each atom, this collector collects possible charges, then it \
    returns their mean.

    Args:
        repo: The Repository to collect charges from.
        rounding_digits: The number of decimals to round charges to.
        nauty: An external Nauty instance to use
    """
    def __init__(
            self,
            repository: Repository,
            rounding_digits: int,
            nauty: Optional[Nauty]=None
            ) -> None:
        """Create a MeanCollector.

        Args:
            repository: The repository to collect charges from
            rounding_digits: Number of digits to round charges to
            nauty: An external Nauty instance to use
        """
        super().__init__(repository, rounding_digits, nauty)

    def _collect(
            self,
            chargeset: EitherChargeSet,
            shell_size: int,
            key: str
            ) -> Tuple[ChargeList, WeightList]:
        """Collect charges for a single atom.

        Queries the given chargeset with the current shell_size (k) and key \
        (hash of the atom's k-neighborhood graph) and returns \
        the mean of the charge values.

        Args:
            chargeset: A dictonary index by shell_size and key that holds the charge values
            shell_size: Shell size k
            key: Hash of the atom's k-neighborhood

        Returns:
            A tuple of the containing the mean charge and an arbitrary weight (1).
        """
        values = chargeset[shell_size][key]
        mean_charge = round(float(np.mean(values)), self._rounding_digits)
        return [mean_charge], [1.0]


class MedianCollector(Collector):
    """A collector that returns the median of all charges found.

    For each atom, this collector collects possible charges, then it \
    returns their median.

    Args:
        repo: The Repository to collect charges from.
        rounding_digits: The number of decimals to round charges to.
        nauty: An external Nauty instance to use
    """
    def __init__(
            self,
            repository: Repository,
            rounding_digits: int,
            nauty: Optional[Nauty]=None
            ) -> None:
        """Create a MedianCollector.

        Args:
            repository: The repository to collect charges from
            rounding_digits: Number of digits to round charges to
            nauty: An external Nauty instance to use
        """
        super().__init__(repository, rounding_digits, nauty)

    def _collect(self,
                 chargeset: EitherChargeSet,
                 shell_size: int,
                 key: str) -> Tuple[ChargeList, WeightList]:
        """Collect charges for a single atom.

        Queries the given chargeset with the current shell_size (k) and key \
        (hash of the atom's k-neighborhood graph) and returns \
        the median of the charge values.

        Args:
            chargeset: A dictonary index by shell_size and key that holds the charge values
            shell_size: Shell size k
            key: Hash of the atom's k-neighborhood

        Returns:
            A tuple of the containing the median charge and an arbitrary weight (1).
        """
        values = chargeset[shell_size][key]
        median_charge = round(median(values), self._rounding_digits)
        return [median_charge], [1.0]


class HistogramCollector(Collector):
    """A collector that returns a histogram of all charges found.

    For each atom, this collector collects possible charges, then it \
    creates a histogram, with each bin becoming a kind of meta-charge \
    in a coarse-grained version of the original list. It returns these \
    meta-charges.

    Args:
        repository: The Repository to collect charges from.
        rounding_digits: The number of decimals to round charges to.
        nauty: An external Nauty instance to use
        scoring: A scoring function for the histogram. See \
             :func:`~charge.collectors.HistogramCollector.score_histogram_count`, \
             :func:`~charge.collectors.HistogramCollector.score_histogram_log`, and \
             :func:`~charge.collectors.HistogramCollector.score_histogram_martin`.
    """
    def __init__(
            self,
            repository: Repository,
            rounding_digits: int,
            nauty: Optional[Nauty]=None,
            scoring: Optional[MethodType]=None,
            max_bins: Optional[int]=MAX_BINS
            ) -> None:
        super().__init__(repository, rounding_digits, nauty)
        self._score_hist = scoring if scoring else HistogramCollector.score_histogram_log
        self._max_bins = max(max_bins, 1)

    def _collect(self,
                 chargeset: EitherChargeSet,
                 shell_size: int,
                 key: str
                 ) -> Tuple[ChargeList, WeightList]:
        """Collect charges for a single atom.

        Queries the given chargeset with the current shell_size (k) and key \
        (hash of the atom's k-neighborhood graph) and returns \
        the histogram of the charge values.

        Args:
            chargeset: A dictonary index by shell_size and key that holds the charge values
            shell_size: Shell size k
            key: Hash of the atom's k-neighborhood

        Returns:
            A tuple of the containing the charge histogram centers and scores.
        """
        values = chargeset[shell_size][key]
        hist = self._calculate_histogram(values, self._max_bins)
        return self._score_hist(hist, float(np.mean(values)))

    def _handle_error(
            self,
            no_vals: List[Atom], shells: List[int]
            ) -> None:
        """Raises an AssignmentError if no_vals is not empty.

        Args:
            no_vals: List of atoms for which no charges could be found.
            shells: List of tried shell sizes.

        Raises:
            AssignmentError: If no_vals is not empty.
        """
        if len(no_vals) > 0:
            err = 'Could not find charges for atoms {0}.'.format(', '.join(map(str, no_vals)))
            if not 0 in shells:
                err += (' Please retry with a smaller "shell" parameter'
                        ' or a SimpleCharger.')
            else:
                err += ' Please retry with a SimpleCharger.'
            raise AssignmentError(err)

    def _calculate_histogram(
            self,
            charges: ChargeList,
            max_bins: int
            ) -> Tuple[ChargeList, WeightList]:
        """Create a histogram from a raw list of charges.

        Histogram bins will be chosen so that:
            - They are all equally wide
            - Their centers fall on n-significant-digit numbers,
                according to self.__rounding_digits
            - Their widths are n-significant-digit numbers,
                according to self.__rounding_digits
            - There are at most max_bins non-zero bins
            - The middle bin center is on the median charge, rounded to
                an n-significant digit number.
            - The bin size is as close to the size chosen by the Freedman-
                Diaconis rule as it can be, given the above constraints.

        Only non-empty bins will be returned.

        Args:
            charges: A list of charges to process into a histogram
            max_bins: The maximum number of bins the histogram should have
        """

        grain = 10**(-self._rounding_digits)

        # calc F-D width
        iqr = third_quartile(charges) - first_quartile(charges)
        fd_width = 2.0 * iqr / (len(charges)**(1./3))
        if fd_width < grain:
            fd_width = grain

        median_charge = round_to(median(charges), grain)

        # subtract one because we start the loop with an increment
        spacing = round_to(fd_width, grain) / grain
        num_bins = max_bins + 1
        while num_bins > max_bins:
            step = grain * spacing

            # align center bin to median and find edge bin centers
            min_charge_bin = (median_charge -
                    round_to(median_charge - min(charges), step))
            max_charge_bin = (median_charge +
                    round_to(max(charges) - median_charge, step))

            # add half a step to include the last value
            bin_centers = np.arange(min_charge_bin, max_charge_bin + 0.5*step, step)

            min_bin_edge = min_charge_bin - 0.5*step
            max_bin_edge = max_charge_bin + 0.5*step
            # add half a step to include the last value
            bin_edges = np.arange(min_bin_edge, max_bin_edge + 0.5*step, step)
            # extend a bit to compensate for round-off error
            bin_edges[-1] += 1e-15

            counts, _ = np.histogram(charges, bins=bin_edges)
            num_bins = np.count_nonzero(counts)
            spacing += 1

        nonzero_bins = np.nonzero(counts)
        counts = counts[nonzero_bins].tolist()
        bin_centers = bin_centers[nonzero_bins].tolist()

        return bin_centers, counts

    @staticmethod
    def score_histogram_count(histogram: Tuple[ChargeList, WeightList],
            *args) -> Tuple[ChargeList, WeightList]:
        """Scores the counts of the charge histgram.

        This function simply returns the original counts as score.

        Args:
            histogram: A list of charges and their counts.

        Returns:
            A lists of charges and their scores.
        """
        return histogram

    @staticmethod
    def score_histogram_log(histogram: Tuple[ChargeList, WeightList],
                            *args) -> Tuple[ChargeList, WeightList]:
        """Scores the counts of the charge histgram.

        This function returns the log counts as score.

        Args:
            histogram: A list of charges and their counts.

        Returns:
            A lists of charges and their scores.
        """
        bin_centers, counts = histogram
        scores = list(map(log, counts))
        return bin_centers, scores

    @staticmethod
    def score_histogram_martin(histogram: Tuple[ChargeList, WeightList],
            mean: float) -> Tuple[ChargeList, WeightList]:
        """Scores the counts of the charge histgram.

        This function returns the score: count / (1 + (charge - mean)**2)

        Args:
            histogram: A list of charges and their counts.
            mean: The mean of the charges.

        Returns:
            A lists of charges and their scores.
        """
        bin_centers, counts = histogram
        scores = []
        for count, center in zip(counts, bin_centers):
            scores.append(count / (1.0 + (center - mean)**2))
        return bin_centers, scores


class ModeCollector(HistogramCollector):
    """A collector that returns the mode of the histogram of all charges found.

    For each atom, this collector collects possible charges, then it \
    creates a histogram, with each bin becoming a kind of meta-charge \
    in a coarse-grained version of the original list. It returns the \
    the mode of these meta-charges. If the histogram is multimodal,
    it returns the mode closest to the median of all possible charges.

    Args:
        repository: The Repository to collect charges from.
        rounding_digits: The number of decimals to round charges to.
        nauty: An external Nauty instance to use
    """
    def __init__(
            self,
            repository: Repository,
            rounding_digits: int,
            nauty: Optional[Nauty] = None,
            max_bins: Optional[int] = MAX_BINS
    ) -> None:
        super().__init__(repository, rounding_digits, nauty, HistogramCollector.score_histogram_count, max_bins)

    def _collect(self,
                 chargeset: EitherChargeSet,
                 shell_size: int,
                 key: str
                 ) -> Tuple[ChargeList, WeightList]:
        """Collect charges for a single atom.

        Queries the given chargeset with the current shell_size (k) and key \
        (hash of the atom's k-neighborhood graph) and returns \
        the mode of the charge values.

        Args:
            chargeset: A dictonary index by shell_size and key that holds the charge values
            shell_size: Shell size k
            key: Hash of the atom's k-neighborhood

        Returns:
            A tuple of the containing the mode of the charges and the number of values that support it.
        """
        values = chargeset[shell_size][key]
        bin_centers, counts = self._calculate_histogram(values, self._max_bins)
        median_charge = median(values)

        mode_charge, mode_count = max(zip(bin_centers, counts),
                                      key=lambda x: (x[1], -abs(median_charge - x[0])))

        return [mode_charge], [mode_count]


class CachingCollector(Collector):
    """A proxy collector that caches return values.

    This collector caches the return values of it another collector. The CachingCollector \
    can detect changes in the Repository of the other collector, if the Repository's \
    versioning attribute is set to True.

    Args:
        basecollector: The collector to cache.
    """
    def __init__(
            self,
            basecollector: Collector):
        super().__init__(basecollector._repository, basecollector._rounding_digits, basecollector._nauty)
        self.__base = basecollector
        self.__cache = dict()

    def _collect(self,
                 chargeset: EitherChargeSet,
                 shell_size: int,
                 key: str
                 ) -> Tuple[ChargeList, WeightList]:
        """Collect charges for a single atom.

        Queries the given chargeset with the current shell_size (k) and key \
        (hash of the atom's k-neighborhood graph) and returns \
        the processed charge values.

        Args:
            chargeset: A dictonary index by shell_size and key that holds the charge values
            shell_size: Shell size k
            key: Hash of the atom's k-neighborhood

        Returns:
            The return value of the base collector.
        """

        cache_key = (id(chargeset), shell_size, key)
        if isinstance(chargeset[shell_size][key], _VersioningList):
            version = chargeset[shell_size][key].version
        else:
            version = None

        if cache_key in self.__cache:
            cache_version, values = self.__cache[cache_key]
            if not version or version == cache_version:
                return values

        values = self.__base._collect(chargeset, shell_size, key)
        self.__cache[cache_key] = (version, values)
        return values
