from abc import ABC
from math import ceil, log
from typing import Any, Dict, List, Optional, Tuple

import networkx as nx
import numpy as np

from charge.util import AssignmentError
from charge.settings import MAX_BINS
from charge.nauty import Nauty
from charge.repository import ChargeSet, Repository

Atom = Any  # TODO: define this for the whole library
ChargeList = List[float]
WeightList = List[float]


class Collector(ABC):
    """Base class for collectors.

    Collectors query the repository for possible charges for each atom \
    in a molecule graph, and return a histogram describing the \
    distribution of the obtained charges for each atom.
    """
    def collect_values(
            self,
            graph: nx.Graph,
            iacm_data_only: bool,
            shells: List[int],
            **kwargs: Any
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
        raise NotImplemented()


class MeanCollector(Collector):
    """A collector that returns the mean of all charges found.

    For each atom, this collector collects possible charges, then it \
    returns their mean.

    Args:
        repo: The Repository to collect charges from.
        rounding_digits: The number of decimals to round charges to.
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
        self.__repository = repository
        self.__rounding_digits = rounding_digits
        self.__nauty = nauty if nauty is not None else Nauty()

    def collect_values(
            self,
            graph: nx.Graph,
            iacm_data_only: bool,
            shells: List[int],
            **kwargs: Any
            ) -> Dict[Atom, Tuple[ChargeList, WeightList]]:
        """Collect charges for a graph's atoms from a Repository.

        For each atom in the graph, this function will determine its \
        neighborhood, then collect all charges for such an atom with \
        such a neighborhood from the given repository. It does this \
        starting with the first shell size in shells, and continues \
        trying with subsequent shell sizes until at least one charge \
        is found (so you probably want to sort shells in descending \
        order). Finally, it calculates the mean of the collected \
        charges and rounds to the given precision.

        The nodes (atoms) in the graph are assumed to have an \
        'atom_type' data attribute associated with them containing the \
        elemental atom type and may have an 'iacm' attribute as well \
        with the IACM atom type.

        Args:
            graph: The graph to collect charges for
            iacm_data_only: If true, only use IACM data in the \
                    repository, and do not fall back to plain element \
                    data.
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

        for atom in graph.nodes():
            for shell_size in shells:
                atom_has_iacm = 'iacm' in graph.node[atom]
                attribute = 'iacm' if atom_has_iacm else 'atom_type'
                charges[atom] = self.__collect(
                        graph, atom, shell_size, attribute,
                        self.__repository.charges_iacm)

                if charges[atom] == [] and not iacm_data_only:
                    charges[atom] = self.__collect(
                            graph, atom, shell_size, 'atom_type',
                            self.__repository.charges_elem)

                if charges[atom] != []:
                    break
            else:
                no_vals.append(atom)

        if len(no_vals) > 0:
            err = 'Could not find charges for atoms {0}.'.format(', '.join(map(str, no_vals)))
            if not 0 in shells:
                err += ' Please retry with a smaller "shell" parameter.'
            raise AssignmentError(err)
        return charges

    def __collect(
            self,
            graph: nx.Graph,
            atom: Atom,
            shell_size: int,
            attribute: str,
            charges: ChargeSet
            ) -> Tuple[float, float]:
        """Collect charges for a particular atom."""
        if shell_size in charges:
            key = self.__nauty.canonize_neighborhood(graph, atom, shell_size, attribute)
            if key in charges[shell_size]:
                values = charges[shell_size][key]
                mean_charge = round(float(np.mean(values)), self.__rounding_digits)
                return [mean_charge], [1.0]
        return []


class HistogramCollector(Collector):
    """A collector that returns a histogram of all charges found.

    For each atom, this collector collects possible charges, then it \
    creates a histogram, with each bin becoming a kind of meta-charge \
    in a coarse-grained version of the original list. It returns these \
    meta-charges.

    Args:
        repository: The Repository to collect charges from.
        rounding_digits: The number of decimals to round charges to.
    """
    def __init__(
            self,
            repository: Repository,
            rounding_digits: int,
            nauty: Optional[Nauty]=None
            ) -> None:
        self.__repository = repository
        self.__rounding_digits = rounding_digits
        self.__nauty = nauty if nauty is not None else Nauty()

    def collect_values(
            self,
            graph: nx.Graph,
            iacm_data_only: bool,
            shells: List[int],
            **kwargs: Any
            ) -> Dict[Atom, Tuple[ChargeList, WeightList]]:
        """Collect charges for a graph's atoms from a Repository.

        For each atom in the graph, this function will determine its \
        neighborhood, then collect all charges for such an atom with \
        such a neighborhood from the given repository. It does this \
        starting with the first shell size in shells, and continues \
        trying with subsequent shell sizes until at least one charge \
        is found (so you probably want to sort shells in descending \
        order). Finally, it summarizes the found charges using a \
        histogram, and returns the course-grained surrogate charges.

        The nodes (atoms) in the graph are assumed to have an \
        'atom_type' data attribute associated with them containing the \
        elemental atom type and may have an 'iacm' attribute as well \
        with the IACM atom type.

        Args:
            graph: The graph to collect charges for
            iacm_data_only: If true, only use IACM data in the \
                    repository, and do not fall back to plain element \
                    data.
            shells: A list of shell sizes to try, in order, until a \
                    match is found.

        Raises:
            AssignmentError: If no charges could be found for at least \
                    one of the atoms in the molecule.

        Returns:
            A dictionary mapping atoms (nodes) in graph to lists of \
                    charges and their weights. The list will be empty \
                    if no charge is found.
        """
        score_hist = self.__score_histogram_log
        histograms = dict()
        no_vals = list()

        max_bins = max(int(kwargs['max_bins']), 1) if 'max_bins' in kwargs else MAX_BINS

        for atom in graph.nodes():
            for shell in shells:
                if shell in self.__repository.charges_iacm:
                    key = self.__nauty.canonize_neighborhood(
                            graph, atom, shell,
                            color_key='iacm' if 'iacm' in graph.node[atom] else 'atom_type')
                    if key in self.__repository.charges_iacm[shell]:
                        charges = self.__repository.charges_iacm[shell][key]
                        hist = self.__calculate_histogram(charges, max_bins)
                        histograms[atom] = score_hist(hist, np.mean(charges))
                        break

                if not iacm_data_only:
                    if shell in self.__repository.charges_elem:
                        key = self.__nauty.canonize_neighborhood(graph, atom, shell)
                        if key in self.__repository.charges_elem[shell]:
                            charges = self.__repository.charges_elem[shell][key]
                            hist = self.__calculate_histogram(charges, max_bins)
                            histograms[atom] = score_hist(hist, np.mean(charges))
                            break
            else:
                no_vals.append(atom)
        # collect charges from first shell size that gives a match, iacm if possible, elem otherwise unless disabled
        # if no matches at any shell size, add to no_vals

        if len(no_vals) > 0:
            err = 'Could not find charges for atoms {0}.'.format(', '.join(map(str, no_vals)))
            if not 0 in shells:
                err += (' Please retry with a smaller "shell" parameter'
                        ' or a SimpleCharger.')
            else:
                err += ' Please retry with a SimpleCharger.'
            raise AssignmentError(err)

        return histograms

    def __calculate_histogram(
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
        def round_to(x, grain):
            # rounds to nearest multiple of grain, with 0.5 rounded down
            return ceil((x / grain) - 0.5) * grain

        def median(values: List[float]):
            n = len(values)
            h = n // 2
            if n % 2 == 0:
                return sum(values[h-1:h+1]) / 2.0
            else:
                return values[h]

        # numpy.percentile isn't quite as accurate for small lists
        def first_quartile(values: List[float]):
            n = len(values)
            if n == 1:
                return values[0]
            if n % 2 == 0:
                return median(values[0:n // 2])
            q = n // 4
            if n % 4 == 1:
                return 0.75 * values[q-1] + 0.25 * values[q]
            if n % 4 == 3:
                return 0.75 * values[q] + 0.25 * values[q+1]

        def third_quartile(values: List[float]):
            n = len(values)
            if n == 1:
                return values[0]
            if n % 2 == 0:
                return median(values[n // 2:])
            q = n // 4
            if n % 4 == 1:
                return 0.25 * values[3 * q] + 0.75 * values[3 * q + 1]
            if n % 4 == 3:
                return 0.25 * values[3 * q + 1] + 0.75 * values[3 * q + 2]


        grain = 10**(-self.__rounding_digits)

        # calc F-D width
        iqr = third_quartile(charges) - first_quartile(charges)
        fd_width = 2.0 * iqr / (len(charges)**(1./3))
        if fd_width == 0.0:
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

    def __score_histogram_count(self, histogram: Tuple[ChargeList, WeightList],
            mean: float) -> Tuple[ChargeList, WeightList]:
        return histogram

    def __score_histogram_log(self, histogram: Tuple[ChargeList, WeightList],
            mean: float) -> Tuple[ChargeList, WeightList]:
        bin_centers, counts = histogram
        scores = list(map(log, counts))
        return bin_centers, scores

    def __score_histogram_martin(self, histogram: Tuple[ChargeList, WeightList],
            mean: float) -> Tuple[ChargeList, WeightList]:
        bin_centers, counts = histogram
        scores = []
        for count, center in zip(counts, bin_centers):
            scores.append(count / (1.0 + (center - mean)**2))
        return bin_centers, scores
