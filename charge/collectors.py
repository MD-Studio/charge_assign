from typing import List

import networkx as nx
import numpy as np

from charge.base import Charger, C
from charge.settings import MAX_BINS


class AtomicMeanCollector(Charger):
    # finds matching charges and returns the mean, rounded to
    # a set number of significant digits

    def _collect_values(self, graph: nx.Graph, iacm_only: bool, shells: List[int], **kwargs) -> C:

        charges = dict()
        no_vals = list()

        for atom in graph.nodes():
            for shell in shells:
                # collect from IACM repo data based on IACM atom types in query molecule
                # if the query molecule doesn't have IACM types, use plain elements, but
                # still match against the IACM repo data.
                key = self._nauty.canonize_neighborhood(graph, atom, shell,
                                                        color_key='iacm' if 'iacm' in graph.node[atom] else 'atom_type')
                if key in self._repo.charges_iacm[shell]:
                    values = self._repo.charges_iacm[shell][key]
                    charges[atom] = (round(float(np.mean(values)), self._rounding_digits), 1.0)
                    break
                elif not iacm_only:
                    # if the IACM didn't match, try with plain elements against plain
                    # element repo data, unless that was disabled
                    key = self._nauty.canonize_neighborhood(graph, atom, shell)
                    if key in self._repo.charges_elem[shell]:
                        values = self._repo.charges_elem[shell][key]
                        charges[atom] = (round(float(np.mean(values)), self._rounding_digits), 1.0)
                        break
            else:
                no_vals.append(atom)

        if len(no_vals) > 0:
            return no_vals
        else:
            return charges


class AtomicHistogramCollector(Charger):
    # Collects charges, builds a coarse-grained surrogate representation
    # using a histogram, then returns that

    def _collect_values(self, graph: nx.Graph, iacm_only: bool, shells: List[int], **kwargs) -> C:

        charges = dict()
        no_vals = list()

        max_bins = max(int(kwargs['max_bins']), 1) if 'max_bins' in kwargs else MAX_BINS

        def process_vals(values, atom):
            # values is a list of charges found for the given atom
            hist, bin_edges = np.histogram(values)
            if len(hist) > max_bins:
                hist, bin_edges = np.histogram(values, bins=max_bins)

            hist = hist.tolist()
            bin_edges = bin_edges.tolist()
            # hist is a list of counts of charges per bin
            # bin_edges is a list of bin edges (one longer than hist, equal width)

            atom_charges = [round(bin_edges[i] + 0.5 * (bin_edges[i + 1] - bin_edges[i]), self._rounding_digits) \
                       for i in range(len(bin_edges) - 1)]
            # atom_charges is the centers of the bin edges, which is the charges corresponding to the bins
            # but rounded, so now we may have bins with the same charge
            filtered_hist, filtered_charges = [hist[0]], [atom_charges[0]]
            for count, charge in zip(hist[1:], atom_charges[1:]):
                if count > 0:
                    if charge == filtered_charges[-1]:
                        filtered_hist[-1] += count
                    else:
                        filtered_charges.append(charge)
                        filtered_hist.append(count)
            # filtered_hist is the list of counts with merged bins, which are still equal-width because they're all 1 or 0.1 or 0.01 wide
            # zero count bins are out though
            # filtered_charges are the corresponding centers

            if 'scale' in kwargs and kwargs['scale']:
                s = sum(filtered_hist)
                # total number of charges
                mean_val = np.mean(values)
                # mean charge
                filtered_hist = [y/s/(1 + (filtered_charges[i] - mean_val)**2) for i, y in enumerate(filtered_hist)]
                # i is the i'th bin, y is the count for that bin
                # for each bin
                    # count = count * (1 + deviation of bin charge from mean squared) / total number of charges
                    # some kind of fudge factor???

            charges[atom] = (filtered_charges, filtered_hist)
            # set charges for this atom to a histogram, charges per bin and count or relative frequency

        for atom in graph.nodes():
            for shell in shells:
                key = self._nauty.canonize_neighborhood(graph, atom, shell,
                                                        color_key='iacm' if 'iacm' in graph.node[
                                                            atom] else 'atom_type')
                if key in self._repo.charges_iacm[shell]:
                    process_vals(self._repo.charges_iacm[shell][key], atom)
                    break
                elif not iacm_only:
                    key = self._nauty.canonize_neighborhood(graph, atom, shell)
                    if key in self._repo.charges_elem[shell]:
                        process_vals(self._repo.charges_elem[shell][key], atom)
                        break
            else:
                no_vals.append(atom)
        # collect charges from first shell size that gives a match, iacm if possible, elem otherwise unless disabled
        # if no matches at any shell size, add to no_vals

        if len(no_vals) > 0:
            return no_vals
        else:
            return charges
