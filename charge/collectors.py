from typing import List

import networkx as nx
import numpy as np

from charge.base import Charger, C
from charge.settings import MAX_BINS


class AtomicMeanCollector(Charger):

    def _collect_values(self, graph: nx.Graph, iacm_only: bool, shells: List[int], **kwargs) -> C:

        charges = dict()
        no_vals = list()

        for atom in graph.nodes():
            for shell in shells:
                key = self._nauty.canonize_neighborhood(graph, atom, shell,
                                                        color_key='iacm' if 'iacm' in graph.node[atom] else 'atom_type')
                if key in self._repo.charges_iacm[shell]:
                    values = self._repo.charges_iacm[shell][key]
                    charges[atom] = (round(float(np.mean(values)), self._rounding_digits), 1.0)
                    break
                elif not iacm_only:
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

    def _collect_values(self, graph: nx.Graph, iacm_only: bool, shells: List[int], **kwargs) -> C:

        charges = dict()
        no_vals = list()

        max_bins = max(int(kwargs['max_bins']), 1) if 'max_bins' in kwargs else MAX_BINS

        def process_vals(values, atom):
            hist, bin_edges = np.histogram(values)
            if len(hist) > max_bins:
                hist, bin_edges = np.histogram(values, bins=max_bins)

            hist = hist.tolist()
            bin_edges = bin_edges.tolist()

            atom_charges = [round(bin_edges[i] + 0.5 * (bin_edges[i + 1] - bin_edges[i]), self._rounding_digits) \
                       for i in range(len(bin_edges) - 1)]
            filtered_hist, filtered_charges = [hist[0]], [atom_charges[0]]
            for count, charge in zip(hist[1:], atom_charges[1:]):
                if count > 0:
                    if charge == filtered_charges[-1]:
                        filtered_hist[-1] += count
                    else:
                        filtered_charges.append(charge)
                        filtered_hist.append(count)
            
            if 'scale' in kwargs and kwargs['scale']:
                s = sum(filtered_hist)
                mean_val = np.mean(values)
                filtered_hist = [y/s/(1 + (filtered_charges[i] - mean_val)**2) for i, y in enumerate(filtered_hist)]

            charges[atom] = (filtered_charges, filtered_hist)

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

        if len(no_vals) > 0:
            return no_vals
        else:
            return charges
