import warnings

import networkx as nx
import numpy

from charge.assign import Charger, AssignmentError


class AtomicCharger(Charger):

    def _set_partial_charges(self, graph: nx.Graph, iacm_only: bool,
                             shell: int, rounding_digits: int,  **kwargs) -> bool:
        shells = sorted(self._repo.charges_iacm.keys(), reverse=True) if shell < 0 else [shell]
        rounding_digits = max(rounding_digits, 0)
        
        def assign(atom):
            for shell in shells:
                key = self._nauty.canonize_neighborhood(graph, atom, shell,
                                                        color_key='iacm' if 'iacm' in graph.node[atom] else 'atom_type')
                if key in self._repo.charges_iacm[shell]:
                    values = self._repo.charges_iacm[shell][key]
                    graph.node[atom]['partial_charge'] = round(float(numpy.mean(values)), rounding_digits)
                    graph.node[atom]['partial_charge_std'] = numpy.std(values)
                    break
                elif not iacm_only:
                    key = self._nauty.canonize_neighborhood(graph, atom, shell)
                    if key in self._repo.charges_elem[shell]:
                        values = self._repo.charges_elem[shell][key]
                        graph.node[atom]['partial_charge'] = round(float(numpy.mean(values)), rounding_digits)
                        graph.node[atom]['partial_charge_std'] = numpy.std(values)
                        break
            else:
                warnings.warn(AssignmentError('Could not assign charge to atom {0}'.format(atom)))
                graph.node[atom]['partial_charge'] = float('nan')
                graph.node[atom]['partial_charge_std'] = float('nan')
                return False

            return True

        return all([assign(atom) for atom in graph.nodes()])
