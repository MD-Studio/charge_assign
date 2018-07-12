import warnings
from typing import List, TypeVar, Dict

import networkx as nx

from charge import util
from charge.nauty import Nauty
from charge.repository import Repository
from charge.settings import ROUNDING_DIGITS, DEFAULT_TOTAL_CHARGE, MAX_ROUNDING_DIGITS

S = TypeVar('S', int, List[int])
C = TypeVar('C', Dict, List)


class AssignmentError(Warning):
    pass


class Charger:

    def __init__(self,
                 repository:Repository=None,
                 nauty:Nauty=None,
                 rounding_digits:int = ROUNDING_DIGITS,
                 **kwargs) -> None:
        self._nauty = nauty or Nauty()
        self._repo = repository or Repository(nauty=self._nauty)
        self._rounding_digits = min(max(rounding_digits, 0), MAX_ROUNDING_DIGITS)

    def _assign_partial_charges(self,
                                graph: nx.Graph,
                                values: Dict,
                                total_charge: int,
                                **kwargs) -> bool:
        pass

    def _collect_values(self,
                        graph: nx.Graph,
                        iacm_only: bool,
                        shells: List[int],
                        **kwargs) -> C:
        pass

    def charge(self,
               graph: nx.Graph,
               iacmize:bool=False,
               iacm_only:bool=False,
               shell:S=None,
               total_charge:int=DEFAULT_TOTAL_CHARGE,
               **kwargs) -> bool:

        if not shell:
            shells = sorted(self._repo.charges_iacm.keys(), reverse=True)
        elif isinstance(shell, int):
            shells = [shell]
        elif isinstance(shell, List[int]):
            shells = shell
        else:
            raise ValueError('shell must be int or List[int]')

        if iacmize:
            graph = util.iacmize(graph)
        
        values = self._collect_values(graph, iacm_only or iacmize, shells, **kwargs)

        if isinstance(values, List):
            err = 'Could find charges for atoms {0}.'.format(', '.join(map(str, values)))
            if not 0 in shells:
                err += ' Please retry with a smaller "shell" parameter.'
            warnings.warn(AssignmentError(err))
            return False
        elif not self._assign_partial_charges(graph, values, total_charge, **kwargs):
            err = 'Could not assign charges. Please retry with a larger "total_charge_diff" parameter ' \
                  'or a SimpleCharger.'
            warnings.warn(AssignmentError(err))
            return False
        else:

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
                min_atom_data['partial_charge_redist'] = round(min_atom_data['partial_charge_redist'] + total_charge - total_charge_redist,
                                                               self._rounding_digits)
                total_charge_redist += total_charge - total_charge_redist
                graph.graph['total_charge_redist'] = round(total_charge_redist, self._rounding_digits)

            return True
