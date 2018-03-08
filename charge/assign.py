
import networkx as nx
from pulp import PULP_CBC_CMD, GLPK_CMD, COIN_CMD, CPLEX_CMD, GUROBI_CMD

from charge.bond_type import BondType
from charge import util
from charge.nauty import Nauty
from charge.repository import Repository
from charge.settings import SOLVER_MAX_SECONDS


class AssignmentError(Warning):
    pass


class Charger:

    def __init__(self, repository:Repository=None, nauty:Nauty=None) -> None:
        self._nauty = nauty or Nauty()
        self._repo = repository or Repository(nauty=self._nauty)

    def _set_partial_charges(self, graph: nx.Graph, iacm_only: bool, shell: int, **kwargs) -> bool:
        pass

    def charge(self, graph: nx.Graph, iacmize:bool=False, iacm_only:bool=False, shell:int=None, **kwargs) -> bool:
        if iacmize:
            graph = util.iacmize(graph)
        return self._set_partial_charges(graph, iacm_only=iacm_only or iacmize,
                                         shell=shell if shell and shell >= 0 else -1, **kwargs)
