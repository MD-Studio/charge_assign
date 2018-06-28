import itertools
import warnings

import networkx as nx
import numpy as np
from pulp import LpProblem, LpMaximize, LpVariable, LpInteger, CPLEX_CMD, GUROBI_CMD, PULP_CBC_CMD, \
    GLPK_CMD, COIN_CMD, LpStatusOptimal

from charge.assign import Charger, AssignmentError
from charge.nauty import Nauty
from charge.repository import Repository
from charge.settings import DEFAULT_TOTAL_CHARGE, DEFAULT_TOTAL_CHARGE_DIFF, SOLVER_MAX_SECONDS, MAX_BINS


class AtomicILPCharger(Charger):

    def __init__(self, repository: Repository = None, nauty: Nauty = None, max_seconds:int = None) -> None:
        super(AtomicILPCharger, self).__init__(repository=repository, nauty=nauty)

        max_seconds = max_seconds and max_seconds > 0 or SOLVER_MAX_SECONDS

        # TODO GUROBI and CPLEX time limits (+ other useful options?)
        if CPLEX_CMD().available():
            self.__solver = CPLEX_CMD(options=[''])
        elif GUROBI_CMD().available():
            self.__solver = GUROBI_CMD(options=[''])
        elif PULP_CBC_CMD().available():
            self.__solver = PULP_CBC_CMD(maxSeconds=max_seconds)
        elif GLPK_CMD().available():
            self.__solver = GLPK_CMD(options=['--tmlim %d' % max_seconds])
        elif COIN_CMD().available():
            self.__solver = COIN_CMD(maxSeconds=max_seconds)
        else:
            raise RuntimeError('No solver found.')

    def _set_partial_charges(self, graph: nx.Graph, iacm_only: bool,
                             shell: int, rounding_digits: int, **kwargs) -> bool:
        shells = sorted(self._repo.charges_iacm.keys(), reverse=True) if shell < 0 else [shell]
        rounding_digits = max(rounding_digits, 0)

        w = dict()
        c = dict()
        idx = dict()

        total_charge = int(kwargs['total_charge']) if 'total_charge' in kwargs else DEFAULT_TOTAL_CHARGE
        total_charge_diff = float(kwargs['total_charge_diff'])\
            if 'total_charge_diff' in kwargs else DEFAULT_TOTAL_CHARGE_DIFF
        max_bins = int(kwargs['max_bins']) if 'max_bins' in kwargs else MAX_BINS

        def process_vals(values, atom):
            hist, bin_edges = np.histogram(values)
            if len(hist) > max_bins:
                hist, bin_edges = np.histogram(values, bins=max_bins)

            charges = [round(bin_edges[i] + 0.5 * (bin_edges[i + 1] - bin_edges[i]), rounding_digits)\
                       for i in range(len(bin_edges) - 1)]
            filtered_hist, filtered_charges = [hist[0]], [charges[0]]
            for count, charge in zip(hist[1:], charges[1:]):
                if count > 0:
                    if charge == filtered_charges[-1]:
                        filtered_hist[-1] += count
                    else:
                        filtered_charges.append(charge)
                        filtered_hist.append(count)

            w[atom] = [count / sum(filtered_hist) for count in filtered_hist]
            c[atom] = filtered_charges
            idx[atom] = list(range(len(filtered_hist)))

        def assign(atom):
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
                warnings.warn(AssignmentError('Could not assign charge to atom {0}'.format(atom)))
                return False

            return True

        if not all([assign(atom) for atom in graph.nodes()]):
            return False

        idx_list = list(itertools.chain.from_iterable(
            [[(atom, i) for i in atom_idx] for atom, atom_idx in idx.items()]
        ))

        x = LpVariable.dicts('x', idx_list, lowBound=0, upBound=1, cat=LpInteger)

        charging_problem = LpProblem("Atomic Charging Problem", LpMaximize)

        # maximize charge frequency
        charging_problem += sum([w[a][i]*x[(a, i)] for a, i in idx_list])

        # select exactly one charge per atom
        for atom, atom_idx in idx.items():
            charging_problem += sum([x[(atom, i)] for i in atom_idx]) == 1

        
        # total charge difference
        charging_problem += sum([c[a][i] * x[(a, i)] for a, i in idx_list]) - total_charge <= total_charge_diff
        charging_problem += sum([c[a][i] * x[(a, i)] for a, i in idx_list]) - total_charge >= -total_charge_diff

        charging_problem.solve(solver=self.__solver)

        if not charging_problem.status == LpStatusOptimal:
            return False

        for a, a_i in idx.items():
            for i in a_i:
                if x[(a, i)].value() == 1.0:
                    graph.nodes[a]['partial_charge'] = c[a][i]

        return True
