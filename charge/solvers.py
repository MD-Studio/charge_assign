import itertools
from time import perf_counter
from typing import Dict

import networkx as nx
from pulp import LpVariable, LpInteger, LpMaximize, LpProblem, LpStatusOptimal, CPLEX_CMD, GUROBI_CMD, PULP_CBC_CMD, \
    GLPK_CMD, COIN_CMD

from charge.base import Charger
from charge.nauty import Nauty
from charge.repository import Repository
from charge.settings import DEFAULT_TOTAL_CHARGE_DIFF, ROUNDING_DIGITS, ILP_SOLVER_MAX_SECONDS


class SimpleSolver(Charger):

    def _assign_partial_charges(self, graph: nx.Graph, values: Dict, total_charge: int, **kwargs) -> bool:

        solutionTime = -perf_counter()

        profit = 0
        charge = 0
        for atom, (pc, score) in values.items():
            graph.node[atom]['partial_charge'] = pc
            graph.node[atom]['score'] = score
            charge += pc
            profit += score

        solutionTime += perf_counter()

        graph.graph['total_charge'] = round(charge, self._rounding_digits)
        graph.graph['score'] = profit
        graph.graph['time'] = solutionTime
        graph.graph['items'] = 'nan'
        graph.graph['scaled_capacity'] = 'nan'

        return True


class ILPSolver(Charger):

    def __init__(self,
                 repository: Repository = None,
                 nauty: Nauty = None,
                 rounding_digits: int = ROUNDING_DIGITS,
                 **kwargs) -> None:
        super().__init__(repository, nauty, rounding_digits, **kwargs)

        max_seconds = int(kwargs['max_seconds']) if 'max_seconds' in kwargs else ILP_SOLVER_MAX_SECONDS

        if CPLEX_CMD().available():
            self.__solver = CPLEX_CMD(timelimit=max_seconds)
        elif GUROBI_CMD().available():
            self.__solver = GUROBI_CMD(options={'timeLimit':max_seconds})
        elif PULP_CBC_CMD().available():
            self.__solver = PULP_CBC_CMD(maxSeconds=max_seconds)
        elif GLPK_CMD().available():
            self.__solver = GLPK_CMD(options=['--tmlim %d' % max_seconds])
        elif COIN_CMD().available():
            self.__solver = COIN_CMD(maxSeconds=max_seconds)
        else:
            raise RuntimeError('No solver found.')

    def _assign_partial_charges(self, graph: nx.Graph, values: Dict, total_charge: int, **kwargs) -> bool:

        total_charge_diff = float(kwargs['total_charge_diff']) \
            if 'total_charge_diff' in kwargs else DEFAULT_TOTAL_CHARGE_DIFF

        atom_idx = dict()
        idx = list()
        # weights = partial charges
        weights = dict()
        # profits = frequencies
        profits = dict()
        pos_total = total_charge
        for k, (atom, (charges, frequencies)) in enumerate(values.items()):
            atom_idx[k] = atom
            idx.append(list(zip(itertools.repeat(k), range(len(charges)))))
            weights[k] = charges
            profits[k] = frequencies
            pos_total -= min(charges)

        x = LpVariable.dicts('x', itertools.chain.from_iterable(idx), lowBound=0, upBound=1, cat=LpInteger)

        charging_problem = LpProblem("Atomic Charging Problem", LpMaximize)

        # maximize profits
        charging_problem += sum([profits[k][i] * x[(k, i)] for k, i in itertools.chain.from_iterable(idx)])

        # select exactly one item per set
        for indices in idx:
            charging_problem += sum([x[(k, i)] for k, i in indices]) == 1

        # total charge difference
        charging_problem +=\
            sum([weights[k][i] * x[(k, i)] for k, i in itertools.chain.from_iterable(idx)]) - total_charge\
            <= total_charge_diff
        charging_problem +=\
            sum([weights[k][i] * x[(k, i)] for k, i in itertools.chain.from_iterable(idx)]) - total_charge\
            >= -total_charge_diff
	
        solutionTime = -perf_counter()
        try:
            charging_problem.solve(solver=self.__solver)
        except:
            return False
        solutionTime += perf_counter()

        if not charging_problem.status == LpStatusOptimal:
            return False

        solution = []
        profit = 0
        charge = 0
        for k, i in itertools.chain.from_iterable(idx):
            if x[(k, i)].value() == 1.0:
                graph.nodes[atom_idx[k]]['partial_charge'] = weights[k][i]
                graph.nodes[atom_idx[k]]['score'] = profits[k][i]
                solution.append(i)
                profit += profits[k][i]
                charge += graph.nodes[atom_idx[k]]['partial_charge']

        graph.graph['total_charge'] = round(charge, self._rounding_digits)
        graph.graph['score'] = profit
        graph.graph['time'] = solutionTime
        graph.graph['items'] = len(x)
        graph.graph['scaled_capacity'] = pos_total + total_charge_diff

        return True


class DPSolver(Charger):

    def _assign_partial_charges(self, graph: nx.Graph, values: Dict, total_charge: int, **kwargs) -> bool:

        total_charge_diff = float(kwargs['total_charge_diff']) \
            if 'total_charge_diff' in kwargs else DEFAULT_TOTAL_CHARGE_DIFF

        blowup = 10 ** self._rounding_digits
        deflate = 10 ** (-self._rounding_digits)

        atom_idx = dict()
        idx = list()
        # item = (index, weight, profit)
        items = list()
        # min weights
        w_min = dict()

        solutionTime = -perf_counter()

        pos_total_charge = total_charge
        for k, (atom, (charges, frequencies)) in enumerate(values.items()):
            atom_idx[k] = atom
            idx.append(zip(itertools.repeat(k), range(len(charges))))
            w_min[k] = min(charges)
            items.append(list(zip(range(len(charges)),
                             [round(blowup * (charge - w_min[k])) for charge in charges],
                             frequencies)))
            pos_total_charge -= w_min[k]

        upper = round(blowup * (pos_total_charge + total_charge_diff))
        lower = round(blowup * (pos_total_charge - total_charge_diff))

        dp = [0] + [-float('inf')] * upper
        tb = [[] for _ in range(upper + 1)]

        for items_l in items:
            for d in range(upper, -1, -1):
                try:
                    idx, dp[d] = max(
                        ((item[0], dp[d - item[1]] + item[2]) for item in items_l if d - item[1] >= 0),
                        key=lambda x: x[1]
                    )
                    if dp[d] >= 0:
                        tb[d] = tb[d - items_l[idx][1]] + [idx]
                except ValueError:
                    dp[d] = -float('inf')

        max_pos, max_val = max(enumerate(dp[lower:upper + 1]), key=lambda x: x[1])

        solutionTime += perf_counter()

        if max_val == -float('inf'):
            return False

        solution = tb[lower + max_pos]

        charge = 0
        for i, j in enumerate(solution):
            graph.node[atom_idx[i]]['partial_charge'] = round((deflate * items[i][j][1]) + w_min[i],
                                                              self._rounding_digits)
            graph.node[atom_idx[i]]['score'] = items[i][j][2]
            charge += graph.node[atom_idx[i]]['partial_charge']

        graph.graph['total_charge'] = round(charge, self._rounding_digits)
        graph.graph['score'] = max_val
        graph.graph['time'] = solutionTime
        graph.graph['items'] = sum(len(i) for i in items)
        graph.graph['scaled_capacity'] = pos_total_charge + total_charge_diff

        return True


class CDPSolver(Charger):

    def _assign_partial_charges(self, graph: nx.Graph, values: Dict, total_charge: int, **kwargs) -> bool:

        import charge.c.dp as dp

        total_charge_diff = float(kwargs['total_charge_diff']) \
            if 'total_charge_diff' in kwargs else DEFAULT_TOTAL_CHARGE_DIFF

        num_sets = len(values)
        num_items = sum(len(charges) for (_, (charges, _)) in values.items())

        atom_idx = dict()

        weights = dp.new_doublea(num_items)
        profits = dp.new_doublea(num_items)
        sets = dp.new_ushorta(num_sets)
        solution = dp.new_ushorta(num_sets)

        offset = 0
        pos_total = total_charge
        for k, (atom, (charges, frequencies)) in enumerate(values.items()):
            atom_idx[k] = atom
            dp.ushorta_setitem(sets, k, len(charges))
            for i, (charge, frequency) in enumerate(zip(charges, frequencies)):
                dp.doublea_setitem(weights, offset + i, charge)
                dp.doublea_setitem(profits, offset + i, frequency)
            pos_total -= min(charges)
            offset += len(charges)

        solutionTime = -perf_counter()

        profit = dp.solve_dp(weights, profits,
                             sets, num_items, num_sets,
                             self._rounding_digits, total_charge, total_charge_diff,
                             solution)

        solutionTime += perf_counter()

        if profit >= 0:
            charge = 0
            offset = 0
            for k in range(num_sets):
                i = dp.ushorta_getitem(solution, k)
                graph.node[atom_idx[k]]['partial_charge'] = dp.doublea_getitem(weights, offset + i)
                graph.node[atom_idx[k]]['score'] = dp.doublea_getitem(profits, offset + i)
                charge += graph.node[atom_idx[k]]['partial_charge']
                offset += dp.ushorta_getitem(sets, k)

            graph.graph['total_charge'] = round(charge, self._rounding_digits)
            graph.graph['score'] = profit
            graph.graph['time'] = solutionTime
            graph.graph['items'] = num_items
            graph.graph['scaled_capacity'] = pos_total + total_charge_diff

        dp.delete_doublea(weights)
        dp.delete_doublea(profits)
        dp.delete_ushorta(sets)
        dp.delete_ushorta(solution)

        return profit >= 0
