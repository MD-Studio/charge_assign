import itertools
from abc import ABC, abstractmethod
from time import perf_counter
from typing import Dict, Tuple, List

import networkx as nx
from pulp import LpVariable, LpInteger, LpMaximize, LpProblem, LpStatusOptimal, CPLEX_CMD, GUROBI_CMD, PULP_CBC_CMD, \
    GLPK_CMD, COIN_CMD

from collections import defaultdict

from charge.charge_types import Atom, ChargeList, WeightList
from charge.settings import DEFAULT_TOTAL_CHARGE_DIFF, ROUNDING_DIGITS, ILP_SOLVER_MAX_SECONDS
from charge.util import AssignmentError


class Solver(ABC):
    """Base class for solvers.

    Solvers assign charges to atoms based on charge distribution \
    histograms obtained by a Collector.
    """
    @abstractmethod
    def solve_partial_charges(
            self,
            graph: nx.Graph,
            charge_dists: Dict[Atom, Tuple[ChargeList, WeightList]],
            total_charge: int,
            keydict: Dict[Atom, str],
            **kwargs
            ) -> None:
        """Assign charges to the atoms in a graph.

        Modify a graph by adding additional attributes describing the \
        atoms' charges and scores. In particular, each atom will get \
        a 'partial_charge' attribute with the partial charge, and a \
        'score' attribute giving a degree of certainty for that charge.

        Args:
            graph: The molecule graph to solve charges for.
            charge_dists: Charge distributions for the atoms, obtained \
                    by a Collector.
            total_charge: The total charge of the molecule.
            keydict: Dictionary containing nauty canonical keys of atom \
                    neighborhoods for each atom of the graph \
                    only used in Symmetric Solvers
        """
        pass


    def compute_atom_neighborhood_classes(self, atom_idx : dict, keydict : dict):
        L = defaultdict(list)

        for atom in atom_idx:
            L[keydict[atom_idx[atom]]].append(atom)

        return [L[l] for l in L]


class SimpleSolver(Solver):
    """A trivial solver that assigns a single statistics of the found charges.

    Use the MeanCollector, MedianCollector or ModeCollector to produce \
    appropriate charge distributions.
    """
    def __init__(self, rounding_digits: int) -> None:
        self.__rounding_digits = rounding_digits

    def solve_partial_charges(
            self,
            graph: nx.Graph,
            charge_dists: Dict[Atom, Tuple[ChargeList, WeightList]],
            total_charge: int,
            keydict: Dict[Atom, str] = None,
            **kwargs
            ) -> None:
        """Assign charges to the atoms in a graph.

        Modify a graph by adding additional attributes describing the \
        atoms' charges and scores. In particular, each atom will get \
        a 'partial_charge' attribute with the partial charge, and a \
        'score' attribute giving a degree of certainty for that charge.

        This solver expects a single charge value for each atom, and \
        assigns it to that atom. With charge distributions from the \
        MeanCollector, it assigns the mean of the found charges.

        Args:
            graph: The molecule graph to solve charges for.
            charge_dists: Charge distributions for the atoms, obtained \
                    by a Collector.
            total_charge: The total charge of the molecule.
        """
        solutionTime = -perf_counter()

        profit = 0
        charge = 0
        for atom, (pcs, scores) in charge_dists.items():
            graph.node[atom]['partial_charge'] = pcs[0]
            graph.node[atom]['score'] = scores[0]
            charge += pcs[0]
            profit += scores[0]

        solutionTime += perf_counter()

        graph.graph['total_charge'] = round(charge, self.__rounding_digits)
        graph.graph['score'] = profit
        graph.graph['time'] = solutionTime
        graph.graph['items'] = 'nan'
        graph.graph['scaled_capacity'] = 'nan'


class ILPSolver(Solver):
    """An optimizing solver using Integer Linear Programming.

    Use the HistogramCollector to produce appropriate charge \
    distributions.
    """

    def __init__(self,
                 rounding_digits: int=ROUNDING_DIGITS,
                 max_seconds: int=ILP_SOLVER_MAX_SECONDS
                 ) -> None:
        """Create an ILPSolver.

        Args:
            rounding_digits: Number of digits to round the charges to.
            max_seconds: Maximum run-time to spend searching for a \
                    solution
        """
        self.__rounding_digits = rounding_digits

        if CPLEX_CMD().available():
            self.__solver = CPLEX_CMD(timelimit=max_seconds)
        elif GUROBI_CMD().available():
            self.__solver = GUROBI_CMD(options=[('timeLimit', max_seconds)])
        elif PULP_CBC_CMD().available():
            self.__solver = PULP_CBC_CMD(maxSeconds=max_seconds)
        elif GLPK_CMD().available():
            self.__solver = GLPK_CMD(options=['--tmlim %d' % max_seconds])
        elif COIN_CMD().available():
            self.__solver = COIN_CMD(maxSeconds=max_seconds)
        else:
            raise RuntimeError('No solver found, there is something'
                    ' wrong with your pulp library setup.')

    def solve_partial_charges(
            self,
            graph: nx.Graph,
            charge_dists: Dict[Atom, Tuple[ChargeList, WeightList]],
            total_charge: int,
            keydict: Dict[Atom, str] = None,
            total_charge_diff: float=DEFAULT_TOTAL_CHARGE_DIFF,
            **kwargs
            ) -> None:
        """Assign charges to the atoms in a graph.

        Modify a graph by adding additional attributes describing the \
        atoms' charges and scores. In particular, each atom will get \
        a 'partial_charge' attribute with the partial charge, and a \
        'score' attribute giving a degree of certainty for that charge.

        This solver formulates the epsilon-Multiple Choice Knapsack \
        Problem as an Integer Linear Programming problem and then uses \
        a generic ILP solver from the pulp library to produce optimised \
        charges.

        Args:
            graph: The molecule graph to solve charges for.
            charge_dists: Charge distributions for the atoms, obtained \
                    by a Collector.
            total_charge: The total charge of the molecule
            total_charge_diff: Maximum allowed deviation from the total charge
        """

        atom_idx = dict()
        idx = list()
        # weights = partial charges
        weights = dict()
        # profits = frequencies
        profits = dict()
        pos_total = total_charge
        for k, (atom, (charges, frequencies)) in enumerate(charge_dists.items()):
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
            raise AssignmentError('Could not solve ILP problem. Please retry'
                    ' with a SimpleCharger')
        solutionTime += perf_counter()

        if not charging_problem.status == LpStatusOptimal:
            raise AssignmentError('Could not solve ILP problem. Please retry'
                    ' with a SimpleCharger')

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

        graph.graph['total_charge'] = round(charge, self.__rounding_digits)
        graph.graph['score'] = profit
        graph.graph['time'] = solutionTime
        graph.graph['items'] = len(x)
        graph.graph['scaled_capacity'] = pos_total + total_charge_diff

class SymmetricILPSolver(Solver):
    """An optimizing solver using Integer Linear Programming.

    Use the HistogramCollector to produce appropriate charge \
    distributions.
    """

    def __init__(self,
                 rounding_digits: int=ROUNDING_DIGITS,
                 max_seconds: int=ILP_SOLVER_MAX_SECONDS
                 ) -> None:
        """Create an ILPSolver.

        Args:
            rounding_digits: Number of digits to round the charges to.
            max_seconds: Maximum run-time to spend searching for a \
                    solution
        """
        self.__rounding_digits = rounding_digits

        if CPLEX_CMD().available():
            self.__solver = CPLEX_CMD(timelimit=max_seconds)
        elif GUROBI_CMD().available():
            self.__solver = GUROBI_CMD(options=[('timeLimit',max_seconds)])
        elif PULP_CBC_CMD().available():
            self.__solver = PULP_CBC_CMD(maxSeconds=max_seconds)
        elif GLPK_CMD().available():
            self.__solver = GLPK_CMD(options=['--tmlim %d' % max_seconds])
        elif COIN_CMD().available():
            self.__solver = COIN_CMD(maxSeconds=max_seconds)
        else:
            raise RuntimeError('No solver found, there is something'
                    ' wrong with your pulp library setup.')

    def solve_partial_charges(
            self,
            graph: nx.Graph,
            charge_dists: Dict[Atom, Tuple[ChargeList, WeightList]],
            total_charge: int,
            keydict: Dict[Atom, str],
            total_charge_diff: float=DEFAULT_TOTAL_CHARGE_DIFF,
            **kwargs
            ) -> None:
        """Assign charges to the atoms in a graph.

        Modify a graph by adding additional attributes describing the \
        atoms' charges and scores. In particular, each atom will get \
        a 'partial_charge' attribute with the partial charge, and a \
        'score' attribute giving a degree of certainty for that charge.

        This solver formulates the epsilon-Multiple Choice Knapsack \
        Problem as an Integer Linear Programming problem and then uses \
        a generic ILP solver from the pulp library to produce optimised \
        charges.

        Args:
            graph: The molecule graph to solve charges for.
            charge_dists: Charge distributions for the atoms, obtained \
                    by a Collector.
            total_charge: The total charge of the molecule.
            keydict: Dictionary containing nauty canonical keys of atom \
                    neighborhoods for each atom of the graph
            total_charge_diff: Maximum allowed deviation from the total charge
        """

        atom_idx = dict()
        idx = list()
        # weights = partial charges
        weights = dict()
        # profits = frequencies
        profits = dict()

        pos_total = total_charge
        for k, (atom, (charges, frequencies)) in enumerate(charge_dists.items()):
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

        #identical neighborhood charge conditions
        neighborhoodclasses = self.compute_atom_neighborhood_classes(atom_idx, keydict)
        for neighborhood_class in neighborhoodclasses:
            i = neighborhood_class[0]
            for j in neighborhood_class[1::]:
                for (_, k) in idx[i]:
                    charging_problem += x[(i, k)] - x[(j, k)] == 0                  # weight k from atom i is selected as partial charge <=> weight k of atom j is selected as partial charge


        solutionTime = -perf_counter()
        try:
            charging_problem.solve(solver=self.__solver)
        except:
            raise AssignmentError('Could not solve ILP problem. Please retry'
                    ' with a SimpleCharger')
        solutionTime += perf_counter()

        if not charging_problem.status == LpStatusOptimal:
            raise AssignmentError('Could not solve ILP problem. Please retry'
                    ' with a SimpleCharger')

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

        graph.graph['total_charge'] = round(charge, self.__rounding_digits)
        graph.graph['score'] = profit
        graph.graph['time'] = solutionTime
        graph.graph['items'] = len(x)
        graph.graph['scaled_capacity'] = pos_total + total_charge_diff
        graph.graph['neighborhoods'] = [[atom_idx[k] for k in i] for i in neighborhoodclasses]


class SymmetricRelaxedILPSolver(Solver):
    """An optimizing solver using Integer Linear Programming.

    Use the HistogramCollector to produce appropriate charge \
    distributions.
    """

    def __init__(self,
                 rounding_digits: int=ROUNDING_DIGITS,
                 max_seconds: int=ILP_SOLVER_MAX_SECONDS
                 ) -> None:
        """Create an ILPSolver.

        Args:
            rounding_digits: Number of digits to round the charges to.
            max_seconds: Maximum run-time to spend searching for a \
                    solution
        """
        self.__rounding_digits = rounding_digits

        if CPLEX_CMD().available():
            self.__solver = CPLEX_CMD(timelimit=max_seconds)
        elif GUROBI_CMD().available():
            self.__solver = GUROBI_CMD(options=[('timeLimit',max_seconds)])
        elif PULP_CBC_CMD().available():
            self.__solver = PULP_CBC_CMD(maxSeconds=max_seconds)
        elif GLPK_CMD().available():
            self.__solver = GLPK_CMD(options=['--tmlim %d' % max_seconds])
        elif COIN_CMD().available():
            self.__solver = COIN_CMD(maxSeconds=max_seconds)
        else:
            raise RuntimeError('No solver found, there is something'
                    ' wrong with your pulp library setup.')

    def solve_partial_charges(
            self,
            graph: nx.Graph,
            charge_dists: Dict[Atom, Tuple[ChargeList, WeightList]],
            total_charge: int,
            keydict: Dict[Atom, str],
            total_charge_diff: float=DEFAULT_TOTAL_CHARGE_DIFF,
            **kwargs
            ) -> None:
        """Assign charges to the atoms in a graph.

        Modify a graph by adding additional attributes describing the \
        atoms' charges and scores. In particular, each atom will get \
        a 'partial_charge' attribute with the partial charge, and a \
        'score' attribute giving a degree of certainty for that charge.

        This solver formulates the epsilon-Multiple Choice Knapsack \
        Problem as an Integer Linear Programming problem and then uses \
        a generic ILP solver from the pulp library to produce optimised \
        charges.

        Args:
            graph: The molecule graph to solve charges for.
            charge_dists: Charge distributions for the atoms, obtained \
                    by a Collector.
            total_charge: The total charge of the molecule.
            keydict: Dictionary containing nauty canonical keys of atom \
                    neighborhoods for each atom of the graph
            total_charge_diff: Maximum allowed deviation from the total charge
        """

        atom_idx = dict()
        idx = list()
        # weights = partial charges
        weights = dict()
        # profits = frequencies
        profits = dict()

        pos_total = total_charge
        for k, (atom, (charges, frequencies)) in enumerate(charge_dists.items()):
            atom_idx[k] = atom
            idx.append(list(zip(itertools.repeat(k), range(len(charges)))))
            weights[k] = charges
            profits[k] = frequencies
            pos_total -= min(charges)

        x = LpVariable.dicts('x', itertools.chain.from_iterable(idx), lowBound=0, upBound=1)

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

        #identical neighborhood charge conditions
        neighborhoodclasses = self.compute_atom_neighborhood_classes(atom_idx, keydict)
        for neighborhood_class in neighborhoodclasses:
            i = neighborhood_class[0]
            for j in neighborhood_class[1::]:
                for (_, k) in idx[i]:
                    charging_problem += x[(i, k)] - x[(j, k)] == 0                  # weight k from atom i is selected as partial charge <=> weight k of atom j is selected as partial charge


        solutionTime = -perf_counter()
        try:
            charging_problem.solve(solver=self.__solver)
        except:
            raise AssignmentError('Could not solve ILP problem. Please retry'
                    ' with a SimpleCharger')
        solutionTime += perf_counter()

        if not charging_problem.status == LpStatusOptimal:
            raise AssignmentError('Could not solve ILP problem. Please retry'
                    ' with a SimpleCharger')

        solution = []
        profit = 0
        charge = 0
        for k, i in itertools.chain.from_iterable(idx):
            if x[(k, i)].value() != 0.0:
                if('partial_charge' in graph.nodes[atom_idx[k]] and 'score' in graph.nodes[atom_idx[k]]):
                    graph.nodes[atom_idx[k]]['partial_charge'] += weights[k][i] * x[(k, i)].value()
                    graph.nodes[atom_idx[k]]['score'] += profits[k][i] * x[(k, i)].value()
                else:
                    graph.nodes[atom_idx[k]]['partial_charge'] = weights[k][i] * x[(k, i)].value()
                    graph.nodes[atom_idx[k]]['score'] = profits[k][i] * x[(k, i)].value()

                solution.append(i)
        for node in graph.nodes:
            profit += graph.nodes[node]['score']
            charge += graph.nodes[node]['partial_charge']

        graph.graph['total_charge'] = round(charge, self.__rounding_digits)
        graph.graph['score'] = profit
        graph.graph['time'] = solutionTime
        graph.graph['items'] = len(x)
        graph.graph['scaled_capacity'] = pos_total + total_charge_diff
        graph.graph['neighborhoods'] = [[atom_idx[k] for k in i] for i in neighborhoodclasses]

class DPSolver(Solver):
    """An optimizing solver using Dynamic Programming.

    Use the HistogramCollector to produce appropriate charge \
    distributions.
    """
    def __init__(self, rounding_digits) -> None:
        """Create a DPSolver.

        Args:
            rounding_digits: How many significant digits to round the \
                    resulting charges to.
        """
        self.__rounding_digits = rounding_digits


    def solve_partial_charges(
            self,
            graph: nx.Graph,
            charge_dists: Dict[Atom, Tuple[ChargeList, WeightList]],
            total_charge: int,
            keydict: Dict[Atom, str] = None,
            total_charge_diff: float=DEFAULT_TOTAL_CHARGE_DIFF,
            **kwargs
            ) -> None:
        """Assign charges to the atoms in a graph.

        Modify a graph by adding additional attributes describing the \
        atoms' charges and scores. In particular, each atom will get \
        a 'partial_charge' attribute with the partial charge, and a \
        'score' attribute giving a degree of certainty for that charge.

        This solver uses Dynamic Programming to solve the \
        epsilon-Multiple Choice Knapsack Problem. This is the Python \
        version of the algorithm, see CDPSolver for a faster \
        implementation.

        Args:
            graph: The molecule graph to solve charges for.
            charge_dists: Charge distributions for the atoms, obtained \
                    by a Collector.
            total_charge: The total charge of the molecule.
            total_charge_diff: Maximum allowed deviation from the total charge
        """

        blowup = 10 ** self.__rounding_digits
        deflate = 10 ** (-self.__rounding_digits)

        atom_idx = dict()
        idx = list()
        # item = (index, weight, profit)
        items = list()
        # min weights
        w_min = dict()

        solutionTime = -perf_counter()

        # transform weights to non-negative integers
        pos_total_charge = total_charge
        max_sum = 0
        for k, (atom, (charges, frequencies)) in enumerate(charge_dists.items()):
            atom_idx[k] = atom
            idx.append(zip(itertools.repeat(k), range(len(charges))))
            w_min[k] = min(charges)
            max_sum += max(charges) - w_min[k]
            items.append(list(zip(range(len(charges)),
                             [round(blowup * (charge - w_min[k])) for charge in charges],
                             frequencies)))
            pos_total_charge -= w_min[k]

        # lower and upper capacity limits
        upper = round(blowup * (pos_total_charge + total_charge_diff))
        lower = max(0, round(blowup * (pos_total_charge - total_charge_diff)))

        # check if feasible solutions may exist
        reachable = round(blowup * max_sum)
        if upper < 0 or lower > reachable:
            # sum of min weights over all sets is larger than the upper bound
            # or sum of max weights over all sets is smaller than the lower bound
            raise AssignmentError('Could not solve DP problem. Please retry'
                                  ' with a SimpleCharger')

        # init DP and traceback tables
        dp = [0] + [-float('inf')] * upper
        tb = [[] for _ in range(upper + 1)]

        # DP
        # iterate over all sets
        for items_l in items:
            # iterate over all capacities
            for d in range(upper, -1, -1):
                try:
                    # find max. profit with capacity i over all items j in set k
                    idx, dp[d] = max(
                        ((item[0], dp[d - item[1]] + item[2]) for item in items_l if d - item[1] >= 0),
                        key=lambda x: x[1]
                    )
                    # copy old traceback indices and add new index to traceback
                    if dp[d] >= 0:
                        tb[d] = tb[d - items_l[idx][1]] + [idx]
                except ValueError:
                    dp[d] = -float('inf')

        # find max profit
        max_pos, max_val = max(enumerate(dp[lower:upper + 1]), key=lambda x: x[1])

        solutionTime += perf_counter()

        if max_val == -float('inf'):
            raise AssignmentError('Could not solve DP problem. Please retry'
                    ' with a SimpleCharger')

        solution = tb[lower + max_pos]

        charge = 0
        score = 0
        for i, j in enumerate(solution):
            graph.node[atom_idx[i]]['partial_charge'] = charge_dists[atom_idx[i]][0][j]
            graph.node[atom_idx[i]]['score'] = charge_dists[atom_idx[i]][1][j]
            charge += graph.node[atom_idx[i]]['partial_charge']
            score += graph.node[atom_idx[i]]['score']

        graph.graph['total_charge'] = round(charge, self.__rounding_digits)
        graph.graph['score'] = score
        graph.graph['time'] = solutionTime
        graph.graph['items'] = sum(len(i) for i in items)
        graph.graph['scaled_capacity'] = pos_total_charge + total_charge_diff

class SymmetricDPSolver(Solver):
    """An optimizing solver using Dynamic Programming.

    Use the HistogramCollector to produce appropriate charge \
    distributions.
    """
    def __init__(self,
                 rounding_digits) -> None:
        """Create a DPSolver.

        Args:
            rounding_digits: How many significant digits to round the \
                    resulting charges to.
        """
        self.__rounding_digits = rounding_digits


    def solve_partial_charges(
            self,
            graph: nx.Graph,
            charge_dists_collector: Dict[Atom, Tuple[ChargeList, WeightList]],
            total_charge: int,
            keydict: Dict[Atom, str],
            total_charge_diff: float=DEFAULT_TOTAL_CHARGE_DIFF,
            **kwargs
            ) -> None:
        """Assign charges to the atoms in a graph.

        Modify a graph by adding additional attributes describing the \
        atoms' charges and scores. In particular, each atom will get \
        a 'partial_charge' attribute with the partial charge, and a \
        'score' attribute giving a degree of certainty for that charge.

        This solver uses Dynamic Programming to solve the \
        epsilon-Multiple Choice Knapsack Problem. This is the Python \
        version of the algorithm, see CDPSolver for a faster \
        implementation.

        Args:
            graph: The molecule graph to solve charges for.
            charge_dists_collector: Charge distributions for the atoms, obtained \
                    by a Collector.
            total_charge: The total charge of the molecule.
            keydict: Dictionary containing nauty canonical keys of atom \
                    neighborhoods for each atom of the graph
            total_charge_diff: Maximum allowed deviation from the total charge
        """

        blowup = 10 ** self.__rounding_digits

        atom_idx = dict()

        for k, (atom, (_, _)) in enumerate(charge_dists_collector.items()):
            atom_idx[k] = atom
        neighborhoodclasses = self.compute_atom_neighborhood_classes(atom_idx, keydict)

        # reduce charge distributions (charges and frequencies of atoms within one neighborhoodclass get combined)
        charge_dists_reduced = self.reduce_charge_distributions(charge_dists_collector, atom_idx, neighborhoodclasses)

        items, pos_total_charge, max_sum = self.transform_weights(charge_dists_reduced, total_charge, blowup)

        solution, solutionTime = self.solve_dp(items, total_charge_diff, pos_total_charge, max_sum, blowup)

        charge = 0
        score = 0
        for i, j in enumerate(solution):
            for k in neighborhoodclasses[i]:
                graph.node[atom_idx[k]]['partial_charge'] = charge_dists_collector[atom_idx[k]][0][j]
                graph.node[atom_idx[k]]['score'] = charge_dists_collector[atom_idx[k]][1][j]
                charge += graph.node[atom_idx[k]]['partial_charge']
                score += graph.node[atom_idx[k]]['score']

        graph.graph['total_charge'] = round(charge, self.__rounding_digits)
        graph.graph['score'] = score
        graph.graph['time'] = solutionTime
        graph.graph['items'] = sum(len(i) for i in items)
        graph.graph['scaled_capacity'] = pos_total_charge + total_charge_diff
        graph.graph['neighborhoods'] = [[atom_idx[k] for k in i] for i in neighborhoodclasses]

    def transform_weights(self, charge_dists, total_charge, blowup):
        atom_idx = dict()
        # item = (index, weight, profit)
        items = list()
        # min weights
        w_min = dict()

        # transform weights to non-negative integers
        pos_total_charge = total_charge
        max_sum = 0
        for k, (atom, (charges, frequencies)) in enumerate(charge_dists.items()):
            atom_idx[k] = atom
            w_min[k] = min(charges)
            max_sum += max(charges) - w_min[k]
            items.append(list(zip(range(len(charges)),
                                  [round(blowup * (charge - w_min[k])) for charge in charges],
                                  frequencies)))
            pos_total_charge -= w_min[k]

        return items, pos_total_charge, max_sum


    def solve_dp(self, items, total_charge_diff, pos_total_charge, max_sum, blowup):


        solutionTime = -perf_counter()
        # lower and upper capacity limits
        upper = round(blowup * (pos_total_charge + total_charge_diff))
        lower = max(0, round(blowup * (pos_total_charge - total_charge_diff)))

        # check if feasible solutions may exist
        reachable = round(blowup * max_sum)
        if upper < 0 or lower > reachable:
            # sum of min weights over all sets is larger than the upper bound
            # or sum of max weights over all sets is smaller than the lower bound
            raise AssignmentError('Could not solve DP problem. Please retry'
                                  ' with a SimpleCharger')

        # init DP and traceback tables
        dp = [0] + [-float('inf')] * upper
        tb = [[] for _ in range(upper + 1)]

        # DP
        # iterate over all sets
        for items_l in items:
            # iterate over all capacities
            for d in range(upper, -1, -1):
                try:
                    # find max. profit with capacity i over all items j in set k
                    idx, dp[d] = max(
                        ((item[0], dp[d - item[1]] + item[2]) for item in items_l if d - item[1] >= 0),
                        key=lambda x: x[1]
                    )
                    # copy old traceback indices and add new index to traceback
                    if dp[d] >= 0:
                        tb[d] = tb[d - items_l[idx][1]] + [idx]
                except ValueError:
                    dp[d] = -float('inf')

        # find max profit
        max_pos, max_val = max(enumerate(dp[lower:upper + 1]), key=lambda x: x[1])

        solutionTime += perf_counter()

        if max_val == -float('inf'):
            raise AssignmentError('Could not solve DP problem. Please retry'
                                  ' with a SimpleCharger')

        solution = tb[lower + max_pos]

        return solution, solutionTime

    def reduce_charge_distributions(self, charge_dists_collector, atom_idx, neighborhoodclasses):
        charge_dists = dict()
        for neighborhoodclass in neighborhoodclasses:
            i = neighborhoodclass[0]
            k = len(neighborhoodclass)
            (charges, frequencies) = charge_dists_collector[atom_idx[i]]
            charge_dists[atom_idx[i]] = ([k * x for x in charges], [k * x for x in frequencies])

        return charge_dists


class CDPSolver(Solver):
    """An optimizing solver using Dynamic Programming, C version.

    Use the HistogramCollector to produce appropriate charge \
    distributions.
    """
    def __init__(self, rounding_digits) -> None:
        """Create a CDPSolver.

        Args:
            rounding_digits: How many significant digits to round the \
                    resulting charges to.
        """
        self.__rounding_digits = rounding_digits

    def solve_partial_charges(
            self,
            graph: nx.Graph,
            charge_dists: Dict[Atom, Tuple[ChargeList, WeightList]],
            total_charge: int,
            keydict: Dict[Atom, str] = None,
            total_charge_diff: float=DEFAULT_TOTAL_CHARGE_DIFF,
            **kwargs
            ) -> None:
        """Assign charges to the atoms in a graph.

        Modify a graph by adding additional attributes describing the \
        atoms' charges and scores. In particular, each atom will get \
        a 'partial_charge' attribute with the partial charge, and a \
        'score' attribute giving a degree of certainty for that charge.

        This solver uses Dynamic Programming to solve the \
        epsilon-Multiple Choice Knapsack Problem. This is the Python \
        version of the algorithm, see DPSolver for the Python \
        implementation.

        Args:
            graph: The molecule graph to solve charges for.
            charge_dists: Charge distributions for the atoms, obtained \
                    by a Collector.
            total_charge: The total charge of the molecule.
            total_charge_diff: Maximum allowed deviation from the total charge
        """

        atom_idx = dict()
        for k, (atom, (_, _)) in enumerate(charge_dists.items()):
            atom_idx[k] = atom

        solution, solutionTime, num_items, scaled_capacity = self.solve_dp_c(charge_dists, total_charge, total_charge_diff)

        charge = 0
        profit = 0
        for (i,j) in solution:
            graph.node[atom_idx[i]]['partial_charge'] = charge_dists[atom_idx[i]][0][j]
            graph.node[atom_idx[i]]['score'] = charge_dists[atom_idx[i]][1][j]
            charge += graph.node[atom_idx[i]]['partial_charge']
            profit += graph.node[atom_idx[i]]['score']

        graph.graph['total_charge'] = round(charge, self.__rounding_digits)
        graph.graph['score'] = profit
        graph.graph['time'] = solutionTime
        graph.graph['items'] = num_items
        graph.graph['scaled_capacity'] = scaled_capacity


    def solve_dp_c(self, charge_dists, total_charge, total_charge_diff):

        import charge.c.dp as dp

        num_sets = len(charge_dists)
        num_items = sum(len(charges) for (_, (charges, _)) in charge_dists.items())

        weights = dp.new_doublea(num_items)
        profits = dp.new_doublea(num_items)
        sets = dp.new_ushorta(num_sets)
        solution = dp.new_ushorta(num_sets)

        offset = 0
        pos_total = total_charge
        for k, (atom, (charges, frequencies)) in enumerate(charge_dists.items()):
            dp.ushorta_setitem(sets, k, len(charges))
            for i, (charge, frequency) in enumerate(zip(charges, frequencies)):
                dp.doublea_setitem(weights, offset + i, charge)
                dp.doublea_setitem(profits, offset + i, frequency)
            pos_total -= min(charges)
            offset += len(charges)

        solutionTime = -perf_counter()

        profit = dp.solve_dp(weights, profits,
                             sets, num_items, num_sets,
                             self.__rounding_digits, total_charge, total_charge_diff,
                             solution)

        solutionTime += perf_counter()

        dp_solution = list()
        if profit >= 0:
            for k in range(num_sets):
                i = dp.ushorta_getitem(solution, k)
                dp_solution.append((k, i))


        dp.delete_doublea(weights)
        dp.delete_doublea(profits)
        dp.delete_ushorta(sets)
        dp.delete_ushorta(solution)

        if profit < 0:
            raise AssignmentError('Could not solve DP problem. Please retry'
                                  ' with a SimpleCharger')

        return dp_solution, solutionTime, num_items, (pos_total + total_charge_diff)


class SymmetricCDPSolver(Solver):
    """An optimizing solver using Dynamic Programming, C version.

    Use the HistogramCollector to produce appropriate charge \
    distributions.
    """
    def __init__(self, rounding_digits) -> None:
        """Create a CDPSolver.

        Args:
            rounding_digits: How many significant digits to round the \
                    resulting charges to.
        """
        self.__rounding_digits = rounding_digits

    def solve_partial_charges(
            self,
            graph: nx.Graph,
            charge_dists_collector: Dict[Atom, Tuple[ChargeList, WeightList]],
            total_charge: int,
            keydict: Dict[Atom, str] = None,
            total_charge_diff: float=DEFAULT_TOTAL_CHARGE_DIFF,
            **kwargs
            ) -> None:
        """Assign charges to the atoms in a graph.

        Modify a graph by adding additional attributes describing the \
        atoms' charges and scores. In particular, each atom will get \
        a 'partial_charge' attribute with the partial charge, and a \
        'score' attribute giving a degree of certainty for that charge.

        This solver uses Dynamic Programming to solve the \
        epsilon-Multiple Choice Knapsack Problem. This is the Python \
        version of the algorithm, see DPSolver for the Python \
        implementation.

        Args:
            graph: The molecule graph to solve charges for.
            charge_dists: Charge distributions for the atoms, obtained \
                    by a Collector.
            total_charge: The total charge of the molecule.
            total_charge_diff: Maximum allowed deviation from the total charge
        """

        atom_idx = dict()
        for k, (atom, (_, _)) in enumerate(charge_dists_collector.items()):
            atom_idx[k] = atom

        neighborhoodclasses = self.compute_atom_neighborhood_classes(atom_idx, keydict)

        # reduce charge distributions (charges and frequencies of atoms within one neighborhoodclass get combined)
        charge_dists_reduced = self.reduce_charge_distributions(charge_dists_collector, atom_idx, neighborhoodclasses)

        solution, solutionTime, num_items, scaled_capacity = self.solve_dp_c(charge_dists_reduced, total_charge, total_charge_diff)

        charge = 0
        profit = 0
        for (i,j) in solution:
            for k in neighborhoodclasses[i]:
                graph.node[atom_idx[k]]['partial_charge'] = charge_dists_collector[atom_idx[k]][0][j]
                graph.node[atom_idx[k]]['score'] = charge_dists_collector[atom_idx[k]][1][j]
                charge += graph.node[atom_idx[k]]['partial_charge']
                profit += graph.node[atom_idx[k]]['score']

        graph.graph['total_charge'] = round(charge, self.__rounding_digits)
        graph.graph['score'] = profit
        graph.graph['time'] = solutionTime
        graph.graph['items'] = num_items
        graph.graph['scaled_capacity'] = scaled_capacity

    def reduce_charge_distributions(self, charge_dists_collector, atom_idx, neighborhoodclasses):
        charge_dists = dict()
        for neighborhoodclass in neighborhoodclasses:
            i = neighborhoodclass[0]
            k = len(neighborhoodclass)
            (charges, frequencies) = charge_dists_collector[atom_idx[i]]
            charge_dists[atom_idx[i]] = ([k * x for x in charges], [k * x for x in frequencies])

        return charge_dists


    def solve_dp_c(self, charge_dists, total_charge, total_charge_diff):

        import charge.c.dp as dp

        num_sets = len(charge_dists)
        num_items = sum(len(charges) for (_, (charges, _)) in charge_dists.items())

        weights = dp.new_doublea(num_items)
        profits = dp.new_doublea(num_items)
        sets = dp.new_ushorta(num_sets)
        solution = dp.new_ushorta(num_sets)

        offset = 0
        pos_total = total_charge
        for k, (atom, (charges, frequencies)) in enumerate(charge_dists.items()):
            dp.ushorta_setitem(sets, k, len(charges))
            for i, (charge, frequency) in enumerate(zip(charges, frequencies)):
                dp.doublea_setitem(weights, offset + i, charge)
                dp.doublea_setitem(profits, offset + i, frequency)
            pos_total -= min(charges)
            offset += len(charges)

        solutionTime = -perf_counter()

        profit = dp.solve_dp(weights, profits,
                             sets, num_items, num_sets,
                             self.__rounding_digits, total_charge, total_charge_diff,
                             solution)

        solutionTime += perf_counter()

        dp_solution = list()
        if profit >= 0:
            for k in range(num_sets):
                i = dp.ushorta_getitem(solution, k)
                dp_solution.append((k, i))


        dp.delete_doublea(weights)
        dp.delete_doublea(profits)
        dp.delete_ushorta(sets)
        dp.delete_ushorta(solution)

        if profit < 0:
            raise AssignmentError('Could not solve DP problem. Please retry'
                                  ' with a SimpleCharger')

        return dp_solution, solutionTime, num_items, (pos_total + total_charge_diff)