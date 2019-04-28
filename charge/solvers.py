import itertools
from abc import ABC, abstractmethod
from time import perf_counter
from typing import Dict, Tuple, Optional, List

import networkx as nx
from pulp import LpVariable, LpInteger, LpMaximize, LpProblem, LpStatusOptimal, CPLEX_CMD, GUROBI_CMD, PULP_CBC_CMD, \
    GLPK_CMD, COIN_CMD

from charge.charge_types import Atom, ChargeList, WeightList
from charge.settings import DEFAULT_TOTAL_CHARGE_DIFF, ROUNDING_DIGITS, ILP_SOLVER_MAX_SECONDS, DEFAULT_SHELL_SIZE
from charge.util import AssignmentError
from charge.nauty import Nauty


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
        """
        pass

    def compute_neighborhood_canonical_key_dict(
            self,
            graph: nx.Graph,
            shells: List[int],
            nauty: Nauty) -> dict():

        """Calculates a nauty canonical key for the neighborhood of each atom.

        Args:
            graph: The molecule graph to solve charges for.
            shells: Shell sizes to use. uses first shellsize
            nauty: A nauty instance

        """

        keydict = dict()

        for atom in graph.nodes():
            shellsize = shells[0]
            atom_has_iacm = 'iacm' in graph.node[atom]

            if atom_has_iacm:
                keydict[atom] = nauty.canonize_neighborhood(graph, atom, shellsize, 'iacm')
            else:
                keydict[atom] = nauty.canonize_neighborhood(graph, atom, shellsize, 'atom_type')

        return keydict

    def compute_atom_neighborhood_classes(self, atom_idx : dict, keydict : dict):
        L = list()
        atoms = list(atom_idx)

        while(len(atoms) > 0):
            i = atoms[0]
            l = list([i])
            for j in atoms[1::]:
                if (keydict[atom_idx[i]] == keydict[atom_idx[j]]):  # if atoms i and j have identical neighborhoods
                    l.append(j)
            L.append(l)
            for j in l:
                atoms.remove(j)
        return L


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
            self.__solver = GUROBI_CMD(options={'timeLimit':max_seconds})
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
                 max_seconds: int=ILP_SOLVER_MAX_SECONDS,
                 nauty: Optional[Nauty]=None
                 ) -> None:
        """Create an ILPSolver.

        Args:
            rounding_digits: Number of digits to round the charges to.
            max_seconds: Maximum run-time to spend searching for a \
                    solution
        """
        self.__rounding_digits = rounding_digits
        self._nauty = nauty if nauty is not None else Nauty()

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
            raise RuntimeError('No solver found, there is something'
                    ' wrong with your pulp library setup.')

    def solve_partial_charges(
            self,
            graph: nx.Graph,
            charge_dists: Dict[Atom, Tuple[ChargeList, WeightList]],
            total_charge: int,
            total_charge_diff: float=DEFAULT_TOTAL_CHARGE_DIFF,
            shells: List[int]=DEFAULT_SHELL_SIZE,
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
        """

        atom_idx = dict()
        idx = list()
        # weights = partial charges
        weights = dict()
        # profits = frequencies
        profits = dict()

        # nauty hashes of k-neighborhoods for each atom
        keydict = dict()

        pos_total = total_charge
        for k, (atom, (charges, frequencies)) in enumerate(charge_dists.items()):
            atom_idx[k] = atom
            idx.append(list(zip(itertools.repeat(k), range(len(charges)))))
            weights[k] = charges
            profits[k] = frequencies
            pos_total -= min(charges)

        # compute k-neighborhood-canonical key for each atom of the molecule
        keydict = self.compute_neighborhood_canonical_key_dict(graph, shells, self._nauty)

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
        for neighborhood_class in self.compute_atom_neighborhood_classes(atom_idx, keydict):
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
        for i, j in enumerate(solution):
            graph.node[atom_idx[i]]['partial_charge'] = round((deflate * items[i][j][1]) + w_min[i],
                                                              self.__rounding_digits)
            graph.node[atom_idx[i]]['score'] = items[i][j][2]
            charge += graph.node[atom_idx[i]]['partial_charge']

        graph.graph['total_charge'] = round(charge, self.__rounding_digits)
        graph.graph['score'] = max_val
        graph.graph['time'] = solutionTime
        graph.graph['items'] = sum(len(i) for i in items)
        graph.graph['scaled_capacity'] = pos_total_charge + total_charge_diff

class SymmetricDPSolver(Solver):
    """An optimizing solver using Dynamic Programming.

    Use the HistogramCollector to produce appropriate charge \
    distributions.
    """
    def __init__(self,
                 rounding_digits,
                 nauty: Optional[Nauty]=None) -> None:
        """Create a DPSolver.

        Args:
            rounding_digits: How many significant digits to round the \
                    resulting charges to.
        """
        self.__rounding_digits = rounding_digits
        self._nauty = nauty if nauty is not None else Nauty()


    def solve_partial_charges(
            self,
            graph: nx.Graph,
            charge_dists_collector: Dict[Atom, Tuple[ChargeList, WeightList]],
            total_charge: int,
            total_charge_diff: float=DEFAULT_TOTAL_CHARGE_DIFF,
            shells: List[int] = DEFAULT_SHELL_SIZE,
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
        """

        blowup = 10 ** self.__rounding_digits

        # reduce charge distributions (charges and frequencies of atoms within one neighborhoodclass get combined)
        charge_dists_reduced, atom_idx, neighborhoodclasses = self.reduce_charge_distributions(graph, charge_dists_collector, shells)

        items, pos_total_charge, max_sum = self.transform_weights(charge_dists_reduced, total_charge, blowup)

        solution, max_val, solutionTime = self.solve_dp(items, total_charge_diff, pos_total_charge, max_sum, blowup)

        charge = 0
        for i, j in enumerate(solution):
            for k in neighborhoodclasses[i]:
                graph.node[atom_idx[k]]['partial_charge'] = charge_dists_collector[atom_idx[k]][0][j]
                graph.node[atom_idx[k]]['score'] = charge_dists_collector[atom_idx[k]][1][j]
                charge += graph.node[atom_idx[k]]['partial_charge']

        graph.graph['total_charge'] = round(charge, self.__rounding_digits)
        graph.graph['score'] = max_val
        graph.graph['time'] = solutionTime
        graph.graph['items'] = sum(len(i) for i in items)
        graph.graph['scaled_capacity'] = pos_total_charge + total_charge_diff

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

        return solution, max_val, solutionTime

    def reduce_charge_distributions(self, graph, charge_dists_collector, shells):
        atom_idx = dict()
        for k, (atom, (charges, frequencies)) in enumerate(charge_dists_collector.items()):
            atom_idx[k] = atom
        keydict = self.compute_neighborhood_canonical_key_dict(graph, shells, self._nauty)
        neighborhoodclasses = self.compute_atom_neighborhood_classes(atom_idx, keydict)

        charge_dists = dict()
        for neighborhoodclass in neighborhoodclasses:
            i = neighborhoodclass[0]
            charge_dists[atom_idx[i]] = charge_dists_collector[atom_idx[i]]
            for j in neighborhoodclass[1::]:
                #Add the charges and frequencies of Atom j to those of atom i
                charges = [round(sum(x),self.__rounding_digits) for x in zip(charge_dists[atom_idx[i]][0], charge_dists_collector[atom_idx[j]][0])]
                frequencies = [round(sum(x),self.__rounding_digits) for x in zip(charge_dists[atom_idx[i]][1], charge_dists_collector[atom_idx[j]][1])]

                charge_dists[atom_idx[i]] = (charges,frequencies)

        return charge_dists, atom_idx, neighborhoodclasses


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
        """

        import charge.c.dp as dp

        num_sets = len(charge_dists)
        num_items = sum(len(charges) for (_, (charges, _)) in charge_dists.items())

        atom_idx = dict()

        weights = dp.new_doublea(num_items)
        profits = dp.new_doublea(num_items)
        sets = dp.new_ushorta(num_sets)
        solution = dp.new_ushorta(num_sets)

        offset = 0
        pos_total = total_charge
        for k, (atom, (charges, frequencies)) in enumerate(charge_dists.items()):
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
                             self.__rounding_digits, total_charge, total_charge_diff,
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

            graph.graph['total_charge'] = round(charge, self.__rounding_digits)
            graph.graph['score'] = profit
            graph.graph['time'] = solutionTime
            graph.graph['items'] = num_items
            graph.graph['scaled_capacity'] = pos_total + total_charge_diff

        dp.delete_doublea(weights)
        dp.delete_doublea(profits)
        dp.delete_ushorta(sets)
        dp.delete_ushorta(solution)

        if profit < 0:
            raise AssignmentError('Could not solve DP problem. Please retry'
                    ' with a SimpleCharger')
