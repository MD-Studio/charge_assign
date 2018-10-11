from collections import deque
from math import ceil
from time import perf_counter
from typing import Any, List

import networkx as nx

from charge.babel import BondType


class AssignmentError(Warning):
    pass


_last_print = 0.0


def print_progress(iteration:int,
                   total:int,
                   prefix:str='',
                   suffix:str='',
                   decimals:int=1,
                   length:int=50,
                   fill:str='#') -> None:
    """Call in a loop to create terminal progress bar

        :param iteration: current iteration
        :type iteration: int
        :param total: total iterations
        :type total: int
        :param prefix: prefix string
        :type prefix: str
        :param suffix: suffix string
        :type suffix: str
        :param decimals: positive number of decimals in percent complete
        :type decimals: int
        :param length: character length of bar
        :type length: int
        :param fill : bar fill character
        :type fill: str
    """
    global _last_print
    cur_time = perf_counter()
    been_a_while = cur_time - _last_print > 0.2

    if iteration == 0 or iteration == total or been_a_while:
        _last_print = cur_time
        percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        print('\r%s |%s| %s%% %s\033[K' % (prefix, bar, percent, suffix), end = '', flush=True)

        if iteration == total:
            print()


def bfs_nodes(G: nx.Graph, source: Any, max_depth:int=0):
    """Iterate over nodes in a breadth-first search.

    The breadth-first search begins at `source` and enqueues the
    neighbors of newly visited nodes specified by the `neighbors`
    function.

    Adapted from `generic_bfs_edges <https://networkx.github.io/documentation/stable/_modules/networkx/algorithms/traversal/breadth_first_search.html>`_.

    :param G: NetworkX graph
    :type G: nx.Graph
    :param source: starting node for the breadth-first search
    :type source: Any
    :param max_depth: maximal depth of the breadth-first search
    :type max_depth: int
    """
    visited = {source}
    queue = deque([(G.neighbors(source), 1)])

    yield source

    while queue:
        children, depth = queue[0]

        if max_depth > 0 and depth > max_depth:
            break

        try:
            child = next(children)
            if child not in visited:
                yield child
                visited.add(child)
                queue.append((G.neighbors(child), depth+1))
        except StopIteration:
            queue.popleft()


def iacmize(graph: nx.Graph) -> nx.Graph:
    def aromatic_neighbors(u) -> list:
        # anything having an aromatic bond
        return list(filter(lambda v: 'bond_type' in graph[u][v] and graph[u][v]['bond_type'] == BondType.AROMATIC,
                           graph.neighbors(u)))

    # See https://pubs.acs.org/doi/full/10.1021/ct200196m
    for atom in graph.nodes():
        element = graph.node[atom]['atom_type']

        bas = list(graph.neighbors(atom))
        if element == 'C':
            bhs = list(filter(lambda a: graph.node[a]['atom_type'] == 'H', bas))
            if len(bas) == 4 and len(bhs) == 0:
                # C atom has four neighbours, none of which are H's
                graph.node[atom]['iacm'] = 'CH0'
            else:
                graph.node[atom]['iacm'] = 'C'
            # Other C IACM types are for united-atom topologies
        elif element == 'H':
            if bas and graph.node[bas[0]]['atom_type'] == 'C':
                # H atom has a C neighbour
                graph.node[atom]['iacm'] = 'HC'
            else:
                graph.node[atom]['iacm'] = 'H'
        elif element == 'O':
            if len(list(filter(lambda a: graph.node[a]['atom_type'] == 'C', bas))) == len(bas) and len(bas) > 1:
                # O between two C's
                graph.node[atom]['iacm'] = 'OE'
            elif len(bas) > 1:
                # O with two neighbours at least one of which is not carbon
                graph.node[atom]['iacm'] = 'OA'
            elif bas and len(list(filter(lambda a: graph.node[a]['atom_type'] == 'O' and \
                            len(list(graph.neighbors(a))) == 1, graph.neighbors(bas[0])))) > 1 and \
                            bas != aromatic_neighbors(atom):
                # O that has a neighbour which has a double-bonded neighbouring O
                graph.node[atom]['iacm'] = 'OM'
            else:
                graph.node[atom]['iacm'] = 'O'
        elif element == 'N':
            if len(bas) > 3:
                # N bound to four other atoms
                graph.node[atom]['iacm'] = 'NL'
            elif len(bas) == 1:
                # Single neighbor
                graph.node[atom]['iacm'] = 'NR'
            elif len(aromatic_neighbors(atom)) > 1:
                # Part of aromatic ring
                graph.node[atom]['iacm'] = 'NR'
            elif len(list(filter(lambda a: graph.node[a]['atom_type'] == 'H', bas))) < 2:
                # Nonaromatic, three neighbors, at most one hydrogen and a carbonyl
                # ! Where's the carbonyl?
                graph.node[atom]['iacm'] = 'N'
            else:
                graph.node[atom]['iacm'] = 'NT'
        elif element == 'S':
            if len(bas) > 2:
                # S with more than two neighbors, modeled as S in DMSO solvent
                graph.node[atom]['iacm'] = 'SDmso'
            else:
                graph.node[atom]['iacm'] = 'S'
        elif element == 'P':
            graph.node[atom]['iacm'] = 'P,SI'
        elif element == 'Si':
            graph.node[atom]['iacm'] = 'AR'
        elif element == 'F':
            graph.node[atom]['iacm'] = 'F'
        elif element == 'Cl':
            graph.node[atom]['iacm'] = 'CL'
        elif element == 'Br':
            graph.node[atom]['iacm'] = 'BR'
        else:
            graph.node[atom]['iacm'] = element

    return graph


def round_to(x, grain):
    # rounds to nearest multiple of grain, with 0.5 rounded down
    return ceil((x / grain) - 0.5) * grain


def median(values: List[float]):
    n = len(values)
    h = n // 2
    if n % 2 == 0:
        return sum(values[h - 1:h + 1]) / 2.0
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
        return 0.75 * values[q - 1] + 0.25 * values[q]
    if n % 4 == 3:
        return 0.75 * values[q] + 0.25 * values[q + 1]


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
