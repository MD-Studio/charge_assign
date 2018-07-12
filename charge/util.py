from collections import deque
from time import perf_counter
from typing import Any

import networkx as nx

from charge.babel import BondType

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
        return list(filter(lambda v: 'bond_type' in graph[u][v] and graph[u][v]['bond_type'] == BondType.AROMATIC,
                           graph.neighbors(u)))

    for atom in graph.nodes():
        element = graph.node[atom]['atom_type']

        bas = list(graph.neighbors(atom))
        if element == 'C':
            bhs = list(filter(lambda a: graph.node[a]['atom_type'] == 'H', bas))
            if len(bas) == 4 and len(bhs) == 0:
                graph.node[atom]['iacm'] = 'CH0'
            else:
                graph.node[atom]['iacm'] = 'C'
        elif element == 'H':
            if bas and graph.node[bas[0]]['atom_type'] == 'C':
                graph.node[atom]['iacm'] = 'HC'
            else:
                graph.node[atom]['iacm'] = 'H'
        elif element == 'O':
            if len(list(filter(lambda a: graph.node[a]['atom_type'] == 'C', bas))) == len(bas) and len(bas) > 1:
                graph.node[atom]['iacm'] = 'OE'
            elif len(bas) > 1:
                graph.node[atom]['iacm'] = 'OA'
            elif bas and len(list(filter(lambda a: graph.node[a]['atom_type'] == 'O' and \
                            len(list(graph.neighbors(a))) == 1, graph.neighbors(bas[0])))) > 1 and \
                            bas != aromatic_neighbors(atom):
                graph.node[atom]['iacm'] = 'OM'
            else:
                graph.node[atom]['iacm'] = 'O'
        elif element == 'N':
            if len(bas) > 3:
                graph.node[atom]['iacm'] = 'NL'
            elif len(bas) == 1:
                graph.node[atom]['iacm'] = 'NR'
            elif len(aromatic_neighbors(atom)) > 1:
                graph.node[atom]['iacm'] = 'NR'
            elif len(list(filter(lambda a: graph.node[a]['atom_type'] == 'H', bas))) < 2:
                graph.node[atom]['iacm'] = 'N'
            else:
                graph.node[atom]['iacm'] = 'NT'
        elif element == 'S':
            if len(bas) > 2:
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
