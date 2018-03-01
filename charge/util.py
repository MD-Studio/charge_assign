from collections import deque
from typing import Any

import networkx as nx

def print_progress(iteration:int,
                   total:int,
                   prefix:str='',
                   suffix:str='',
                   decimals:int=1,
                   length:int=80,
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
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '', flush=True)

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
