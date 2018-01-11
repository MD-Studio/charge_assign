import warnings
from collections import defaultdict, deque
from typing import Any

import networkx as nx

from charge.babel import BondType
from charge.nauty import Nauty
from charge.repository import Repository


class AssignmentError(Warning):
    pass


class AtomicCharger:

    def __init__(self, repository: Repository=None, nauty: Nauty=None) -> None:
        if not nauty:
            self.__nauty = Nauty()
        if not repository:
            self.__repo = Repository(nauty=self.__nauty)

    def __get_canon_key(self, graph: nx.Graph, atom: Any, shell: int) -> str:

        fragment = nx.Graph()

        visited = set()
        depth = defaultdict(int)
        queue = deque()

        queue.append(atom)
        visited.add(atom)

        while len(queue) > 0:
            current_atom = queue.popleft()
            fragment.add_node(current_atom, attr_dict=graph.node[current_atom].copy())

            if depth[current_atom] < shell:
                for other_atom in graph.neighbors(current_atom):
                    if not other_atom in visited:
                        visited.add(other_atom)
                        depth[other_atom] = depth[current_atom] + 1
                        queue.append(other_atom)

        nodes = fragment.nodes()
        for i, u in enumerate(nodes):
            for v in nodes[i+1:]:
                if graph.has_edge(u, v):
                    fragment.add_edge(u, v)

        return self.__nauty.canonize(fragment)

    def __set_partial_charges(self, graph: nx.Graph, iacm_only: bool, shell: int) -> nx.Graph:
        shells = sorted(self.__repo.charges_iacm.keys(), reverse=True) if shell < 1 else [shell]
        for atom, data in graph.nodes_iter(data=True):
            # TODO optional: fixed shell size
            for sh in shells:
                key = self.__get_canon_key(graph, atom, sh)
                if key in self.__repo.charges_iacm[sh]:
                    # TODO get charge and uncertainty from the histogram
                    charge, uncertainty = self.__repo.charges_iacm[sh][key]
                    graph.node[atom]['partial_charge'] = charge
                    graph.node[atom]['uncertainty'] = uncertainty
                    break
                elif not iacm_only:
                    key = self.__get_canon_key(graph, atom, sh)
                    if key in self.__repo.charges_elem[sh]:
                        # TODO get charge and uncertainty from the histogram
                        charge, uncertainty = self.__repo.charges_elem[sh][key]
                        graph.node[atom]['partial_charge'] = charge
                        graph.node[atom]['uncertainty'] = uncertainty
                        break
            else:
                warnings.warn(AssignmentError('Could not assign charge_assign to atom {0}'.format(atom)))

        return graph

    def __iacmize(self, graph: nx.Graph) -> nx.Graph:
        def aromatic_neighbors(u) -> list:
            return list(filter(lambda v: 'bond_type' in graph[u][v] and graph[u][v]['bond_type'] == BondType.AROMATIC,
                               graph.neighbors(u)))

        for atom in graph.nodes_iter():
            element = graph.node[atom]['atom_type']

            bas = graph.neighbors(atom)
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
                                len(graph.neighbors(a)) == 1, graph.neighbors(bas[0])))) > 1 and \
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
                elif len(list(filter(lambda a: a.GetSymbol() == 'H', bas))) < 2:
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

    def charge(self, graph: nx.Graph, iacmize:bool=True, iacm_only:bool=False, shell:int=None) -> nx.Graph:
        if iacmize:
            graph = self.__iacmize(graph)
        if not iacmize:
            iacm_only=True
        return self.__set_partial_charges(graph, iacm_only=iacm_only, shell=shell if shell and shell > 0 else -1)
