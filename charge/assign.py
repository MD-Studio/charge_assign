import warnings

import networkx as nx
import numpy

from charge.bond_type import BondType
from charge.nauty import Nauty
from charge.repository import Repository


class AssignmentError(Warning):
    pass


class AtomicCharger:

    def __init__(self, repository: Repository=None, nauty: Nauty=None) -> None:
        self.__nauty = nauty or Nauty()
        self.__repo = repository or Repository(nauty=self.__nauty)

    def __set_partial_charges(self, graph: nx.Graph, iacm_only: bool, shell: int) -> nx.Graph:
        shells = sorted(self.__repo.charges_iacm.keys(), reverse=True) if shell < 0 else [shell]
        for atom in graph.nodes():
            for shell in shells:
                key = self.__nauty.canonize_neighborhood(graph, atom, shell,
                                                         color_key='iacm' if 'iacm' in graph.node[atom] else 'atom_type')
                if key in self.__repo.charges_iacm[shell]:
                    values = self.__repo.charges_iacm[shell][key]
                    graph.node[atom]['partial_charge'] = numpy.mean(values)
                    graph.node[atom]['partial_charge_std'] = numpy.std(values)
                    break
                elif not iacm_only:
                    key = self.__nauty.canonize_neighborhood(graph, atom, shell)
                    if key in self.__repo.charges_elem[shell]:
                        values = self.__repo.charges_elem[shell][key]
                        graph.node[atom]['partial_charge'] = numpy.mean(values)
                        graph.node[atom]['partial_charge_std'] = numpy.std(values)
                        break
            else:
                warnings.warn(AssignmentError('Could not assign charge to atom {0}'.format(atom)))
                graph.node[atom]['partial_charge'] = float('nan')
                graph.node[atom]['partial_charge_std'] = float('nan')

        return graph

    def __iacmize(self, graph: nx.Graph) -> nx.Graph:
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

    def charge(self, graph: nx.Graph, iacmize:bool=False, iacm_only:bool=False, shell:int=None) -> nx.Graph:
        if iacmize:
            graph = self.__iacmize(graph)
        return self.__set_partial_charges(graph, iacm_only=iacm_only or iacmize,
                                          shell=shell if shell and shell >= 0 else -1)
