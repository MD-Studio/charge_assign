import bisect
import os
from collections import defaultdict
from itertools import groupby
from typing import Callable

import networkx as nx

from charge.babel import convert_from, IOType
from charge.nauty import Nauty
from charge.settings import REPO_LOCATION, IACM_MAP


class LoadError(Exception):
    pass


class Repository:

    def __init__(self, location: str= REPO_LOCATION, data_location: str=None, nauty: Nauty=None,
                 max_shell: int=7) -> None:

        self.__nauty = nauty or Nauty()
        self.__max_shell = max_shell

        self.charges_iacm = defaultdict(lambda: defaultdict(list))
        self.charges_elem = defaultdict(lambda: defaultdict(list))

        if data_location:
            self.__create(data_location)
            self.__create(data_location, iacm_to_elements=True)

        # TODO else read from input file (location)

    def __create(self, data_location: str, iacm_to_elements: bool=False) -> None:
        # TODO add support for ITF files
        molids = [int(fn.replace('.lgf', ''))
                  for fn in os.listdir(data_location) if fn.endswith('.lgf')]

        canons = dict()
        for molid in molids:
            with open(os.path.join(data_location, '%d.lgf' % molid), 'r') as f:
                print(molid)
                graph = convert_from(f.read(), IOType.LGF)
                if iacm_to_elements:
                    for v, data in graph.nodes(data=True):
                        graph.node[v]['atom_type'] = IACM_MAP[data['atom_type']]
                canons[molid] = self.__nauty.canonize(graph, with_core=False)

                for shell in range(1, self.__max_shell + 1):
                    print(molid, shell)
                    if not iacm_to_elements:
                        for key, partial_charge in self.__iter_atomic_fragments(graph, shell):
                            self.charges_iacm[shell][key].append(partial_charge)
                    else:
                        for key, partial_charge in self.__iter_atomic_fragments(graph, shell):
                            self.charges_elem[shell][key].append(partial_charge)

        if not iacm_to_elements:
            self.__iso_iacm = defaultdict(list)
        else:
            self.__iso_elem = defaultdict(list)
        for _, group in groupby(molids, key=lambda molid: canons[molid]):
            isomorphics = list(group)
            if len(isomorphics) > 1:
                for molid in isomorphics:
                    self.__iso_iacm[molid] = isomorphics
        for _, group in groupby(molids, key=lambda molid: canons[molid]):
            isomorphics = list(group)
            if len(isomorphics) > 1:
                for molid in isomorphics:
                    self.__iso_elem[molid] = isomorphics

        for shell in range(1, self.__max_shell + 1):
            if not iacm_to_elements:
                for key, values in self.charges_iacm[shell].items():
                    self.charges_iacm[shell][key] = sorted(values)
            else:
                for key, values in self.charges_elem[shell].items():
                    self.charges_elem[shell][key] = sorted(values)

    def __iter_atomic_fragments(self, graph: nx.Graph, shell: int):
        for atom in graph.nodes():
            partial_charge = float(graph.node[atom]['partial_charge'])
            yield self.__nauty.canonize_neighborhood(graph, atom, shell), partial_charge

    def __iterate(self, data_location: str, molid: int,
                  callable_iacm: Callable[[int, str, float], None],
                  callable_elem: Callable[[int, str, float], None]):
        ids_iacm = set(molid)
        ids_elem = set(molid)
        if molid in self.__iso_iacm:
            ids_iacm.union(set(self.__iso_iacm[molid]))
        if molid in self.__iso_elem:
            ids_iacm.union(set(self.__iso_elem[molid]))

        for molid in ids_iacm:
            with open(os.path.join(data_location, '%d.lgf' % molid), 'r') as f:
                graph = convert_from(f.read(), IOType.LGF)
                for shell in range(1, self.__max_shell + 1):
                    for key, partial_charge in self.__iter_atomic_fragments(graph, shell):
                        callable_iacm(key, partial_charge)

        for molid in ids_elem:
            with open(os.path.join(data_location, '%d.lgf' % molid), 'r') as f:
                graph = convert_from(f.read(), IOType.LGF)
                for v, data in graph.nodes(data=True):
                    graph.node[v]['atom_type'] = IACM_MAP[data['atom_type']]
                for shell in range(1, self.__max_shell + 1):
                    for key, partial_charge in self.__iter_atomic_fragments(graph, shell):
                        callable_elem(key, partial_charge)

    def add(self, data_location: str, molid: int):
        self.__iterate(data_location, molid,
            lambda shell, key, partial_charge: bisect.insort_left(self.charges_iacm[shell][key], partial_charge),
            lambda shell, key, partial_charge: bisect.insort_left(self.charges_iacm[shell][key], partial_charge)
        )

    def subtract(self,  data_location: str, molid: int):
        self.__iterate(data_location, molid,
            lambda shell, key, partial_charge: \
                self.charges_iacm[shell][key].pop(bisect.bisect_left(self.charges_iacm[shell][key], partial_charge)),
            lambda shell, key, partial_charge: \
               self.charges_elem[shell][key].pop(bisect.bisect_left(self.charges_elem[shell][key], partial_charge))
                       )

    def write(self, out: str):
        # TODO write to output file
        pass
