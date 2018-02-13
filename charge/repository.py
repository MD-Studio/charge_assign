import bisect
import os
from collections import defaultdict
from itertools import groupby
from multiprocessing.pool import Pool
from typing import Callable, List

import math
from zipfile import ZipFile

import itertools
import msgpack
import networkx as nx

from charge.babel import convert_from, IOType
from charge.nauty import Nauty
from charge.settings import REPO_LOCATION, IACM_MAP


class LoadError(Exception):
    pass


class Repository:

    def __init__(self,
                 location: str= REPO_LOCATION,
                 data_location: str=None,
                 nauty: Nauty=None,
                 min_shell: int=1,
                 max_shell: int=7,
                 processes: int=None) -> None:

        self.__nauty = nauty or Nauty()
        self.__min_shell = max(min_shell, 0)
        self.__max_shell = max_shell

        self.charges_iacm = defaultdict(lambda: defaultdict(list))
        self.charges_elem = defaultdict(lambda: defaultdict(list))

        if data_location:
            if processes:
                processes = max(processes, 1)
            self.__create(data_location, processes)
            self.__create(data_location, processes, iacm_to_elements=True)

        else:
            with ZipFile(location, mode='r') as zf:
                self.__min_shell, self.__max_shell = msgpack.unpackb(zf.read('meta'), encoding='utf-8')
                self.charges_iacm = msgpack.unpackb(zf.read('charges_iacm'), encoding='utf-8')
                self.charges_elem = msgpack.unpackb(zf.read('charges_elem'), encoding='utf-8')
                self.__iso_iacm = msgpack.unpackb(zf.read('iso_iacm'), encoding='utf-8')
                self.__iso_elem = msgpack.unpackb(zf.read('iso_elem'), encoding='utf-8')

    def __create(self, data_location: str, processes:int, iacm_to_elements: bool=False) -> None:
        # TODO add support for ITF files
        molids = [int(fn.replace('.lgf', ''))
                  for fn in os.listdir(data_location) if fn.endswith('.lgf')]

        if iacm_to_elements:
            print('IACM to element mode...')
        print('reading files', end='', flush=True)

        chunksize = int(math.ceil(len(molids) / 100))
        chunks = map(lambda i:
                     (molids[i:i+chunksize], data_location, iacm_to_elements, self.__nauty),
                     range(0, len(molids), chunksize))
        with Pool(processes) as pool:
            graphs = pool.starmap(read, chunks)
        print(flush=True)
        graphs, canons = zip(*itertools.chain.from_iterable(graphs))
        canons = dict(zip(molids, canons))

        for shell in range(self.__min_shell, self.__max_shell + 1):
            print('shell %d' % shell, end='', flush=True)
            chunks = map(lambda i:
                         (graphs[i:i + chunksize], self.__nauty, shell, iacm_to_elements),
                         range(0, len(graphs), chunksize))
            with Pool(processes) as pool:
                charges = pool.starmap(get_charges, chunks)
            print(flush=True)
            if not iacm_to_elements:
                for c in charges:
                    for key, values in c.items():
                        self.charges_iacm[shell][key] += values
            else:
                for c in charges:
                    for key, values in c.items():
                        self.charges_elem[shell][key] += values

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

    def __iterate(self, data_location: str, molid: int,
                  callable_iacm: Callable[[int, str, float], None],
                  callable_elem: Callable[[int, str, float], None]):
        ids_iacm = {molid}
        ids_elem = {molid}
        if molid in self.__iso_iacm:
            ids_iacm.union(set(self.__iso_iacm[molid]))
        if molid in self.__iso_elem:
            ids_iacm.union(set(self.__iso_elem[molid]))

        for molid in ids_iacm:
            with open(os.path.join(data_location, '%d.lgf' % molid), 'r') as f:
                graph = convert_from(f.read(), IOType.LGF)
                for shell in range(1, self.__max_shell + 1):
                    for key, partial_charge in iter_atomic_fragments(graph, self.__nauty, shell):
                        callable_iacm(shell, key, partial_charge)

        for molid in ids_elem:
            with open(os.path.join(data_location, '%d.lgf' % molid), 'r') as f:
                graph = convert_from(f.read(), IOType.LGF)
                for v, data in graph.nodes(data=True):
                    graph.node[v]['atom_type'] = IACM_MAP[data['atom_type']]
                for shell in range(1, self.__max_shell + 1):
                    for key, partial_charge in iter_atomic_fragments(graph, self.__nauty, shell):
                        callable_elem(shell, key, partial_charge)

    def add(self, data_location: str, molid: int):
        def a(shell, key, partial_charge, repo):
            if not shell in repo:
                repo[shell] = dict()
            if not key in repo[shell]:
                repo[shell][key] = []
            bisect.insort_left(repo[shell][key], partial_charge)

        self.__iterate(data_location, molid,
            lambda shell, key, partial_charge: a(shell, key, partial_charge, self.charges_iacm),
            lambda shell, key, partial_charge: a(shell, key, partial_charge, self.charges_elem)
        )

    def subtract(self,  data_location: str, molid: int):
        def s(shell, key, partial_charge, repo):
            repo[shell][key].pop(bisect.bisect_left(repo[shell][key], partial_charge))
            if len(repo[shell][key]) == 0:
                del repo[shell][key]
            if len(repo[shell]) == 0:
                del repo[shell]

        self.__iterate(data_location, molid,
            lambda shell, key, partial_charge: s(shell, key, partial_charge, self.charges_iacm),
            lambda shell, key, partial_charge: s(shell, key, partial_charge, self.charges_elem))

    def write(self, out: str):
        with ZipFile(out, mode='w') as zf:
            zf.writestr('meta', msgpack.packb((self.__min_shell, self.__max_shell)))
            zf.writestr('charges_iacm', msgpack.packb(self.charges_iacm))
            zf.writestr('charges_elem', msgpack.packb(self.charges_elem))
            zf.writestr('iso_iacm', msgpack.packb(self.__iso_iacm))
            zf.writestr('iso_elem', msgpack.packb(self.__iso_elem))


def read(molids: List[int], data_location: str, iacm_to_elements: bool, nauty: Nauty):
    res = []
    for molid in molids:
        with open(os.path.join(data_location, '%d.lgf' % molid), 'r') as f:
            graph = convert_from(f.read(), IOType.LGF)
            if iacm_to_elements:
                for v, data in graph.nodes(data=True):
                    graph.node[v]['atom_type'] = IACM_MAP[data['atom_type']]
            canon = nauty.canonize(graph, with_core=False)
            res.append((graph, canon))
    print('.', end='', flush=True)
    return res

def get_charges(graphs: List[nx.Graph], nauty: Nauty, shell: int, iacm_to_elements: bool):
    charges = defaultdict(list)
    for graph in graphs:
        if not iacm_to_elements:
            for key, partial_charge in iter_atomic_fragments(graph, nauty, shell):
                charges[key].append(partial_charge)
        else:
            for key, partial_charge in iter_atomic_fragments(graph, nauty, shell):
                charges[key].append(partial_charge)
    print('.', end='', flush=True)
    return charges

def iter_atomic_fragments(graph: nx.Graph, nauty: Nauty, shell: int):
    for atom in graph.nodes():
        partial_charge = float(graph.node[atom]['partial_charge'])
        yield nauty.canonize_neighborhood(graph, atom, shell), partial_charge
