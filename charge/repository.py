import bisect
import os
from collections import defaultdict
from itertools import groupby
from multiprocessing import Value, Process, JoinableQueue
from queue import Queue, Empty
from typing import Callable, List
from zipfile import ZipFile

import math
import msgpack
import networkx as nx
import time

from charge.babel import convert_from, IOType
from charge.nauty import Nauty
from charge.settings import REPO_LOCATION, IACM_MAP, NAUTY_EXC
from charge.multiprocessor import MultiProcessor
from charge.util import print_progress


class Repository:

    def __init__(self,
                 location: str= REPO_LOCATION,
                 data_location: str=None,
                 data_type: IOType=IOType.LGF,
                 min_shell: int=1,
                 max_shell: int=7,
                 processes: int=None) -> None:

        self.__nauty = Nauty()
        self.__min_shell = max(min_shell, 0)
        self.__max_shell = max_shell

        self.charges_iacm = defaultdict(lambda: defaultdict(list))
        self.charges_elem = defaultdict(lambda: defaultdict(list))

        if data_type == IOType.LGF:
            self.__ext = '.lgf'
        elif data_type == IOType.GML:
            self.__ext = '.gml'
        elif data_type == IOType.ITP:
            self.__ext = '.itp'
        else:
            raise ValueError('Unsupported file type: {}'.format(data_type.name))
        self.__data_type = data_type

        if data_location:
            cpus = os.cpu_count()
            if processes:
                processes = min(max(processes-1, 1), cpus-1)
            else:
                processes = cpus-1
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
        molids = [int(fn.replace(self.__ext, ''))
                  for fn in os.listdir(data_location) if fn.endswith(self.__ext)]

        if iacm_to_elements:
            print('IACM to element mode...')

        canons = dict()
        graphs = []

        with MultiProcessor(ReadWorker, (data_location, iacm_to_elements, self.__ext, self.__data_type)) as mp:
            for molid, graph, canon in mp.processed(molids):
                graphs.append(graph)
                canons[molid] = canon

                if len(graphs) % 20 == 0 or len(graphs) == len(molids):
                    print_progress(len(graphs), len(molids), 'reading files:')

        for shell in range(self.__min_shell, self.__max_shell + 1):
            progress = 0
            with MultiProcessor(ChargeWorker, shell) as mp:
                for c in mp.processed(graphs):
                    for key, values in c.items():
                        if not iacm_to_elements:
                            self.charges_iacm[shell][key] += values
                        else:
                            self.charges_elem[shell][key] += values
                    progress += 1
                    if progress % 20 == 0 or progress == len(graphs):
                        print_progress(progress, len(graphs), prefix='shell %d:' % shell)


        if not iacm_to_elements:
            self.__iso_iacm = defaultdict(list)
        else:
            self.__iso_elem = defaultdict(list)
        for _, group in groupby(molids, key=lambda molid: canons[molid]):
            isomorphics = list(group)
            if len(isomorphics) > 1:
                for molid in isomorphics:
                    if not iacm_to_elements:
                        self.__iso_iacm[molid] = isomorphics
                    else:
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
            with open(os.path.join(data_location, '%d%s' % (molid, self.__ext)), 'r') as f:
                graph = convert_from(f.read(), self.__data_type)
                for shell in range(1, self.__max_shell + 1):
                    for key, partial_charge in iter_atomic_fragments(graph, self.__nauty, shell):
                        callable_iacm(shell, key, partial_charge)

        for molid in ids_elem:
            with open(os.path.join(data_location, '%d%s' % (molid, self.__ext)), 'r') as f:
                graph = convert_from(f.read(), self.__data_type)
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


class ReadWorker:
    def __init__(self, data_location: str, iacm_to_elements: bool, ext: str, data_type: IOType):
        self.__data_location = data_location
        self.__iacm_to_elements = iacm_to_elements
        self.__ext = ext
        self.__data_type = data_type

        self.__nauty = Nauty()


    def process(self, molid: int) -> None:
        with open(os.path.join(self.__data_location, '%d%s' % (molid, self.__ext)), 'r') as f:
            graph = convert_from(f.read(), self.__data_type)
            if self.__iacm_to_elements:
                for v, data in graph.nodes(data=True):
                    graph.node[v]['atom_type'] = IACM_MAP[data['atom_type']]
            canon = self.__nauty.canonize(graph)
            return molid, graph, canon


class ChargeWorker:
    def __init__(self, shell: int):
        self.__shell = shell
        self.__nauty = Nauty()


    def process(self, graph: nx.Graph) -> None:
        charges = defaultdict(list)

        for key, partial_charge in iter_atomic_fragments(graph, self.__nauty, self.__shell):
            charges[key].append(partial_charge)

        return charges


def iter_atomic_fragments(graph: nx.Graph, nauty: Nauty, shell: int):
    for atom in graph.nodes():
        if not 'partial_charge' in graph.node[atom]:
            raise KeyError('Missing property "partial_charge" for atom {}'.format(atom))
        partial_charge = float(graph.node[atom]['partial_charge'])
        yield nauty.canonize_neighborhood(graph, atom, shell), partial_charge
