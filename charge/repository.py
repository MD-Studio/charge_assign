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

        out_q = JoinableQueue()

        canons = dict()
        graphs = []
        pool = []

        chunksize = int(math.ceil(len(molids) / processes))
        chunks = map(lambda i: molids[i:i + chunksize], range(0, len(molids), chunksize))

        for chunk in chunks:
            process = Process(target=read_worker,
                              args=(chunk, data_location, iacm_to_elements, self.__ext, self.__data_type,
                                    out_q))
            pool.append(process)
            process.start()

        progress = 0
        for molid, graph, canon in iter_queue(pool, out_q):
            graphs.append(graph)
            canons[molid] = canon

            progress += 1
            print_progress(progress, len(molids), 'reading files:')

        for worker in pool:
            worker.join()
            if worker.exitcode is None:
                worker.terminate()

        for shell in range(self.__min_shell, self.__max_shell + 1):
            out_q = JoinableQueue()

            chunks = map(lambda i: graphs[i:i + chunksize], range(0, len(graphs), chunksize))

            for chunk in chunks:
                process = Process(target=charge_worker,
                                  args=(chunk, shell, iacm_to_elements, out_q))
                pool.append(process)
                process.start()

            progress = 0
            for c in iter_queue(pool, out_q):
                for key, values in c.items():
                    if not iacm_to_elements:
                        self.charges_iacm[shell][key] += values
                    else:
                        self.charges_elem[shell][key] += values
                progress += 1
                print_progress(progress, len(graphs), prefix='shell %d:' % shell)

            for worker in pool:
                worker.join()
                if worker.exitcode is None:
                    worker.terminate()

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


def iter_queue(pool: List[Process], queue: Queue):
    live_workers = list(pool)
    while live_workers:
        try:
            while True:
                yield queue.get(block=False)
        except Empty:
            pass

        time.sleep(0.1)
        if not queue.empty():
            continue
        live_workers = [p for p in live_workers if p.is_alive()]

    while not queue.empty():
        yield queue.get()


def read_worker(molids: List[int], data_location: str, iacm_to_elements: bool, ext: str, data_type: IOType,
                out_q:Queue):

    nauty = Nauty()

    for molid in molids:
        with open(os.path.join(data_location, '%d%s' % (molid, ext)), 'r') as f:
            graph = convert_from(f.read(), data_type)
            if iacm_to_elements:
                for v, data in graph.nodes(data=True):
                    graph.node[v]['atom_type'] = IACM_MAP[data['atom_type']]
            canon = nauty.canonize(graph)
            out_q.put((molid, graph, canon))


def charge_worker(graphs: List[nx.Graph], shell: int, iacm_to_elements: bool,
                  out_q:Queue):

    nauty = Nauty()

    for graph in graphs:
        charges = defaultdict(list)

        if not iacm_to_elements:
            for key, partial_charge in iter_atomic_fragments(graph, nauty, shell):
                charges[key].append(partial_charge)
        else:
            for key, partial_charge in iter_atomic_fragments(graph, nauty, shell):
                charges[key].append(partial_charge)
        out_q.put(charges)


def iter_atomic_fragments(graph: nx.Graph, nauty: Nauty, shell: int):
    for atom in graph.nodes():
        if not 'partial_charge' in graph.node[atom]:
            raise KeyError('Missing property "partial_charge" for atom {}'.format(atom))
        partial_charge = float(graph.node[atom]['partial_charge'])
        yield nauty.canonize_neighborhood(graph, atom, shell), partial_charge
