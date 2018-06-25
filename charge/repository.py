import bisect
import math
import os
import time
from collections import defaultdict
from itertools import groupby
from multiprocessing import Value, Process, JoinableQueue
from queue import Queue, Empty
from typing import Callable, List
from zipfile import ZipFile, ZIP_DEFLATED

import msgpack
import networkx as nx

from charge.babel import convert_from, IOType
from charge.nauty import Nauty
from charge.settings import REPO_LOCATION, IACM_MAP
from charge.util import print_progress


class Repository:

    def __init__(self,
                 location: str= REPO_LOCATION,
                 data_location: str=None,
                 data_type: IOType=IOType.LGF,
                 nauty: Nauty = None,
                 min_shell: int=1,
                 max_shell: int=7,
                 processes: int=None) -> None:

        self.__nauty = nauty or Nauty()
        self.__nauty_exe = self.__nauty.exe
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
                self.iso_iacm = msgpack.unpackb(zf.read('iso_iacm'), encoding='utf-8')
                self.iso_elem = msgpack.unpackb(zf.read('iso_elem'), encoding='utf-8')

    def __create(self, data_location: str, processes:int, iacm_to_elements: bool=False) -> None:
        molids = [fn.replace(self.__ext, '')
                  for fn in os.listdir(data_location) if fn.endswith(self.__ext)]

        if iacm_to_elements:
            print('IACM to element mode...')

        progress = Value('i', 0)
        total = Value('i', 100)
        out_q = JoinableQueue()

        progress.value = 0
        total.value = len(molids)

        canons = dict()
        graphs = []
        pool = []

        chunksize = int(math.ceil(len(molids) / processes))
        chunks = map(lambda i: molids[i:i + chunksize], range(0, len(molids), chunksize))

        for chunk in chunks:
            process = Process(target=read_worker,
                              args=(chunk, data_location, iacm_to_elements, self.__ext, self.__data_type,
                                    self.__nauty_exe, out_q, progress, total))
            pool.append(process)
            process.start()

        for molid, graph, canon in iter_queue(pool, out_q):
            graphs.append(graph)
            canons[molid] = canon

        if progress.value != total.value:
            print_progress(total.value, total.value, prefix='reading files:')

        for worker in pool:
            worker.join()
            if worker.exitcode is None:
                worker.terminate()

        for shell in range(self.__min_shell, self.__max_shell + 1):
            progress.value = 0
            total.value = len(graphs)
            out_q = JoinableQueue()

            chunks = map(lambda i: graphs[i:i + chunksize], range(0, len(graphs), chunksize))

            for chunk in chunks:
                process = Process(target=charge_worker,
                                  args=(chunk, shell, iacm_to_elements, self.__nauty_exe, out_q, progress, total))
                pool.append(process)
                process.start()

            for c in iter_queue(pool, out_q):
                for key, values in c.items():
                    if not iacm_to_elements:
                        self.charges_iacm[shell][key] += values
                    else:
                        self.charges_elem[shell][key] += values

            if progress.value != total.value:
                print_progress(total.value, total.value, prefix='shell %d:' % shell)

            for worker in pool:
                worker.join()
                if worker.exitcode is None:
                    worker.terminate()

        if not iacm_to_elements:
            self.iso_iacm = defaultdict(list)
        else:
            self.iso_elem = defaultdict(list)

        molids.sort(key=lambda molid: canons[molid])
        for _, group in groupby(molids, key=lambda molid: canons[molid]):
            isomorphics = list(group)
            if len(isomorphics) > 1:
                for molid in isomorphics:
                    if not iacm_to_elements:
                        self.iso_iacm[molid] = isomorphics
                    else:
                        self.iso_elem[molid] = isomorphics

        for shell in range(self.__min_shell, self.__max_shell + 1):
            if not iacm_to_elements:
                for key, values in self.charges_iacm[shell].items():
                    self.charges_iacm[shell][key] = sorted(values)
            else:
                for key, values in self.charges_elem[shell].items():
                    self.charges_elem[shell][key] = sorted(values)

    def __iterate(self, data_location: str, molid: int, with_iso: bool,
                  callable_iacm: Callable[[int, str, float], None],
                  callable_elem: Callable[[int, str, float], None]):
        ids_iacm = {molid}
        ids_elem = {molid}
        if with_iso:
            if molid in self.iso_iacm:
                ids_iacm.union(set(self.iso_iacm[molid]))
            if molid in self.iso_elem:
                ids_iacm.union(set(self.iso_elem[molid]))

        for molid in ids_iacm:
            with open(os.path.join(data_location, '%d%s' % (molid, self.__ext)), 'r') as f:
                graph = convert_from(f.read(), self.__data_type)
                for shell in range(self.__min_shell, self.__max_shell + 1):
                    for key, partial_charge in iter_atomic_fragments(graph, self.__nauty, shell):
                        callable_iacm(shell, key, partial_charge)

        for molid in ids_elem:
            with open(os.path.join(data_location, '%d%s' % (molid, self.__ext)), 'r') as f:
                graph = convert_from(f.read(), self.__data_type)
                for v, data in graph.nodes(data=True):
                    graph.node[v]['atom_type'] = IACM_MAP[data['atom_type']]
                for shell in range(self.__min_shell, self.__max_shell + 1):
                    for key, partial_charge in iter_atomic_fragments(graph, self.__nauty, shell):
                        callable_elem(shell, key, partial_charge)

    # TODO optional: add/subtract isomorphic molids
    def add(self, data_location: str, molid: int, with_iso: bool=False):
        def a(shell, key, partial_charge, repo):
            if not shell in repo:
                repo[shell] = dict()
            if not key in repo[shell]:
                repo[shell][key] = []
            bisect.insort_left(repo[shell][key], partial_charge)

        self.__iterate(data_location, molid, with_iso,
            lambda shell, key, partial_charge: a(shell, key, partial_charge, self.charges_iacm),
            lambda shell, key, partial_charge: a(shell, key, partial_charge, self.charges_elem)
        )

    def subtract(self,  data_location: str, molid: int, with_iso: bool=False):
        def s(shell, key, partial_charge, repo):
            if shell in repo and key in repo[shell]:
                repo[shell][key].pop(bisect.bisect_left(repo[shell][key], partial_charge))
                if len(repo[shell][key]) == 0:
                    del repo[shell][key]
                if len(repo[shell]) == 0:
                    del repo[shell]

        self.__iterate(data_location, molid, with_iso,
            lambda shell, key, partial_charge: s(shell, key, partial_charge, self.charges_iacm),
            lambda shell, key, partial_charge: s(shell, key, partial_charge, self.charges_elem))

    def write(self, out: str):
        with ZipFile(out, mode='w', compression=ZIP_DEFLATED) as zf:
            zf.writestr('meta', msgpack.packb((self.__min_shell, self.__max_shell)))
            zf.writestr('charges_iacm', msgpack.packb(self.charges_iacm))
            zf.writestr('charges_elem', msgpack.packb(self.charges_elem))
            zf.writestr('iso_iacm', msgpack.packb(self.iso_iacm))
            zf.writestr('iso_elem', msgpack.packb(self.iso_elem))


def iter_queue(pool: List[Process], queue: Queue, sleep: float=0.1):
    live_workers = list(pool)
    while live_workers:
        try:
            while True:
                yield queue.get(block=False)
        except Empty:
            pass

        time.sleep(sleep)
        if not queue.empty():
            continue
        live_workers = [p for p in live_workers if p.is_alive()]

    while not queue.empty():
        yield queue.get()


def read_worker(molids: List, data_location: str, iacm_to_elements: bool, ext: str, data_type: IOType,
                nauty_exe:str, out_q:Queue, progress: Value, total: Value):

    nauty = Nauty(nauty_exe)

    for molid in molids:
        with open(os.path.join(data_location, '{}{}'.format(molid, ext)), 'r') as f:
            graph = convert_from(f.read(), data_type)
            if iacm_to_elements:
                for v, data in graph.nodes(data=True):
                    graph.node[v]['atom_type'] = IACM_MAP[data['atom_type']]

            canon = nauty.canonize(graph)
            out_q.put((molid, graph, canon))
        progress.value += 1
        print_progress(progress.value, total.value, prefix='reading files:')


def charge_worker(graphs: List[nx.Graph], shell: int, iacm_to_elements: bool, nauty_exe:str,
                  out_q:Queue, progress: Value, total: Value):

    nauty = Nauty(nauty_exe)

    for graph in graphs:
        charges = defaultdict(list)

        if not iacm_to_elements:
            for key, partial_charge in iter_atomic_fragments(graph, nauty, shell):
                charges[key].append(partial_charge)
        else:
            for key, partial_charge in iter_atomic_fragments(graph, nauty, shell):
                charges[key].append(partial_charge)
        out_q.put(charges)

        progress.value += 1
        print_progress(progress.value, total.value, prefix='shell %d:' % shell)


def iter_atomic_fragments(graph: nx.Graph, nauty: Nauty, shell: int):
    for atom in graph.nodes():
        if not 'partial_charge' in graph.node[atom]:
            raise KeyError('Missing property "partial_charge" for atom {}'.format(atom))
        partial_charge = float(graph.node[atom]['partial_charge'])
        yield nauty.canonize_neighborhood(graph, atom, shell), partial_charge
