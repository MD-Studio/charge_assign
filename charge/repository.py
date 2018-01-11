import json
import os
import re
import zipfile
from collections import defaultdict
from itertools import groupby, product
from typing import Dict, Tuple, List, Callable
from zipfile import ZipFile

import jsonschema
import networkx as nx
import numpy

from charge.babel import convert_from, IOType
from charge.nauty import Nauty
from charge.settings import ATOMIC_SCHEMA, DB_SCHEMA, REPO_LOCATION, IACM_MAP


class LoadError(Exception):
    pass


class Repository:

    def __init__(self, location: str= REPO_LOCATION,
                 atomic_location: str=None, data_location: str=None, nauty: Nauty=None) -> None:

        self.__nauty = nauty or Nauty()

        self.charges_iacm = defaultdict(dict)
        self.charges_elem = defaultdict(dict)

        if atomic_location and data_location:
            self.__create(atomic_location, data_location, nauty)
            return

        if not zipfile.is_zipfile(location):
            raise LoadError('%s is not a valid zip file.' % location)

        with ZipFile(location, 'r') as zf:
            for name in zf.namelist():
                m = re.match('iacm_shell_(\d+)\.json$', name)
                if not m:
                    iacm = False
                    m = re.match('elements_shell_(\d+)\.json$', name)
                    if not m:
                        continue
                else:
                    iacm = True
                shell = int(m.group(1))

                data = json.loads(zf.read(name).decode())
                if not data or len(data) == 0:
                    raise LoadError('%s is invalid.' % name)
                # TODO charges/uncertainties -> histogram (hist, bin_edges)
                jsonschema.validate(data, DB_SCHEMA)

                for _, obj in data.items():
                    graph = self.__db_json_to_nx(obj)
                    key = nauty.canonize(graph)
                    if iacm:
                        self.charges_iacm[shell][key] = (obj['charges'][0], obj['uncertainties'][0])
                    else:
                        self.charges_elem[shell][key] = (obj['charges'][0], obj['uncertainties'][0])

    def __create(self, atomic_loc: str, data_loc: str) -> None:
        self.__atomic_loc = atomic_loc
        self.__data_loc = data_loc

        molids = [int(fn.replace('.lgf', ''))
                  for fn in os.listdir(data_loc) if fn.endswith('.lgf')]

        charges = dict()
        graphs = dict()
        canons = dict()
        for id in molids:
            charges[id], graphs[id], canons[id] = self.__read_lgf(os.path.join(data_loc, '%d.lgf' % id), self.__nauty)

        self.__iso = defaultdict(set)
        for iacm_key, group in groupby(molids, key=lambda id: canons[id]):
            isomorphics = list(group)
            if len(isomorphics) > 1:
                for id1, id2 in product(isomorphics):
                    self.__iso[id1].add(id2)
                    self.__iso[id2].add(id1)

        # TODO optimal linspace params (41, 401, 4001) -> scipy sparse matrices (needs special I/O ops)?
        def create_iacm(values: List[int], shell: int, key: str):
            histogram = numpy.histogram(values, numpy.linspace(-2, 2, 41))
            self.charges_iacm[shell][key] = histogram

        def create_elem(values: List[int], shell: int, key: str):
            histogram = numpy.histogram(values, numpy.linspace(-2, 2, 41))
            self.charges_elem[shell][key] = histogram

        self.__iterate(charges, create_iacm, create_elem)


    def __iterate(self, charges: Dict[int, Dict[int, float]],
                  func_iacm: Callable[[List[int], int, str], None],
                  func_elem: Callable[[List[int], int, str], None]):
        subdirs = dict((int(re.match('(\d+)$', fn).group(1)), fn) for fn in os.listdir(self.__data_loc) if os.path.isdir(fn))
        for shell in sorted(subdirs.keys()):
            # TODO check/change atomicfragments cmd to also export elementary atomic fragments
            for fn in filter(lambda fn: fn.endswith('.iacm.json'), os.listdir(os.path.join(self.__atomic_loc, subdirs[shell]))):
                with open(os.path.join(self.__atomic_loc, fn)) as f:
                    fragment = json.loads(f.read())
                    jsonschema.validate(fragment, ATOMIC_SCHEMA)
                    key = self.__nauty.canonize(self.__db_json_to_nx(fragment))
                    values = [charges[int(molid)][fragment['atom_mappings'][int(molid)][fragment['core_ids'][0]]] \
                              for molid in fragment['atom_mappings'].keys()]
                    func_iacm(values, shell, key)
            for fn in filter(lambda fn: fn.endswith('.elem.json'), os.listdir(os.path.join(self.__atomic_loc, subdirs[shell]))):
                with open(os.path.join(self.__atomic_loc, fn)) as f:
                    fragment = json.loads(f.read())
                    jsonschema.validate(fragment, ATOMIC_SCHEMA)
                    key = self.__nauty.canonize(self.__db_json_to_nx(fragment))
                    values = [charges[int(molid)][fragment['atom_mappings'][int(molid)][fragment['core_ids'][0]]] \
                              for molid in fragment['atom_mappings'].keys()]
                    func_elem(values, shell, key)

    def add(self, molid: int):
        charges = dict()
        for id in set(molid).union(self.__iso[molid]):
            charges[id], _, _ = self.__read_lgf(os.path.join(self.__data_loc, '%d.lgf' % id), self.__nauty)

        def add_val(values: List[int], hist: List[int], bin_edges: List[float]):
            for val in values:
                for k, right_edge in enumerate(bin_edges[1:]):
                    if bin_edges[k] <= val <= right_edge:
                        hist[k] += 1

        def add_iacm(values: List[int], shell: int, key: str):
            hist, bin_edges = self.charges_iacm[shell][key]
            add_val(values, hist, bin_edges)

        def add_elem(values: List[int], shell: int, key: str):
            hist, bin_edges = self.charges_elem[shell][key]
            add_val(values, hist, bin_edges)

        self.__iterate(charges, add_iacm, add_elem)


    def subtract(self, molid: int):
        charges = dict()
        for id in set(molid).union(self.__iso[molid]):
            charges[id], _, _ = self.__read_lgf(os.path.join(self.__data_loc, '%d.lgf' % id), self.__nauty)

        def sub_val(values: List[int], hist: List[int], bin_edges: List[float]):
            for val in values:
                for k, right_edge in enumerate(bin_edges[1:]):
                    if bin_edges[k] <= val <= right_edge:
                        hist[k] -= 1

        def sub_iacm(values: List[int], shell: int, key: str):
            hist, bin_edges = self.charges_iacm[shell][key]
            sub_val(values, hist, bin_edges)

        def sub_elem(values: List[int], shell: int, key: str):
            hist, bin_edges = self.charges_elem[shell][key]
            sub_val(values, hist, bin_edges)

        self.__iterate(charges, sub_iacm, sub_elem)

    def write(self, out: str):
        # TODO write to output file
        pass

    def __db_json_to_nx(self, obj: dict) -> nx.Graph:
        graph = nx.Graph()
        for i, (idx, type) in enumerate(obj['atoms'].items()):
            idx = int(idx)
            graph.add_node(idx, attr_dict={
                'atom_type': type,
                'core': idx in obj['core_ids'],
                'idx': i
            })
        for bond in obj['bonds']:
            graph.add_edge(bond[0], bond[1])
        return graph

    def __read_lgf(self, lgf: str, nauty: Nauty) -> Tuple[Dict[int, float], nx.Graph, str]:
        charges = {}
        graph = convert_from(lgf, IOType.LGF)
        for node, data in graph.nodes_iter(data=True):
            if not 'partial_charge' in data:
                raise ValueError('Missing attribute "partial_charge".')
            charges[int(node)] = float(data['partial_charge'])

        return charges, graph, nauty.canonize(graph, with_core=False)
