import hashlib
import os
import subprocess
import tempfile
from collections import defaultdict, deque
from itertools import groupby
from typing import Any, Dict, Tuple

import msgpack
import networkx as nx

from charge.settings import NAUTY_EXC


class Nauty:

    def __init__(self, executable: str=NAUTY_EXC) -> None:
        if not os.path.isfile(executable) or not os.access(executable, os.X_OK):
            raise ValueError('Could not find dreadnaut executable at: "%s". Did you install nauty (http://users.cecs.'
                             'anu.edu.au/~bdm/nauty/)?' % executable)
        self.__exe = executable

    def canonize_neighborhood(self, graph: nx.Graph, atom: Any, shell: int, color_key='atom_type') -> str:
        fragment = nx.Graph()

        visited = {atom}
        depth = defaultdict(int)
        queue = deque([atom])

        while len(queue) > 0:
            v = queue.popleft()

            attr = graph.node[v].copy()
            attr['core'] = v == atom
            fragment.add_node(v, **attr)

            if depth[v] < shell:
                for child in graph.neighbors(v):
                    if not child in visited:
                        depth[child] = depth[v] + 1
                        visited.add(child)
                        queue.append(child)

        nodes = list(fragment.nodes())
        for i, u in enumerate(nodes):
            for v in nodes[i+1:]:
                if graph.has_edge(u, v):
                    fragment.add_edge(u, v)

        return self.canonize(fragment, color_key=color_key)

    def canonize(self, graph: nx.Graph, color_key='atom_type', with_core:bool=True) -> str:

        with tempfile.TemporaryFile(buffering=0) as tmp:
            cmd, idx = self.__nauty_input(graph, color_key=color_key)
            tmp.write(cmd.encode())
            tmp.seek(0)

            p = subprocess.Popen(
                [self.__exe],
                shell=True,
                stdin=tmp,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )

            out, err = p.communicate()

            if len(err) > 0:
                raise Exception(err.decode())
            if len(out) == 0:
                raise Exception()

        return self.__nauty_output(out.strip().decode(), graph, idx, color_key=color_key, with_core=with_core)

    def __nauty_input(self, graph: nx.Graph, color_key='atom_type') -> Tuple[str, Dict[int, Any]]:
        idx = {v:i for i,v in enumerate(graph.nodes())}
        input_str = ' n={num_atoms} g {edges}. f=[{node_partition}] cxb'.format(
            num_atoms=graph.number_of_nodes(),
            edges=self.__nauty_edges(graph, idx),
            node_partition=self.__nauty_node_partition(graph, idx, color_key=color_key),
        )

        return input_str, idx

    def __nauty_edges(self, graph: nx.Graph, idx: Dict[int, Any]) -> str:
        bonds = sorted(map(lambda e: sorted((idx[e[0]], idx[e[1]])), graph.edges()))

        return ';'.join(map(lambda bond: '{0}:{1}'.format(bond[0], bond[1]), bonds))

    def __nauty_node_partition(self, graph: nx.Graph, idx: Dict[int, Any], color_key:str='atom_type') -> str:
        def to_nx_node(v):
            if not color_key in v[1]:
                raise ValueError('Missing attribute: %s' % color_key)
            return (idx[v[0]], v[1][color_key])

        atoms = sorted(map(to_nx_node, graph.nodes(data=True)))

        return '|'.join(map(lambda group: ','.join(map(lambda atom: str(atom[0]), group[1])),
                            groupby(atoms, key=lambda atom: atom[1])))

    def __nauty_output(self, nautstr: str, graph: nx.Graph, idx: Dict[int, Any], color_key:str='atom_type',
                       with_core:bool=True) -> str:
        nodes = {v: k for k, v in idx.items()}
        lines = [line.strip() for line in nautstr.split('seconds')[-1].strip().split('\n')]

        adj = [[int(val) for val in line[:-1].split(':')[-1].split()] for line in lines[1:]]

        if with_core:
            colors, core = list(zip(*map(
                lambda idx: (graph.node[nodes[idx]][color_key], graph.node[nodes[idx]]['core']),
                map(int, lines[0].split()))))
            return hashlib.md5(msgpack.packb([adj, colors, core])).hexdigest()
        else:
            colors = [graph.node[nodes[idx]][color_key] for idx in map(int, lines[0].split())]
            return hashlib.md5(msgpack.packb([adj, colors])).hexdigest()
