import hashlib
import os
import subprocess
from itertools import groupby
from typing import Any, Dict, Tuple

import msgpack
import networkx as nx

from charge.settings import NAUTY_EXC
from charge.util import bfs_nodes


class Nauty:

    def __init__(self, executable: str=NAUTY_EXC) -> None:
        if not os.path.isfile(executable) or not os.access(executable, os.X_OK):
            raise ValueError('Could not find dreadnaut executable at: "%s". Did you install nauty (http://users.cecs.'
                             'anu.edu.au/~bdm/nauty/)?' % executable)
        self.exe = executable
        self.__process = subprocess.Popen(
            [self.exe],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            bufsize=0,
            close_fds=True
        )

    def __del__(self):
        if hasattr(self, '__process'):
            if not self.__process.poll():
                self.__process.terminate()

    def canonize_neighborhood(self, graph: nx.Graph, atom: Any, shell: int, color_key='atom_type') -> str:

        if shell > 0:
            fragment = graph.subgraph(bfs_nodes(graph, atom, max_depth=shell))
        else:
            fragment = graph.subgraph(atom)

        return self.canonize(fragment, color_key=color_key, core=atom)

    def canonize(self, graph: nx.Graph, color_key='atom_type', core:Any=None) -> str:

        cmd, idx = self.__nauty_input(graph, color_key, core)

        if self.__process.poll():
            self.__process = subprocess.Popen(
                [self.exe],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                bufsize=0,
                close_fds=True
            )

        self.__process.stdin.write(cmd.encode())
        self.__process.stdin.flush()

        out = self.__process.stdout.read(1000)
        while not b'END' in out:
            out += self.__process.stdout.read(1000)

        return self.__nauty_output(out.strip().decode(), graph, idx, color_key, core)

    def __nauty_input(self, graph: nx.Graph, color_key:str, core: Any) -> Tuple[str, Dict[int, Any]]:
        idx = {v: i for i, v in enumerate(graph.nodes())}
        input_str = ' n={num_atoms} g {edges}. f=[{node_partition}] cxb"END\n"->>\n'.format(
            num_atoms=graph.number_of_nodes(),
            edges=self.__nauty_edges(graph, idx),
            node_partition=self.__nauty_node_partition(graph, idx, color_key, core),
        )
        return input_str, idx

    def __nauty_edges(self, graph: nx.Graph, idx: Dict[int, Any]) -> str:
        bonds = sorted(sorted((idx[u], idx[v])) for u, v in  graph.edges())

        return ';'.join('{0}:{1}'.format(u, v) for u, v in bonds)

    def __nauty_node_partition(self, graph: nx.Graph, idx: Dict[int, Any], color_key:str, core:Any) -> str:
        if core:
            nodes = sorted((v == core, data[color_key], idx[v]) for v, data in graph.nodes(data=True))
            nodes = groupby(nodes, key=lambda x: (x[0], x[1]))
            return '|'.join(','.join(str(v) for _, _, v in group) for _, group in nodes)
        else:
            nodes = sorted((data[color_key], idx[v]) for v, data in graph.nodes(data=True))
            nodes = groupby(nodes, key=lambda x: x[0])
            return '|'.join(','.join(str(v) for _, v in group) for _, group in nodes)

    def __nauty_output(self, nautstr: str, graph: nx.Graph, idx: Dict[int, Any], color_key:str, core:Any) -> str:
        nodes = {v: k for k, v in idx.items()}
        lines = [line.strip() for line in nautstr.split('seconds')[-1].strip().split('\n')]

        adj = [[int(val) for val in line[:-1].split(':')[-1].split()] for line in lines if ':' in line]

        if core:
            colors, core = list(zip(*(
                (graph.node[nodes[idx]][color_key], nodes[idx] == core) for idx in map(int, lines[0].split()))))
            print([adj, colors, core])
            return hashlib.md5(msgpack.packb([adj, colors, core])).hexdigest()
        else:
            colors = [graph.node[nodes[idx]][color_key] for idx in map(int, lines[0].split())]
            return hashlib.md5(msgpack.packb([adj, colors])).hexdigest()
