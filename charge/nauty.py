import hashlib
import os
import subprocess
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

        cmd, idx = self.__nauty_input(graph, color_key=color_key)

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

    def __nauty_input(self, graph: nx.Graph, color_key:str) -> Tuple[str, Dict[int, Any]]:
        idx = {v: i for i, v in enumerate(graph.nodes())}
        input_str = ' n={num_atoms} g {edges}. f=[{node_partition}] cxb"END\n"->>\n'.format(
            num_atoms=graph.number_of_nodes(),
            edges=self.__nauty_edges(graph, idx),
            node_partition=self.__nauty_node_partition(graph, idx, color_key),
        )

        return input_str, idx

    def __nauty_edges(self, graph: nx.Graph, idx: Dict[int, Any]) -> str:
        bonds = sorted(map(lambda e: sorted((idx[e[0]], idx[e[1]])), graph.edges()))

        return ';'.join(map(lambda bond: '{0}:{1}'.format(bond[0], bond[1]), bonds))

    def __nauty_node_partition(self, graph: nx.Graph, idx: Dict[int, Any], color_key:str) -> str:
        colors = sorted(set(data[color_key] for _, data in graph.nodes(data=True)))
        return '|'.join(
                ','.join(map(str, sorted(idx[v] for v, data in graph.nodes(data=True) if data[color_key] == color)))
            for color in colors)

    def __nauty_output(self, nautstr: str, graph: nx.Graph, idx: Dict[int, Any], color_key:str, core:Any) -> str:
        nodes = {v: k for k, v in idx.items()}
        lines = [line.strip() for line in nautstr.split('seconds')[-1].strip().split('\n')]

        adj = [[int(val) for val in line[:-1].split(':')[-1].split()] for line in lines if ':' in line]

        if core:
            colors, core = list(zip(*map(
                lambda idx: (graph.node[nodes[idx]][color_key], nodes[idx] == core),
                map(int, lines[0].split()))))
            return hashlib.md5(msgpack.packb([adj, colors, core])).hexdigest()
        else:
            colors = [graph.node[nodes[idx]][color_key] for idx in map(int, lines[0].split())]
            return hashlib.md5(msgpack.packb([adj, colors])).hexdigest()
