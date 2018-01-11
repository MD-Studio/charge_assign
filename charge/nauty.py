import hashlib
import os
import subprocess
import tempfile
from itertools import groupby

import msgpack
import networkx as nx

from charge.settings import NAUTY_EXC


class Nauty:

    def __init__(self, executable: str=NAUTY_EXC) -> None:
        if not os.path.isfile(executable) or not os.access(executable, os.X_OK):
            raise ValueError('Could not find dreadnaut executable at: "%s". Did you install nauty (http://users.cecs.'
                             'anu.edu.au/~bdm/nauty/)?' % executable)
        self.__exe = executable

    def canonize(self, graph: nx.Graph, color_key='atom_type', with_core=True) -> str:

        with tempfile.TemporaryFile(buffering=0) as tmp:
            tmp.write(self.__nauty_input(graph, color_key=color_key).encode())
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

        return self.__nauty_output(out.strip().decode(), graph, color_key=color_key, with_core=with_core)

    def __nauty_input(self, graph: nx.Graph, color_key='atom_type') -> str:
        input_str = ' n={num_atoms} g {edges}. f=[{node_partition}] cxb'.format(
            num_atoms=graph.number_of_nodes(),
            edges=self.__nauty_edges(graph),
            node_partition=self.__nauty_node_partition(graph, color_key=color_key),
        )

        return input_str

    def __nauty_edges(self, graph: nx.Graph) -> str:
        bonds = sorted(map(lambda e: sorted((graph.node[e[0]]['idx'], graph.node[e[1]]['idx'])), graph.edges_iter()))

        return ';'.join(map(lambda bond: '{0}:{1}'.format(bond[0], bond[1]), bonds))

    def __nauty_node_partition(self, graph: nx.Graph, color_key='atom_type') -> str:
        atoms = sorted(map(lambda v: (v[1]['idx'], v[1][color_key]), graph.nodes_iter(data=True)))

        return '|'.join(map(lambda group: ','.join(map(lambda atom: str(atom[0]), group[1])),
                            groupby(atoms, key=lambda atom: atom[1])))

    def __nauty_output(self, nautstr: str, graph: nx.Graph, color_key='atom_type', with_core=True) -> str:
        idxmap = {v: k for k, v in nx.get_node_attributes(graph, 'idx').items()}
        lines = [line.strip() for line in nautstr.split('seconds')[-1].strip().split('\n')]

        adj = [[int(val) for val in line[:-1].split(':')[-1].split()] for line in lines[1:]]

        if with_core:
            colors, core = zip(*map(
                lambda idx: (graph.node[idxmap[idx]][color_key], graph.node[idxmap[idx]]['core']), map(int, lines[0].split())))
            return hashlib.md5(msgpack.packb([adj, colors, core])).hexdigest()
        else:
            colors = map(lambda idx: graph.node[idxmap[idx]][color_key], map(int, lines[0].split()))
            return hashlib.md5(msgpack.packb([adj, colors])).hexdigest()
