import hashlib
import os
import subprocess
from itertools import groupby
from typing import Any, Dict, Tuple, List

import msgpack
import networkx as nx

from charge.settings import NAUTY_EXC
from charge.util import bfs_nodes


Color = Tuple[bool, str]
NautyEdges = List[Tuple[int, int]]


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
        try:
            if not self.__process.poll():
                self.__process.stdin.write('q'.encode('utf-8'))
                self.__process.stdin.close()
                self.__process.stdout.close()
                self.__process.stderr.close()
                self.__process.wait(timeout=1)
        except ValueError:
            pass

    def canonize_neighborhood(self, graph: nx.Graph, atom: Any, shell: int, color_key='atom_type') -> str:
        if shell > 0:
            fragment = graph.subgraph(bfs_nodes(graph, atom, max_depth=shell))
        else:
            fragment = graph.subgraph(atom)

        result = self.canonize(fragment, color_key=color_key, core=atom)
        return result

    def canonize(self, graph: nx.Graph, color_key='atom_type', core: Any=None) -> str:
        node_colors = list()
        for node, color_str in graph.nodes(data=color_key):
            node_colors.append((node == core, color_str))

        input_str, colors = self.__make_nauty_input(graph, node_colors)

        if self.__process.poll():
            self.__process = subprocess.Popen(
                [self.exe],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                bufsize=0,
                close_fds=True
            )

        self.__process.stdin.write(input_str.encode())
        self.__process.stdin.flush()

        out = self.__process.stdout.read(1000)
        while not b'END' in out:
            out += self.__process.stdout.read(1000)

        output_str = out.strip().decode()
        canonical_nodes, adjacency_list = self.__parse_nauty_output(output_str, node_colors)
        key = self.__make_hash(canonical_nodes, adjacency_list, colors)
        return key

    def __make_nauty_input(
            self,
            graph: nx.Graph,
            node_colors: List[Color]
            ) -> Tuple[str, List[Color]]:

        to_nauty_id = { v: i for i, v in enumerate(graph.nodes()) }

        nauty_edges = self.__make_nauty_edges(graph.edges(), to_nauty_id)
        colors, partition = self.__make_partition(list(graph.nodes()), node_colors, to_nauty_id)

        edges_str = self.__format_edges(nauty_edges)
        partition_str = self.__format_partition(partition)

        input_str = ' n={num_atoms} g {edges}. f=[{partition}] cxb"END\n"->>\n'.format(
                num_atoms=graph.number_of_nodes(),
                edges=edges_str,
                partition=partition_str)

        return input_str, colors

    def __make_nauty_edges(
            self,
            edges: List[Tuple[Any, Any]],
            to_nauty_id: Dict[Any, int]
            ) -> NautyEdges:

        nauty_edges = list()
        for u, v in edges:
            nauty_edges.append((to_nauty_id[u], to_nauty_id[v]))
        nauty_edges.sort()
        return nauty_edges

    def __make_partition(
            self,
            nodes: List[Any],
            node_colors: List[Tuple[bool, str]],
            to_nauty_id: Dict[Any, int]
            ) -> Tuple[List[Color], Any]:

        def by_color(node_and_color: Tuple[int, Color]) -> Color:
            return node_and_color[1]

        def get_node(node_and_color: Tuple[int, Color]) -> int:
            return node_and_color[0]

        colored_nauty_nodes = list()
        for node_id, color in enumerate(node_colors):
            colored_nauty_nodes.append((to_nauty_id[nodes[node_id]], color))

        colored_nauty_nodes.sort(key=by_color)

        colors = list()
        partition = list()
        for color, node_and_colors in groupby(colored_nauty_nodes, key=by_color):
            colors.append(color)
            nauty_ids = sorted(map(get_node, node_and_colors))
            partition.append((color, nauty_ids))

        return colors, partition

    def __format_edges(self, nauty_edges: NautyEdges) -> str:
        nauty_edge_strings = list()
        for u, v in nauty_edges:
            nauty_edge_strings.append('{}:{}'.format(u, v))

        return ';'.join(nauty_edge_strings)

    def __format_partition(self, partition: Any) -> str:
        nauty_cell_strings = list()
        for _, nauty_ids in partition:
            nauty_id_strs = map(str, nauty_ids)
            nauty_cell_strings.append(','.join(nauty_id_strs))

        return '|'.join(nauty_cell_strings)

    def __parse_nauty_output(
            self,
            nauty_output: str,
            node_colors: List[Color]
            ) -> Tuple[List[int], NautyEdges]:

        def get_color(nauty_id_str: str) -> Color:
            return node_colors[int(nauty_id_str)]

        data = nauty_output.split('seconds')[-1].strip()
        lines = data.split('\n')

        canonical_nodes_str = ''
        i = 0
        while ':' not in lines[i]:
            canonical_nodes_str += lines[i][0:-1]
            i += 1

        canonical_nodes_strs = canonical_nodes_str.split()
        canonical_nodes = list(map(get_color, canonical_nodes_strs))

        adjacency_list_lines = lines[i:-1]
        adjacency_list = list()
        for line in adjacency_list_lines:
            parts = line.split(':')
            node_id = int(parts[0])
            neighbors_str = parts[1][0:-1].strip()
            neighbors_strs = neighbors_str.split(' ')
            neighbors = list(map(int, neighbors_strs))
            adjacency_list.append((node_id, neighbors))

        adjacency_list.sort()

        return canonical_nodes, adjacency_list

    def __make_hash(
            self,
            canonical_nodes: List[int],
            adjacency_list: List[Tuple[int, List[int]]],
            colors: List[Color]
            ) -> str:

        canonical_signature = [canonical_nodes, adjacency_list, colors]
        canonical_bytes = msgpack.packb(canonical_signature)
        return hashlib.md5(canonical_bytes).hexdigest()
