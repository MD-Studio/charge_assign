import bisect
import math
import os
import time
from collections import defaultdict
from itertools import groupby
from typing import Callable, Dict, List, Tuple, Union
from zipfile import ZipFile, ZIP_DEFLATED

import msgpack
import networkx as nx

from charge.babel import convert_from, IOType
from charge.nauty import Nauty
from charge.settings import REPO_LOCATION, IACM_MAP
from charge.charge_types import Atom
from charge.multiprocessor import MultiProcessor


ChargeSet = Dict[int, Dict[str, List[float]]]
"""A collection of possible charges, indexed by shell size and \
        neighborhood canonical key.
"""


TraceableChargeSet = Dict[int, Dict[str, List[Tuple[float, int, Atom]]]]
"""A collection of possible charges with the molid of the molecule \
        they came from and the core atom of the neighborhood,
"""


EitherChargeSet = Union[ChargeSet, TraceableChargeSet]


class Repository:
    """A collection of atom charges by neighborhood.

    Args:
        min_shell: Minimum shell size in this repository.
        max_shell: Maximum shell size in this repository.

    Attributes:
        charges_iacm: A dictionary, keyed by shell size, of \
                dictionaries, keyed by neighborhood hash, of lists of \
                charges (floats) for the atom at the center of the \
                neighborhood. Atoms use IACM types. Optionally, may \
                contain tuples of (charge, molid, atom) if the \
                repository is traceable.
        charges_elem: A dictionary, keyed by shell size, of \
                dictionaries, keyed by neighborhood hash, of lists of \
                charges (floats) for the atom at the center of the \
                neighborhood. Atoms use plain elements. Optionally, may \
                contain tuples of (charge, molid, atom) if the \
                repository is traceable.
    """
    def __init__(self,
                 min_shell: int=1,
                 max_shell: int=7) -> None:

        self.__nauty = Nauty()
        self.__min_shell = max(min_shell, 0)
        self.__max_shell = max_shell

        self.charges_iacm = defaultdict(lambda: defaultdict(list))  # type: EitherChargeSet
        self.charges_elem = defaultdict(lambda: defaultdict(list))  # type: EitherChargeSet
        self.iso_iacm = defaultdict(list)
        self.iso_elem = defaultdict(list)

    @staticmethod
    def create_from(
            data_location: str,
            data_type: IOType=IOType.LGF,
            min_shell: int=1,
            max_shell: int=7,
            traceable: bool=False
            ) -> 'Repository':

        """Creates a new Repository from a directory of files.

        Args:
            data_location: Path to the data directory.
            data_type: Type of the files to read.
            min_shell: Minimum shell size to compute.
            max_shell: Maximum shell size to compute.

        Returns:
            A new Repository with data read and processed.
        """
        repo = Repository(min_shell, max_shell)
        extension = data_type.get_extension()

        molids = [int(fn.replace(extension, ''))
                  for fn in os.listdir(data_location)
                  if fn.endswith(extension)]

        # load graphs
        graphs = repo.__read_graphs(
                molids, data_location, extension, data_type)

        # process with iacm atom types
        repo.charges_iacm = repo.__generate_charges(graphs, 'iacm', traceable)
        canons = repo.__make_canons(graphs)
        repo.iso_iacm = repo.__make_isomorphics(molids, canons)

        # process as plain elements
        repo.charges_elem = repo.__generate_charges(graphs, 'atom_type', traceable)
        canons = repo.__make_canons(graphs)
        repo.iso_elem = repo.__make_isomorphics(molids, canons)

        return repo

    @staticmethod
    def read(location: str = REPO_LOCATION) -> 'Repository':
        """Create a Repository by loading from a zip file.

        The zip file must have been created by a call to write().

        Args:
            location: Path to the zip file to be read.

        Returns:
            A new Repository.
        """
        repo = Repository()
        with ZipFile(location, mode='r') as zf:
            repo.__min_shell, repo.__max_shell = msgpack.unpackb(
                    zf.read('meta'), encoding='utf-8')
            repo.charges_iacm = msgpack.unpackb(
                    zf.read('charges_iacm'), encoding='utf-8')
            repo.charges_elem = msgpack.unpackb(
                    zf.read('charges_elem'), encoding='utf-8')
            repo.iso_iacm = msgpack.unpackb(
                    zf.read('iso_iacm'), encoding='utf-8')
            repo.iso_elem = msgpack.unpackb(
                    zf.read('iso_elem'), encoding='utf-8')
        return repo

    def write(self, out: str) -> None:
        """Write the repository to disk as a zip file.

        Args:
            out: Path to the zip file to be written.
        """
        with ZipFile(out, mode='w') as zf:
            zf.writestr('meta', msgpack.packb(
                (self.__min_shell, self.__max_shell)))
            zf.writestr('charges_iacm', msgpack.packb(self.charges_iacm))
            zf.writestr('charges_elem', msgpack.packb(self.charges_elem))
            zf.writestr('iso_iacm', msgpack.packb(self.iso_iacm))
            zf.writestr('iso_elem', msgpack.packb(self.iso_elem))

    # TODO optional: add/subtract isomorphic molids
    def add(self, data_location: str, molid: int, data_type: IOType) -> None:
        """Add a new molecule to the Repository.

        Args:
            data_location: Path to the data directory.
            molid: Molecule id to load.
            data_type: Type of the file to load.
        """
        def a(shell, key, partial_charge, repo):
            if shell not in repo:
                repo[shell] = dict()
            if key not in repo[shell]:
                repo[shell][key] = []
            bisect.insort_left(repo[shell][key], partial_charge)

        self.__iterate(
                data_location, molid, data_type,
                lambda shell, key, partial_charge: a(
                    shell, key, partial_charge, self.charges_iacm),
                lambda shell, key, partial_charge: a(
                    shell, key, partial_charge, self.charges_elem))

    def subtract(
            self,
            data_location: str,
            molid: int,
            data_type: IOType
            ) -> None:
        """Remove a molecule from the Repository.

        Args:
            data_location: Path to the data directory.
            molid: Molecule id to remove.
            data_type Type of file to load.
        """
        def s(shell, key, partial_charge, repo):
            if shell in repo and key in repo[shell]:
                repo[shell][key].pop(
                        bisect.bisect_left(repo[shell][key], partial_charge))
                if len(repo[shell][key]) == 0:
                    del repo[shell][key]
                if len(repo[shell]) == 0:
                    del repo[shell]

        self.__iterate(
                data_location, molid, data_type,
                lambda shell, key, partial_charge: s(
                    shell, key, partial_charge, self.charges_iacm),
                lambda shell, key, partial_charge: s(
                    shell, key, partial_charge, self.charges_elem))

    def __read_graphs(
            self,
            molids: List[int],
            data_location: str,
            ext: str,
            data_type: IOType
            ) -> List[Tuple[int, nx.Graph]]:
        """Read graphs from a directory of input files."""

        graphs = []
        with MultiProcessor(
                _ReadWorker, (data_location, ext, data_type)) as mp:
            for molid, graph in mp.processed(molids, 'reading files'):
                graphs.append((molid, graph))

        return graphs

    def __generate_charges(
            self,
            graphs: List[Tuple[int, nx.Graph]],
            color_key: str,
            traceable: bool=False
            ) -> Dict[int, Dict[str, List[float]]]:
        """Generate charges for all shell sizes and neighborhoods."""
        charges = defaultdict(lambda: defaultdict(list))

        if traceable:
            Worker = _TraceableChargeWorker
        else:
            Worker = _ChargeWorker

        for shell in range(self.__min_shell, self.__max_shell + 1):
            with MultiProcessor(Worker, (shell, color_key)) as mp:
                for c in mp.processed(graphs, 'shell %d' % shell):
                    for key, values in c.items():
                        charges[shell][key] += values

        for shell in range(self.__min_shell, self.__max_shell + 1):
            for key, values in charges[shell].items():
                charges[shell][key] = sorted(values)

        return charges

    def __make_isomorphics(
            self,
            molids: List[int],
            canons: Dict[int, str]
            ) -> Dict[int, List[int]]:
        """Find isomorphic molids and create map of them."""
        isomorphics = defaultdict(list)
        molids_by_key = sorted(molids, key=lambda molid: canons[molid])
        for _, group in groupby(molids_by_key, key=lambda molid: canons[molid]):
            isogroup = list(group)
            if len(isogroup) > 1:
                for molid in isogroup:
                    isomorphics[molid] = isogroup
        return isomorphics

    def __make_canons(
            self,
            graphs: List[Tuple[int, nx.Graph]]
            ) -> Dict[int, str]:
        """Canonicalize the given graphs using Nauty."""
        canons = dict()
        with MultiProcessor(_CanonicalizationWorker) as mp:
            for molid, canon in mp.processed(graphs):
                canons[molid] = canon
        return canons

    def __iterate(
            self,
            data_location: str, molid: int, data_type: IOType,
            callable_iacm: Callable[[int, str, float], None],
            callable_elem: Callable[[int, str, float], None]):
        ids_iacm = {molid}
        ids_elem = {molid}
        if molid in self.iso_iacm:
            ids_iacm.union(set(self.iso_iacm[molid]))
        if molid in self.iso_elem:
            ids_iacm.union(set(self.iso_elem[molid]))

        extension = data_type.get_extension()

        for molid in ids_iacm:
            path = os.path.join(data_location, '%d%s' % (molid, extension))
            with open(path, 'r') as f:
                graph = convert_from(f.read(), data_type)
                for shell in range(1, self.__max_shell + 1):
                    for key, partial_charge, _ in _iter_atomic_fragments(
                            graph, self.__nauty, shell):
                        callable_iacm(shell, key, partial_charge)

        for molid in ids_elem:
            path = os.path.join(data_location, '%d%s' % (molid, extension))
            with open(path, 'r') as f:
                graph = convert_from(f.read(), data_type)
                for v, data in graph.nodes(data=True):
                    graph.node[v]['atom_type'] = IACM_MAP[data['atom_type']]
                for shell in range(1, self.__max_shell + 1):
                    for key, partial_charge, _ in _iter_atomic_fragments(
                            graph, self.__nauty, shell):
                        callable_elem(shell, key, partial_charge)


class _ReadWorker:
    """Reads a graph from a file."""
    def __init__(self, data_location: str, extension: str, data_type: IOType):
        self.__data_location = data_location
        self.__extension = extension
        self.__data_type = data_type

    def process(self, molid: int) -> None:
        filename = os.path.join(
                self.__data_location, '%d%s' % (molid, self.__extension))
        with open(filename, 'r') as f:
            graph = convert_from(f.read(), self.__data_type)
            return molid, graph


class _CanonicalizationWorker:
    """Returns a canonical hash of a graph.

    Isomorphic graphs return the same hash (key).
    """
    def __init__(self):
        self.__nauty = Nauty()

    def process(self, molid: int, graph: nx.Graph) -> str:
        return molid, self.__nauty.canonize(graph)


class _ChargeWorker:
    """Collects charges per neighborhood from the given graph."""
    def __init__(self, shell: int, color_key: str):
        self.__shell = shell
        self.__color_key = color_key
        self.__nauty = Nauty()

    def process(self, molid: int, graph: nx.Graph) -> defaultdict(list):
        charges = defaultdict(list)

        for key, partial_charge, _ in _iter_atomic_fragments(
                graph, self.__nauty, self.__shell, self.__color_key):
            charges[key].append(partial_charge)

        return charges


class _TraceableChargeWorker:
    """Collects charges per neighborhood from the given graph."""
    def __init__(self, shell: int, color_key: str):
        self.__shell = shell
        self.__color_key = color_key
        self.__nauty = Nauty()

    def process(self, molid: int, graph: nx.Graph) -> defaultdict(list):
        charges = defaultdict(list)

        for key, partial_charge, atom in _iter_atomic_fragments(
                graph, self.__nauty, self.__shell, self.__color_key):
            charges[key].append((partial_charge, molid, atom))

        return charges


def _iter_atomic_fragments(graph: nx.Graph, nauty: Nauty, shell: int, color_key: str):
    """Yields all atomic neighborhoods in the graph of the given shell size."""
    for atom in graph.nodes():
        if 'partial_charge' not in graph.node[atom]:
            raise KeyError(
                'Missing property "partial_charge" for atom {}'.format(atom))
        partial_charge = float(graph.node[atom]['partial_charge'])
        yield nauty.canonize_neighborhood(graph, atom, shell, color_key), partial_charge, atom
