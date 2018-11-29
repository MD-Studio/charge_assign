import os
from collections import defaultdict
from io import IOBase
from itertools import groupby
from typing import Dict, List, Tuple, Union, Optional, AnyStr
from uuid import uuid4
from zipfile import ZipFile

import msgpack
import networkx as nx
from typing.io import IO

from charge.babel import convert_from, IOType
from charge.charge_types import Atom
from charge.molecule import atoms_neighborhoods_charges
from charge.multiprocessor import MultiProcessor
from charge.nauty import Nauty
from charge.settings import REPO_LOCATION

ChargeSet = Dict[int, Dict[str, List[float]]]
"""A collection of possible charges, indexed by shell size and \
        neighborhood canonical key.
"""


TraceableChargeSet = Dict[int, Dict[str, List[Tuple[float, int, Atom]]]]
"""A collection of possible charges with the molid of the molecule \
        they came from and the core atom of the neighborhood, \
        indexed by shell size and neighborhood canonical key.
"""


EitherChargeSet = Union[ChargeSet, TraceableChargeSet]

FileOrFileLike = Union[str, os.PathLike, IOBase, IO[AnyStr]]


class Repository:
    """A collection of atom charges by neighborhood.

    Args:
        min_shell: Minimum shell size in this repository.
        max_shell: Maximum shell size in this repository.
        nauty: Nauty instance.
        versioning: If True, assigns a unique int to the lists of charges on each change.

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
        iso_iacm: A dictionary mapping molids to lists of isomorphic \
                molids. Atoms use IACM types.
        iso_elem: A dictionary mapping molids to lists of isomorphic \
                molids. Atoms use plain elements.
    """
    def __init__(self,
                 min_shell: int=1,
                 max_shell: int=7,
                 nauty: Optional[Nauty]=None,
                 versioning: Optional[bool]=False) -> None:

        self.__nauty = nauty if nauty else Nauty()
        self.__min_shell = max(min_shell, 0)
        self.__max_shell = max_shell

        if not versioning:
            self.charges_iacm = defaultdict(lambda: defaultdict(list))  # type: EitherChargeSet
            self.charges_elem = defaultdict(lambda: defaultdict(list))  # type: EitherChargeSet
        else:
            self.charges_iacm = defaultdict(lambda: defaultdict(_VersioningList))  # type: EitherChargeSet
            self.charges_elem = defaultdict(lambda: defaultdict(_VersioningList))  # type: EitherChargeSet

        self.iso_iacm = defaultdict(list)
        self.iso_elem = defaultdict(list)

    @staticmethod
    def create_from(
            data_location: str,
            data_type: IOType=IOType.LGF,
            min_shell: int=1,
            max_shell: int=7,
            traceable: Optional[bool]=False,
            versioning: Optional[bool]=False,
            nauty: Optional[Nauty]=None
            ) -> 'Repository':

        """Creates a new Repository from a directory of files.

        Args:
            data_location: Path to the data directory.
            data_type: Type of the files to read.
            min_shell: Minimum shell size to compute.
            max_shell: Maximum shell size to compute.
            nauty: Nauty instance.
            versioning: If True, assigns a unique int to the lists of charges on each change.

        Returns:
            A new Repository with data read and processed.
        """
        repo = Repository(min_shell, max_shell, nauty, versioning)
        extension = data_type.get_extension()

        molids = [int(fn.replace(extension, ''))
                  for fn in os.listdir(data_location)
                  if fn.endswith(extension)]

        # load graphs
        graphs = repo.__read_graphs(
                molids, data_location, extension, data_type)

        # process with iacm atom types
        repo.charges_iacm = repo.__generate_charges(graphs, 'iacm', traceable, versioning)
        canons = repo.__make_canons(graphs)
        repo.iso_iacm = repo.__make_isomorphics(molids, canons)

        # process as plain elements
        repo.charges_elem = repo.__generate_charges(graphs, 'atom_type', traceable, versioning)
        canons = repo.__make_canons(graphs)
        repo.iso_elem = repo.__make_isomorphics(molids, canons)

        return repo

    @staticmethod
    def read(
            location: Optional[FileOrFileLike] = REPO_LOCATION,
            versioning: Optional[bool] = False,
            nauty: Optional[Nauty]=None
            ) -> 'Repository':
        """Create a Repository by loading from a zip file.

        The zip file must have been created by a call to write().

        Args:
            location: Path to the zip file (or file like object) to be read.
            nauty: Nauty instance.
            versioning: If True, assigns a unique int to the lists of charges on each change.

        Returns:
            A new Repository.
        """
        repo = Repository(nauty=nauty, versioning=versioning)
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

            if versioning:
                for _, chdct in repo.charges_iacm.items():
                    for key, charges in chdct.items():
                        chdct[key] = _VersioningList(charges)
                for _, chdct in repo.charges_elem.items():
                    for key, charges in chdct.items():
                        chdct[key] = _VersioningList(charges)
        return repo

    def write(self, out: Optional[FileOrFileLike] = REPO_LOCATION) -> None:
        """Write the repository to disk as a zip file.

        Args:
            out: Path to the zip file (or file like object) to be written.
        """
        with ZipFile(out, mode='w') as zf:
            zf.writestr('meta', msgpack.packb(
                (self.__min_shell, self.__max_shell)))
            zf.writestr('charges_iacm', msgpack.packb(self.charges_iacm))
            zf.writestr('charges_elem', msgpack.packb(self.charges_elem))
            zf.writestr('iso_iacm', msgpack.packb(self.iso_iacm))
            zf.writestr('iso_elem', msgpack.packb(self.iso_elem))

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
            traceable: bool=False,
            versioing: bool=False
            ) -> Dict[int, Dict[str, List[float]]]:
        """Generate charges for all shell sizes and neighborhoods."""
        if not versioing:
            charges = defaultdict(lambda: defaultdict(list))
        else:
            charges = defaultdict(lambda: defaultdict(_VersioningList))

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
                charges[shell][key].sort()

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


class _ReadWorker:
    """Reads a graph from a file."""
    def __init__(self, data_location: str, extension: str, data_type: IOType):
        self.__data_location = data_location
        self.__extension = extension
        self.__data_type = data_type

    def process(self, molid: int) -> Tuple[int, nx.Graph]:
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

    def process(self, molid: int, graph: nx.Graph) -> Tuple[int, str]:
        return molid, self.__nauty.canonize(graph)


class _ChargeWorker:
    """Collects charges per neighborhood from the given graph."""
    def __init__(self, shell: int, color_key: str):
        self.__shell = shell
        self.__color_key = color_key
        self.__nauty = Nauty()

    def process(self, molid: int, graph: nx.Graph) -> Dict[str, List]:
        charges = defaultdict(list)

        for _, key, partial_charge in atoms_neighborhoods_charges(
                graph, self.__nauty, self.__shell, self.__color_key):
            charges[key].append(partial_charge)

        return charges


class _TraceableChargeWorker:
    """Collects charges per neighborhood from the given graph.

    Charges come with the molid and atom they came from, so you get \
    lists of triples in the repository, rather than lists of floats.
    """
    def __init__(self, shell: int, color_key: str):
        self.__shell = shell
        self.__color_key = color_key
        self.__nauty = Nauty()

    def process(self, molid: int, graph: nx.Graph) -> Dict[str, List]:
        charges = defaultdict(list)

        for atom, key, partial_charge in atoms_neighborhoods_charges(
                graph, self.__nauty, self.__shell, self.__color_key):
            charges[key].append((partial_charge, molid, atom))

        return charges


def proxy():

    changer_methods = {'append', 'clear', 'extend', 'insert', 'pop', 'remove',
                       '__add__', '__delitem__', '__iadd__', '__init__', '__imul__', '__mul__', '__rmul__', '__setitem__'}

    def proxy_decorator(func):
        def wrapper(*args, **kw):
            # self.version = uuid4().int
            args[0].version = uuid4().int
            return func(*args, **kw)

        wrapper.__name__ = func.__name__
        return wrapper

    dct = list.__dict__.copy()
    for key, value in dct.items():
        if key in changer_methods:
            dct[key] = proxy_decorator(value)
    return type("proxy_list", (list,), dct)


_VersioningList = proxy()
