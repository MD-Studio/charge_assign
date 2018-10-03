import networkx as nx
from typing import Generator, Tuple

from charge.nauty import  Nauty
from charge.types import Atom


def atoms_neighborhoods_charges(
        graph: nx.Graph,
        nauty: Nauty,
        shell: int,
        atom_type_key: str
        ) -> Generator[Tuple[Atom, str, float], None, None]:
    """Yields neighborhood hash and partial charge for each atom.

    Args:
        nauty: The Nauty instance to use to canonize the neighborhoods.
        shell: The shell size to use to make the neighborhood
        atom_type_key: The name of the atom type attribute to use

    Yields:
        Tuples containing an atom, the neighborhood hash, and the \
                partial charge of the atom.
    """
    for atom in graph.nodes():
        if 'partial_charge' not in graph.node[atom]:
            raise KeyError(
                'Missing property "partial_charge" for atom {}'.format(atom))
        partial_charge = graph.node[atom]['partial_charge']
        key = nauty.canonize_neighborhood(graph, atom, shell, atom_type_key)
        yield atom, key, partial_charge
