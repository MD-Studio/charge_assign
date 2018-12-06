import argparse
import networkx as nx
from typing import Optional

from charge.repository import Repository
from charge.chargers import Charger, CDPCharger
from charge.collectors import HistogramCollector


_charger = None  # type: Optional[Charger]


def parse_arguments() -> str:
    """Parses command line arguments.

    Returns:
        The location of the repository to use, as specified using
                the -r command line option.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--repository', type=str,
                        help='Path to Repository (zip file) to use for'
                        ' reference charges.')
    args = parser.parse_args()
    return args.repository


def init(repo_location: Optional[str] = None) -> None:
    """Create the charger from a repository.

    Gets the repository location from the command line if none is
    given.

    Args:
        repo_location: The path to the repository to use.

    """
    global _charger
    if not repo_location:
        repo_location = parse_arguments()
    repo = Repository.read(repo_location)
    _charger = CDPCharger(repo,
                          rounding_digits=3)


def charge(graph: nx.Graph, total_charge: int) -> None:
    """Adds charges to a molecule.

    Updates the graph object.

    Args:
        graph: The molecule to add charges to.
        total_charge: The total charge to solve for.
    """
    global _charger
    _charger.charge(graph, total_charge, True, False, range(20, -1, -1))
