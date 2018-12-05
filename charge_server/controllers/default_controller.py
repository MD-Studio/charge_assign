import connexion
import six

from charge_server import util
from charge_server import charge_server

from charge.babel import convert_from, convert_to, IOType


def charge_molecule(molecule: bytes, total_charge: int) -> str:
    """Submit a molecule for charging

    Accepts input and produces output in Lemon Graph Format, for which
    there is no MIME type, so this specifies text/plain (which it is).
    See http://lemon.cs.elte.hu/pub/doc/1.2.3/a00002.html.

    Args:
        molecule: Description of the input molecule.
        total_charge: Desired total charge of the molecule.
    """
    graph = convert_from(molecule.decode('utf-8'), IOType.LGF)
    charge_server.charge(graph, total_charge)
    lgf_output = convert_to(graph, IOType.LGF)
    return lgf_output.encode()
