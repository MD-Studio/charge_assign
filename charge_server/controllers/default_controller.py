import connexion
import six

from charge_server import util
from charge_server import charge_server

from charge.babel import convert_from, convert_to, IOType
from charge.util import AssignmentError


def charge_molecule(molecule: bytes, total_charge: int) -> str:
    """Submit a molecule for charging

    Accepts input and produces output in Lemon Graph Format, for which
    there is no MIME type, so this specifies text/plain (which it is).
    See http://lemon.cs.elte.hu/pub/doc/1.2.3/a00002.html.

    Args:
        molecule: Description of the input molecule.
        total_charge: Desired total charge of the molecule.
    """
    try:
        graph = convert_from(molecule.decode('utf-8'), IOType.LGF)
    except (ValueError, AttributeError):
        return ('Error decoding input, is it valid LGF, and sent as'
                ' "Content-Type: text/plain" ?'), 400

    try:
        charge_server.charge(graph, total_charge)
    except AssignmentError:
        return ('Charges could not be assigned due to lack of data or because'
                ' the total charge was too far off from our reference charges.'
                ), 404

    lgf_output = convert_to(graph, IOType.LGF)
    return lgf_output.encode()
