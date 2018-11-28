import connexion
import six

from charge_server import util


def charge_molecule(body):  # noqa: E501
    """Submit a molecule for charging

    Accepts input and produces output in Lemon Graph Format, for which
    there is no MIME type, so this specifies text/plain (which it is).
    See http://lemon.cs.elte.hu/pub/doc/1.2.3/a00002.html.

    :param body: Description of the input molecule
    :type body: str

    :rtype: None
    """
    return 'do some magic!'
