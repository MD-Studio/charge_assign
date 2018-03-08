"""
Definitions for bond types.

.. moduleauthor:: Martin S. Engler <martin.engler@cwi.nl>

"""

from enum import Enum

class BondType(Enum):
    """Supported chemical bond types."""

    SINGLE = 'SINGLE'
    DOUBLE = 'DOUBLE'
    TRIPLE = 'TRIPLE'
    AROMATIC = 'AROMATIC'
    UNKNOWN = 'UNKNOWN'
