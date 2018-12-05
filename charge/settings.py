

import os

NAUTY_EXC = None
"""Location of the dreadnaut executable."""

if 'NAUTY_EXC' in os.environ:
    fpath = os.path.join(os.environ['NAUTY_EXC'], 'dreadnaut')
    if os.path.isfile(fpath) and os.access(fpath, os.X_OK):
        NAUTY_EXC = fpath

if not NAUTY_EXC:
    for path in os.environ["PATH"].split(os.pathsep):
        path = path.strip('"')
        fpath = os.path.join(path, 'dreadnaut')
        if os.path.isfile(fpath) and os.access(fpath, os.X_OK):
            NAUTY_EXC = fpath
            break
    else:
        raise Exception('Could not find nauty executable.')

ILP_SOLVER_MAX_SECONDS = 60
"""Time limit for the ILP solver in seconds."""

DEFAULT_TOTAL_CHARGE = 0
"""Default target total charge."""

DEFAULT_TOTAL_CHARGE_DIFF = 0.01
"""Default allowed deviation from target total charge."""

MAX_BINS = 25
"""Maximal number of bins for the histogram calculations."""

ROUNDING_DIGITS = 3
"""Default number of digits after the decimal point for the assigned charges."""

MAX_ROUNDING_DIGITS = 9
"""Maximal number of digits after the decimal point for the assigned charges."""

REPO_LOCATION=None
"""Default repository location."""

if 'REPO_LOCATION' in os.environ:
    fpath = os.environ['REPO_LOCATION']
    if os.path.isfile(fpath) and os.access(fpath, os.R_OK):
        REPO_LOCATION = fpath

IACM_MAP = {
        'O': 'O', 'OM': 'O', 'OA': 'O', 'OE': 'O', 'OW': 'O', 'N': 'N',
        'NT': 'N', 'NL': 'N', 'NR': 'N', 'NZ': 'N', 'NE': 'N',
        'C': 'C', 'CH0': 'C', 'CH1': 'C', 'CH2': 'C', 'CH3': 'C', 'CH4': 'C', 'CH2r': 'C', 'CR1': 'C', 'CC14': 'C', 'HC': 'H', 'H': 'H',
        'DUM': None, 'S': 'S', 'CU1+': 'Cu', 'CU2+': 'Cu', 'FE': 'Fe', 'ZN2+': 'Zn', 'MG2+': 'Mg', 'CA2+': 'Ca', 'P,SI': 'P',
        'AR': 'Ar', 'F': 'F', 'CL': 'Cl', 'BR': 'Br', 'CMet': 'C', 'OMet': 'O', 'NA+': 'Na', 'CL-': 'Cl', 'CChl': 'C',
        'CLChl': 'Cl', 'HChl': 'H', 'SDmso': 'S', 'CDmso': 'C', 'ODmso': 'O', 'CCl4': 'C', 'CLCl4': 'Cl', 'FTFE': 'F',
        'CTFE': 'C', 'CHTFE': 'C', 'OTFE': 'O', 'CUrea': 'C', 'OUrea': 'O', 'NUrea': 'N', 'CH3p': 'C', 'I': 'I', 'CLOpt': 'Cl',
        'B': 'B', 'SE': 'Se', 'HS14': 'H', 'CLAro': 'Cl', 'BROpt': 'Br', 'OEOpt': 'O', 'NOpt': 'N', 'CAro': 'C', 'CPos': 'C',
        'NPri': 'N', 'NTer': 'N', 'OAlc': 'O', 'P': 'P', 'SI': 'Si',
    }

IACM_ELEMENTS = [
    'O', 'OM', 'OA', 'OE', 'OW',
    'N', 'NT', 'NL', 'NR', 'NZ', 'NE',
    'C', 'CH0', 'CH1', 'CH2', 'CH3', 'CH4', 'CH2r', 'CR1',
    'HC', 'H',
    'DUM',
    'S',
    'CU1+', 'CU2+',
    'FE', 'ZN2+', 'MG2+', 'CA2+', 'P,SI', 'AR', 'F', 'CL', 'BR',
    'CMet', 'OMet',
    'NA+', 'CL-',
    'CChl', 'CLChl', 'HChl',
    'SDmso', 'CDmso', 'ODmso',
    'CCl4', 'CLCl4',
    'FTFE', 'CTFE', 'CHTFE', 'OTFE',
    'CUrea', 'OUrea', 'NUrea',
    'CH3p', 'I', 'CLOpt', 'B', 'SE', 'HS14', 'CLAro',
    'BROpt', 'OEOpt', 'NOpt', 'CAro', 'CPos',
    'NPri', 'NTer',
    'OAlc',
    'SI', 'P']
"""GROMOS IACM elements"""
