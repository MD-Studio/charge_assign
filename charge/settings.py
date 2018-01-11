

import os

for path in os.environ["PATH"].split(os.pathsep):
    path = path.strip('"')
    fpath = os.path.join(path, 'dreadnaut')
    if os.path.isfile(fpath) and os.access(fpath, os.X_OK):
        NAUTY_EXC = fpath
        break
else:
    NAUTY_EXC = os.path.expanduser('~/workspace/nauty26r7/dreadnaut')

REPO_LOCATION='atb.zip'

IACM_MAP = dict([('O', 'O'), ('OM', 'O'), ('OA', 'O'), ('OE', 'O'), ('OW', 'O'), ('OMet', 'O'), ('ODmso', 'O'),
                 ('OTFE', 'O'), ('OUrea', 'O'),
                 ('N', 'N'), ('NT', 'N'), ('NL', 'N'), ('NR', 'N'), ('NZ', 'N'), ('NE', 'N'), ('NUrea', 'N'),
                 ('C', 'C'), ('CH0', 'C'), ('CH1', 'C'), ('CH2', 'C'), ('CH3', 'C'), ('CH4', 'C'), ('CH2r', 'C'),
                 ('CR1', 'C'), ('CMet', 'C'), ('CChl', 'C'), ('CDmso', 'C'), ('CCl4', 'C'), ('CTFE', 'C'),
                 ('CHTFE', 'C'), ('CUrea', 'C'), ('CH3p', 'C'),
                 ('HC', 'H'), ('H', 'H'), ('HChl', 'H'), ('HS14', 'H'),
                 ('S', 'S'), ('SDmso', 'S'),
                 ('CU1+', 'Cu'), ('CU2+', 'Cu'),
                 ('F', 'F'), ('FTFE', 'F'),
                 ('CL', 'Cl'), ('CL-', 'Cl'), ('CLChl', 'Cl'), ('CLCl4', 'Cl'), ('CLOpt', 'Cl'), ('CLAro', 'Cl'),
                 ('BR', 'Br'), ('BROpt', 'Br'),
                 ('NA+', 'Na'), ('I', 'I'), ('B', 'B'), ('SE', 'Se'), ('FE', 'Fe'), ('ZN2+', 'Zn'), ('MG2+', 'Mg'),
                 ('CA2+', 'Ca'), ('P,SI', 'P'), ('AR', 'Si')])

IACM_ELEMENTS = ['O', 'OM', 'OA',  'OE',  'OW',  'N', 'NT',
                 'NL', 'NR', 'NZ', 'NE', 'C', 'CH0',
                 'CH1', 'CH2', 'CH3', 'CH4', 'CH2r',
                 'CR1', 'HC', 'H', 'DUM', 'S', 'CU1+',
                 'CU2+', 'FE', 'ZN2+', 'MG2+', 'CA2+',
                 'P,SI', 'AR', 'F', 'CL', 'BR', 'CMet',
                 'OMet', 'NA+', 'CL-', 'CChl', 'CLChl', 'HChl',
                 'SDmso', 'CDmso', 'ODmso', 'CCl4', 'CLCl4',
                 'FTFE', 'CTFE', 'CHTFE', 'OTFE', 'CUrea',
                 'OUrea', 'NUrea', 'CH3p', 'I', 'CLOpt', 'B',
                 'SE', 'HS14', 'CLAro', 'BROpt']

DB_SCHEMA = {
    "type": "object",
    "additionalProperties": {
        "type": "object",
        "properties": {
            "core_ids": { "type": "array", "items": { "type": "number" } },
            "atoms": { "type": "object", "additionalProperties": { "type": "string" } },
            "bonds": { "type": "array", "items": { "type": "array", "items": { "type": "number" } } },
            "charges": { "type": "array", "items": { "type": "number" } },
            "uncertainties": { "type": "array", "items": { "type": "number" } }
        }
    }
}

ATOMIC_SCHEMA = {
    "type": "object",
    "properties": {
        "core_atoms": { "type": "array", "items": { "type": "number" } },
        "atom_mappings": {
           "type": "object",
            "additionalProperties": {
                "type": "array",
                "items": { "type": "object", "additionalProperties": { "type": "number" } }
            }
        }
    }
}
