# Subcomplex_dict
SUBCOMPLEX_DICT = {
    'M': ['NdhA', 'NdhB', 'NdhC', 'NdhD', 'NdhE', 'NdhF', 'NdhG'],
    'A': ['NdhH', 'NdhI', 'NdhJ', 'NdhK', 'NdhL', 'NdhM', 'NdhN', 'NdhO'],
    'EDB': ['NdhT', 'NdhS', 'NdhU', 'NdhV'],
    'B': ['PnsB1', 'PnsB2', 'PnsB3', 'PnsB4', 'PnsB5'],
    'L': ['PnsL1', 'PnsL2', 'PnsL3', 'PnsL4', 'PnsL5'],
    'CRR': ['CRR1', 'CRR2', 'CRR21', 'CRR27', 'CRR3', 'CRR4', 'CRR41', 'CRR42', 'CRR6', 'CRR7'],
    'FKBP': ['FKBP12', 'FKBP13', 'FKBP15-1', 'FKBP15-2', 'FKBP15-3', 'FKBP16-3', 'FKBP20-1'],
    'LHCA': ['LHCA1', 'LHCA2', 'LHCA3', 'LHCA4', 'LHCA5', 'LHCA6'],
    'PPD': ['PPD1', 'PPD2', 'PPD3', 'PPD4', 'PPD5', 'PPD6', 'PPD7', 'PPD8'],
    'PSB': ['PSBP-1', 'PSBP-2', 'PSBQ-1', 'PSBQ-2', 'PPL1', 'PQL3'],
    'PGR': ['PGR5', 'PGRL1A', 'PGRL1B', 'HCF101', 'HCF136', 'NDF5']
}

# Create reverse lookup dictionary for subcomplexes
SUBUNIT_TO_SUBCOMPLEX = {}
for subcomplex, subunits in SUBCOMPLEX_DICT.items():
    for subunit in subunits:
        SUBUNIT_TO_SUBCOMPLEX[subunit] = subcomplex
