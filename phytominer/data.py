import pandas as pd
import os
import string
from phytominer.config import TSV_DIR

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

# Create reverse lookup dictionary for sub-complexes
SUBUNIT_TO_SUBCOMPLEX = {}
for subcomplex, subunits in SUBCOMPLEX_DICT.items():
    for subunit in subunits:
        SUBUNIT_TO_SUBCOMPLEX[subunit] = subcomplex

def read_all_tsv_files(tsv_dir_path: str) -> pd.DataFrame:
    """
    Reads TSV files (ndhA-O.tsv), extracts gene IDs, and assigns subunit from filename.
    """
    all_tsv_data = []
    letters = string.ascii_uppercase[:15]  # A to O

    for letter in letters:
        file_basename = f"ndh{letter}"
        file_path = os.path.join(tsv_dir_path, f"{file_basename}.tsv")

        if not os.path.exists(file_path):
            print(f"Warning: TSV file not found: {file_path}. Skipping.")
            continue

        try:
            if file_basename == "ndhO":  # ndhO.tsv has a header
                temp_df = pd.read_csv(file_path, header=0, usecols=[0], engine='python', skipinitialspace=True)
                if not temp_df.empty:
                    current_df = pd.DataFrame({'Gene_ID_from_TSV': temp_df.iloc[:, 0]})
            else:  # Other files have no header
                current_df = pd.read_csv(file_path, header=None, names=['Gene_ID_from_TSV'], engine='python', skipinitialspace=True)

            if current_df is not None and not current_df.empty:
                current_df.dropna(subset=['Gene_ID_from_TSV'], inplace=True)
                if not current_df.empty:
                    current_df['Validated_Subunit_from_file'] = file_basename
                    current_df['Gene_ID_from_TSV'] = current_df['Gene_ID_from_TSV'].astype(str).str.strip()
                    all_tsv_data.append(current_df)
                else:
                    print(f"Info: TSV file {file_path} resulted in an empty DataFrame after dropping NA gene IDs. Skipping.")
            elif current_df is None:
                print(f"Warning: TSV file {file_path} (ndhO) might be empty or header issue. Skipping.")
            else:
                print(f"Info: TSV file {file_path} is empty. Skipping.")

        except pd.errors.EmptyDataError:
            print(f"Warning: TSV file {file_path} is empty. Skipping.")
        except Exception as e:
            print(f"Error reading TSV file {file_path}: {e}")

    if all_tsv_data:
        return pd.concat(all_tsv_data, ignore_index=True)
    else:
        return pd.DataFrame()
