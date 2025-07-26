import pandas as pd
from pathlib import Path
import logging
logger = logging.getLogger(__name__)
from .config import TSV_DIR

# --- Subcomplex and Subunit Definitions ---
SUBCOMPLEX_DICT = {
    'M': ['NDHA', 'NDHB', 'NDHC', 'NDHD', 'NDHE', 'NDHF', 'NDHG'],
    'A': ['NDHH', 'NDHI', 'NDHJ', 'NDHK', 'NDHL', 'NDHM', 'NDHN'],
    'EDB': ['NDHT', 'NDHS', 'NDHU', 'NDHV'],
    'B': ['PNSB1', 'PNSB2', 'PNSB3', 'PNSB4', 'PNSB5'],
    'L': ['PNSL1', 'PNSL2', 'PNSL3', 'PNSL4', 'PNSL5'],
    'CRR': ['CRR1', 'CRR2', 'CRR21', 'CRR27', 'CRR3', 'CRR4', 'CRR41', 'CRR42', 'CRR6', 'CRR7'],
    'FKBP': ['FKBP12', 'FKBP13', 'FKBP15-1', 'FKBP15-2', 'FKBP15-3', 'FKBP16-3', 'FKBP20-1'],
    'LHCA': ['LHCA1', 'LHCA2', 'LHCA3', 'LHCA4', 'LHCA5', 'LHCA6'],
    'PPD': ['PPD1', 'PPD2', 'PPD3', 'PPD4', 'PPD5', 'PPD6', 'PPD7', 'PPD8'],
    'PSB': ['PSBP-1', 'PSBP-2', 'PSBQ-1', 'PSBQ-2', 'PPL1', 'PQL3'],
    'PGR': ['PGR5', 'PGRL1A', 'PGRL1B', 'HCF101', 'HCF136', 'NDF5']
}

# Reverse lookup dictionary for finding a subcomplex by its subunit.
SUBUNIT_TO_SUBCOMPLEX = {
    subunit: subcomplex
    for subcomplex, subunits in SUBCOMPLEX_DICT.items()
    for subunit in subunits
}

def read_tsv_files(tsv_dir_path: str) -> pd.DataFrame:
    """
    Reads 'ndh*.tsv' files from a directory, extracts gene IDs, and assigns
    a subunit name based on the filename.

    This function dynamically finds all files matching the pattern, making it
    robust to the addition or removal of TSV files. It handles files with
    and without headers, standardizes column names, and cleans the data.

    Args:
        tsv_dir_path: The absolute path to the directory containing the TSV files.

    Returns:
        A pandas DataFrame containing the aggregated data from all valid TSV files.
        The DataFrame will have 'Gene_ID_from_TSV' and 'subunit2' columns.
        Returns an empty DataFrame if the directory does not exist or no valid data can be loaded.
    """
    all_tsv_data = []
    tsv_dir = Path(tsv_dir_path)

    if not tsv_dir.is_dir():
        logger.error(f"TSV directory not found: {tsv_dir_path}")
        return pd.DataFrame()

    # Use glob to dynamically find all files starting with 'ndh'
    for file_path in tsv_dir.glob("ndh*.tsv"):
        file_basename = file_path.stem 

        try:
            # The 'ndhO.tsv' file has a header, others do not.
            header = 0 if file_basename == "ndhO" else None

            # Read only the first column from the TSV file.
            df = pd.read_csv(
                file_path,
                header=header,
                usecols=[0],
                engine='python',
                skipinitialspace=True
            )

            if df.empty:
                logger.info(f"TSV file is empty or has no data rows: {file_path}. Skipping.")
                continue

            # Standardize the column name and clean the data.
            df.columns = ['Gene_ID_from_TSV']
            df.dropna(subset=['Gene_ID_from_TSV'], inplace=True)
            df['Gene_ID_from_TSV'] = df['Gene_ID_from_TSV'].astype(str).str.strip()

            # If the DataFrame is empty after cleaning, skip it.
            if df.empty:
                logger.info(f"TSV file is empty after cleaning NA values: {file_path}. Skipping.")
                continue

            df['subunit2'] = file_basename.upper()
            all_tsv_data.append(df)

        except pd.errors.EmptyDataError:
            logger.warning(f"TSV file is empty and could not be read: {file_path}. Skipping.")
        except Exception as e:
            logger.error(f"Failed to read or process TSV file {file_path}: {e}")

    if not all_tsv_data:
        logger.warning("No valid TSV data was loaded from the directory.")
        return pd.DataFrame()

    return pd.concat(all_tsv_data, ignore_index=True)