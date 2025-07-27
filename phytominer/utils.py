import logging
logger = logging.getLogger(__name__)
from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from .config import HEATMAP_DIR

def pivotmap(dataframe, index='organism.shortName', columns='subunit1', values='primaryIdentifier'):
    """
    Creates a pivot table and visualises it with a heatmap.

    This function is a convenient wrapper that combines pandas pivot_table
    and seaborn heatmap functionalities for quick analysis.

    Parameters:
        dataframe (pd.DataFrame): The input dataframe.
        index (str): Column to use for the pivot table index.
        columns (str): Column to use for the pivot table columns.
        values (str): Column to aggregate for the pivot table values.
    Returns:
        pd.DataFrame: The generated pivot table. Returns an empty DataFrame
                      if the required columns are not found.
    """
    required_cols = [index, columns, values]
    if not all(col in dataframe.columns for col in required_cols):
        logger.error(f"DataFrame is missing one or more required columns for pivotmap: {required_cols}")
        return pd.DataFrame()

    pivot_homolog = dataframe.pivot_table(index=index, columns=columns, values=values, aggfunc='count')

    try:
        plt.figure(figsize=(15, 10))
        sns.heatmap(pivot_homolog, cmap='viridis', annot=True, fmt='g')
        plt.title(f'Heatmap of {values} Counts by {index} of {columns}')
        
        # Ensure the directory exists before saving
        Path(HEATMAP_DIR).parent.mkdir(parents=True, exist_ok=True)
        
        plt.savefig(HEATMAP_DIR, bbox_inches='tight')
        logger.info(f"Heatmap saved to {HEATMAP_DIR}")
        plt.close()  # Close the figure to free up memory
    except Exception as e:
        logger.error(f"Failed to generate or save heatmap: {e}", exc_info=True)
    return pivot_homolog

def log_summary(df: pd.DataFrame, stage_message: str = "DataFrame Summary"):
    """
    Logs a concise summary of a DataFrame's properties.

    Args:
        df: The DataFrame to summarize.
        stage_message: A message to provide context for the summary.
    """
    if not isinstance(df, pd.DataFrame):
        logger.warning(f"{stage_message}: Input is not a pandas DataFrame.")
        return

    logger.info(f"--- {stage_message} ---")
    logger.info(f"Shape: {df.shape}")

    if df.empty:
        logger.info("DataFrame is empty.")
        return

    logger.info(f"Columns: {df.columns.tolist()}")
    mem_usage = df.memory_usage(deep=True).sum() / (1024 * 1024)  # in MiB
    logger.info(f"Memory Usage: {mem_usage:.2f} MiB")

    if 'organism.shortName' in df.columns:
        logger.info(f"Unique Homolog Organisms: {df['organism.shortName'].nunique()}")
    if 'subunit1' in df.columns:
        logger.info(f"Unique Subunits processed: {df['subunit1'].nunique()}")
    logger.info("--- End of Summary ---")

def log_message(message: str):
    """
    Logs a message to the console.
    """
    logger.info(message)