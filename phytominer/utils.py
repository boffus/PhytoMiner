import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def pivotmap(dataframe, index='organism.shortName', columns='subunit1', values='primaryIdentifier'):
    """
    Creates a pivot table and visualises it with a heatmap.

    Parameters:
        dataframe (pd.DataFrame): The input dataframe.
        index (str): Column to use for the pivot table index.
        columns (str): Column to use for the pivot table columns.
        values (str): Column to aggregate for the pivot table values.
    Returns:
        pd.DataFrame: The pivot table.
    """
    if not all(col in dataframe.columns for col in [index, columns, values]):
        print("Error: DataFrame is missing one or more required columns for pivotmap.")
        return pd.DataFrame()

    pivot_homolog = dataframe.pivot_table(index=index, columns=columns, values=values, aggfunc='count')

    plt.figure(figsize=(15, 10))
    sns.heatmap(pivot_homolog, cmap='viridis', annot=True, fmt='g')
    plt.title(f'Heatmap of {values} Counts by {index} and {columns}')
    plt.show()

    return pivot_homolog

def print_summary(df, stage_message="DataFrame Summary"):
    """
    Prints a concise summary of the DataFrame.
    """
    print(f"\n{stage_message}:")
    print(f"\n  - Shape: {df.shape}")
    print(f"\n  - Info: {df.info()}")
    if 'organism.shortName' in df.columns:
        print(f"  - Unique Homolog Organisms: {df['organism.shortName'].nunique()}")
    if 'Subunit' in df.columns:
        print(f"  - Unique Subunits processed: {df['subunit1'].nunique()}")

def log_message(message: str):
    """
    Logs a message to the console.
    """
    print(f"[INFO] {message}")

