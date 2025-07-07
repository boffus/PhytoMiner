import pandas as pd
import os
import time
from .api import initial_fetch, subsequent_fetch
from .processing import process_homolog_data, merge_homolog_and_tsv_data
from .utils import print_summary, pivotmap, log_message
from .config import (
    DEFAULT_MAX_WORKERS,
    HOMOLOGS_OUTPUT_FILE,
    TSV_DIR,
    JOIN2_OUTPUT_FILE,
)
from .data import read_all_tsv_files

def run_homologs_pipeline(
    initial_organism,
    initial_genes_dict,
    subsequent_organisms,
    max_workers=DEFAULT_MAX_WORKERS,
    checkpoint_dir="homolog_checkpoints"
):
    """
    Orchestrates a multi-step homolog search pipeline.

    This function performs an initial fetch from a source organism, then iteratively
    fetches homologs for a list of subsequent organisms, using checkpoints to
    save and resume progress.

    Parameters:
        initial_organism (str): The short name of the initial organism (e.g., "athaliana").
        initial_genes_dict (dict): A dictionary mapping transcript IDs to subunit names for the initial fetch.
        subsequent_organisms (list): A list of organism short names to process iteratively.
        max_workers (int): The maximum number of parallel threads for fetching.
        checkpoint_dir (str): The directory to save/load intermediate results. If None, checkpointing is disabled.

    Returns:
        pd.DataFrame: The final, processed master DataFrame containing all homolog data.
    """
    homolog_df = pd.DataFrame()

    if checkpoint_dir:
        os.makedirs(checkpoint_dir, exist_ok=True)
        initial_checkpoint_file = os.path.join(checkpoint_dir, "step1_initial_fetch.csv")

        if os.path.exists(initial_checkpoint_file):
            print(f"Loading initial data from checkpoint: {initial_checkpoint_file}")
            homolog_df = pd.read_csv(initial_checkpoint_file)
            print_summary(homolog_df, "Loaded Initial Data Summary")
        else:
            print(f"--- Starting Fetch for {initial_organism} ---")
            homolog_df = initial_fetch(
                source_organism_name=initial_organism,
                transcript_names=list(initial_genes_dict.keys()),
                subunit_dict=initial_genes_dict,
                max_workers=DEFAULT_MAX_WORKERS
            )
            time.sleep(1)

            if not homolog_df.empty:
                homolog_df = process_homolog_data(homolog_df)
                homolog_df.to_csv(initial_checkpoint_file, index=False)
                print(f"Saved initial data to checkpoint: {initial_checkpoint_file}")
                print_summary(homolog_df, "Initial Fetch and Processing Complete")
            else:
                print(f"Initial fetch for {initial_organism} yielded no results. Exiting.")
                return pd.DataFrame()
    else: # No checkpointing
        print(f"--- Starting Fetch for {initial_organism} (no checkpoints) ---")
        homolog_df = initial_fetch(initial_organism, list(initial_genes_dict.keys()), initial_genes_dict, max_workers)
        if not homolog_df.empty:
            homolog_df = process_homolog_data(homolog_df)
        else:
            return pd.DataFrame()

    # Process subsequent organisms
    for organism_name in subsequent_organisms:
        new_homologs_df = pd.DataFrame()
        if checkpoint_dir:
            safe_org_name = "".join(c if c.isalnum() else "_" for c in organism_name)
            organism_checkpoint_file = os.path.join(checkpoint_dir, f"{safe_org_name}.csv")

            if os.path.exists(organism_checkpoint_file):
                print(f"\nLoading data for {organism_name} from checkpoint: {organism_checkpoint_file}")
                new_homologs_df = pd.read_csv(organism_checkpoint_file)
            else:
                print(f"\n--- Fetching homologs for subsequent organism: {organism_name} ---")
                new_homologs_df = subsequent_fetch(homolog_df, organism_name, max_workers)
                if not new_homologs_df.empty:
                    new_homologs_df = process_homolog_data(new_homologs_df)
                    new_homologs_df.to_csv(organism_checkpoint_file, index=False)
        else: # No checkpointing
            new_homologs_df = subsequent_fetch(homolog_df, organism_name, max_workers)
            if not new_homologs_df.empty:
                new_homologs_df = process_homolog_data(new_homologs_df)

        if not new_homologs_df.empty:
            homolog_df = pd.concat([homolog_df, new_homologs_df], ignore_index=True)
            homolog_df = process_homolog_data(homolog_df)
            print_summary(homolog_df, f"Updated Master DataFrame after adding {organism_name}")

    pivotmap(homolog_df, index='organism.shortName', columns='subunit1', values='primaryIdentifier')

    print("\n--- Homolog Pipeline Finished ---")
    return homolog_df

def run_workflow2():
    """
    Runs the complete data processing workflow.
    """
    log_message("Reading homolog data...")
    homolog_df = pd.read_csv(HOMOLOGS_OUTPUT_FILE)

    log_message("Reading TSV files...")
    tsv_df = read_all_tsv_files(TSV_DIR)

    log_message("Merging data...")
    merged_df = merge_homolog_and_tsv_data(homolog_df, tsv_df)

    log_message(f"Saving merged data to {JOIN2_OUTPUT_FILE}...")
    merged_df.to_csv(JOIN2_OUTPUT_FILE, index=False)
    log_message("Workflow completed successfully.")