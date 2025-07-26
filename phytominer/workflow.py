import logging
import os
from pathlib import Path

import pandas as pd

from .config import (
    DEFAULT_CHUNK,
    DEFAULT_WORKERS,
    GENE_CHECKPOINT_DIR,
    GENE_OUTPUT_FILE,
    HOMOLOGS_CHECKPOINT_DIR,
    HOMOLOGS_OUTPUT_FILE,
    STEP2_OUTPUT_FILE,
    TSV_DIR
)
from .data import read_tsv_files
from .processing import (
    initial_fetch,
    genes_fetch,
    process_homolog_data,
    step2_merge,
    subsequent_fetch,
    load_master_df
)
from .utils import log_summary, pivotmap

logger = logging.getLogger(__name__)

__all__ = ["homologs_pipe", "join_tsvs", "genes_pipe"]

def homologs_pipe(
    initial_organism: str,
    initial_genes_dict: dict,
    subsequent_organisms: list,
    max_workers: int = DEFAULT_WORKERS,
    chunk_size: int = DEFAULT_CHUNK,
    checkpoint_dir: str = HOMOLOGS_CHECKPOINT_DIR,
    output_file: str = HOMOLOGS_OUTPUT_FILE,
) -> pd.DataFrame:
    """
    Orchestrates a multi-step homolog search pipeline with checkpointing.
    """
    os.makedirs(checkpoint_dir, exist_ok=True)
    all_homolog_dfs = []

    # Step 1: Initial Fetch
    initial_checkpoint_file = Path(checkpoint_dir) / "step1_initial_fetch.csv"
    if initial_checkpoint_file.exists():
        logger.info(f"Loading initial data from checkpoint: {initial_checkpoint_file}")
        homolog_df = pd.read_csv(initial_checkpoint_file)
        log_summary(homolog_df, "Loaded Initial Data Summary")
    else:
        logger.info(f"--- Starting Initial Fetch for {initial_organism} ---")
        homolog_df = initial_fetch(
            source_organism_name=initial_organism,
            transcript_names=list(initial_genes_dict.keys()),
            subunit_dict=initial_genes_dict,
            max_workers=max_workers,
            chunk_size=chunk_size,
        )
        if not homolog_df.empty:
            homolog_df = process_homolog_data(homolog_df)
            homolog_df.to_csv(initial_checkpoint_file, index=False)
            logger.info(f"Saved initial data to checkpoint: {initial_checkpoint_file}")
            log_summary(homolog_df, "Initial Fetch and Processing Complete")
    all_homolog_dfs.append(homolog_df)

    # Step 2: Subsequent Fetches
    for organism_name in subsequent_organisms:
        safe_org_name = "".join(c if c.isalnum() else "_" for c in organism_name)
        organism_checkpoint_file = Path(checkpoint_dir) / f"{safe_org_name}.csv"

        if organism_checkpoint_file.exists():
            logger.info(f"Loading data for {organism_name} from checkpoint: {organism_checkpoint_file}")
            new_homologs_df = pd.read_csv(organism_checkpoint_file)
        else:
            logger.info(f"--- Fetching homologs for subsequent organism: {organism_name} ---")
            current_master_df = pd.concat(all_homolog_dfs, ignore_index=True)
            new_homologs_df = subsequent_fetch(
                current_master_df, organism_name, max_workers, chunk_size
            )
            if not new_homologs_df.empty:
                new_homologs_df = process_homolog_data(new_homologs_df)
                new_homologs_df.to_csv(organism_checkpoint_file, index=False)
                logger.info(f"Saved checkpoint for {organism_name}")

        if not new_homologs_df.empty:
            all_homolog_dfs.append(new_homologs_df)

    # Step 3: Final Combination and Processing
    master_df = pd.concat(all_homolog_dfs, ignore_index=True)
    master_df = process_homolog_data(master_df) 

    log_summary(master_df, "Final DataFrame")
    pivotmap(master_df)
    master_df.to_csv(output_file, index=False)
    logger.info(f"--- Homolog Pipeline Finished. Output saved to {output_file} ---")
    return master_df

def join_tsvs(
    homolog_file: str = HOMOLOGS_OUTPUT_FILE,
    tsv_dir: str = TSV_DIR,
    output_file: str = STEP2_OUTPUT_FILE,
) -> pd.DataFrame:
    """
    Reads homolog data and TSV files, merges them, and saves the result.
    """
    logger.info("--- Starting Step 2: Merging with TSV data ---")
    try:
        homolog_df = pd.read_csv(homolog_file)
        logger.info(f"Loaded homolog data from {homolog_file}")
    except FileNotFoundError:
        logger.error(f"Homolog file not found: {homolog_file}. Aborting merge.")
        return pd.DataFrame()

    tsv_df = read_tsv_files(tsv_dir)
    if tsv_df.empty:
        logger.warning("TSV data is empty. Merge will not add any new columns.")

    merged_df = step2_merge(homolog_df, tsv_df)
    log_summary(merged_df, "Merged DataFrame Summary")
    merged_df.to_csv(output_file, index=False)
    logger.info(f"--- Merge complete. Output saved to {output_file} ---")
    return merged_df

def genes_pipe(
    master_file: str = STEP2_OUTPUT_FILE,
    checkpoint_dir: str = GENE_CHECKPOINT_DIR,
    output_file: str = GENE_OUTPUT_FILE,
    max_workers: int = DEFAULT_WORKERS,
    chunk_size: int = DEFAULT_CHUNK
) -> pd.DataFrame:
    """
    Fetches detailed gene data for homologs found in the master file.
    """
    logger.info("--- Starting Step 3: Fetching detailed gene data ---")
    master_df = load_master_df(master_file)
    if master_df is None or master_df.empty:
        logger.error(f"Master file at {master_file} could not be loaded or is empty. Exiting.")
        return pd.DataFrame()

    os.makedirs(checkpoint_dir, exist_ok=True)
    all_gene_data_parts = []

    # Group by the initial subunit classification
    for subunit, sub_df in master_df.groupby('subunit1'):
        safe_subunit_name = "".join(c if c.isalnum() else "_" for c in str(subunit))
        checkpoint_file = Path(checkpoint_dir) / f"genes_{safe_subunit_name}.csv"

        if checkpoint_file.exists():
            logger.info(f"Loading checkpoint for subunit '{subunit}'")
            df_chunk = pd.read_csv(checkpoint_file)
            all_gene_data_parts.append(df_chunk)
            continue

        gene_ids = sub_df['primaryIdentifier'].dropna().unique().tolist()
        if not gene_ids:
            logger.info(f"No unique gene IDs for subunit '{subunit}', skipping.")
            continue

        logger.info(f"Fetching {len(gene_ids)} genes for subunit '{subunit}'...")
        df_chunk = genes_fetch(gene_ids, max_workers, chunk_size)

        if not df_chunk.empty:
            df_chunk['subunit1'] = subunit
            df_chunk.to_csv(checkpoint_file, index=False)
            logger.info(f"Saved checkpoint for subunit '{subunit}' ({len(df_chunk)} rows)")
            all_gene_data_parts.append(df_chunk)

    if not all_gene_data_parts:
        logger.warning("No gene data was fetched across all subunits.")
        return pd.DataFrame()

    genes_df = pd.concat(all_gene_data_parts, ignore_index=True)
    log_summary(genes_df, "Final Gene Data Summary")

    genes_df.to_csv(output_file, index=False)
    logger.info(f"--- Gene data fetching finished. Output saved to {output_file} ---")
    return genes_df