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
    DEFAULT_SLEEP_SECONDS
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

    Performs an initial fetch from source organism, then iteratively
    fetches homologs from a list of subsequent organisms using checkpoints to
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
            time.sleep(DEFAULT_SLEEP_SECONDS)
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

def run_expressions_workflow(
    master_file,
    checkpoint_dir,
    output_file,
    fetch_expression_data_for_gene_chunk,
    load_master_df,
    chunk_size=24,
    max_workers=8
):
    """
    Orchestrates the fetching of expression data for all subunits and merges with master homolog data.
    """
    import os
    import pandas as pd
    import concurrent.futures
    import time

    master_df = load_master_df(master_file)
    if master_df is None:
        print("Master DataFrame could not be loaded. Exiting workflow.")
        return

    os.makedirs(checkpoint_dir, exist_ok=True)
    all_expression_results = []
    if 'Subunit' in master_df.columns:
        master_df['Subunit'] = master_df['Subunit'].astype(str).str.upper()
        unique_subunits = master_df['Subunit'].dropna().unique()
        print(f"\nFound {len(unique_subunits)} unique subunits to process.")

        for subunit_val in unique_subunits:
            print(f"\n--- Preparing tasks for Subunit Group: {subunit_val} ---")
            subunit_specific_df = master_df[master_df['Subunit'] == subunit_val].copy()
            if subunit_specific_df.empty:
                print(f"  No data for subunit {subunit_val}. Skipping.")
                continue

            safe_subunit_filename = "".join(c if c.isalnum() else "_" for c in str(subunit_val))
            subunit_checkpoint_file = os.path.join(checkpoint_dir, f"expression_data_{safe_subunit_filename}.csv")

            if os.path.exists(subunit_checkpoint_file):
                print(f"  Checkpoint found for subunit {subunit_val}. Loading...")
                try:
                    expression_data_for_this_subunit = pd.read_csv(subunit_checkpoint_file)
                    all_expression_results.append(expression_data_for_this_subunit)
                    print(f"  Loaded {len(expression_data_for_this_subunit)} records from checkpoint for {subunit_val}.")
                    continue 
                except Exception as e:
                    print(f"  Error loading checkpoint for {subunit_val}: {e}. Will fetch.")
                    if os.path.exists(subunit_checkpoint_file): os.remove(subunit_checkpoint_file)
            
            gene_ids_for_subunit = subunit_specific_df['Homolog_Gene'].dropna().unique().tolist()
            if not gene_ids_for_subunit:
                print(f"  No unique 'Homolog_Gene' values for subunit {subunit_val}. Skipping.")
                continue

            gene_id_chunks = [gene_ids_for_subunit[i:i + chunk_size] for i in range(0, len(gene_ids_for_subunit), chunk_size)]
            print(f"  Split {len(gene_ids_for_subunit)} genes for {subunit_val} into {len(gene_id_chunks)} chunks.")

            subunit_chunk_results = []
            tasks = [{'chunk': chunk, 'name': subunit_val, 'chunk_num': i} for i, chunk in enumerate(gene_id_chunks)]

            if tasks:
                print(f"  Submitting {len(tasks)} chunks for {subunit_val} for parallel fetching...")
                with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
                    future_to_task = {
                        executor.submit(fetch_expression_data_for_gene_chunk, task['chunk'], task['name'], task['chunk_num'], len(gene_id_chunks)): task
                        for task in tasks
                    }
                    print(f"    All {len(future_to_task)} futures submitted for {subunit_val}.")
                    for future in concurrent.futures.as_completed(future_to_task):
                        task_info = future_to_task[future]
                        print(f"      Processing future for chunk {task_info['chunk_num'] + 1} of {subunit_val}...")
                        try:
                            result_df_chunk = future.result()
                            if not result_df_chunk.empty:
                                subunit_chunk_results.append(result_df_chunk)
                                print(f"        Appended result from chunk {task_info['chunk_num'] + 1} of {subunit_val}. Length of subunit_chunk_results: {len(subunit_chunk_results)}")
                            else:
                                print(f"        Result from chunk {task_info['chunk_num'] + 1} of {subunit_val} was empty.")
                        except Exception as exc:
                            print(f"    Chunk {task_info['chunk_num'] + 1} for {task_info['name']} generated an exception: {exc}")
                        time.sleep(0.1)
                    print(f"    Finished processing all futures for {subunit_val}.")
            
            if subunit_chunk_results:
                print(f"  Attempting to concatenate {len(subunit_chunk_results)} chunk results for {subunit_val}...")
                combined_data_for_subunit = pd.concat(subunit_chunk_results, ignore_index=True)
                print(f"    Concatenation successful for {subunit_val}. Shape: {combined_data_for_subunit.shape}")
                combined_data_for_subunit.to_csv(subunit_checkpoint_file, index=False)
                all_expression_results.append(combined_data_for_subunit)
                print(f"  All chunks for {subunit_val} processed. Checkpoint SAVED.")
            else:
                print(f"  No data fetched for any chunk of {subunit_val}.")
            print(f"--- Finished processing for Subunit Group: {subunit_val} ---")
    else:
        print("\n'Subunit' column not found. Cannot fetch by subunit.")

    if all_expression_results:
        fetched_expression_df = pd.concat(all_expression_results, ignore_index=True)
        print("\n--- Combined All Fetched Subunit Expression Data ---")
        print(f"Shape: {fetched_expression_df.shape}")

        master_df_for_join = master_df.copy()
        if 'Subunit' in master_df_for_join.columns:
            master_df_for_join['Subunit'] = master_df_for_join['Subunit'].astype(str).str.upper()
        
        print(f"Merging original homolog data into fetched expression data...")
        final_intermediate_df = pd.merge(
            fetched_expression_df,
            master_df_for_join,
            left_on=['Input_PrimaryID', 'Gene.organism.shortName', 'Step2_Subunit'],
            right_on=['Homolog_Gene', 'Organism_ShortName', 'Subunit'],
            how='left',
            suffixes=('', '_homologs_master')
        )
        print(f"  Shape after merging with homolog data: {final_intermediate_df.shape}")

        try:
            final_intermediate_df.to_csv(output_file, index=False)
            print(f"\nFetched expression data with homolog context saved to: {output_file}")
        except Exception as e:
            print(f"\nError saving intermediate fetched expression data: {e}")
    else:
        print("\nNo expression data was fetched for any subunit group.")

    print("\nFetching script finished.")