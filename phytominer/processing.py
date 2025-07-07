import pandas as pd

def process_homolog_data(df_combined):
    """
    Processes a DataFrame of homolog data.
    - Calculates occurrence count for homologs.
    - Deduplicates entries and keeps most relevant homologs.

    Parameters:
        df_combined (pd.DataFrame): DataFrame containing combined homolog data.

    Returns:
        pd.DataFrame: Processed - deduplicated DataFrame.
    """
    if df_combined.empty:
        print("Input DataFrame is empty. No processing done.")
        return df_combined

    print(f"Processing DataFrame with {len(df_combined)} rows...")
    processed_df = df_combined.copy()

    # 1. Categorize Relationship
    relationship_categories = ['one-to-one', 'one-to-many', 'many-to-one', 'many-to-many']
    processed_df['relationship'] = pd.Categorical(
        processed_df['relationship'],
        categories=relationship_categories,
        ordered=True
    )

    # 2. Define key for grouping and calculate homolog occurrences
    origin_key_cols = ['primaryIdentifier', 'organism.shortName']
    processed_df['homolog.occurrences'] = processed_df.groupby(origin_key_cols, observed=False)['source.gene'] \
        .transform('size') # Counts how many source genes point to this homolog

    # 3. Deduplication
    sort_by_cols = ['subunit1', 'relationship', 'homolog.occurrences', 'organism.shortName', 'primaryIdentifier', 'source.organism']
    sort_by_cols = [col for col in sort_by_cols if col in processed_df.columns]

    ascending_map = {'relationship': True, 'homolog.occurrences': False}
    ascending_order = [ascending_map.get(col, True) for col in sort_by_cols]

    dedup_subset_cols = ['subunit1', 'primaryIdentifier', 'organism.shortName']
    dedup_subset_cols = [col for col in dedup_subset_cols if col in processed_df.columns]

    if dedup_subset_cols:
        processed_df = processed_df.sort_values(by=sort_by_cols, ascending=ascending_order)
        processed_df = processed_df.drop_duplicates(subset=dedup_subset_cols, keep='first')

    # 4. Define final column order
    final_columns = [
        'source.organism', 'source.gene', 'subunit1', 'relationship', 'primaryIdentifier',
        'secondaryIdentifier', 'organism.commonName', 'organism.shortName', 'organism.proteomeId',
        'gene.length', 'sequence.length', 'sequence.residues', 'homolog.occurrences'
    ]

    # Reorder columns to a standard format, keeping any extra columns at the end
    existing_final_columns = [col for col in final_columns if col in processed_df.columns]
    other_columns = [col for col in processed_df.columns if col not in existing_final_columns]
    processed_df = processed_df[existing_final_columns + other_columns]

    print(f"Processing complete. DataFrame shape after processing: {processed_df.shape}")
    return processed_df

def merge_homolog_and_tsv_data(homolog_df: pd.DataFrame, tsv_df: pd.DataFrame) -> pd.DataFrame:
    """
    Merges homolog data with TSV data on gene IDs.
    """
    merged_df = pd.merge(homolog_df, tsv_df, left_on='gene_id', right_on='Gene_ID_from_TSV', how='left')
    return merged_df

def load_master_df(filepath):
    """
    Loads the master homolog DataFrame.
    """
    print(f"Loading master homolog data from: {filepath}")
    try:
        df = pd.read_csv(filepath)
        print(f"Successfully loaded master homolog data. Shape: {df.shape}")
        if not all(col in df.columns for col in MASTER_JOIN_COL):
            missing_cols = [col for col in MASTER_JOIN_COL if col not in df.columns]
            print(f"Error: One or more master join columns were not found in {filepath}. Missing: {missing_cols}")
            return None
        if 'Subunit' not in df.columns:
            print(f"Error: 'Subunit' column not found in {filepath}. This column is required for splitting.")
            return None
        return df
    except FileNotFoundError:
        print(f"Error: Master homolog file not found at {filepath}.")
        return None
    except Exception as e:
        print(f"An unexpected error occurred while loading {filepath}: {e}")
        return None
    
def fetch_expression_data(gene_id_chunk, subunit_name_for_context, chunk_num, total_chunks):
    """
    Fetches Phytozome data for a list (chunk) of gene IDs.
    """
    if not gene_id_chunk:
        print(f"  NOTE: Empty gene ID chunk provided for subunit {subunit_name_for_context}. Skipping fetch.")
        return pd.DataFrame()

    all_results = []
    num_genes_to_fetch = len(gene_id_chunk)

    for idx, gene_id in enumerate(gene_id_chunk):
        primary_gene_id = str(gene_id) if pd.notna(gene_id) else None
        if not primary_gene_id:
            print(f"    Skipping missing or invalid gene ID.")
            continue

        print(f"    Fetching gene {idx + 1}/{num_genes_to_fetch} (PrimaryID: {primary_gene_id}) for subunit {subunit_name_for_context}...")
        fetched_data = fetch_genes(primary_gene_id)
        time.sleep(1) 

        if not fetched_data.empty and not fetched_data.isna().all().all():
            if 'secondaryIdentifier' in fetched_data.columns:
                fetched_data.rename(columns={'secondaryIdentifier': 'Gene.secondaryIdentifier'}, inplace=True)
            
            fetched_data['Step2_Subunit'] = subunit_name_for_context 
            fetched_data['Input_PrimaryID'] = primary_gene_id 

            if 'Gene.secondaryIdentifier' in fetched_data.columns and not fetched_data['Gene.secondaryIdentifier'].empty:
                sec_id_val = fetched_data['Gene.secondaryIdentifier'].iloc[0]
                fetched_data['Input_ID'] = str(sec_id_val) if pd.notna(sec_id_val) and str(sec_id_val).strip() else primary_gene_id
            else:
                fetched_data['Input_ID'] = primary_gene_id

            all_results.append(fetched_data)
            print(f"      SUCCESS: Data fetched for PrimaryID {primary_gene_id}")
        else:
            print(f"      NOTE: No usable data fetched for PrimaryID {primary_gene_id}")
            all_results.append(pd.DataFrame([{
                'Input_PrimaryID': primary_gene_id,
                'Input_ID': primary_gene_id,
                'Step2_Subunit': subunit_name_for_context
            }]))

    if all_results:
        processed_results_for_concat = [df_item.dropna(axis=1, how='all') for df_item in all_results if not df_item.empty]
        processed_results_for_concat = [df_item for df_item in processed_results_for_concat if not df_item.empty]

        if processed_results_for_concat:
            combined_subunit_df = pd.concat(processed_results_for_concat, ignore_index=True)
            print(f"  Completed fetching {subunit_name_for_context} chunk {chunk_num + 1}/{total_chunks}. Total records: {len(combined_subunit_df)}")
            return combined_subunit_df
    print(f"  No valid data records to process for gene chunk (subunit {subunit_name_for_context}).")
    return pd.DataFrame()