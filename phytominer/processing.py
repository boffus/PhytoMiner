import pandas as pd

def process_homolog_data(df_combined):
    """
    Processes a combined DataFrame of homolog data.
    - Aggregates origin source organisms for each homolog.
    - Calculates occurrence count for each homolog.
    - Deduplicates entries to keep the most relevant homolog.

    Parameters:
        df_combined (pd.DataFrame): DataFrame containing combined homolog data.

    Returns:
        pd.DataFrame: Processed and deduplicated DataFrame.
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

    # 2. Aggregate Origin Source Organisms
    origin_key_cols = ['subunit1', 'primaryIdentifier']
    if 'source.organism' in processed_df.columns:
        origin_map = processed_df.groupby(origin_key_cols, observed=False)['source.organism'] \
            .agg(lambda x: tuple(sorted(x.dropna().unique()))).reset_index()
        origin_map.rename(columns={'source.organism': 'origin.source.organisms'}, inplace=True)
        processed_df = pd.merge(processed_df.drop(columns=['origin.source.organisms'], errors='ignore'),
                                origin_map, on=origin_key_cols, how='left')
    else:
        processed_df['origin.source.organisms'] = None

    # 3. Calculate Homolog Occurrences
    processed_df['homolog.occurrences'] = processed_df.groupby(origin_key_cols, observed=False)['source.gene'] \
        .transform('size')

    # 4. Deduplication
    sort_by_cols = ['subunit1', 'relationship', 'homolog.occurrences', 'organism.shortName', 'primaryIdentifier', 'source.organism']
    sort_by_cols = [col for col in sort_by_cols if col in processed_df.columns]

    ascending_map = {'relationship': True, 'homolog.occurrences': False}
    ascending_order = [ascending_map.get(col, True) for col in sort_by_cols]

    dedup_subset_cols = ['subunit1', 'primaryIdentifier', 'organism.shortName']
    dedup_subset_cols = [col for col in dedup_subset_cols if col in processed_df.columns]

    if dedup_subset_cols:
        processed_df = processed_df.sort_values(by=sort_by_cols, ascending=ascending_order)
        processed_df = processed_df.drop_duplicates(subset=dedup_subset_cols, keep='first')

    # 5. Define final column order
    final_columns = [
        'subunit1', 'source.organism', 'source.gene', 'relationship',
        'primaryIdentifier', 'secondaryIdentifier', 'organism.shortName', 'organism.commonName',
        'organism.proteomeId', 'gene.length', 'sequence.length', 'sequence.residues',
        'homolog.occurrences', 'origin.source.organisms'
    ]

    # Reorder columns to a standard format, keeping any extra columns at the end
    existing_final_columns = [col for col in final_columns if col in processed_df.columns]
    other_columns = [col for col in processed_df.columns if col not in existing_final_columns]
    processed_df = processed_df[existing_final_columns + other_columns]

    print(f"Processing complete. DataFrame shape after processing: {processed_df.shape}")
    return processed_df
