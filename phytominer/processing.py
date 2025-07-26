import pandas as pd
import logging
logger = logging.getLogger(__name__)
import concurrent.futures
from typing import Union
from .api import phytozome_homologs, phytozome_genes
from .config import (
    DEFAULT_WORKERS,
    DEFAULT_CHUNK,
    JOIN_COL,
    HOMOLOG_GENE_ID,
    TSV_GENE_ID
)

def parallel_chunk_fetch(
    items: list,
    fetch_fn: callable,
    fetch_args: tuple = None,
    max_workers: int = DEFAULT_WORKERS,
    chunk_size: int = DEFAULT_CHUNK,
    context_msg: str = ""
) -> pd.DataFrame:
    """
    Process items in parallel chunks using ThreadPoolExecutor.
    Ideal for I/O-bound tasks like API calls.
    """ 
    if not items:
        logger.info(f"{context_msg}No items to fetch.")
        return pd.DataFrame()

    chunks = [items[i:i + chunk_size] for i in range(0, len(items), chunk_size)]
    logger.info(f"{context_msg}Split {len(items)} items into {len(chunks)} chunks (size={chunk_size}).")

    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_chunk = {
            executor.submit(fetch_fn, chunk, *(fetch_args or [])): chunk
            for chunk in chunks
        }

        for future in concurrent.futures.as_completed(future_to_chunk):
            chunk = future_to_chunk[future]
            try:
                df = future.result()
                if isinstance(df, pd.DataFrame) and not df.empty:
                    results.append(df)
            except Exception as exc:
                logger.error(f"{context_msg}Chunk starting with '{chunk[:2]}...' generated an exception: {exc}", exc_info=True)

    if not results:
        logger.warning(f"{context_msg}All chunks resulted in empty data or errors.")
        return pd.DataFrame()

    df_all = pd.concat(results, ignore_index=True)
    logger.info(f"{context_msg}Successfully fetched data from {len(results)} chunks, totaling {len(df_all)} rows.")
    return df_all

def process_homolog_data(df_combined: pd.DataFrame) -> pd.DataFrame:
    """
    Processes a DataFrame of homolog data by calculating occurrence counts and deduplicating.
    """
    if df_combined.empty:
        logger.info("Input DataFrame is empty. No processing done.")
        return df_combined

    logger.info(f"Processing DataFrame with {len(df_combined)} rows...")
    processed_df = df_combined.copy()
    relationship_categories = ['one-to-one', 'one-to-many', 'many-to-one', 'many-to-many']
    processed_df['relationship'] = pd.Categorical(
        processed_df['relationship'],
        categories=relationship_categories,
        ordered=True
    )

    origin_key_cols = JOIN_COL
    processed_df['homolog.occurrences'] = processed_df.groupby(origin_key_cols, observed=False)['source.gene'].transform('size')
    sort_by_cols = ['subunit1', 'relationship', 'homolog.occurrences', 'organism.shortName', 'primaryIdentifier', 'source.organism']
    ascending_map = {'relationship': True, 'homolog.occurrences': False}
    ascending_order = [ascending_map.get(col, True) for col in sort_by_cols]
    dedup_subset_cols = ['subunit1', 'primaryIdentifier', 'organism.shortName']
    processed_df = processed_df.sort_values(by=sort_by_cols, ascending=ascending_order)
    processed_df = processed_df.drop_duplicates(subset=dedup_subset_cols, keep='first')
    return processed_df

def initial_fetch(
    source_organism_name: str,
    transcript_names: list,
    subunit_dict: dict,
    max_workers: int = DEFAULT_WORKERS,
    chunk_size: int = DEFAULT_CHUNK
) -> pd.DataFrame:
    """
    Fetches homologs from Phytozome for a list of transcripts.
    """
    subunit_map = {tid: subunit_dict.get(tid) for tid in transcript_names}

    def fetch_fn_wrapper(chunk, *args):
        return phytozome_homologs(args[0], chunk, args[1])

    return parallel_chunk_fetch(
        items=transcript_names,
        fetch_fn=fetch_fn_wrapper,
        fetch_args=(source_organism_name, subunit_map),
        max_workers=max_workers,
        chunk_size=chunk_size,
        context_msg=f"[{source_organism_name}] "
    )

def subsequent_fetch(
    current_df: pd.DataFrame,
    target_organism_name: str,
    max_workers: int = DEFAULT_WORKERS,
    chunk_size: int = DEFAULT_CHUNK
) -> pd.DataFrame:
    """
    Identifies genes from a target organism in a homolog dataset and fetches their homologs.
    """
    logger.info(f"Preparing subsequent fetch for organism: {target_organism_name}")
    next_query_df = current_df[current_df['organism.shortName'] == target_organism_name]
    if next_query_df.empty:
        return pd.DataFrame()

    next_transcript_ids = next_query_df['primaryIdentifier'].unique().tolist()
    subunit_map_df = next_query_df[['primaryIdentifier', 'subunit1']].drop_duplicates()
    next_subunit_map = dict(zip(subunit_map_df.primaryIdentifier, subunit_map_df.subunit1))

    return initial_fetch(
        source_organism_name=target_organism_name,
        transcript_names=next_transcript_ids,
        subunit_dict=next_subunit_map,
        max_workers=max_workers,
        chunk_size=chunk_size
    )

def step2_merge(
    homolog_df: pd.DataFrame,
    tsv_df: pd.DataFrame,
    homolog_gene_id: str = HOMOLOG_GENE_ID,
    tsv_gene_id: str = TSV_GENE_ID
) -> pd.DataFrame:
    """
    Merges homolog data with TSV data on gene IDs.
    """
    return pd.merge(
        homolog_df,
        tsv_df,
        left_on=homolog_gene_id,
        right_on=tsv_gene_id,
        how='left'
    )

def load_master_df(filepath: str) -> Union[pd.DataFrame, None]:
    """
    Loads the master homolog DataFrame from a CSV file and validates required columns.
    """
    logger.info(f"Loading master homolog data from: {filepath}")
    try:
        df = pd.read_csv(filepath)
        logger.info(f"Successfully loaded master homolog data. Shape: {df.shape}")
        required_cols = JOIN_COL + ['subunit1']
        if not all(col in df.columns for col in required_cols):
            return None
        return df
    except FileNotFoundError:
        return None

def genes_fetch(
    gene_ids: list,
    max_workers: int = DEFAULT_WORKERS,
    chunk_size: int = DEFAULT_CHUNK
) -> pd.DataFrame:
    """
    Fetches gene data for a list of gene IDs in parallel.
    """
    logger.info(f"Fetching details for {len(gene_ids)} gene IDs.")
    return parallel_chunk_fetch(
        items=gene_ids,
        fetch_fn=phytozome_genes,
        max_workers=max_workers,
        chunk_size=chunk_size,
        context_msg="[Gene Fetch] "
    )