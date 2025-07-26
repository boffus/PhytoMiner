import pandas as pd
import logging
logger = logging.getLogger(__name__)
import time
from typing import Dict, List
from intermine.webservice import Service
from .config import (
    PHYTOZOME_SERVICE_URL,
    DEFAULT_WORKERS,
    DEFAULT_CHUNK,
    DEFAULT_SLEEP,
    HOMOLOGS_QUERY_VIEW,
    HOMOLOGS_COLUMN_ORDER, 
    GENES_QUERY_VIEW    
    )

def phytozome_homologs(
    source_organism: str,
    transcript_chunk: List[str],
    subunit_map_for_transcripts: Dict[str, str]
) -> pd.DataFrame:
    """
    Fetches homologs from Phytozome in chunk's.

    Parameters:
        source_organism_name (str): The shortName of the organism being queried.
        transcript_chunk (list): A CHUNK of gene primaryIdentifiers.
        subunit_map_for_transcripts (dict): A dictionary mapping transcript_identifiers to Subunit names.
    Returns:
        pd.DataFrame: DataFrame containing fetched homolog data for the chunk.
                      Returns an empty DataFrame if no homologs are found or an error occurs.
    """
    if not transcript_chunk:
        return pd.DataFrame()

    service = Service(PHYTOZOME_SERVICE_URL)
    logger.info(f"Processing {len(transcript_chunk)} transcripts in chunk from {source_organism}...")

    query = service.new_query("Homolog")
    query.add_view(*HOMOLOGS_QUERY_VIEW)
    query.add_sort_order("Homolog.relationship", "DESC")
    query.add_constraint("gene.primaryIdentifier", "ONE OF", transcript_chunk, code="A")
    query.add_constraint("organism.shortName", "=", source_organism, code="B")
    query.add_constraint("ortholog_organism.shortName", "!=", source_organism, code="C")
    query.set_logic("A and B and C")

    all_rows = []
    max_retries = 5
    backoff = 0.5
    for attempt in range(1, max_retries + 1):
        try:
            all_rows = [{
                    "source.organism": source_organism,
                    "source.gene": row["gene.primaryIdentifier"],
                    "primaryIdentifier": row["ortholog_gene.primaryIdentifier"],
                    "secondaryIdentifier": row["ortholog_gene.secondaryIdentifier"],
                    "gene.length": row["ortholog_gene.length"],
                    "gene.genomicOrder": row["ortholog_gene.genomicOrder"],
                    "sequence.length": row["ortholog_gene.sequence.length"],
                    "sequence.residues": row["ortholog_gene.sequence.residues"],
                    "organism.shortName": row["ortholog_organism.shortName"],
                    "organism.commonName": row["ortholog_organism.commonName"],
                    "organism.proteomeId": row["ortholog_organism.proteomeId"],
                    "relationship": row["relationship"]
                }
                for row in query.rows()
            ]
            logger.info(f"[{source_organism}] chunk fetched ({len(transcript_chunk)} items), "
                        f"{len(all_rows)} homologs collected.")
            break
        except Exception as exc:
            logger.warning(f"[{source_organism}] Attempt {attempt} failed: {exc}")
            if attempt == max_retries:
                logger.error(f"[{source_organism}] Max retries reached—skipping chunk.")
                return pd.DataFrame()
            sleep_time = backoff * (2 ** (attempt - 1))
            time.sleep(sleep_time)

    if not all_rows:
        return pd.DataFrame()

    df = pd.DataFrame(all_rows)
    df["subunit1"] = df["source.gene"].map(subunit_map_for_transcripts)

    available_cols = [c for c in HOMOLOGS_COLUMN_ORDER if c in df.columns]
    other_cols = [c for c in df.columns if c not in available_cols]
    df = df[available_cols + other_cols]
    return df

def phytozome_genes(gene_primary_ids: List[str]) -> pd.DataFrame:
    """
    Batch query Gene objects in Phytozome via InterMine, returns a DataFrame.

    Parameters:
        gene_primary_ids (List[str]): A list of gene primaryIdentifier strings.

    Returns:
        pd.DataFrame: A DataFrame containing the fetched gene data.
                      Returns an empty DataFrame if the input is invalid, no genes are found,
                      or an error occurs.
    """
    if not gene_primary_ids or not isinstance(gene_primary_ids, (list, tuple)):
        logger.warning(f"Invalid gene list: {gene_primary_ids}")
        return pd.DataFrame()
    
    service = Service(PHYTOZOME_SERVICE_URL)
    query = service.new_query("Gene")
    query.add_constraint("primaryIdentifier", "ONE OF", gene_primary_ids, code="G")
    query.outerjoin("proteins")
    query.outerjoin("rnaSeqExpressions")
    query.outerjoin("coexpressions")
    query.outerjoin("CDSs")
        
    query.add_view(*GENES_QUERY_VIEW)
        
    try:
        rows = list(query.rows())
        if not rows:
            logger.info("No gene records found for the provided IDs.")
            return pd.DataFrame()
        df = pd.DataFrame(rows)
        logger.info(f"Fetched {len(df)} gene records.")
        return df
    except Exception as exc:
        logger.error(f"Service error fetching genes: {exc}")
        return pd.DataFrame()