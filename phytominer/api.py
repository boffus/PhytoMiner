import pandas as pd
import concurrent.futures
import time
import math
from intermine.webservice import Service
from phytominer.config import PHYTOZOME_SERVICE_URL

def phytozome_homologs(source_organism_name, transcript_chunk, subunit_map_for_transcripts):
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
    all_results_for_chunk = []
    print(f"  Processing {len(transcript_chunk)} transcripts in chunk from {source_organism_name}...")

    query = service.new_query("Homolog")
    query.add_view(
        "gene.primaryIdentifier", "relationship", "ortholog_organism.commonName",
        "ortholog_organism.shortName", "ortholog_organism.proteomeId",
        "ortholog_gene.primaryIdentifier", "ortholog_gene.length",
        "ortholog_gene.secondaryIdentifier", "ortholog_gene.genomicOrder",
        "ortholog_gene.sequence.length", "ortholog_gene.sequence.residues")   
    query.add_sort_order("Homolog.relationship", "DESC")   
    query.add_constraint("gene.primaryIdentifier", "ONE OF", transcript_chunk, code="A")
    query.add_constraint("organism.shortName", "=", source_organism_name, code="B")
    query.add_constraint("ortholog_organism.shortName", "!=", source_organism_name, code="C")
    query.set_logic("A and B and C")

    # Retry logic for connection timeouts/failures
    max_retries = 5
    backoff_factor = 0.5  # seconds
    for attempt in range(1, max_retries + 1):
        try:
            for row in query.rows():
                all_results_for_chunk.append({
                    "source.organism": source_organism_name,
                    "source.gene": row["gene.primaryIdentifier"],
                    "primaryIdentifier": row["ortholog_gene.primaryIdentifier"],
                    "secondaryIdentifier": row["ortholog_gene.secondaryIdentifier"],
                    "gene.length": row["ortholog_gene.length"],
                    "sequence.length": row["ortholog_gene.sequence.length"],
                    "sequence.residues": row["ortholog_gene.sequence.residues"],
                    "organism.shortName": row["ortholog_organism.shortName"],
                    "organism.commonName": row["ortholog_organism.commonName"],
                    "organism.proteomeId": row["ortholog_organism.proteomeId"],
                    "relationship": row["relationship"],
                })
            print(f"  Retrieved {len(all_results_for_chunk)} homologs for chunk starting with {transcript_chunk[0]}")
            break  # Success, exit retry loop
        except Exception as e:
            print(f"  Error querying chunk for {source_organism_name} (attempt {attempt}/{max_retries}): {e}")
            if attempt == max_retries:
                print(f"  Max retries reached for chunk from {source_organism_name}. Skipping this chunk.")
                return pd.DataFrame()
            sleep_time = backoff_factor * (2 ** (attempt - 1))
            print(f"  Retrying in {sleep_time:.1f} seconds...")
            time.sleep(sleep_time)

    if not all_results_for_chunk:
        return pd.DataFrame()

    chunk_df = pd.DataFrame(all_results_for_chunk)
    chunk_df['subunit1'] = chunk_df['source.gene'].map(subunit_map_for_transcripts)

    # Reorder columns for clarity
    ordered_columns = ['source.organism', 'source.gene', 'relationship', 'subunit1', 'primaryIdentifier',
        'secondaryIdentifier', 'organism.commonName', 'organism.shortName', 'organism.proteomeId',
        'gene.length', 'sequence.length', 'sequence.residues']
    # Ensure all ordered columns exist, and add any others at the end
    existing_ordered_columns = [col for col in ordered_columns if col in chunk_df.columns]
    other_columns = [col for col in chunk_df.columns if col not in existing_ordered_columns]
    chunk_df = chunk_df[existing_ordered_columns + other_columns]
    return chunk_df

def initial_fetch(source_organism_name, transcript_names, subunit_dict, max_workers=8):
    """
    Fetches homologs from Phytozome for an initial list of transcripts from a source organism.
    Handles chunking and parallel execution using ThreadPoolExecutor.

    Parameters:
        source_organism_name (str): The shortName of the organism whose transcripts are being queried.
        transcript_names (list): A list of gene primaryIdentifiers to query for homologs.
        subunit_dict (dict): A dictionary mapping transcript_identifiers to their Subunit names.
        max_workers (int): The maximum number of parallel threads to use.

    Returns:
        pd.DataFrame: DataFrame containing fetched homolog data for all transcripts.
    """
    if not transcript_names:
        print(f"No transcript names provided for initial fetch from {source_organism_name}. Returning empty DataFrame.")
        return pd.DataFrame()

    print(f"Initiating initial homolog fetch for {len(transcript_names)} transcripts from {source_organism_name}...")

    subunit_map_for_transcripts = {tid: subunit_dict.get(tid) for tid in transcript_names}

    chunk_size = math.ceil(len(transcript_names) / max_workers) if transcript_names else 1
    chunk_size = max(1, chunk_size)
    transcript_id_chunks = [transcript_names[i:i + chunk_size] for i in range(0, len(transcript_names), chunk_size)]
    print(f"  Split {len(transcript_names)} transcripts into {len(transcript_id_chunks)} chunks.")

    all_fetched_data = []
    if transcript_id_chunks:
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_chunk = {
                executor.submit(pythozome_homologs, source_organism_name, chunk, subunit_map_for_transcripts): chunk
                for chunk in transcript_id_chunks
            }
            for future in concurrent.futures.as_completed(future_to_chunk):
                try:
                    result_df_chunk = future.result()
                    if not result_df_chunk.empty:
                        all_fetched_data.append(result_df_chunk)
                except Exception as exc:
                    print(f"  A chunk for {source_organism_name} generated an exception: {exc}")

    return pd.concat(all_fetched_data, ignore_index=True) if all_fetched_data else pd.DataFrame()

def subsequent_fetch(current_df, target_organism_name, max_workers=8):
    """
    Identifies genes from a target organism within an existing homolog dataset,
    and then fetches their homologs from Phytozome.
    """
    print(f"Preparing to fetch homologs for subsequent organism: {target_organism_name}")

    transcripts_for_next_query_df = current_df[current_df['organism.shortName'] == target_organism_name]

    if transcripts_for_next_query_df.empty:
        print(f"  No existing homologs found in master list for {target_organism_name} to use as query.")
        return pd.DataFrame()

    next_transcript_ids = transcripts_for_next_query_df['primaryIdentifier'].unique().tolist()
    if not next_transcript_ids:
        print(f"  No unique transcript IDs derived for {target_organism_name}. Skipping fetch.")
        return pd.DataFrame()

    next_subunit_map = pd.Series(
        transcripts_for_next_query_df.subunit1.values,
        index=transcripts_for_next_query_df.primaryIdentifier
    ).to_dict()

    # Re-use parallel fetching
    return initial_fetch(target_organism_name, next_transcript_ids, next_subunit_map, max_workers)

def fetch_genes(gene_primary_id):
    """
    Fetches detailed gene information from Phytozome for a given gene primary identifier.
    """
    if not isinstance(gene_primary_id, str) or not gene_primary_id.strip():
        print(f"NOTE: Invalid gene identifier: {gene_primary_id}. Skipping.")
        return pd.DataFrame()
    try:
        service = Service(PHYTOZOME_SERVICE_URL)
        query = service.new_query("Gene")
        query.outerjoin("proteins")
        query.outerjoin("rnaSeqExpressions")
        query.outerjoin("coexpressions")
        
        query.add_view(
            "length", "primaryIdentifier", "secondaryIdentifier", "organism.commonName",
            "organism.proteomeId", "organism.shortName", "organism.species", "coexpressions.highRange", 
            "coexpressions.lowRange", "coexpressions.JSON", "proteins.length","proteins.primaryIdentifier",
            "rnaSeqExpressions.abundance", "rnaSeqExpressions.confhi", "rnaSeqExpressions.conflo", 
            "rnaSeqExpressions.count", "rnaSeqExpressions.countdispersionvar", "rnaSeqExpressions.countuncertaintyvar", 
            "rnaSeqExpressions.countvariance", "rnaSeqExpressions.libraryExpressionLevel", "rnaSeqExpressions.experiment.name", 
            "rnaSeqExpressions.experiment.experimentGroup", "sequence.length", "sequence.residues",
            "CDSs.gene.CDSs.name", "CDSs.primaryIdentifier", "CDSs.organism.name"
        )
        query.add_constraint("primaryIdentifier", "=", gene_primary_id) 
        
        records = [dict(row) for row in query.rows()]
        return pd.DataFrame(records)
        
    except Exception as e:
        print(f"    Service Error: Could not fetch info for {gene_primary_id}: {e}")
        return pd.DataFrame()


