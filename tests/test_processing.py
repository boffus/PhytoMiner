import logging
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal
from phytominer.processing import process_homolog_data, genes_fetch

def make_homolog_df(data):
    columns = [
        'subunit1', 'source.organism', 'source.gene', 'relationship',
        'primaryIdentifier', 'organism.shortName'
    ]
    return pd.DataFrame(data, columns=columns)

@pytest.mark.parametrize(
    "data,expected_len,expected_occurrences,expected_relationship",
    [
        # Complex case: deduplication and occurrence counting
         (
            [
                {'subunit1': 'NDHA', 'source.organism': 'A. thaliana TAIR10', 'source.gene': 'AT_NDHA', 'relationship': 'one-to-one', 'primaryIdentifier': 'sbicolor_gene1', 'organism.shortName': 'S. bicolor v3.1.1'},
                {'subunit1': 'NDHB', 'source.organism': 'A. thaliana TAIR10', 'source.gene': 'AT_NDHB_1', 'relationship': 'one-to-one', 'primaryIdentifier': 'sbicolor_gene2', 'organism.shortName': 'S. bicolor v3.1.1'},
                {'subunit1': 'NDHB', 'source.organism': 'A. thaliana TAIR10', 'source.gene': 'AT_NDHB_2', 'relationship': 'one-to-many', 'primaryIdentifier': 'sbicolor_gene2', 'organism.shortName': 'S. bicolor v3.1.1'},
            ],
            2,
            {'sbicolor_gene1': 1, 'sbicolor_gene2': 2},
            {'sbicolor_gene2': 'one-to-one'}
        ),
    ]
)
def test_process_homolog_data_cases(data, expected_len, expected_occurrences, expected_relationship):
    df = make_homolog_df(data)
    processed_df = process_homolog_data(df)
    assert len(processed_df) == expected_len
    assert 'homolog.occurrences' in processed_df.columns
    processed_df = processed_df.set_index('primaryIdentifier')
    for gene, occ in expected_occurrences.items():
        assert processed_df.loc[gene]['homolog.occurrences'] == occ
    for gene, rel in expected_relationship.items():
        assert processed_df.loc[gene]['relationship'] == rel

def test_process_homolog_data_with_empty_input():
    empty_df = make_homolog_df([])
    processed_df = process_homolog_data(empty_df)
    assert processed_df.empty
    assert_frame_equal(processed_df, empty_df)

@pytest.fixture
def mock_phytozome_genes(mocker):
    """
    Mocks the phytozome_genes API call where it is used: in the processing module.
    """
    return mocker.patch("phytominer.processing.phytozome_genes")

def test_genes_fetch_all_genes_found(mock_phytozome_genes):
    gene_ids = ['geneA', 'geneB', 'geneC']
    all_data = {
        'geneA': {'value': 10},
        'geneB': {'value': 20},
        'geneC': {'value': 30}
    }
    def side_effect_func(chunk):
        data = [{'primaryIdentifier': gid, **all_data[gid]} for gid in chunk if gid in all_data]
        return pd.DataFrame(data)

    mock_phytozome_genes.side_effect = side_effect_func

    df = genes_fetch(gene_ids, chunk_size=2)
    assert not df.empty
    assert set(df['primaryIdentifier']) == set(gene_ids)
    assert set(df['value']) == {10, 20, 30}

def test_genes_fetch_some_genes_missing(mock_phytozome_genes):
    gene_ids = ['geneX', 'geneY', 'geneZ']
    all_data = {'geneX': {'value': 1}, 'geneZ': {'value': 3}}

    def side_effect_func(chunk):
        data = [{'primaryIdentifier': gid, **all_data[gid]} for gid in chunk if gid in all_data]
        return pd.DataFrame(data)

    mock_phytozome_genes.side_effect = side_effect_func

    df = genes_fetch(gene_ids, chunk_size=3)
    assert len(df) == 2
    assert set(df['primaryIdentifier']) == {'geneX', 'geneZ'}

def test_genes_fetch_no_genes_found(mock_phytozome_genes, caplog):
    gene_ids = ['geneA', 'geneB']
    mock_phytozome_genes.return_value = pd.DataFrame()

    with caplog.at_level(logging.WARNING):
        df = genes_fetch(gene_ids, chunk_size=1)

    assert df.empty
    assert "All chunks resulted in empty data or errors." in caplog.text