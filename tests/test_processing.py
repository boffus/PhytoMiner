import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal
from phytominer import processing

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
                {'subunit1': 'NDHC', 'source.organism': 'A. thaliana TAIR10', 'source.gene': 'AT_NDHC_1', 'relationship': 'one-to-one', 'primaryIdentifier': 'sbicolor_gene3', 'organism.shortName': 'S. bicolor v3.1.1'},
                {'subunit1': 'NDHC', 'source.organism': 'A. thaliana TAIR10', 'source.gene': 'AT_NDHC_2', 'relationship': 'one-to-one', 'primaryIdentifier': 'sbicolor_gene3', 'organism.shortName': 'S. bicolor v3.1.1'},
            ],
            3,
            {'sbicolor_gene1': 1, 'sbicolor_gene2': 2, 'sbicolor_gene3': 2},
            {'sbicolor_gene2': 'one-to-one'}
        ),
        # No duplicates
        (
            [
                {'subunit1': 'NDHA', 'source.organism': 'A. thaliana TAIR10', 'source.gene': 'AT_NDHA', 'relationship': 'one-to-one', 'primaryIdentifier': 'sbicolor_gene1', 'organism.shortName': 'S. bicolor v3.1.1'},
                {'subunit1': 'NDHB', 'source.organism': 'A. thaliana TAIR10', 'source.gene': 'AT_NDHB', 'relationship': 'one-to-one', 'primaryIdentifier': 'osativa_gene2', 'organism.shortName': 'O. sativa Kitaake v3.1'}
            ],
            2,
            {'sbicolor_gene1': 1, 'osativa_gene2': 1},
            {}
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
