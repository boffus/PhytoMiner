import pandas as pd
import pytest
from pandas.testing import assert_frame_equal
from phytominer.processing import process_homolog_data


@pytest.fixture
def sample_homolog_df():
    """Provides a complex DataFrame to test processing logic."""
    data = [
        # Case 1: Simple, unique homolog that should be kept as is.
        {'Subunit': 'NDHA', 'source.organism': 'A. thaliana TAIR10', 'source.gene': 'AT_NDHA', 'relationship': 'one-to-one', 'primaryIdentifier': 'sbicolor_gene1', 'organism.shortName': 'S. bicolor v3.1.1'},

        # Case 2: Deduplication based on relationship. 'one-to-one' should be preferred.
        {'Subunit': 'NDHB', 'source.organism': 'A. thaliana TAIR10', 'source.gene': 'AT_NDHB_1', 'relationship': 'one-to-one', 'primaryIdentifier': 'sbicolor_gene2', 'organism.shortName': 'S. bicolor v3.1.1'},
        {'Subunit': 'NDHB', 'source.organism': 'A. thaliana TAIR10', 'source.gene': 'AT_NDHB_2', 'relationship': 'one-to-many', 'primaryIdentifier': 'sbicolor_gene2', 'organism.shortName': 'S. bicolor v3.1.1'},

        # Case 3: Tests homolog.occurrences. 'gene3' is found from two source genes. Count should be 2.
        {'Subunit': 'NDHC', 'source.organism': 'A. thaliana TAIR10', 'source.gene': 'AT_NDHC_1', 'relationship': 'one-to-one', 'primaryIdentifier': 'sbicolor_gene3', 'organism.shortName': 'S. bicolor v3.1.1'},
        {'Subunit': 'NDHC', 'source.organism': 'A. thaliana TAIR10', 'source.gene': 'AT_NDHC_2', 'relationship': 'one-to-one', 'primaryIdentifier': 'sbicolor_gene3', 'organism.shortName': 'S. bicolor v3.1.1'},

        # Case 4: Tests origin.source.organisms aggregation. 'gene4' is from two source organisms.
        {'Subunit': 'NDHD', 'source.organism': 'A. thaliana TAIR10', 'source.gene': 'AT_NDHD', 'relationship': 'one-to-one', 'primaryIdentifier': 'sbicolor_gene4', 'organism.shortName': 'S. bicolor v3.1.1'},
        {'Subunit': 'NDHD', 'source.organism': 'O. sativa Kitaake v3.1', 'source.gene': 'OS_NDHD', 'relationship': 'one-to-one', 'primaryIdentifier': 'sbicolor_gene4', 'organism.shortName': 'S. bicolor v3.1.1'},
    ]
    columns = ['Subunit', 'source.organism', 'source.gene', 'relationship', 'primaryIdentifier', 'organism.shortName']
    return pd.DataFrame(data, columns=columns)


def test_process_homolog_data_with_complex_cases(sample_homolog_df):
    """
    Tests deduplication, occurrence counting, and source aggregation.
    """
    processed_df = process_homolog_data(sample_homolog_df)

    # Expected number of rows after deduplication is 4 (one for each unique homolog gene)
    assert len(processed_df) == 4
    assert 'homolog.occurrences' in processed_df.columns
    assert 'origin.source.organisms' in processed_df.columns

    # Set index for easy lookup
    processed_df = processed_df.set_index('primaryIdentifier')

    # Test Case 1: Simple homolog
    assert processed_df.loc['sbicolor_gene1']['homolog.occurrences'] == 1
    assert processed_df.loc['sbicolor_gene1']['origin.source.organisms'] == ('A. thaliana TAIR10')

    # Test Case 2: Deduplication by relationship
    # The 'one-to-one' relationship should have been kept.
    assert processed_df.loc['sbicolor_gene2']['relationship'] == 'one-to-one'
    # The occurrence count should be 2, calculated before deduplication.
    assert processed_df.loc['sbicolor_gene2']['homolog.occurrences'] == 2

    # Test Case 3: Occurrence counting
    assert processed_df.loc['sbicolor_gene3']['homolog.occurrences'] == 2

    # Test Case 4: Source organism aggregation
    # The tuple should contain both source organisms, sorted alphabetically.
    assert processed_df.loc['sbicolor_gene4']['origin.source.organisms'] == ('A. thaliana TAIR10', 'O. sativa Kitaake v3.1')
    assert processed_df.loc['sbicolor_gene4']['homolog.occurrences'] == 2


def test_process_homolog_data_with_empty_input():
    """
    Tests that the function handles an empty DataFrame gracefully.
    """
    empty_df = pd.DataFrame(columns=[
        'subunit1', 'source.organism', 'source.gene', 'relationship',
        'primaryIdentifier', 'organism.shortName'
    ])
    processed_df = process_homolog_data(empty_df)
    assert processed_df.empty
    assert_frame_equal(processed_df, empty_df)


def test_process_homolog_data_with_no_duplicates():
    """
    Tests that the function works correctly when there is nothing to deduplicate.
    """
    data = [
        {'Subunit': 'NDHA', 'source.organism': 'A. thaliana TAIR10', 'source.gene': 'AT_NDHA', 'relationship': 'one-to-one', 'primaryIdentifier': 'sbicolor_gene1', 'organism.shortName': 'S. bicolor v3.1.1'},
        {'Subunit': 'NDHB', 'source.organism': 'A. thaliana TAIR10', 'source.gene': 'AT_NDHB', 'relationship': 'one-to-one', 'primaryIdentifier': 'osativa_gene2', 'organism.shortName': 'O. sativa Kitaake v3.1'}
    ]
    df = pd.DataFrame(data)
    processed_df = process_homolog_data(df)

    # No rows should be dropped
    assert len(processed_df) == 2

    # Check occurrence and origin for the first row
    row1 = processed_df[processed_df['primaryIdentifier'] == 'sbicolor_gene1'].iloc[0]
    assert row1['homolog.occurrences'] == 1
    assert row1['origin.source.organisms'] == ('A. thaliana TAIR10')


