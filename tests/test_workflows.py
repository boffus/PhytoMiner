import pandas as pd
import pytest
from pandas.testing import assert_frame_equal
from phytominer.workflow import (
    homologs_pipe,
    join_tsvs,
    genes_pipe,
)

@pytest.fixture
def sample_data():
    """Provides sample DataFrames for testing."""
    return {
        "initial_fetch": pd.DataFrame({
            'source.gene': ['AT1G01090'],
            'primaryIdentifier': ['LOC_Os01g01010'],
            'organism.shortName': ['osativa'],
            'subunit1': ['NDHA']
        }),
        "subsequent_fetch": pd.DataFrame({
            'source.gene': ['LOC_Os01g01010'],
            'primaryIdentifier': ['Solyc01g005000'],
            'organism.shortName': ['slycopersicum'],
            'subunit1': ['NDHA']
        }),
        "processed": pd.DataFrame({'gene_id': ['processed_gene']}),
        "tsv": pd.DataFrame({
            'Gene_ID_from_TSV': ['AT1G01090'],
            'subunit2': ['NDHA']
        }),
        "master_for_step3": pd.DataFrame({
            'subunit1': ['NDHA', 'NDHA', 'NDHB'],
            'primaryIdentifier': ['gene1', 'gene2', 'gene3'],
            'organism.shortName': ['org1', 'org2', 'org1']
        }),
        "fetched_genes": pd.DataFrame({
            'primaryIdentifier': ['gene1', 'gene2', 'gene3'],
            'sequence.length': [100, 200, 300]
        })
    }

@pytest.fixture
def mock_dependencies(mocker, sample_data):
    """Mocks all external dependencies for the workflow functions."""
    mocks = {
        "initial_fetch": mocker.patch("phytominer.workflow.initial_fetch", return_value=sample_data["initial_fetch"]),
        "subsequent_fetch": mocker.patch("phytominer.workflow.subsequent_fetch", return_value=sample_data["subsequent_fetch"]),
        "process_homolog_data": mocker.patch("phytominer.workflow.process_homolog_data", side_effect=lambda df: df),
        "read_tsv_files": mocker.patch("phytominer.workflow.read_tsv_files", return_value=sample_data["tsv"]),
        "load_master_df": mocker.patch("phytominer.workflow.load_master_df", return_value=sample_data["master_for_step3"]),
        "genes_fetch": mocker.patch("phytominer.workflow.genes_fetch", return_value=sample_data["fetched_genes"]),
        "pivotmap": mocker.patch("phytominer.workflow.pivotmap"),
        "to_csv": mocker.patch("pandas.DataFrame.to_csv"),
        "read_csv": mocker.patch("pandas.read_csv"),
        "path_exists": mocker.patch("pathlib.Path.exists")
    }
    return mocks

def test_step1_homologs_pipe_no_checkpoints(mock_dependencies):
    """Test the happy path for step1_homolog_pipe where no checkpoints exist."""
    mock_dependencies["path_exists"].return_value = False

    homologs_pipe(
        initial_organism='athaliana',
        initial_genes_dict={'AT1G01090': 'NDHA'},
        subsequent_organisms=['osativa']
    )

    mock_dependencies["initial_fetch"].assert_called_once()
    mock_dependencies["subsequent_fetch"].assert_called_once()


def test_step1_homologs_pipe_with_checkpoints(mock_dependencies, sample_data):
    """Test that step1_homolog_pipe loads from checkpoints instead of fetching."""
    mock_dependencies["path_exists"].return_value = True
    mock_dependencies["read_csv"].return_value = sample_data["initial_fetch"]
    homologs_pipe(
        initial_organism='athaliana',
        initial_genes_dict={'AT1G01090': 'NDHA'},
        subsequent_organisms=['osativa']
    )

    mock_dependencies["initial_fetch"].assert_not_called()
    mock_dependencies["subsequent_fetch"].assert_not_called()
    assert mock_dependencies["read_csv"].call_count == 2 # Initial + subsequent

def test_step2_join_tsvs(mock_dependencies, sample_data):
    """Test the step2 merge pipeline."""
    mock_dependencies["read_csv"].return_value = sample_data["initial_fetch"]

    result_df = join_tsvs()

    mock_dependencies["read_tsv_files"].assert_called_once()
    assert 'subunit2' in result_df.columns

def test_step3_genes_pipe(mock_dependencies):
    """Test the step3 genes fetching pipeline."""
    mock_dependencies["path_exists"].return_value = False

    genes_pipe()

    mock_dependencies["load_master_df"].assert_called_once()
    assert mock_dependencies["genes_fetch"].call_count == 2