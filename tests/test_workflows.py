import unittest
from unittest.mock import patch
import pandas as pd
from phytominer.workflow import run_homologs_pipeline
from phytominer.config import DEFAULT_MAX_WORKERS, JOIN2_OUTPUT_FILE

class TestWorkflow(unittest.TestCase):

    def setUp(self):
        """Set up test data and configurations."""
        self.initial_genes = {'AT1G01090': 'NdhA', 'AT1G01120': 'NdhB'}
        self.subsequent_orgs = ['osativa', 'slycopersicum']
        self.initial_org = 'athaliana'
        self.checkpoint_dir = "test_checkpoints"

        # Sample DataFrames to be returned by mocked functions
        self.sample_initial_df = pd.DataFrame({
            'source.gene': ['AT1G01090'],
            'primaryIdentifier': ['LOC_Os01g01010'],
            'organism.shortName': ['osativa'],
            'subunit1': ['NdhA']
        })

        self.sample_subsequent_df = pd.DataFrame({
            'source.gene': ['LOC_Os01g01010'],
            'primaryIdentifier': ['Solyc01g005000'],
            'organism.shortName': ['slycopersicum'],
            'subunit1': ['NdhA']
        })

        self.sample_homolog_df = pd.DataFrame({'gene_id': ['AT1G01090', 'AT1G01120']})
        self.sample_tsv_df = pd.DataFrame({
            'Gene_ID_from_TSV': ['AT1G01090'],
            'Validated_Subunit_from_file': ['ndhA']
        })

    @patch('phytominer.workflow.pivotmap')
    @patch('phytominer.workflow.print_summary')
    @patch('phytominer.workflow.process_homolog_data')
    @patch('phytominer.workflow.subsequent_fetch')
    @patch('phytominer.workflow.initial_fetch')
    @patch('pandas.DataFrame.to_csv')
    @patch('pandas.read_csv')
    @patch('os.path.exists')
    @patch('os.makedirs')
    def test_run_homologs_pipeline_happy_path_with_checkpoints(
        self, mock_makedirs, mock_exists, mock_read_csv, mock_to_csv,
        mock_initial_fetch, mock_subsequent_fetch, mock_process_homolog_data,
        mock_print_summary, mock_pivotmap
    ):
        # Setup mocks
        mock_exists.return_value = False
        mock_initial_fetch.return_value = self.sample_initial_df
        mock_process_homolog_data.return_value = self.sample_homolog_df

        # Run the pipeline with correct argument order
        run_homologs_pipeline(
            self.initial_org,
            self.initial_genes,
            self.subsequent_orgs,
            max_workers=DEFAULT_MAX_WORKERS,
            checkpoint_dir=self.checkpoint_dir
        )
        self.assertEqual(mock_subsequent_fetch.call_count, 2)
        calls = mock_subsequent_fetch.call_args_list
        self.assertEqual(calls[0][0][1], 'osativa')
        self.assertEqual(calls[1][0][1], 'slycopersicum')
        self.assertEqual(calls[0][0][2], DEFAULT_MAX_WORKERS)
        self.assertEqual(calls[1][0][2], DEFAULT_MAX_WORKERS)

if __name__ == '__main__':
    unittest.main()