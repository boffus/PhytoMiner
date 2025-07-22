"""
PhytoMiner: A toolkit to query genomic data from Phytozome.
"""

__version__ = "0.1.4"

__all__ = [
    'config',
    'data',
    'processing',
    'utils',
    'workflow',
    'run_homologs_pipeline',
    'run_workflow2',
    'run_expressions_workflow',
    'load_master_df',
    'fetch_expression_data'
]

# Import essential components
from . import config
from . import data
from . import processing
from . import utils
from . import workflow

from .workflow import run_homologs_pipeline, run_workflow2, run_expressions_workflow
from .processing import load_master_df, fetch_expression_data