"""
PhytoMiner: A toolkit for analysing plant genomic data from Phytozome.
"""

__version__ = "0.1.2"

__all__ = [
    'config',
    'data',
    'processing',
    'utils',
    'workflow',
    'run_homologs_pipeline',
    'run_workflow2'
]

# Import essential components
from . import config
from . import data
from . import processing
from . import utils
from . import workflow
from .workflow import run_homologs_pipeline, run_workflow2

