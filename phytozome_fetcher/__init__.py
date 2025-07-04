"""
A Python library to fetch and process homolog data from Phytozome.
"""

__version__ = "0.1.0"

from .api import initial_fetch, subsequent_fetch
from .processing import process_homolog_data
from .utils import pivotmap, print_summary
from .data import SUBCOMPLEX_DICT, SUBUNIT_TO_SUBCOMPLEX
from .workflow import run_homolog_pipeline
