"""
PhytoMiner: A toolkit for fetching and processing gene data from the Phytozome database.
"""
import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())

__version__ = "0.2.0"

__all__ = [
    # Main workflow functions
    "homologs_pipe",
    "join_tsvs",
    "genes_pipe",
    # Key utility functions
    "pivotmap",
    "log_summary"
]

from .workflow import homologs_pipe, join_tsvs, genes_pipe
from .utils import pivotmap, log_summary