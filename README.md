# PhytoMiner
This is a package for fetching Phytozome data

[![CI](https://github.com/boffus/PhytoMiner/actions/workflows/python-publish.yml/badge.svg)](https://github.com/boffus/PhytoMiner/actions/workflows/python-publish.yml)
[![PyPI version](https://badge.fury.io/py/phytominer.svg)](https://badge.fury.io/py/phytominer)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Python library for efficiently fetching and processing gene homolog data from the [Phytozome](https://phytozome-next.jgi.doe.gov/) database via its InterMine API.

This library is designed to simplify complex, iterative bioinformatic queries, allowing researchers to trace gene homology across multiple species with ease.

## Features

- **Three-Step Pipeline**: A clear, sequential workflow for fetching homologs, merging local data, and retrieving detailed gene information.
- **Iterative Search**: Automatically performs chained searches using homologs found in previous steps to build a comprehensive dataset.
- **Parallel Processing**: Utilizes multithreading for efficient, parallel data fetching, significantly speeding up large queries.
- **Checkpointing**: Automatically saves and loads intermediate results to prevent losing progress and allow for easy resumption of long-running jobs.
- **Data Processing & Visualization**: Includes functions to clean, de-duplicate, and enrich data, plus a utility to quickly generate a heatmap of homolog distribution.

## Installation

You can install the latest `PhytoMiner` release directly from PyPI:

```bash
pip install phytominer
```

## Usage

Here is a complete example of the three-step workflow:

Define a set of known genes in a source organism (e.g., A. thaliana).
Run homologs_pipe to find homologs in other species.
Run join_tsvs to combine the homolog data with local metadata from TSV files.
Run genes_pipe to fetch detailed gene data for the final homolog set.

```python
import logging
from phytominer.workflow import step1_homolog_pipe, step2_merge_pipe, step3_gene_pipe

# It's highly recommended to configure logging to see the progress
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# 1. Define initial search parameters
# The initial organism to start the search from
initial_organism = "A. thaliana TAIR10"

# A dictionary of initial gene IDs and their corresponding subunit names
initial_genes = {
    "AT1G01090": "NDHA",
    "AT1G01120": "NDHB",
    "ATCG00520": "NDHC",
}

# A list of other organisms to find homologs in
subsequent_organisms = [
    "S. bicolor v3.1.1",
    "O. sativa Kitaake v3.1",
    "S. viridis v2.1"
]

# 2. Run the three-step pipeline
# Step 1: Fetch all homolog data, starting with the initial organism
# and iterating through the subsequent ones.
step1_df = homologs_pipe(
    initial_organism=initial_organism,
    initial_genes_dict=initial_genes,
    subsequent_organisms=subsequent_organisms
)

# Step 2: Merge the homolog data with local TSV files containing additional metadata.
# This step assumes you have a directory with TSV files (e.g., 'data/tsv/').
step2_df = join_tsvs()

# Step 3: Fetch detailed gene data (e.g., expression, sequence) for the homologs found.
step3_df = genes_pipe()

# The final DataFrames are saved to CSV files at each step (e.g., step1output.csv).
print("PhytoMiner workflow complete!")

```

## API Overview

The phytominer library is structured around a sequential, three-step workflow.

Workflow Functions
These are the main functions you'll use, found in phytominer.workflow.

homologs_pipe(...): Orchestrates the entire homolog search. It starts with an initial set of genes, finds their homologs, and then iteratively searches for homologs of the results in other specified organisms. It handles checkpointing and produces a final, processed DataFrame of homolog data.
join_tsvs(...): Takes the output from Step 1 and merges it with local TSV files containing supplementary data (e.g., subunit validation).
genes_pipe(...): Takes the output from Step 2 and fetches detailed gene information (sequences, expression data, etc.) for all unique homologs identified in the pipeline.

### Utility Functions

These helper functions are available in phytominer.utils.

pivotmap(dataframe, ...): Generates a pivot table and a corresponding heatmap to visualize the count of homologs across different species and subunits.
log_summary(df, ...): Logs a concise summary of a DataFrame's shape, columns, memory usage, and other key statistics.

## Continuous Integration & Deployment

This project uses [GitHub Actions](https://github.com/features/actions) for automated testing and publishing.

- **Automated Testing:**  
  Every push to the `main` branch triggers the test suite using Python 3.9.
- **Automated Publishing:**  
  When a new release is published on GitHub, the package is automatically built and uploaded to PyPI.

You can find the workflow configuration in [`.github/workflows/python-publish.yml`](.github/workflows/python-publish.yml).

## Contributing

Contributions are welcome! If you have a suggestion or find a bug, please open an issue. Pull requests are also encouraged.

1.  Fork the repository.
2.  Create your feature branch (`git checkout -b feature/AmazingFeature`).
3.  Commit your changes (`git commit -m 'Add some AmazingFeature'`).
4.  Push to the branch (`git push origin feature/AmazingFeature`).
5.  Open a Pull Request.

### Running Tests Locally

To run the test suite locally:

```bash
pip install -e .[dev]
pytest
```

## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Contact

Author: Kris Kari
Email: toffe.kari@gmail.com