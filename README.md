# PhytoMiner
This is a package for fetching Phytozome data

![CI](https://github.com/boffus/PhytoMiner/actions/workflows/python-publish.yml/badge.svg)
[![PyPI version](https://badge.fury.io/py/phytominer.svg)](https://badge.fury.io/py/phytominer)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Python library for efficiently fetching and processing gene homolog data from the [Phytozome](https://phytozome-next.jgi.doe.gov/) database via its InterMine API.

This library is designed to simplify complex, iterative bioinformatic queries, allowing researchers to trace gene homology across multiple species with ease.

## Features

- **Initial Fetch**: Start a search with a list of genes from a source organism.
- **Iterative Search**: Perform chained or subsequent searches using homologs found in previous steps.
- **Parallel Processing**: Utilizes multithreading for efficient, parallel data fetching, significantly speeding up large queries.
- **Data Processing**: Includes functions to clean, de-duplicate, and enrich the fetched data by calculating occurrence counts and aggregating source information.
- **Visualisation**: Comes with a utility function to quickly generate a heatmap of homolog distribution across species and subunits.

## Installation

You can install the latest `PhytoMiner` release directly from PyPI:

```bash
pip install phytominer
```

## Usage

Here is a complete example of a common workflow:
1.  Start with a set of known genes in a source organism (e.g. `A. thaliana TAIR10`).
2.  Perform an `initial_fetch` to find homologs in other species.
3.  Use the results to perform a `subsequent_fetch` for a specific target organism (`S. bicolor v3.1.1`).
4.  Combine and process the data.
5.  Visualize the homolog distribution with `pivotmap`.

```python
import pandas as pd
from phytominer import (
    run_homolog_pipeline,
    initial_fetch,
    subsequent_fetch,
    process_homolog_data,
    pivotmap,
    print_summary
)

# 1. Define initial query genes for Arabidopsis thaliana
# (Using a small subset for this example)
athaliana_genes = {
    'AT5G52100': 'CRR1',
    'AT3G46790': 'CRR2',
    'AT2G01590': 'CRR3',
}

# 2. Run the whole pipeline
from phytominer.workflow import run_homologs_pipeline
from phytominer.config import DEFAULT_MAX_WORKERS

results = run_homologs_pipeline(
    initial_organism="athaliana",
    initial_genes_dict=athaliana_genes,
    subsequent_organisms=["osativa", "slycopersicum"],
    max_workers=DEFAULT_MAX_WORKERS,
    checkpoint_dir="homolog_checkpoints"
)

# 3. Expression Data Fetch Workflow
from phytominer.workflow import run_expressions_workflow
from phytominer.config import (
    JOIN2_OUTPUT_FILE,
    EXPRESSION_CHECKPOINT_DIR,
    EXPRESSIONS_OUTPUT_FILE
)
from phytominer.processing import load_master_df, fetch_expression_data

run_expressions_workflow(
    master_file=JOIN2_OUTPUT_FILE,
    checkpoint_dir=EXPRESSION_CHECKPOINT_DIR,
    output_file=EXPRESSIONS_OUTPUT_FILE,
    fetch_expression_data=fetch_expression_data,
    load_master_df=load_master_df
)

# 4. Alternatively perform the initial fetch from Arabidopsis thaliana
print("--- Starting Initial Fetch ---")
initial_df = initial_fetch(
    source_organism_name="A. thaliana TAIR10",
    transcript_names=list(athaliana_genes.keys()),
    subunit_dict=athaliana_genes,
    max_workers=4
)
print_summary(initial_df, "Initial Fetch Results")

# 5. Perform a subsequent fetch using homologs found in Sorghum bicolor
print("\n--- Starting Subsequent Fetch for Sorghum bicolor ---")
subsequent_df = subsequent_fetch(
    current_master_df=initial_df,
    target_organism_name="S. bicolor v3.1.1",
    max_workers=4
)
print_summary(s_df, "Subsequent Fetch Results for Sorghum")

# 6. Combine and process the data
print("\n--- Combining and Processing Data ---")
master_df = pd.concat([initial_df, subsequent_df], ignore_index=True)
processed_df = process_homolog_data(master_df)
print_summary(processed_df, "Final Processed DataFrame")

# 7. Visualize missing genes
print("\n--- Generating Heatmap ---")
# For a cleaner plot, let's display the top 15 organisms by homolog count
top_organisms = processed_df['organism.shortName'].value_counts().nlargest(15).index
filtered_df = processed_df[processed_df['organism.shortName'].isin(top_organisms)]

pivot_table = pivotmap(filtered_df)
print("\nPivot Table Head:")
print(pivot_table.head())

```

## API Overview

### Core Functions

- `initial_fetch(source_organism_name, transcript_names, subunit_dict, max_workers)`: Kicks off the homolog search with a defined set of genes.
- `subsequent_fetch(current_master_df, target_organism_name, max_workers)`: Expands the search by using the results from a previous fetch as input for a new target organism.
- run_homologs_pipeline(initial_organism, initial_genes_dict, subsequent_organisms, max_workers, checkpoint_dir): Run the full homolog search pipeline.
- run_expressions_workflow(master_file, checkpoint_dir, output_file, fetch_expression_data_for_gene_chunk, load_master_df, ...): Fetch and merge expression data for all subunits.

### Utility Functions

- `pivotmap(dataframe, index, columns, values)`: Generates a pivot table and a corresponding heatmap to visualize the count of homologs.
- `print_summary(df, stage_message)`: Prints a quick summary of a DataFrame's shape and contents.
- load_master_df(filepath): Load and validate the master homolog DataFrame.
- fetch_expression_data(gene_id_chunk, subunit_name_for_context, chunk_num, total_chunks): Fetch expression data for chunks of gene IDs.

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