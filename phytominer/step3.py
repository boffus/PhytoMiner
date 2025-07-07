from phytominer.workflow import run_expressions_workflow
from phytominer.processing import load_master_df, fetch_expression_data
from phytominer.config import (
    JOIN2_OUTPUT_FILE,
    EXPRESSION_CHECKPOINT_DIR,
    EXPRESSIONS_OUTPUT_FILE
)

if __name__ == '__main__':
    run_expressions_workflow(
        master_file=JOIN2_OUTPUT_FILE,
        checkpoint_dir=EXPRESSION_CHECKPOINT_DIR,
        output_file=EXPRESSIONS_OUTPUT_FILE,
        fetch_expression_data_for_gene_chunk=fetch_expression_data,
        load_master_df=load_master_df
    )