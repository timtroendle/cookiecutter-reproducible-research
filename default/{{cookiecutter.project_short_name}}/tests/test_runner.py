import sys
from pathlib import Path

import pandas as pd
import pytest


def run_test(snakemake):
    exit_code = pytest.main(
        [
            snakemake.input.test_dir,
            f"--html={snakemake.log[0]}",
            "--self-contained-html",
            "--verbose"
        ],
        plugins=[
            _create_config_plugin(snakemake=snakemake)
        ]
    )
    if exit_code == 0:
        Path(snakemake.output[0]).touch()
    sys.exit(exit_code)


def _create_config_plugin(snakemake):
    """Creates fixtures from Snakemake configuration."""

    class SnakemakeConfigPlugin():

        @pytest.fixture()
        def model_results(self):
            return pd.read_pickle(snakemake.input.model_results)

    return SnakemakeConfigPlugin()


if __name__ == "__main__":
    run_test(snakemake)
