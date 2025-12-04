"""
This runs the burden association tests using the LoadModule class. Note that if you are running REGENIE step1 for
these tests, they will take a while. Once you have run step1 once, you can save the output files to save time
on running step2. Step1 files are also provided in case you want to copy them into the test directory and avoid
running step1.

Note that BOLT will always take a long time to run.

For access to the test data, please contact one of the authors of the burden package.
"""

import shutil
from pathlib import Path

import pandas as pd
import pytest

from extract.loader import LoadModule

# Set this flag to True if you want to keep (copy) the temporary output files
KEEP_TEMP = True

# test data directory
test_data_dir = Path(__file__).parent / 'test_data'


@pytest.fixture
def temporary_path(tmp_path, monkeypatch):
    """
    Prepare a temporary working directory that contains a copy of the test_data
    directory, then change the working directory to it.

    If KEEP_TEMP is True, after the test the entire temporary directory will be copied
    to a folder 'temp_test_outputs' in the project root.
    """
    # Determine where the original test_data directory is located.
    # (Assumes it is at <project_root>/test_data)
    test_data_source = Path(__file__).parent / "test_data"

    # Create the destination folder inside the tmp_path.
    destination = tmp_path / "test_data"
    destination.parent.mkdir(parents=True, exist_ok=True)

    # Copy the entire test_data directory into the temporary directory.
    shutil.copytree(test_data_source, destination)

    # Change the current working directory to the temporary directory.
    monkeypatch.chdir(tmp_path)

    # Yield the temporary directory to the test.
    yield tmp_path

    # After the test, if KEEP_TEMP is True, copy the temporary directory to a persistent location.
    if KEEP_TEMP:
        persistent_dir = Path(__file__).parent / "temp_test_outputs" / tmp_path.name
        persistent_dir.parent.mkdir(exist_ok=True)
        shutil.copytree(tmp_path, persistent_dir, dirs_exist_ok=True)
        print(f"Temporary output files have been copied to: {persistent_dir}")


@pytest.mark.parametrize("input_args, filename, expected_output", [
    (
            (
                    f"--association_tarballs {test_data_dir}/HC_PTV-MAF_001.tar.gz "
                    f"--bgen_index {test_data_dir}/bgen_locs.tsv "
                    f"--sparse_grm {test_data_dir}/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx "
                    f"--sparse_grm_sample {test_data_dir}/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx.sampleIDs.txt "
                    f"--gene_ids OR4F5 MATN1 "
                    f"--phenofile {test_data_dir}/phenotype.tsv "
                    f"--transcript_index {test_data_dir}/transcripts.tsv.gz "
                    f"--base_covariates {test_data_dir}/base_covariates.covariates "
            ),
            "test.genes.STAAR_glm.stats.tsv.gz",
            Path(__file__).parent / "expected_results/test.genes.STAAR_glm.stats.tsv.gz"
    ),
])
def test_load_module_run(input_args, filename, expected_output, temporary_path):
    """
    Test the LoadModule class: initialization, running start_module, and outputs.
    """

    loader = LoadModule(output_prefix="test", input_args=input_args)

    # Run the module
    loader.start_module()

    # Get outputs set by the extract tool
    outputs = loader.get_outputs()

    # Basic checks on outputs
    assert outputs is not None, "Outputs should not be None"

    # assert the output
    df_results = pd.read_csv(f'{filename}', sep='\t', compression='gzip')
    df_expected = pd.read_csv(expected_output, sep='\t', compression='gzip')

    # Basic shape/column checks (optional but nice)
    assert df_results.shape == df_expected.shape
    assert list(df_results.columns) == list(df_expected.columns)

    # Sort rows so row order doesn't matter
    if "variant" in df_results.columns:
        sort_cols = ["variant"]
    else:
        # fall back to sorting by all columns if you want it general
        sort_cols = list(df_results.columns)

    df_results_sorted = df_results.sort_values(sort_cols).reset_index(drop=True)
    df_expected_sorted = df_expected.sort_values(sort_cols).reset_index(drop=True)

    pd.testing.assert_frame_equal(df_results_sorted, df_expected_sorted)