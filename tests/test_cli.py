import os

import tszip
from click.testing import CliRunner

from . import test_preprocess
from tsbrowse import __main__ as main


def test_preprocess_cli(tmpdir):
    tszip_path = os.path.join(tmpdir, "test_input.tszip")
    default_output_path = os.path.join(tmpdir, "test_input.tsbrowse")
    custom_output_path = os.path.join(tmpdir, "custom_input.tsbrowse")

    ts = test_preprocess.single_tree_example_ts()
    tszip.compress(ts, tszip_path)

    runner = CliRunner()
    result = runner.invoke(main.cli, ["preprocess", tszip_path])
    assert result.exit_code == 0
    assert os.path.exists(default_output_path)
    tszip.load(default_output_path).tables.assert_equals(ts.tables)

    result = runner.invoke(
        main.cli, ["preprocess", tszip_path, "--output", custom_output_path]
    )
    assert result.exit_code == 0
    assert os.path.exists(custom_output_path)
    tszip.load(custom_output_path).tables.assert_equals(ts.tables)

    # TODO Load into model and check that the model is correct
