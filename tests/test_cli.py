import os

import tszip
from click.testing import CliRunner
from PIL import Image

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


def test_screenshot_cli(tmpdir):
    tszip_path = os.path.join(tmpdir, "test_input.tszip")
    tsbrowse_path = os.path.join(tmpdir, "test_input.tsbrowse")
    ts = test_preprocess.single_tree_example_ts()
    tszip.compress(ts, tszip_path)
    runner = CliRunner()
    result = runner.invoke(main.cli, ["preprocess", tszip_path])
    assert result.exit_code == 0

    # Test with default output
    result = runner.invoke(main.cli, ["screenshot", tsbrowse_path, "overview"])
    assert result.exit_code == 0
    default_output = os.path.join(tmpdir, "test_input_overview.png")
    assert os.path.exists(default_output)
    with Image.open(default_output) as img:
        width, height = img.size
        assert width > 1000
        assert height == 347

    # Test with path
    custom_output = os.path.join(tmpdir, "custom_screenshot.png")
    result = runner.invoke(
        main.cli,
        [
            "screenshot",
            tsbrowse_path,
            "overview",
            "--output",
            custom_output,
        ],
    )
    assert result.exit_code == 0
    assert os.path.exists(custom_output)

    # Test with invalid page
    result = runner.invoke(main.cli, ["screenshot", tsbrowse_path, "InvalidPage"])
    assert result.exit_code != 0
    assert "Invalid value" in result.output
