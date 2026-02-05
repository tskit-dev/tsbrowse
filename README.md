# tsbrowse
[![License](https://img.shields.io/github/license/tskit-dev/tsbrowse)](https://github.com/tskit-dev/tsbrowse/blob/main/LICENSE) [![PyPI version](https://img.shields.io/pypi/v/tsbrowse.svg)](https://pypi.org/project/tsbrowse/) [![Supported Python Versions](https://img.shields.io/pypi/pyversions/tsbrowse.svg)](https://pypi.org/project/tsbrowse/) [![Docs Build](https://github.com/tskit-dev/tsbrowse/actions/workflows/docs.yml/badge.svg)](https://github.com/tskit-dev/tsbrowse/actions/workflows/docs.yml) [![Tests](https://github.com/tskit-dev/tsbrowse/actions/workflows/tests.yml/badge.svg)](https://github.com/tskit-dev/tsbrowse/actions/workflows/tests.yml)
Inspect large genetic genealogies (e.g. ARGs) stored in the [tskit](https://tskit.dev) "succinct tree sequence" format, via a web app.
_Tsbrowse_ can scale to ARGs with millions of samples.

Documentation: https://tskit.dev/tsbrowse/docs/stable/ (latest: https://tskit.dev/tsbrowse/docs/latest/)

It is particularly useful to help evaluate ARGs that have been inferred using tools such as
[tsinfer](https://github.com/tskit-dev/tsinfer),
[sc2ts](https://github.com/jeromekelleher/sc2ts),
[Relate](https://github.com/MyersGroup/relate),
[KwARG](https://github.com/a-ignatieva/kwarg),
[Threads](https://pypi.org/project/threads-arg/), etc.

## Quickstart

## Installation

Install from PyPI:

```
python -m pip install tsbrowse
```

A tskit tree sequence must first be preprocessed to a `.tsbrowse` file:

`python -m tsbrowse preprocess /path/to/trees-file`

This will write a `.tsbrowse` file that can then be viewed using the web app:

`python -m tsbrowse serve /path/to/tsbrowse-file`

This command will launch a web server that can be accessed at `http://localhost:8080` in a web browser.

To display the genes track use:

`python -m tsbrowse serve /path/to/tsbrowse-file --annotations-file genes.csv`

(where `genes.csv` is a semicolon-separated text file containing a header and information about one gene per row in the order: `chr;start;end;strand;ensembl ID;gene name`)

An example tree sequence file can be found here: [example.trees](https://raw.githubusercontent.com/tskit-dev/tsbrowse/refs/heads/main/example/example.trees.tsz).

## Tips

If you are using Windows Subsystem for Linux (WSL) you may need to disable Numba's CUDA support:

`NUMBA_DISABLE_CUDA=1 python -m tsbrowse serve /path/to/tsbrowse-file`

A PNG of a specific page in tsbrowse can be generated using the `screenshot` command:

`python -m tsbrowse screenshot /path/to/tsbrowse-file mutations`


## Development

Tests are run with pytest:

`python -m pytest`

To run the UI tests so you can see what the browser is doing use

`python -m pytest --headed --slowmo 1000 tests/test_ui.py`

`playwright codegen` is also useful for writing UI test code.

## Citing tsbrowse
If you use `tsbrowse` in your work, please cite [this paper](<https://doi.org/10.1093/bioinformatics/btaf393>)
