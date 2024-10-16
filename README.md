# tsbrowse
Inspect large genetic genealogies (e.g. ARGs) stored in the [tskit](https://tskit.dev) "succinct tree sequence" format,
via a genome browser style app. _Tsbrowse_ can deal with ARGs of thousands or potentially millions of samples.
It is particularly useful to help evaluate ARGs that have been inferred using tools such as
[tsinfer](https://github.com/tskit-dev/tsinfer),
[sc2ts](https://github.com/tskit-dev/sc2ts),
[Relate](https://github.com/MyersGroup/relate),
[KwARG](https://github.com/a-ignatieva/kwarg),
[Threads](https://pypi.org/project/threads-arg/), etc.

To view a tskit tree sequence or tszip file first pre-process it:

`python -m tsbrowse preprocess /path/to/trees-file`

This will write a `.tsbrowse` file

To launch the app use:

`python -m tsbrowse serve /path/to/tsbrowse-file`

Or to generate a PNG of a specific page use, e.g:

`python -m tsbrowse screenshot /path/to/tsbrowse-file mutations`

On WSL, it may be necessary to disable Numba's CUDA support:

`NUMBA_DISABLE_CUDA=1 python -m tsbrowse serve /path/to/tsbrowse-file`

## Installation

tsbrowse is currently in development. To install the latest dev version from github, try

```
python -m pip install git+https://github.com/tskit-dev/tsbrowse
```

## Development

Test are run with pytest:

`python -m pytest`

To run the UI tests so you can see what the browser is doing use

`python -m pytest --headed --slowmo 1000 tests/test_ui.py`

`playwright codegen` is also useful for writing UI test code.
