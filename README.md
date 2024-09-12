# tsbrowse

Inspect large genetic genealogies (e.g. ARGs) stored in the [tskit](https://tskit.dev) "succinct tree sequence" format,
via a genome browser style app. This can deal with ARGs of thousands or potentially millions of samples.
It is particularly useful to help evaluating ARGs that have been inferred using inference tools such as
`tsinfer`, `sc2ts`, `Relate`, `KwARG`, `Threads`, etc.


To launch the app use:

`python -m tsbrowse /path/to/trees-file`

On WSL, it may be necessary to disable Numba's CUDA support:

`NUMBA_DISABLE_CUDA=1 python -m tsbrowse /path/to/trees-file`

## Installation

tsbrowse is currently in development. To install the latest dev version from github, try

```
python -m pip install git+https://github.com/tskit-dev/tsbrowse
```
