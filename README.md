# tsbrowse
Inspect large genetic genealogies (e.g. ARGs) stored in the [tskit](https://tskit.dev) "succinct tree sequence" format,
via a genome browser style app. _Tsbrowse_ can deal with ARGs of thousands or potentially millions of samples.
It is particularly useful to help evaluate ARGs that have been inferred using tools such as
[tsinfer](https://github.com/tskit-dev/tsinfer),
[sc2ts](https://github.com/tskit-dev/sc2ts),
[Relate](https://github.com/MyersGroup/relate),
[KwARG](https://github.com/a-ignatieva/kwarg),
[Threads](https://pypi.org/project/threads-arg/), etc.

To launch the app use:

`python -m tsbrowse /path/to/trees-file`

On WSL, it may be necessary to disable Numba's CUDA support:

`NUMBA_DISABLE_CUDA=1 python -m tsbrowse /path/to/trees-file`

## Installation

tsbrowse is currently in development. To install the latest dev version from github, try

```
python -m pip install git+https://github.com/tskit-dev/tsbrowse
```
