# ts-qc
Utilities for evaluating large genetic genealogies (e.g. ARGs) stored in the [tskit](https://tskit.dev) "succinct tree sequence" format.
This is targetted at ARGs that have been inferred using a variety of inference tools such as `tsinfer`, `sc2ts`, `Relate`, `kwARG`, `ARGneedle`, etc. 


To launch the app use:

`python -m tsqc /path/to/trees-file`

On WSL, it may be necessary to disable Numba's CUDA support:

`NUMBA_DISABLE_CUDA=1 python -m tsqc /path/to/trees-file`

## Installation

Tsqc is currently in development. To install the latest dev version from github, try

```
python -m pip install git+https://github.com/tskit-dev/tsqc
```
