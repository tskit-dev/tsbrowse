(contributing)=

# Contributing to tsbrowse

All contributions, bug reports, documentation improvements and ideas are welcome. If you think
there is anything missing, please open an [issue](https://github.com/tskit-dev/tsbrowse/issues)
or [pull request](https://github.com/tskit-dev/tsbrowse/pulls) on Github.

## Quick start
### Using a virtual environment
First, set up a virtual environment:
```
$ python3 -m venv tsbrowse_env
```
Then, activate the virtual environment:
```
$ source tsbrowse_env/bin/activate
```
Within the virtual environment, install requirements:
```
$ python3 -m pip install -r requirements.txt
```

### Github workflow
1.    Fork the tsbrowse github repository at https://github.com/tskit-dev/tsbrowse/ to create your own copy of the project.

2.    Clone a local copy:
```
$ git clone https://github.com/your-user-name/tsbrowse.git
$ cd tsbrowse
$ git remote add upstream https://github.com/tskit-dev/tsbrowse.git
```
This creates the directory tsbrowse and connects your repository to the upstream (main project) tsbrowse repository.
3.    Install pre-commit hooks. These will be executed automatically on every commit, to ensure the code meets development requirements.
```
$ pre-commit install
```
4.    Run the test suite and make sure tests pass:
```
$ python3 -m pytest
```
5.    Create a topic branch and add commits:
```
$ git fetch upstream
$ git checkout upstream/main
$ git checkout -b topic_branch_name
```
6.    When you are ready to share your contribution, push changes to your fork:
```
$ git push origin branch-name
```
7.    Create a [pull request](https://help.github.com/articles/about-pull-requests/) on GitHub.

### Documentation
Tsbrowse documentation is written in markdown. These files are contained in the docs directory. To build the documentation, run `make` in the `docs` directory.
