# ImmunoProbs

[![Build Status](https://img.shields.io/travis/penuts7644/ImmunoProbs.svg?branch=master&longCache=true&style=for-the-badge)](https://travis-ci.org/penuts7644/ImmunoProbs)
[![PyPI version shields.io](https://img.shields.io/pypi/v/immuno-probs.svg?longCache=true&style=for-the-badge)](https://pypi.python.org/pypi/immuno-probs/)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/immuno-probs.svg?longCache=true&style=for-the-badge)](https://pypi.python.org/pypi/immuno-probs/)
[![PyPI license](https://img.shields.io/pypi/l/immuno-probs.svg?longCache=true&style=for-the-badge)](https://pypi.python.org/pypi/immuno-probs/)

ImmunoProbs Python package uses a simplified manner for calculating the generation probability of V(D)J and CDR3 sequences.

### How to install

ImmunoProbs is installable via PyPI using the following terminal command `pip install immuno-probs`. To see all available tools for ImmunoProbs use the `immuno-probs -h` command after installation. Each ImmunoProbs tool has its own specific arguments. The options for a specific tools are available via `immuno-probs <TOOL NAME> -h`.

### Requirements

All Python dependencies that are used by this package are installed through pip upon installation of ImmunoProbs. However, some software (not available via pip) needs to installed manually when planning on using certain commandline tools from ImmunoProbs:

| Function | Requirement |
| -------- | ----------- |
| create-cdr3-anchors | This requires [MUSCLE](http://www.drive5.com/muscle/) to be installed on the computer. For macOS users, muscle can also be install via HomeBrew by tapping into `brewsci/bio` and installing MUSCLE via `brew install muscle` |
| create-igor-model | This will use Python's subprocess package to execute commands for the [IGoR](https://github.com/qmarcou/IGoR) package. For this tool to work properly, make sure that you have at least installed IGoR 1.3.0 using the installation steps from [IGoR's documentation](https://qmarcou.github.io/IGoR/#install). |

### Development

Development of ImmunoProbs is currently active. If you would like to see new features, please open a new [issue](https://github.com/penuts7644/ImmunoProbs/issues/new).

It is also possible to help out by developing your own new features by forking this project and creating a [pull request](https://github.com/penuts7644/ImmunoProbs/compare).

When using a forked copy ImmunoProbs, make sure to have the correct Python version installed. Install the local Python development requirements via `make setup` command.

### Package structure

```
immuno_probs
├── __init__.py
├── __main__.py
├── alignment
│   ├── __init__.py
│   ├── align_read.py
│   └── muscle_aligner.py
├── cdr3
│   ├── __init__.py
│   └── anchor_locator.py
├── cli
│   ├── __init__.py
│   └── create_cdr3_anchors.py
└── util
    ├── __init__.py
    ├── cli.py
    ├── constant.py
    ├── conversion.py
    ├── exception.py
    ├── io.py
    └── processing.py
```
