# ImmunoProbs

[![Build Status](https://img.shields.io/travis/penuts7644/ImmunoProbs.svg?branch=master&longCache=true&style=for-the-badge)](https://travis-ci.org/penuts7644/ImmunoProbs)
[![PyPI version shields.io](https://img.shields.io/pypi/v/immuno-probs.svg?longCache=true&style=for-the-badge)](https://pypi.python.org/pypi/immuno-probs/)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/immuno-probs.svg?longCache=true&style=for-the-badge)](https://pypi.python.org/pypi/immuno-probs/)
[![PyPI license](https://img.shields.io/pypi/l/immuno-probs.svg?longCache=true&style=for-the-badge)](https://pypi.python.org/pypi/immuno-probs/)

ImmunoProbs Python package uses simplified manner for calculating the generation probability of V(D)J and CDR3 sequences.

### How to install

ImmunoProbs is installable via PyPI using the following terminal command `pip install immuno-probs`.

### Requirements

All Python dependencies that are used by this package are installed when installing ImmunoProbs via the pip command. However, some software needs to installed manually when planning on using specific functionalities within ImmunoProbs:

| Function | Requirement |
| -------- | ----------- |
| Creating CDR3 anchors (`cdr3` module) | This requires [MUSCLE](http://www.drive5.com/muscle/) to be installed on the computer. For macOS users, muscle can also be install via HomeBrew by tapping into `brewsci/bio` and installing MUSCLE via `brew install muscle` |

### Development

Development of ImmunoProbs is active. If you would like to see new features, please open a new [issue](https://github.com/penuts7644/ImmunoProbs/issues/new).

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
