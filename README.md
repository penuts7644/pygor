# pygor

[![Build Status](https://img.shields.io/travis/penuts7644/pygor.svg?branch=master&longCache=true&style=for-the-badge)](https://travis-ci.org/penuts7644/pygor)
[![PyPI version shields.io](https://img.shields.io/pypi/v/pygor.svg?longCache=true&style=for-the-badge)](https://pypi.python.org/pypi/pygor/)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/pygor.svg?longCache=true&style=for-the-badge)](https://pypi.python.org/pypi/pygor/)
[![PyPI license](https://img.shields.io/pypi/l/pygor.svg?longCache=true&style=for-the-badge)](https://pypi.python.org/pypi/pygor/)

Pygor is part of the IGoR (Inference and Generation of Repertoires) software, which can be accessed via Github: [https://github.com/qmarcou/IGoR](https://github.com/qmarcou/IGoR).

Pygor Python package can be used to post process files generated by IGoR.

### How to install

Pygor is installable via PyPI using the following terminal command `pip install pygor`.

### Requirements

All Python dependencies that are used by this package are installed when installing pygor via the pip command. However, some software needs to installed manually when planning on using specific functionalities within pygor:

| Function | Requirement |
| -------- | ----------- |
| Creating CDR3 anchors (`cdr3` module) | This requires [MUSCLE](http://www.drive5.com/muscle/) to be installed on the computer. For macOS users, muscle can also be install via HomeBrew by tapping into `brewsci/bio` and installing MUSCLE via `brew install muscle` |

### Development

Development of pygor is active. If you would like to see new features, please open a new [issue](https://github.com/penuts7644/pygor/issues/new).

It is also possible to help out by developing your own new features by forking this project and creating a [pull request](https://github.com/penuts7644/pygor/compare).

When using a forked copy pygor, make sure to have the correct Python version installed. Install the local Python development requirements via `make setup` command.

### Package structure

```
pygor
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
├── counter
│   ├── __init__.py
│   ├── bestscenarios
│   │   ├── __init__.py
│   │   └── bestscenarios.py
│   └── coverage
│       ├── __init__.py
│       └── coverage.py
├── model
│   ├── __init__.py
│   ├── entropy.py
│   ├── genmodel.py
│   └── hypermglobal.py
└── util
    ├── __init__.py
    ├── cli.py
    ├── constant.py
    ├── conversion.py
    ├── exception.py
    ├── io.py
    └── processing.py
```
