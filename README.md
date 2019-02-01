# ImmunoProbs

[![Build Status](https://img.shields.io/travis/penuts7644/ImmunoProbs.svg?branch=master&longCache=true&style=for-the-badge)](https://travis-ci.org/penuts7644/ImmunoProbs)
[![PyPI version shields.io](https://img.shields.io/pypi/v/immuno-probs.svg?longCache=true&style=for-the-badge)](https://pypi.python.org/pypi/immuno-probs/)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/immuno-probs.svg?longCache=true&style=for-the-badge)](https://pypi.python.org/pypi/immuno-probs/)
[![PyPI license](https://img.shields.io/pypi/l/immuno-probs.svg?longCache=true&style=for-the-badge)](https://pypi.python.org/pypi/immuno-probs/)

ImmunoProbs Python package able to calculate the generation probability of V(D)J and CDR3 sequences.

### How to install

ImmunoProbs is installable via PyPI using the following terminal command:

```
pip install immuno-probs
```

In order to see all available tools/functionalities for ImmunoProbs use the help command after installation:

```
immuno-probs -h
```

Each ImmunoProbs tool has its own specific arguments (accompanying the top-level ones). The options for a specific tools are available via the tool's help command:

```
immuno-probs <TOOL NAME> -h
```

Make sure to install the necessary requirements (see the requirements section).

##### Using ImmunoProbs via docker image

It is also possible to use a docker image of ImmunoProbs with all necessary requirement pre-installed in an ubuntu environment. Make sure to install [docker](https://www.docker.com) first and pull the most recent version of the image with the following command:

```
docker pull penuts7644/immuno-probs:0.1.4
```

Now you can use the ImmunoProbs image however you like.

You can also execute the ImmunoProbs image while using your local machine's working directory:

```
docker run --rm --volume "$PWD":/tmp penuts7644/immuno-probs:0.1.4 <TOOL NAME>
```

Note that the container is removed after execution and that you'll need to link to the `/tmp` directory of the container since that is the location where output files are created.

### Requirements

All Python dependencies that are used by this package are installed through pip upon installation of ImmunoProbs. However, some software (not available via pip) needs to installed manually when planning on using certain commandline tools from ImmunoProbs:

| Command | Requirement |
| ------- | ----------- |
| `locate-cdr3-anchors` | This requires [MUSCLE](http://www.drive5.com/muscle/) to be installed. For linux, MUSCLE can be installed via `apt-get install muscle`. For macOS with HomeBrew, tap into `brewsci/bio` and install MUSCLE via `brew install muscle` |
| `build-igor-model` `generate_seqs` `evaluate_seqs` | This will use Python's subprocess package to pass the user arguments to [IGoR](https://github.com/qmarcou/IGoR). For these tools to work properly, make sure that you have at least compiled and installed IGoR version 1.3.0 using the guide in [IGoR's documentation](https://qmarcou.github.io/IGoR/#install). |

### Development

Development of ImmunoProbs is in an active state. If you would like to see new features, please open a new [issue](https://github.com/penuts7644/ImmunoProbs/issues/new).

It is also possible to help out ImmunoProbs development by forking this project and creating a [pull request](https://github.com/penuts7644/ImmunoProbs/compare) to submit new features.

When using a forked copy ImmunoProbs, make sure to have the correct Python version installed. Install the local Python development requirements using the make command:

```
make setup
```

### Package structure

```
immuno_probs
├── __init__.py
├── alignment
│   ├── __init__.py
│   ├── align_read.py
│   └── muscle_aligner.py
├── cdr3
│   ├── __init__.py
│   ├── anchor_locator.py
│   └── olga_container.py
├── cli
│   ├── __init__.py
│   ├── __main__.py
│   ├── build_igor_model.py
│   ├── evaluate_seqs.py
│   ├── generate_seqs.py
│   └── locate_cdr3_anchors.py
├── model
│   ├── __init__.py
│   ├── igor_interface.py
│   └── igor_loader.py
└── util
    ├── __init__.py
    ├── cli.py
    ├── constant.py
    ├── conversion.py
    ├── exception.py
    ├── io.py
    └── processing.py
```
