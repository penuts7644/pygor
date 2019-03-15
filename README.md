# ImmunoProbs

[![Build Status](https://img.shields.io/travis/penuts7644/ImmunoProbs.svg?branch=master&longCache=true&style=for-the-badge)](https://travis-ci.org/penuts7644/ImmunoProbs)
[![PyPI version shields.io](https://img.shields.io/pypi/v/immuno-probs.svg?longCache=true&style=for-the-badge)](https://pypi.python.org/pypi/immuno-probs/)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/immuno-probs.svg?longCache=true&style=for-the-badge)](https://pypi.python.org/pypi/immuno-probs/)
[![PyPI license](https://img.shields.io/pypi/l/immuno-probs.svg?longCache=true&style=for-the-badge)](https://pypi.python.org/pypi/immuno-probs/)

Create IGoR models and calculate the generation probability of V(D)J and CDR3 sequences.

For **installation**, **use cases** and a **tutorial**, please have a look at the [ImmunoProbs wiki page](https://github.com/penuts7644/ImmunoProbs/wiki).

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
├── data
│   ├── human_b_heavy
│   │   ├── J_gene_CDR3_anchors.csv
│   │   ├── V_gene_CDR3_anchors.csv
│   │   ├── genomic_D.fasta
│   │   ├── genomic_J.fasta
│   │   ├── genomic_V.fasta
│   │   ├── model_marginals.txt
│   │   └── model_params.txt
│   ├── human_t_alpha
│   │   ├── J_gene_CDR3_anchors.csv
│   │   ├── V_gene_CDR3_anchors.csv
│   │   ├── genomic_J.fasta
│   │   ├── genomic_V.fasta
│   │   ├── model_marginals.txt
│   │   └── model_params.txt
│   ├── human_t_beta
│   │   ├── J_gene_CDR3_anchors.csv
│   │   ├── V_gene_CDR3_anchors.csv
│   │   ├── genomic_D.fasta
│   │   ├── genomic_J.fasta
│   │   ├── genomic_V.fasta
│   │   ├── model_marginals.txt
│   │   └── model_params.txt
│   ├── mouse_t_beta
│   │   ├── J_gene_CDR3_anchors.csv
│   │   ├── V_gene_CDR3_anchors.csv
│   │   ├── genomic_D.fasta
│   │   ├── genomic_J.fasta
│   │   ├── genomic_V.fasta
│   │   ├── model_marginals.txt
│   │   └── model_params.txt
│   └── tutorial_model
│       ├── CDR3_sequences.csv
│       └── VDJ_sequences.csv
├── model
│   ├── __init__.py
│   ├── default_models.py
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
