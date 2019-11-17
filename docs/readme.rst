
Readme
======

|Build Status| |PyPI version| |PyPI pyversions| |PyPI license|

Create IGoR models and calculate the generation probability of V(D)J and CDR3 sequences.

Getting started
^^^^^^^^^^^^^^^

ImmunoProbs Python package can be used to create IGoR models and calculate the generation probability of V(D)J and CDR3 sequences. In order to get setup, please follow the steps below:

1. Install ImmunoProbs via :ref:`installation:PyPI`. Make sure to install all the additional :ref:`installation:Requirements`.


Additional information
~~~~~~~~~~~~~~~~~~~~~~

For more information you can have a look at the following wiki pages. These pages explain more in-depth on how to use ImmunoProbs and which features are supported:

-  Detailed look at the use cases for ImmunoProbs, containing the required and additional input parameters as well as the expected output. :ref:`use_cases:Use cases`

-  The command options that are available in ImmunoProbs for each of the tools and how to use them. :ref:`usage:Usage`

-  A description of the included pre-trained IGoR models in ImmunoProbs. :ref:`models:Models`

-  How to use ImmunoProbs docker image in combination with a galaxy server setup. :ref:`installation:Galaxy server`


Development
^^^^^^^^^^^

Development of ImmunoProbs is in an active state. If you would like to see new features, please open a new `issue <https://github.com/penuts7644/ImmunoProbs/issues/new>`__.

It is also possible to help out ImmunoProbs development by forking this project and creating a `pull request <https://github.com/penuts7644/ImmunoProbs/compare>`__ to submit new features.

When using a forked copy ImmunoProbs, make sure to have the correct Python version installed. Install the local Python development requirements using the make command:

.. code-block:: none

    make setup

.. |Build Status| image:: https://img.shields.io/travis/penuts7644/ImmunoProbs/master?style=for-the-badge
   :target: https://github.com/penuts7644/ImmunoProbs
.. |PyPI version| image:: https://img.shields.io/pypi/v/immuno-probs?style=for-the-badge
   :target: https://pypi.python.org/pypi/immuno-probs/
.. |PyPI pyversions| image:: https://img.shields.io/pypi/pyversions/immuno-probs?style=for-the-badge
   :target: https://pypi.python.org/pypi/immuno-probs/
.. |PyPI license| image:: https://img.shields.io/pypi/l/immuno-probs?style=for-the-badge
   :target: https://pypi.python.org/pypi/immuno-probs/
