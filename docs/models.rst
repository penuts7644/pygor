
Models
======

Pre-trained
^^^^^^^^^^^

A few pre-trained models come pre-installed with the ImmunoProbs package. These models can be selected by providing the ``-model`` parameter in combination with the model name.

Currently the following pre-trained IGoR models can be selected with the ``-model`` parameter:

-  Human TCR alpha model (``human-t-alpha``)
-  Human TCR beta model (``human-t-beta``)
-  Human IG heavy model (``human-b-heavy``)
-  Mouse TCR beta model (``mouse-t-beta``)

Custom
^^^^^^

You can also build your own IGoR models using ImmunoProbs. Instead of using the ``-model`` parameter, make sure to use ``-custom-model``.

Supported types
~~~~~~~~~~~~~~~

Currently there are 4 four supported model types (``-type``) for which models can be build: ``alpha``, ``beta``, ``light`` and ``heavy`` chain.
