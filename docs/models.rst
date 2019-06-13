
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

For tutorial purposes, there is the ``tutorial-model`` value. This model makes use of the **human TCR beta** model. Have a look at the :ref:`tutorial:Using pre-trained models` tutorial to see how to use one of the pre-trained models included with ImmunoProbs.

Custom
^^^^^^

You can build your own IGoR models using ImmunoProbs by following the :ref:`tutorial:Building your own model` tutorial. Instead of using the ``-model`` parameter, make sure to use ``-custom-model``.

Supported types
~~~~~~~~~~~~~~~

Currently there are 4 four supported model types (``-type``) for which models can be build: ``alpha``, ``beta``, ``light`` and ``heavy`` chain.
