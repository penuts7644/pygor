
Tutorial
========

Getting started
^^^^^^^^^^^^^^^

Make sure to have `Docker <https://www.docker.com>`__ installed and have pulled the latest version of the ImmunoProbs docker image (see the :ref:`installation:Docker` installation guide). All commands below will be run inside the docker image container. Log in to the container's terminal with the following command:

.. code-block:: none

    docker run --rm -it penuts7644/immuno-probs bash

You can exit the terminal inside the container with ``exit``. The container is removed after exiting because of the ``--rm`` flag.

Using pre-trained models
^^^^^^^^^^^^^^^^^^^^^^^^

In the first section of the tutorial, we will perform the pre-trained model :ref:`use_cases:Use cases` use case for a pre-trained IGoR model, included with ImmunoProbs. For this tutorial we will use the `tutorial and default models <https://github.com/penuts7644/ImmunoProbs/tree/master/immuno_probs/data>`__ located in the GitHub project. These files come pre-installed with ImmunoProbs package.

The included IGoR models can be used to generate and evaluate V(D)J or CDR3 sequences. A list of supported models can be found in the :ref:`models:Pre-trained` models section. For the rest of this tutorial, we’ll be using the special tutorial model (``tutorial-model``) based on the **human TCR beta** model.

Generate V(D)J sequences
~~~~~~~~~~~~~~~~~~~~~~~~

Generating sequences with a predefined IGoR model can be done by specifying the model you would like to use (``-model``) and the number of sequences to generate (``-n-gen``).

.. code-block:: none

    immuno-probs \
      generate \
        -model tutorial-model \
        -n-gen 100

Calculate the generation probabilities for V(D)J sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We can calculate the generation probability of our VDJ sequences by specifying the input data through the ``-seqs`` command.

.. code-block:: none

    immuno-probs \
      evaluate \
        -model tutorial-model \
        -seqs VDJ_sequences.tsv

Generate CDR3 sequences
~~~~~~~~~~~~~~~~~~~~~~~

Since the CDR3 sequences generation and evaluation requires additional V and J anchor index positions to locate the CDR3 sequences, they require two extra files. However, these files are also included with ImmunoProbs's models, making the command similar to the one in section **a**. We just add the ``-cdr3`` flag to get CDR3 sequences instead of VDJ sequences.

.. code-block:: none

    immuno-probs \
      generate \
        -model tutorial-model \
        -n-gen 100 \
        -cdr3

Calculate the generation probabilities for CDR3 sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Similar to step **b**, but now we are using CDR3 input sequences so we’ll have to set the CDR3 flag.

.. code-block:: none

    immuno-probs \
      evaluate \
        -model tutorial-model \
        -seqs CDR3_sequences.tsv \
        -cdr3

Building your own model
^^^^^^^^^^^^^^^^^^^^^^^

In the second part of the ImmunoProbs tutorial, we will perform the general and custom model :ref:`use_cases:Use cases` to create our own IGoR model. For this tutorial we will use the `zipped data files <https://github.com/penuts7644/ImmunoProbs/tree/master/tutorial_data.zip>`__ located in the GitHub project to create a human VDJ beta chain model. Finally, we are going to generate and evaluate sequences using our created model. The zipped data files are included in the ``tutorial_data`` directory in the root of the ImmunoProbs docker image.

Building a model
~~~~~~~~~~~~~~~~

We'll start by specifying the reference genomic template FASTA files (``-ref``) for the V, D and J gene as well as some input sequences (``-seqs``). We specify the number of training rounds (``-n-iter``) and the desired type of the model we would like to build (``-type``). This might take a while depending on your system configuration.

.. code-block:: none

    immuno-probs \
      build \
        -ref V /tutorial_data/TRBV.fasta \
        -ref D /tutorial_data/TRBD.fasta \
        -ref J /tutorial_data/TRBJ.fasta \
        -seqs /tutorial_data/1000_sample_seqs.tsv \
        -n-iter 10 \
        -type beta

Locate CDR3 anchors positions for CDR3 sequence generation and evaluation steps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CDR3 anchor positions are required in order to accurately generate and evaluate CDR3 sequences. Specify the V and J genomic reference files with the ``-ref`` option.

.. code-block:: none

    immuno-probs \
      locate \
        -ref V /tutorial_data/TRBV.fasta \
        -ref J /tutorial_data/TRBJ.fasta

Generate VJ, VDJ or CDR3 sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We need to specify our model marginals and parameters files as well as the model type (``-type``).

.. code-block:: none

    immuno-probs \
      generate \
        -custom-model /tutorial_data/model_params.txt /tutorial_data/model_marginals.txt \
        -n-gen 100 \
        -type beta

To generate some CDR3 sequences, we'll add the ``-cdr3`` flag at the end of the command and specify the anchor position files created in section **b** with ``-anchor``.

.. code-block:: none

    immuno-probs \
      generate \
        -custom-model /tutorial_data/model_params.txt /tutorial_data/model_marginals.txt \
        -n-gen 100 \
        -type beta \
        -cdr3 \
        -anchor V /tutorial_data/V_gene_CDR3_anchors.tsv \
        -anchor J /tutorial_data/J_gene_CDR3_anchors.tsv

Calculate the generation probabilities for VJ, VDJ or CDR3 sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We are selecting the sequences generated in the previous step (``-seqs``), the model parameters and marginals (``-custom-model``), the type of the input model and the genomic templates (``-ref``)

.. code-block:: none

    immuno-probs \
      evaluate \
        -custom-model /tutorial_data/model_params.txt /tutorial_data/model_marginals.txt \
        -seqs /tutorial_data/generated_seqs_beta.tsv \
        -type beta \
        -ref V /tutorial_data/TRBV.fasta \
        -ref D /tutorial_data/TRBD.fasta \
        -ref J /tutorial_data/TRBJ.fasta

To evaluate CDR3 sequences generated in the previous section, we'll add the ``-cdr3`` flag at the end the command and replace the ``-ref`` options with ``-anchor``. For CDR3 we don't need genomic templates.

.. code-block:: none

    immuno-probs \
      evaluate \
        -custom-model /tutorial_data/model_params.txt /tutorial_data/model_marginals.txt \
        -seqs /tutorial_data/generated_seqs_beta_CDR3.tsv \
        -type beta \
        -cdr3 \
        -anchor V /tutorial_data/V_gene_CDR3_anchors.tsv \
        -anchor J /tutorial_data/J_gene_CDR3_anchors.tsv
