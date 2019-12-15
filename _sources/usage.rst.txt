
Usage
=====

Tools
^^^^^

ImmunoProbs has four tools: for converting adaptive sequence files, building IGoR models, locating the CDR3 anchor positions in sequences, generating and evaluating sequences. For the generate and evaluate options, ImmunoProbs can be used with one of the four pre-trained IGoR models that are included in the package. These models are described in the :ref:`models:Pre-trained` models section.

ImmunoProbs has a number of global options that are used throughout the other ImmunoProbs tools. The command is followed by the tool name and its respective options.

.. code-block:: none

    immuno-probs \
      -separator <SEPARATOR CHARS> \
      -threads <NUM THREADS> \
      -set-wd <DIRECTORY> \
      -out-name <OUTPUT NAME>
      -config-file <CONFIG FILE>
        [TOOL NAME] \
          <TOOL OPTIONS>

Converting adaptive file format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can convert adaptive sequence data to ImmunoProbs's standard by specifying the reference germline template FASTA files (``ref``) for the V and J gene as well as the adaptive input sequence data to be converted (``seqs``). The template files can be downloaded from `IMGT <http://www.imgt.org/vquest/refseqh.html>`__.
There are 4 output data files for 1 run of the tool: full length productive VDJ sequences, full length unproductive VDJ sequences, the total full length VDJ sequences and finally one with CDR3 sequences. The full length sequences can be used for training of new models and the CDR3 extract can be evaluated.
Some assumptions are made when converting the resolved genes. If ``use-allele`` is True, the allele from the input value is used, if this can't be found, the default allele in the config is used instead. The resolved column always requires a family to be found. If there is no gene, only the allele is inserted. If there is a gene, the family, gene and allele values are recombined. In addition if that gene value equals 1, an extra pattern is attached with only the family and allele.

.. code-block:: none

    immuno-probs \
      convert \
        -ref <GENE> <FASTA> \
        -seqs <SEPARATED>

You can convert a random subset of sequences from the given input file. Just specify the ``n-random`` flag, followed by a number. You can also use ``use-allele`` flag to use the allele information from the input data file during conversion. You might want to have a look at the :ref:`usage:Configuration file setup` for ImmunoProbs in order to specify you own file column names.

Building a model
~~~~~~~~~~~~~~~~

You can create your own IGoR model by specifying the reference IMGT germline template files (``ref``) for the V, D and J gene as well as some input sequences for the training of the model (``seqs``). The template files can be downloaded from `IMGT <http://www.imgt.org/vquest/refseqh.html>`__. Finally, specify the number of training rounds (``n-iter``) and the desired type of the model you would like to build (``type``).

.. code-block:: none

    immuno-probs \
      build \
        -ref <GENE> <FASTA> \
        -seqs <SEPARATED/FASTA> \
        -n-iter <NUM ITERATIONS> \
        -type <MODEL TYPE>

Locate CDR3 anchors positions for CDR3 sequence generation and evaluation steps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CDR3 anchor positions are required in order to accurately generate and evaluate CDR3 sequences. Specify the V and J germline reference files with the ``ref`` option. You can download the templates from `IMGT <http://www.imgt.org/vquest/refseqh.html>`__. Optionally, you could specify multiple motif parameters (``motif``).

.. code-block:: none

    immuno-probs \
      locate \
        -ref <GENE> <FASTA>
        -motif <MOTIF>

Generate VJ, VDJ or CDR3 sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Generation of sequences can be done with either an included model (``model``) or by selecting your own model marginals and parameters files (``custom-model``).

Generating sequences with a predefined IGoR model can be done by specifying the model you would like to use in combination with the number of sequences to generate (``n-gen``).

.. code-block:: none

    immuno-probs \
      generate \
        -model <MODEL NAME> \
        -n-gen <NUM SEQUENCES>

When using one of your now IGoR models, you'll have to specify the model with parameters and marginals and the type of the input model (``type``).

.. code-block:: none

    immuno-probs \
      generate \
        -custom-model <PARAMETERS> <MARGINALS> \
        -n-gen <NUM SEQUENCES> \
        -type <MODEL TYPE>

Both of the scenarios above will generate VJ or VDJ sequences. If you rather want CDR3 sequences, you'll need to add the ``cdr3`` flag at the end of either of the commands. When using a custom model, you also want to specify the anchor position files created in section **b** by adding: ``anchor <GENE> <SEPARATED>``.

Calculate the generation probabilities for VJ, VDJ or CDR3 sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The evaluation of sequences can be done with either an included model (``model``) or by selecting your own model marginals and parameters files (``custom-model``).

With the included models, we calculate the generation probability by specifying the sequences through the ``seqs`` command as well as the model files.

.. code-block:: none

    immuno-probs \
      evaluate \
        -model <MODEL NAME> \
        -seqs <SEPARATED/FASTA>

With a custom model: select the sequences (``seqs``), the model parameters and marginals (``custom-model``), the type of the input model and the germline templates (``ref``) from `IMGT <http://www.imgt.org/vquest/refseqh.html>`__.

.. code-block:: none

    immuno-probs \
      evaluate \
        -custom-model <PARAMETERS> <MARGINALS> \
        -seqs <SEPARATED/FASTA> \
        -ref <GENE> <FASTA> \
        -type <MODEL TYPE>

Both of the scenarios above can be used for evaluating VJ or VDJ sequences. If your input data consists of CDR3 sequences, you'll need to add the ``cdr3`` flag at the end of either of the commands. You can also use ``use-allele`` flag to use allele information from the input data to calculate the generation probability. When using a custom model, you also want to replace the ``ref`` command with ``anchor <GENE> <SEPARATED>``. Note that for CDR3, we don't need germline templates.

Parameters
^^^^^^^^^^

+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
| Tool         | Option                | Description                                                                                                                                                                       | Default                                                                                  | Required                                         |
+==============+=======================+===================================================================================================================================================================================+==========================================================================================+==================================================+
|              | ``separator``         | The separator character used for input files and for writing new files.                                                                                                           | Tab character (``\t``)                                                                   |                                                  |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
|              | ``threads``           | The number of threads the program is allowed to use.                                                                                                                              | Max available threads in system                                                          |                                                  |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
|              | ``set-wd``            | An optional location for writing files.                                                                                                                                           | The current working directory                                                            |                                                  |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
|              | ``out-name``          | An optional output file name. If multiple files are created, the value is used as a prefix for the file.                                                                          |                                                                                          |                                                  |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
|              | ``config-file``       | An optional configuration file path for ImmunoProbs. This file is combined with the default configuration to make up missing values.                                              |                                                                                          |                                                  |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
| ``convert``  | ``ref``               | A gene (V or J) followed by a reference genome FASTA file. Note: the FASTA reference genome files needs to conform to IGMT annotation (separated by vertical bar character).      |                                                                                          | Yes                                              |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
| ``convert``  | ``seqs``              | An input FASTA or separated data file with sequences for training the model.                                                                                                      |                                                                                          | Yes                                              |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
| ``convert``  | ``n-random``          | The number of random sequences to convert from the input adaptive data file (only if higher than 0).                                                                              | 0                                                                                        |                                                  |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
| ``convert``  | ``use-allele``        | If specified, the allele information from the gene resolved fields is used into the converted output file.                                                                        | Allele ``01`` is used for each gene.                                                     |                                                  |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
| ``build``    | ``ref``               | A gene (V, D or J) followed by a reference genome FASTA file. Note: the FASTA reference genome files needs to conform to IGMT annotation (separated by vertical bar character).   |                                                                                          | Yes                                              |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
| ``build``    | ``seqs``              | An input FASTA or separated data file with sequences for training the model.                                                                                                      |                                                                                          | Yes                                              |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
| ``build``    | ``n-iter``            | The number of inference iterations to perform when creating the model.                                                                                                            | 1                                                                                        |                                                  |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
| ``build``    | ``type``              | The type of model to create. (select one: ``alpha``, ``beta``, ``light`` or ``heavy``.                                                                                            |                                                                                          | Yes                                              |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
| ``locate``   | ``ref``               | A gene (V or J) followed by a reference genome FASTA file. Note: the FASTA reference genome files needs to conform to IGMT annotation (separated by vertical bar character).      |                                                                                          | Yes                                              |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
| ``locate``   | ``motif``             | The motif to look for. Can be used multiple times.                                                                                                                                | ``V`` (Cystein - TGT and TGC) or ``J`` (Tryptophan - TGG, Phenylalanine - TTC and TTT)   |                                                  |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
| ``generate`` | ``model``             | Specify a pre-installed model for generation. (select one: ``human-t-alpha``, ``human-t-beta``, ``human-b-heavy`` or ``mouse-t-beta``).                                           |                                                                                          | If ``custom-model`` NOT specified                |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
| ``generate`` | ``custom-model``      | A IGoR parameters file followed by an IGoR marginals file.                                                                                                                        |                                                                                          |                                                  |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
| ``generate`` | ``n-gen``             | The number of sequences to generate.                                                                                                                                              | 1                                                                                        |                                                  |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
| ``generate`` | ``type``              | The type of model to create. (select one: ``alpha``, ``beta``, ``light`` or ``heavy``.                                                                                            |                                                                                          | If ``custom-model`` specified                    |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
| ``generate`` | ``cdr3``              | Generate CDR3 sequences instead.                                                                                                                                                  | Generate V(D)J full length sequences.                                                    |                                                  |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
| ``generate`` | ``anchor``            | A gene (V or J) followed by a CDR3 anchor separated data file. Note: need to contain gene in the first column, anchor index in the second and gene function in the third.         |                                                                                          | If ``cdr3`` and ``custom-model`` specified       |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
| ``evaluate`` | ``model``             | Specify a pre-installed model for generation. (select one: ``human-t-alpha``, ``human-t-beta``, ``human-b-heavy`` or ``mouse-t-beta``).                                           |                                                                                          | If ``custom-model`` NOT specified                |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
| ``evaluate`` | ``custom-model``      | A IGoR parameters file followed by an IGoR marginals file.                                                                                                                        |                                                                                          |                                                  |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
| ``evaluate`` | ``seqs``              | An input FASTA or separated data file with sequences for training the model.                                                                                                      |                                                                                          | Yes                                              |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
| ``evaluate`` | ``ref``               | A gene (V, D or J) followed by a reference genome FASTA file. Note: the FASTA reference genome files needs to conform to IGMT annotation (separated by vertical bar character).   |                                                                                          | If ``custom-model`` without ``cdr3`` specified   |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
| ``evaluate`` | ``type``              | The type of model to create. (select one: ``alpha``, ``beta``, ``light`` or ``heavy``.                                                                                            |                                                                                          | If ``custom-model`` specified                    |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
| ``evaluate`` | ``cdr3``              | Generate CDR3 sequences instead.                                                                                                                                                  | Generate V(D)J full length sequences.                                                    |                                                  |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
| ``evaluate`` | ``anchor``            | A gene (V or J) followed by a CDR3 anchor separated data file. Note: need to contain gene in the first column, anchor index in the second and gene function in the third.         |                                                                                          | If ``cdr3`` and ``custom-model`` specified       |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+
| ``evaluate`` | ``use-allele``        | If specified in combination with the ``cdr3`` flag, the allele information from the gene resolved fields is used to calculate the generation probability.                         | Allele ``01`` is used for each gene.                                                     |                                                  |
+--------------+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+--------------------------------------------------+

Configuration file setup
^^^^^^^^^^^^^^^^^^^^^^^^

ImmunoProbs supports user specified run configurations to modify additional settings that are not available to the user via the commandline tools. When a user specifies a configuration file to ImmunoProbs is will be merged with ImmunoProbs default configuration to make sure that all variables are set. The configuration is separated into a number of general sections:

* ``COMMON`` - Parameters that are common to all ImmunoProbs tools. Note that the flags given in the commandline will overwrite the ones in the configuration file (priority: ``default ImmunoProbs configuration < user specified configuration < commandline parameters``).
* ``EXPERT`` - Parameters that will likely never get modified. These could solve some system depending (e.g a compute cluster) issues when executing ImmunoProbs.

Additionally to the general sections, there are sections for each tool (e.g ``LOCATE``). These contain variables that are only used within that specific tool. The complete default configuration file of ImmunoProbs is shown in the code block below. Remember that the user does not have to specify each section and variable in their own configuration file. Only the variables with corresponding section that are of interest.

.. code-block:: ini

    ; Contains common parameters used throughout ImmunoProbs and its tools.
    [COMMON]
    ; The number of threads the system can use. By default max threads to system.
    NUM_THREADS
    ; The separator character for file in/out. Default tab character.
    SEPARATOR
    ; The directory for ImmunoProbs for writing files to. Default current directory.
    WORKING_DIR
    ; The output filename (or prefix value) that should be used for any given ImmunoProbs tool. Default None
    OUT_NAME
    ; The name of the column to use that identifies the each row in the input file.
    ROW_ID_COL = row_id
    ; The column name to use for the sequence filename idetifier.
    FILE_NAME_ID_COL = file_name_id
    ; Name of the column containing the sequence indices.
    I_COL = seq_index
    ; Name of the column containing the nucleotide sequences.
    NT_COL = nt_sequence
    ; Name of the column containing the nucleotide pgen scores.
    NT_P_COL = nt_pgen_estimate
    ; Name of the column containing the amino acid sequences.
    AA_COL = aa_sequence
    ; Name of the column containing the amino acid pgen scores.
    AA_P_COL = aa_pgen_estimate
    ; Name of the column specifying the frame type of the sequences.
    FRAME_TYPE_COL = frame_type
    ; Name of the column specifying the length of the CDR3 sequences.
    CDR3_LENGTH_COL = cdr3_length
    ; Name of the column containing the resolved V gene name string.
    V_RESOLVED_COL = v_resolved
    ; Name of the column containing the resolved J gene name string.
    J_RESOLVED_COL = j_resolved
    ; Name of the column containing the V gene choice string.
    V_GENE_CHOICE_COL = v_gene_choice
    ; Name of the column containing the D gene choice string.
    D_GENE_CHOICE_COL = d_gene_choice
    ; Name of the column containing the J gene choice string.
    J_GENE_CHOICE_COL = j_gene_choice

    ; Parameters specific for the 'convert' tool.
    [CONVERT]
    ; The number of random sequences to convert.
    NUM_RANDOM = 0
    ; The default allele value to use when USE_ALLELE = false.
    DEFAULT_ALLELE = 01
    ; If true, use the the allele information from the input file
    USE_ALLELE = false

    ; Parameters specific for the 'build' tool.
    [BUILD]
    ; The default number of inference iteration for training a model.
    NUM_ITERATIONS = 1

    ; Parameters specific for the 'locate' tool.
    [LOCATE]
    ; The default search motifs for the V gene.
    V_MOTIFS = TGT,TGC
    ; The default search motifs for the J gene.
    J_MOTIFS = TGG,TTC,TTT

    ; Parameters specific for the 'generate' tool.
    [GENERATE]
    ; The number of sequences to generate.
    NUM_GENERATE = 1
    ; If the tool should evaluate CDR3 sequnces instead of VDJ ones.
    EVAL_CDR3 = false

    ; Parameters specific for the 'evaluate' tool.
    [EVALUATE]
    ; If the tool should evaluate CDR3 sequnces instead of VDJ ones.
    EVAL_CDR3 = false
    ; The default allele value to use when USE_ALLELE = false.
    DEFAULT_ALLELE = 01
    ; If true, use the the allele information from the input file
    USE_ALLELE = false

    ; Contains expert parameters that should never have to be modified with normal usage of ImmunoProbs.
    [EXPERT]
    ; Should ImmunoProbs use the system's temporary directory (default) or use the WORKING_DIR location?
    USE_SYSTEM_TEMP = true
    ; Should ImmunoProbs use remove the temporary working directory after execution?
    REMOVE_TEMP_DIR = true
    ; The name of the temporary directory used by ImmunoProbs.
    TEMP_DIR = immuno_probs_tmp
