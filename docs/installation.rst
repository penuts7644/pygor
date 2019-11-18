
Installation
============

PyPI
^^^^

ImmunoProbs is installable via PyPI using the following terminal command:

.. code-block:: none

    pip install immuno-probs

In order to see all available tools/functionalities for ImmunoProbs use the help command after installation:

.. code-block:: none

    immuno-probs -h

Each ImmunoProbs tool has its own specific arguments (accompanying the top-level ones). The options for a specific tools are available via the tool's help command:

.. code-block:: none

    immuno-probs <TOOL NAME> -h

Make sure to install the necessary requirements (see the requirements below).

Requirements
~~~~~~~~~~~~

All Python dependencies that are used by this package are installed through pip upon installation of ImmunoProbs. However, some software (not available via pip) needs to installed manually when planning on using certain commandline tools from ImmunoProbs:

+-----------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Command                                       | Requirement                                                                                                                                                                                                                                                                                                               |
+===============================================+===========================================================================================================================================================================================================================================================================================================================+
| ``locate``                                    | This requires `MUSCLE <http://www.drive5.com/muscle/>`__ to be installed. For linux, MUSCLE can be installed via ``apt-get install muscle``. For macOS with HomeBrew, tap into ``brewsci/bio`` and install MUSCLE via ``brew install muscle``                                                                             |
+-----------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``build`` ``generate`` ``evaluate``           | This will use Python's subprocess package to pass the user arguments to `IGoR <https://github.com/qmarcou/IGoR>`__. For these tools to work properly, make sure that you have at least compiled and installed IGoR version 1.3.0 using the guide in `IGoR's documentation <https://qmarcou.github.io/IGoR/#install>`__.   |
+-----------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
