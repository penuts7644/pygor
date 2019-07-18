
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

Docker
^^^^^^

It is also possible to use a docker image of ImmunoProbs with all necessary requirement pre-installed in an ubuntu environment. Make sure to install `Docker <https://www.docker.com>`__ first and pull the most recent version of the image with the following command.

.. code-block:: none

    docker pull penuts7644/immuno-probs

Now you can use the ImmunoProbs image however you like (through calling ImmunoProbs directly or by opening an interactive session in the bash shell).

You can also execute the ImmunoProbs image while using files in your local machine's working directory as a mounted volume:

.. code-block:: none

    docker run --rm -v "$PWD":/tmp penuts7644/immuno-probs \
      immuno-probs \
        [TOOL NAME] \
          <TOOL OPTIONS>

Note that the container is removed after execution and that you'll need to link to the ``/tmp`` directory of the container since that is the location where output files are created.

Some additional information on running docker images can be found below:

-  `Copy files <https://docs.docker.com/engine/reference/commandline/cp/>`__ to/from the container (``docker cp``).

-  Run an image in a `detached <https://docs.docker.com/engine/reference/commandline/run/>`__ container (``docker run -d``) and `attach <https://docs.docker.com/engine/reference/commandline/attach/>`__ it later on (``docker attach``).

Galaxy server
~~~~~~~~~~~~~

The ImmunoProbs docker image can be integrated as a galaxy tool by modifying the tool configuration XML file (``tool_conf.xml``). For each of the build in ImmunoProbs tools, there is a wrapper file located in the ``galaxy`` directory. The tools can be installed in a galaxy server by copying each of them over to your tools directory in the galaxy. Replace ``<LOCATION>`` to the location of your galaxy tool directory.

.. code-block:: none

    wget -P <LOCATION> "https://raw.githubusercontent.com/penuts7644/ImmunoProbs/master/galaxy/build_igor_model.xml" & \
    wget -P <LOCATION> "https://raw.githubusercontent.com/penuts7644/ImmunoProbs/master/galaxy/locate_cdr3_anchors.xml" & \
    wget -P <LOCATION> "https://raw.githubusercontent.com/penuts7644/ImmunoProbs/master/galaxy/generate_sequences.xml" & \
    wget -P <LOCATION> "https://raw.githubusercontent.com/penuts7644/ImmunoProbs/master/galaxy/evaluate_sequences.xml"

Finally, add the section with the ImmunoProbs tools to the ``tool_conf.xml``. Replace ``<LOCATION>`` to the location of each ImmunoProbs tool within you galaxy tools directory.

.. code:: xml

    <section id="immuno_probs" name="ImmunoProbs">
        <tool file="<LOCATION>/build_igor_model.xml" />
        <tool file="<LOCATION>/locate_cdr3_anchors.xml" />
        <tool file="<LOCATION>/generate_sequences.xml" />
        <tool file="<LOCATION>/evaluate_sequences.xml" />
    </section>

Make sure to have setup you galaxy server to be able to use docker images. This can be done inside the ``job_conf.xml`` file by adding the following:

.. code:: xml

    <destinations default="docker_local">
        <destination id="local" runner="local"/>
            <destination id="docker_local" runner="local">
                <param id="docker_enabled">true</param>
            </destination>
    </destinations>
