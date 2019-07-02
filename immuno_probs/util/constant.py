# Create IGoR models and calculate the generation probability of V(D)J and
# CDR3 sequences. Copyright (C) 2019 Wout van Helvoirt

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


"""Contains a collection of global constant variables."""


import os
import re
from ConfigParser import RawConfigParser
from pkg_resources import resource_filename

import pathos.helpers as ph

from immuno_probs.util.exception import SeparatorNotValidException, \
NumThreadsValueException, DirectoryNonExistingException


CONFIG_DATA = None


def set_config_data(value=None):
    """Sets and updates the global CONFIG_DATA variable by parsing config files.

    Parameters
    ----------
    value : str, optional
        An optional ImmunoProbs configuration file path to parse besides the
        default file.

    """
    # Parse default configuration file.
    pkg_name = __name__.split('.')[0]
    config_file_path = resource_filename(
        pkg_name, os.path.join('config', 'default.ini'))
    conf_parser = RawConfigParser(allow_no_value=True)
    conf_parser.read(config_file_path)

    # If given parse additional configuration.
    if value:
        conf_parser.read(value)

    # Set the global config data object.
    globals().update(CONFIG_DATA=conf_parser)

    # Overwrite default values if not given.
    if not conf_parser.get('ImmunoProbs', 'NUM_THREADS'):
        set_num_threads()
    if not conf_parser.get('ImmunoProbs', 'SEPARATOR'):
        set_separator()
    if conf_parser.get('ImmunoProbs', 'SEPARATOR') == '|':
        set_separator('|')
    if not conf_parser.get('ImmunoProbs', 'WORKING_DIR'):
        set_working_dir()
    if not conf_parser.get('ImmunoProbs', 'OUT_NAME'):
        set_out_name()

def get_config_data(value):
    """Collects and returns the global CONFIG_DATA variable.

    Parameters
    ----------
    value : str
        The option to return its value from.

    Returns
    -------
    str
        The value of the option within the configuration file.

    """
    return CONFIG_DATA.get('ImmunoProbs', value)

def set_num_threads(value=ph.cpu_count()):
    """Sets and updates the global NUM_THREADS variable.

    Parameters
    ----------
    value : int, optional
        The number of threads the program is allowed to use (default: max
        available threads).

    Raises
    ------
    NumThreadsValueException
        When the NUM_THREADS global variable is not an integer or is smaller
        then 1.

    """
    if not isinstance(value, int) or value < 1:
        raise NumThreadsValueException(
            "The NUM_THREADS variable needs to be of type integer and higher " \
            "than zero", value)
    else:
        CONFIG_DATA.set('ImmunoProbs', 'NUM_THREADS', str(value))

def set_separator(value='\t'):
    """Sets and updates the global SEPARATOR variable.

    Parameters
    ----------
    value : str, optional
        The separator character to be used when writing files (default: tab
        character).

    Raises
    ------
    SeparatorNotValidException
        When the SEPARATOR global variable is not of type string or is a single
        '|' character.

    """
    if not isinstance(value, str) or value == '|':
        raise SeparatorNotValidException(
            "The SEPARATOR variable needs to be of type string and cannot " \
            "be '|' character", value)
    else:
        CONFIG_DATA.set('ImmunoProbs', 'SEPARATOR', value)

def set_working_dir(value=os.getcwd()):
    """Sets and updates the global WORKING_DIR variable.

    Parameters
    ----------
    value : str, optional
        The directory path to use when writing output files (default: the
        current working directory).

    Raises
    ------
    DirectoryNonExistingException
        When the WORKING_DIR global variable is not of type string and does not
        exist.

    """
    if not isinstance(value, str) or not os.path.isdir(value):
        raise DirectoryNonExistingException(
            "The WORKING_DIR variable needs to be of type string and exist " \
            "on the system", value)
    else:
        CONFIG_DATA.set('ImmunoProbs', 'WORKING_DIR', value)

def set_out_name(value=''):
    """Sets and updates the global OUT_NAME variable.

    Parameters
    ----------
    value : str, optional
        The output file name string to use when writing output files or when
        prefixing output files (default: none).

    """
    updated_value = re.sub(r'\s+', '', value)
    CONFIG_DATA.set('ImmunoProbs', 'OUT_NAME', updated_value)

# Set the default config data object.
set_config_data()
