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


CONFIG_DATA = None


def set_config_data(value=None):
    """Sets and updates the global CONFIG_DATA variable by parsing config files.

    Parameters
    ----------
    value : str, optional
        An optional ImmunoProbs configuration file path to parse besides the default file.

    """
    # Parse default configuration file.
    pkg_name = __name__.split('.')[0]
    config_file_path = resource_filename(pkg_name, os.path.join('config', 'default.ini'))
    conf_parser = RawConfigParser(allow_no_value=True)
    conf_parser.read(config_file_path)

    # If given parse additional configuration.
    if value:
        conf_parser.read(value)

    # Set the global config data object.
    globals().update(CONFIG_DATA=conf_parser)

    # Overwrite default values if not given.
    if not conf_parser.get('COMMON', 'NUM_THREADS'):
        set_num_threads()
    if not conf_parser.get('COMMON', 'SEPARATOR'):
        set_separator()
    if not conf_parser.get('COMMON', 'WORKING_DIR'):
        set_working_dir()
    if not conf_parser.get('COMMON', 'OUT_NAME'):
        set_out_name()


def get_config_data(section, value, option_type=None):
    """Collects and returns the global CONFIG_DATA variable.

    Parameters
    ----------
    section : str
        The section where the given option is located.
    value : str
        The option to return its value from.
    option_type : str, optional
        The type of the option to return its value from, by default returns a string. Currently supported values are 'bool'
        for boolean, 'int' for integer and 'float' for float.

    Returns
    -------
    str
        The value of the option within the configuration file.

    """
    if CONFIG_DATA.has_option(section, value):
        if option_type == 'bool':
            return CONFIG_DATA.getboolean(section, value)
        if option_type == 'int':
            return CONFIG_DATA.getint(section, value)
        if option_type == 'float':
            return CONFIG_DATA.getfloat(section, value)
        return CONFIG_DATA.get(section, value)


def set_num_threads(value=ph.cpu_count()):
    """Sets and updates the global NUM_THREADS variable.

    Parameters
    ----------
    value : int, optional
        The number of threads the program is allowed to use (default: max available threads).

    Raises
    ------
    TypeError
        When the NUM_THREADS global variable is not an integer.
    ValueError
        When the NUM_THREADS global variable is smaller then 1.

    """
    if not isinstance(value, int):
        raise TypeError("The NUM_THREADS variable needs to be of type integer", value)
    if value < 1:
        raise ValueError("The NUM_THREADS variable needs to be higher than zero", value)
    else:
        CONFIG_DATA.set('COMMON', 'NUM_THREADS', str(value))


def set_separator(value='tab'):
    """Sets and updates the global SEPARATOR variable.

    Parameters
    ----------
    value : str, optional
        The separator character to be used when writing files (default: tab character).

    Raises
    ------
    TypeError
        When the SEPARATOR global variable is not of type string.

    """
    separators = {'tab': '\t', 'semi-colon': ';', 'comma': ','}
    if not isinstance(value, str):
        raise TypeError("The SEPARATOR variable needs to be of type string", value)
    else:
        CONFIG_DATA.set('COMMON', 'SEPARATOR', separators[value])


def set_working_dir(value=os.getcwd()):
    """Sets and updates the global WORKING_DIR variable.

    Parameters
    ----------
    value : str, optional
        The directory path to use when writing output files (default: the current working directory).

    Raises
    ------
    TypeError
        When the WORKING_DIR global variable is not of type string.
    IOError
        When the WORKING_DIR global variable directory does not exist on the system.

    """
    if not isinstance(value, str):
        raise TypeError("The WORKING_DIR variable needs to be of type string", value)
    if not os.path.isdir(value):
        raise IOError("The WORKING_DIR variable needs to be an existing directory", value)
    else:
        CONFIG_DATA.set('COMMON', 'WORKING_DIR', value)


def set_out_name(value=None):
    """Sets and updates the global OUT_NAME variable.

    Parameters
    ----------
    value : str, optional
        The output file name string to use when writing output files or when prefixing output files (default: None).

    """
    if value:
        value = re.sub(r'\s+', '', value)
    CONFIG_DATA.set('COMMON', 'OUT_NAME', value)


# Set the default config data object.
set_config_data()
