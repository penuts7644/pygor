# ImmunoProbs Python package able to calculate the generation probability of
# V(D)J and CDR3 sequences. Copyright (C) 2019 Wout van Helvoirt

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


"""Contains global constant variables used in immuno_probs."""


import os

import pathos.helpers as ph

from immuno_probs.util.exception import SeparatorNotValidException, NumThreadsValueException, DirectoryNonExistingException


NUM_THREADS = ph.cpu_count()
SEPARATOR = ';'
WORKING_DIR = os.getcwd()


def set_num_threads(value):
    """Updates the global NUM_THREADS variable.

    Parameters
    ----------
    value : int
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
        globals().update(NUM_THREADS=value)

def get_num_threads():
    """Returns the global NUM_THREADS variable.

    Returns
    -------
    str
        The globally set value for the number of threads.

    """
    return NUM_THREADS

def set_separator(value):
    """Updates the global SEPARATOR variable.

    Parameters
    ----------
    value : str
        The separator character to be used when writing files (default:
        semicolon character).

    Raises
    ------
    SeparatorNotValidException
        When the SEPARATOR global variable is not of type string.

    """
    if not isinstance(value, str):
        raise SeparatorNotValidException(
            "The SEPARATOR variable needs to be of type string", value)
    else:
        globals().update(SEPARATOR=value)

def get_separator():
    """Returns the global SEPARATOR variable.

    Returns
    -------
    str
        The globally set separator value.

    """
    return SEPARATOR

def set_working_dir(value):
    """Updates the global WORKING_DIR variable.

    Parameters
    ----------
    value : str
        The directory path to use when writing output files (default:
        the current working directory).

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
        globals().update(WORKING_DIR=value)

def get_working_dir():
    """Returns the global WORKING_DIR variable.

    Returns
    -------
    str
        The globally set working directory value.

    """
    return WORKING_DIR
