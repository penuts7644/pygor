# ImmunoProbs Python package uses a simplified manner for calculating the
# generation probability of V(D)J and CDR3 sequences.
# Copyright (C) 2018 Wout van Helvoirt

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


import pathos.helpers as ph


NUM_THREADS = ph.cpu_count()
SEPARATOR = ','


def set_num_threads(value):
    """Updates the global NUM_THREADS variable.

    Parameters
    ----------
    value : int
        The number of threads the program is allowed to use (default: max
        available threads).

    """
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
        comma character).

    """
    globals().update(SEPARATOR=value)

def get_separator():
    """Returns the global SEPARATOR variable.

    Returns
    -------
    str
        The globally set separator value.

    """
    return SEPARATOR
