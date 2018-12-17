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


def set_num_threads(num_threads):
    """Sets the global NUM_THREADS variable."""
    globals().update(NUM_THREADS=num_threads)

def get_num_threads():
    """Returns the global NUM_THREADS variable."""
    return NUM_THREADS

def set_separator(separator):
    """Sets the global SEPARATOR variable."""
    globals().update(SEPARATOR=separator)

def get_separator():
    """Returns the global SEPARATOR variable."""
    return SEPARATOR
