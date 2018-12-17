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


"""Contains processing functions used in immuno_probs."""


import numpy
import pathos.pools as pp
import pathos.helpers as help

from immuno_probs.util.constant import get_num_threads
from immuno_probs.util.exception import NumThreadsValueException


def multiprocess_array(ary, func, **kwargs):
    """Applies multiprocessing on a multi array using the given function.

    Parameters
    ----------
    ary : list
        List 'like' object to be split for multiple workers.
    func : Object
        A function object that the workers should apply.
    **kwargs
        The remaining arguments to be given to the input function.

    Returns
    -------
    list
        Containing the results from each of the workers.

    Raises
    ------
    NumThreadsValueException
        When the NUM_THREADS global variable is not an integer or is smaller
        then 1.

    Notes
    -----
        This function uses the NUM_THREADS constant from immuno_probs.util.constant and
        will limit the number of workers to create based on this value. Overwrite
        NUM_THREADS constant to increase the number of workers. By default uses the
        cpu count from pathos package.

    """
    # Check out available worker count and adjust accordingly.
    num_workers = get_num_threads()
    if not isinstance(num_workers, int) or num_workers < 1:
        raise NumThreadsValueException("The NUM_THREADS variable needs to be of " \
                                       "type integer and higher than zero", num_workers)
    if len(ary) < num_workers:
        num_workers = len(ary)

    # Divide the array into chucks for the workers.
    pool = pp.ProcessPool(nodes=num_workers)
    result = pool.map(func, [(d, kwargs)
                             for d in numpy.array_split(ary, num_workers)])
    return result
