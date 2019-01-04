# ImmunoProbs Python package able to calculate the generation probability of
# V(D)J and CDR3 sequences. Copyright (C) 2018 Wout van Helvoirt

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


def multiprocess_array(ary, func, num_workers, **kwargs):
    """Applies multiprocessing on a multi array using the given function.

    Parameters
    ----------
    ary : list
        List 'like' object to be split for multiple workers.
    func : Object
        A function object that the workers should apply.
    num_workers : int
        The number of threads the program is allowed to use.
    **kwargs
        The remaining arguments to be given to the input function.

    Returns
    -------
    list
        Containing the results from each of the workers.

    """
    # Check out available worker count and adjust accordingly.
    if len(ary) < num_workers:
        num_workers = len(ary)

    # Divide the array into chucks for the workers.
    pool = pp.ProcessPool(nodes=num_workers)
    result = pool.amap(func, [(d, kwargs)
                              for d in numpy.array_split(ary, num_workers)])
    return result.get()
