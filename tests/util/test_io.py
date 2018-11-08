# Pygor is part of the IGoR (Inference and Generation of Repertoires) software.
# Pygor Python package can be used to post process files generated by IGoR.
# Copyright (C) 2018 Quentin Marcou & Wout van Helvoirt

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


"""Test file for testing pygor.util.io file."""


import pandas
import pytest

from pygor.util.io import read_fasta_as_dataframe


@pytest.mark.parametrize('infile, expected', [
    ('tests/test_data/IGL_J_mouse.fasta', pandas.DataFrame(
        [['V00813|IGLJ1*01|Mus', 'CTGGGTGTTCGGTGGAGGAACCAAACTGACTGTCCTAG'],
         ['J00593|IGLJ2*01|Mus', 'TTATGTTTTCGGCGGTGGAACCAAGGTCACTGTCCTAG'],
         ['J00583|IGLJ3*01|Mus', 'GTTTATTTTCGGCAGTGGAACCAAGGTCACTGTCCTAG'],
         ['J00584|IGLJ3P*01|Mus', 'AGGTTCTTTTTCCTCAAATGGCCTATTGTATGCAGGAG'],
         ['J00596|IGLJ4*01|Mus', 'TTGGGTGTTCGGAGGTGGAACCAGATTGACTGTCCTAGATGA'],
         ['M16555|IGLJ4*01|Mus', 'TTGGGTGTTCGGAGGTGGAACCAGATTGACTGTCCTAG'],
         ['AF357974|IGLJ5*01|Mus', 'TTGGGTGTTCGGAGGTGGAACCAGATTGACTGTCCTAG']],
        columns=['name', 'sequence']))
])
def test_read_fasta_as_dataframe(infile, expected):
    """Test if a FASTA file can be read as pandas.DataFrame.

    The dataframe contains label name and sequence columns containing the
    corresponding FASTA data.

    Parameters
    ----------
    filename : string
        Location of the FASTA file to be read in.

    Raises
    -------
    AssertionError
        If the performed test failed.

    """
    result = read_fasta_as_dataframe(infile=infile)
    assert result.equals(expected)
