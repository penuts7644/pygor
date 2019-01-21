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


"""Test file for testing immuno_probs.util.io file."""


import pandas
import pytest

from immuno_probs.util.io import read_fasta_as_dataframe


@pytest.mark.parametrize('infile, expected', [
    ('tests/test_data/IGH_mus_musculus/ref_genomes/genomicJs.fasta',
     pandas.DataFrame(
         [['V00762|IGHJ1*01|Mus musculus_BALB/c|F|J-REGION|444..496|53 nt|2| | | | |53+0=53| | |', 'CTACTGGTACTTCGATGTCTGGGGCGCAGGGACCACGGTCACCGTCTCCTCAG'],
          ['V00770|IGHJ1*02|Mus musculus|F|J-REGION|65..117|53 nt|2| | | | |53+0=53| | |', 'CTACTGGTACTTCGATGTCTGGGGCGCAGGGACCACGGTCACCGTTTCCTCAG'],
          ['X63164|IGHJ1*03|Mus musculus_A/J|F|J-REGION|12..64|53 nt|2| | | | |53+0=53| | |', 'CTACTGGTACTTCGATGTCTGGGGCACAGGGACCACGGTCACCGTCTCCTCAG'],
          ['V00770|IGHJ2*01|Mus musculus|F|J-REGION|383..430|48 nt|3| | | | |48+0=48| | |', 'ACTACTTTGACTACTGGGGCCAAGGCACCACTCTCACAGTCTCCTCAG'],
          ['S73821|IGHJ2*02|Mus musculus|F|J-REGION|267..314|48 nt|3| | | | |48+0=48| | |', 'ACTACTTTGACTACTGGGGCCAAGGCACCTCTCTCACAGTCTCCTCAG']],
         columns=['header', 'sequence']))
])
def test_read_fasta_as_dataframe(infile, expected):
    """Test if a FASTA file can be read as pandas.DataFrame.

    The dataframe contains header name and sequence columns containing the
    corresponding FASTA data.

    Parameters
    ----------
    filename : string
        Location of the FASTA file to be read in.
    expected : pandas.DataFrame
        The expected output pandas dataframe.

    Raises
    -------
    AssertionError
        If the performed test failed.

    """
    result = read_fasta_as_dataframe(infile=infile)
    assert result.head().equals(expected)
