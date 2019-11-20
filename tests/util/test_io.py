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


"""Test file for testing immuno_probs.util.io file."""


import pandas
import pytest

from immuno_probs.util.io import read_fasta_as_dataframe


@pytest.mark.parametrize(
    'file, expected',
    [
        (
            'tests/data/human_t_beta/ref_genomes/TRBJ.fasta',
            pandas.DataFrame(
                [
                    ['TGAACACTGAAGCTTTCTTTGGACAAGGCACCAGACTCACAGTTGTAG'],
                    ['CTAACTATGGCTACACCTTCGGTTCGGGGACCAGGTTAACCGTTGTAG'],
                    ['CTCTGGAAACACCATATATTTTGGAGAGGGAAGTTGGCTCACTGTTGTAG'],
                    ['CAACTAATGAAAAACTGTTTTTTGGCAGTGGAACCCAGCTCTCTGTCTTGG'],
                    ['TAGCAATCAGCCCCAGCATTTTGGTGATGGGACTCGACTCTCCATCCTAG']
                ],
                columns=['nt_sequence']
            )
        )
    ]
)
def test_read_fasta_as_dataframe(file, expected):
    """Test if a FASTA file can be read as pandas.DataFrame.

    The dataframe contains header name and sequence columns containing the
    corresponding FASTA data.

    Parameters
    ----------
    filename : str
        Location of the FASTA file to be read in.
    expected : pandas.DataFrame
        The expected output pandas dataframe.

    Raises
    -------
    AssertionError
        If the performed test failed.

    """
    result = read_fasta_as_dataframe(file=file, col='nt_sequence')
    assert (result.head() == expected).all().all()
