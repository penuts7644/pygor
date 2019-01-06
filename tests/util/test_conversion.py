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


"""Test file for testing immuno_probs.util.conversion file."""


import pytest

from immuno_probs.util.conversion import nucleotides_to_integers
from immuno_probs.util.conversion import integers_to_nucleotides
from immuno_probs.util.conversion import reverse_complement
from immuno_probs.util.conversion import string_array_to_list


@pytest.mark.parametrize('seq, expected', [
    ('ACGT', '0123'),
    pytest.param('ACGT', '1302', marks=pytest.mark.xfail)
])
def test_nucleotides_to_integers(seq, expected):
    """Test if nucleotide sequence can be converted to interger representation.

    Parameters
    ----------
    seq : string
        A nucleotide sequence string.
    expected : str
        The expected output string.

    Raises
    -------
    AssertionError
        If the performed test failed.

    """
    out = nucleotides_to_integers(seq=seq)
    assert out == expected


@pytest.mark.parametrize('int_seq, expected', [
    ('0123', 'ACGT'),
    pytest.param('1302', 'ACGT', marks=pytest.mark.xfail)
])
def test_integers_to_nucleotides(int_seq, expected):
    """Test if integer sequence can be converted to nucleotide representation.

    Parameters
    ----------
    int_seq : string
        A integer sequence string.
    expected : str
        The expected output string.

    Raises
    -------
    AssertionError
        If the performed test failed.

    """
    out = integers_to_nucleotides(int_seq=int_seq)
    assert out == expected


@pytest.mark.parametrize('seq, expected', [
    ('ACGT', 'TGCA'),
    pytest.param('CTAG', 'ACGT', marks=pytest.mark.xfail)
])
def test_reverse_complement(seq, expected):
    """Test if nucleotide sequence can be converted to reverse complement.

    Parameters
    ----------
    seq : string
        A nucleotide sequence string.
    expected : str
        The expected output string.

    Raises
    -------
    AssertionError
        If the performed test failed.

    """
    out = reverse_complement(seq=seq)
    assert out == expected


@pytest.mark.parametrize('in_str, dtype, l_bound, r_bound, sep, expected', [
    ('(1, 2, 3, 4)', int, '(', ')', ',', [1, 2, 3, 4]),
    ('(1, 2, 3, 4)', float, '(', ')', ',', [1.0, 2.0, 3.0, 4.0]),
    pytest.param('((1, 2, 3, 4))', int, '(', ')', ',', [1, 2, 3, 4], marks=pytest.mark.xfail),
    pytest.param('(1, 2, 3, 4)', int, '[', ']', ',', [1, 2, 3, 4], marks=pytest.mark.xfail),
    pytest.param('[1, 2, 3, 4]', int, '[', ']', '.', [1, 2, 3, 4], marks=pytest.mark.xfail)
])
def test_string_array_to_list(in_str, dtype, l_bound, r_bound, sep, expected):
    """Test if integer sequence can be converted to nucleotide representation.

    Parameters
    ----------
    in_str : string
        A array representated as string.
    dtype : type
        The dtype to used for converting the individual the list elements.
    l_bound : string
        A string specifying the left boundary character(s).
    r_bound : string
        A string specifying the right boundary character(s).
    sep : string
        The separator character used in the input string.
    expected : list
        The expected output list.

    Raises
    -------
    AssertionError
        If the performed test failed.

    """
    out = string_array_to_list(in_str=in_str, dtype=dtype, l_bound=l_bound,
                               r_bound=r_bound, sep=sep)
    assert out == expected
