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


"""Contains conversion functions used in immuno_probs."""


from immuno_probs.util.exception import CharacterNotFoundException


def nucleotides_to_integers(seq):
    """Converts a nucleotide sequence to an interger representation.

    The base characters in the nucleotide string (A, C, G and T) are converted
    to the following: A -> 0, C -> 1, G -> 2 and T -> 3. The combined uppercase
    string is returned.

    Parameters
    ----------
    seq : string
        A nucleotide sequence string.

    Returns
    -------
    string
        The interger representation string for the given nucleotide sequence.

    """
    int_sequence = []
    for i in seq.upper():
        if i == 'A':
            int_sequence.append(str(0))
        elif i == 'C':
            int_sequence.append(str(1))
        elif i == 'G':
            int_sequence.append(str(2))
        elif i == 'T':
            int_sequence.append(str(3))
    return ''.join(int_sequence)


def integers_to_nucleotides(int_seq):
    """Converts a integer sequence to an nucleotide representation.

    The base characters in the integer string (0, 1, 2 and 3) are converted
    to the following: 0 -> A, 1 -> C, 2 -> G and 3 -> T. The combined string
    is returned.

    Parameters
    ----------
    int_seq : string
        A integer sequence string.

    Returns
    -------
    string
        The nucleotide representation string for the given integer sequence.

    """
    nuc_sequence = []
    for i in int_seq:
        if int(i) == 0:
            nuc_sequence.append('A')
        elif int(i) == 1:
            nuc_sequence.append('C')
        elif int(i) == 2:
            nuc_sequence.append('G')
        elif int(i) == 3:
            nuc_sequence.append('T')
    return ''.join(nuc_sequence)


def reverse_complement(seq):
    """Converts a nucleotide sequence to reverse complement.

    The base characters in the nucleotide string (A, C, G and T) are converted
    to the following: A <-> T and C <-> G. The combined uppercase string is
    returned.

    Parameters
    ----------
    seq : string
        A nucleotide sequence string.

    Returns
    -------
    string
        The reverse complemented nucleotide sequence.

    """
    reverse_complement_seq = []
    for i in seq.upper():
        if i == 'A':
            reverse_complement_seq.append('T')
        elif i == 'C':
            reverse_complement_seq.append('G')
        elif i == 'G':
            reverse_complement_seq.append('C')
        elif i == 'T':
            reverse_complement_seq.append('A')
    return ''.join(reverse_complement_seq)


def string_array_to_list(in_str, dtype=float, l_bound='(', r_bound=')', sep=','):
    """Converts a string representation of an array to a python list.

    Removes the given boundary characters from the string and separates the
    individual items on the given seperator character. Each item is converted to
    the given dtype. The python list is returned.

    Parameters
    ----------
    in_str : string
        A array representated as string.
    dtype : type, optional
        The dtype to used for converting the individual the list elements. By
        default uses float.
    l_bound : string, optional
        A string specifying the left boundary character(s). By default '('.
    r_bound : string, optional
        A string specifying the right boundary character(s). By default ')'.
    sep : string, optional
        The separator character used in the input string. By default ','.

    Returns
    -------
    list
        The converted input string as python list.

    Raises
    -------
    CharacterNotFoundException
        When the given seperator or L/R bound characters are not found.

    """
    if len(in_str) > (len(l_bound) + len(r_bound)):

        # Check if start and end of the string match the boundary characters.
        if in_str[: len(l_bound)] != l_bound:
            raise CharacterNotFoundException('Start character not found', l_bound)
        elif in_str[len(in_str) - len(r_bound) :] != r_bound:
            raise CharacterNotFoundException('End character not found', r_bound)
        elif in_str.find(sep) == -1:
            raise CharacterNotFoundException('Seperator character not found', sep)

        # Strip the boundary characters, split on seperator and small cleanup.
        converted_str = [dtype(i.strip(' \"\''))
                         for i in in_str[len(l_bound) : len(in_str) - len(r_bound)]
                         .split(sep)]
    return converted_str