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


"""Test file for testing pygor.cdr3.anchor_locator file."""


import pandas
import pytest

from pygor.alignment.muscle_aligner import MuscleAligner
from pygor.cdr3.anchor_locator import AnchorLocator


def create_alignment():
    """Create an alignment to use for testing."""
    aligner = MuscleAligner(infile='tests/test_data/IGL_J_mouse.fasta')
    return aligner.get_muscle_alignment()


@pytest.mark.parametrize('gene, motif, expected', [
    ('J', None, pandas.DataFrame(
        [['TGG', 'J00593|IGLJ2*01|Mus', 15.0],
         ['TGG', 'J00583|IGLJ3*01|Mus', 15.0],
         ['TGG', 'J00596|IGLJ4*01|Mus', 15.0],
         ['TGG', 'M16555|IGLJ4*01|Mus', 15.0],
         ['TGG', 'AF357974|IGLJ5*01|Mus', 15.0],
         ['TGG', 'J00584|IGLJ3P*01|Mus', 18.0],
         ['TTT', 'J00593|IGLJ2*01|Mus', 5.0],
         ['TTT', 'J00583|IGLJ3*01|Mus', 5.0],
         ['TTT', 'J00584|IGLJ3P*01|Mus', 8.0],
         ['TTC', 'J00593|IGLJ2*01|Mus', 7.0],
         ['TTC', 'J00583|IGLJ3*01|Mus', 7.0],
         ['TTC', 'J00596|IGLJ4*01|Mus', 7.0],
         ['TTC', 'M16555|IGLJ4*01|Mus', 7.0],
         ['TTC', 'AF357974|IGLJ5*01|Mus', 7.0],
         ['TTC', 'V00813|IGLJ1*01|Mus', 7.0]],
        columns=['motif', 'seq_id', 'start_index'])
    ),
    pytest.param('J', 'TGG', pandas.DataFrame(
        [['TGG', 'J00593|IGLJ2*01|Mus', 15.0],
         ['TGG', 'J00583|IGLJ3*01|Mus', 15.0],
         ['TGG', 'J00596|IGLJ4*01|Mus', 15.0],
         ['TGG', 'M16555|IGLJ4*01|Mus', 15.0],
         ['TGG', 'AF357974|IGLJ5*01|Mus', 15.0],
         ['TGG', 'J00584|IGLJ3P*01|Mus', 18.0]],
        columns=['motif', 'seq_id', 'start_index'])),
    pytest.param('X', None, None, marks=pytest.mark.xfail)
])
def test_anchor_locator(gene, motif, expected):
    """Test if correct indices of conserved motif regions are returned.

    Parameters
    ----------
    gene : string
        A gene identifier, either V or J, specifying the alignment's origin.
    motif : string
        A custom motif string to use for the search.
    expected : pandas.DataFrame
        The expected output pandas.Dataframe with coreect columns and values.

    Raises
    -------
    AssertionError
        If the performed test failed.

    """
    locator = AnchorLocator(alignment=create_alignment(), gene=gene)
    if motif is not None:
        result = locator.get_indices_motifs(motif)
    else:
        result = locator.get_indices_motifs()
    assert result.equals(expected)
