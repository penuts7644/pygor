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


"""Test file for testing immuno_probs.cdr3.anchor_locator file."""


import pandas
import pytest

from immuno_probs.alignment.muscle_aligner import MuscleAligner
from immuno_probs.cdr3.anchor_locator import AnchorLocator


def create_alignment():
    """Create an alignment to use for testing."""
    filename = 'tests/data/mouse_B_heavy/ref_genomes/genomicJs.fasta'
    aligner = MuscleAligner(infile=filename)
    return aligner.get_muscle_alignment()


@pytest.mark.parametrize('gene, motif, expected', [
    ('J', 'TTT', pandas.DataFrame(
        [['V00770|IGHJ1*02|Mus musculus|F|J-REGION|65..117|53 nt|2| | | | |53+0=53| | |', 44, 'TTT']],
        columns=['name', 'anchor_index', 'motif'])
    ),
    pytest.param('J', 'TGG', pandas.DataFrame(
        [['V00770|IGHJ3*01|Mus musculus|F|J-REGION|766..813|48 nt|3| | | | |48+0=48| | |', 14, 'TGG'],
         ['S73821|IGHJ3*02|Mus musculus|P|J-REGION|650..697|48 nt|3| | | | |48+0=48| | |', 14, 'TGG'],
         ['V00770|IGHJ1*02|Mus musculus|F|J-REGION|65..117|53 nt|2| | | | |53+0=53| | |', 19, 'TGG'],
         ['X63164|IGHJ1*03|Mus musculus_A/J|F|J-REGION|12..64|53 nt|2| | | | |53+0=53| | |', 19, 'TGG'],
         ['V00762|IGHJ1*01|Mus musculus_BALB/c|F|J-REGION|444..496|53 nt|2| | | | |53+0=53| | |', 19, 'TGG']],
        columns=['name', 'anchor_index', 'motif'])),
    pytest.param('X', None, None, marks=pytest.mark.xfail)
])
def test_anchor_locator(gene, motif, expected):
    """Test if correct indices of conserved motif regions are returned.

    Parameters
    ----------
    gene : str
        A gene identifier, either V or J, specifying the alignment's origin.
    motif : str
        A custom motif string to use for the search.
    expected : pandas.DataFrame
        The expected output pandas.Dataframe with correct columns and values.

    Raises
    -------
    AssertionError
        If the performed test failed.

    """
    locator = AnchorLocator(alignment=create_alignment(), gene=gene)
    if motif is not None:
        result = locator.get_indices_motifs(motif).head()
    else:
        result = locator.get_indices_motifs().head()
    assert (result == expected).all().all()
