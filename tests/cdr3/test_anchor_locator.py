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
    filename = 'tests/data/human_t_beta/ref_genomes/TRBJ.fasta'
    aligner = MuscleAligner(infile=filename)
    return aligner.get_muscle_alignment()


@pytest.mark.parametrize('gene, motif, expected', [
    ('J', 'TTT', pandas.DataFrame(
        [['M14158|TRBJ1-3*01|Homo sapiens|F|J-REGION|1499..1548|50 nt|2| | | | |50+0=50| | |', 19, 'TTT'],
         ['M14158|TRBJ1-4*01|Homo sapiens|F|J-REGION|2095..2145|51 nt|3| | | | |51+0=51| | |', 20, 'TTT'],
         ['M14158|TRBJ1-5*01|Homo sapiens|F|J-REGION|2368..2417|50 nt|2| | | | |50+0=50| | |', 19, 'TTT'],
         ['X02987|TRBJ2-2*01|Homo sapiens|F|J-REGION|995..1045|51 nt|3| | | | |51+0=51| | |', 20, 'TTT'],
         ['K02545|TRBJ1-1*01|Homo sapiens|F|J-REGION|749..796|48 nt|3| | | | |48+0=48| | |', 17, 'TTT']],
        columns=['name', 'anchor_index', 'motif'])
    ),
    pytest.param('J', 'TGG', pandas.DataFrame(
        [['M14158|TRBJ1-3*01|Homo sapiens|F|J-REGION|1499..1548|50 nt|2| | | | |50+0=50| | |', 21, 'TGG'],
         ['M14158|TRBJ1-4*01|Homo sapiens|F|J-REGION|2095..2145|51 nt|3| | | | |51+0=51| | |', 22, 'TGG'],
         ['M14158|TRBJ1-5*01|Homo sapiens|F|J-REGION|2368..2417|50 nt|2| | | | |50+0=50| | |', 21, 'TGG'],
         ['X02987|TRBJ2-2*01|Homo sapiens|F|J-REGION|995..1045|51 nt|3| | | | |51+0=51| | |', 22, 'TGG'],
         ['K02545|TRBJ1-1*01|Homo sapiens|F|J-REGION|749..796|48 nt|3| | | | |48+0=48| | |', 19, 'TGG']],
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
        result = locator.get_indices_motifs(1, motif).head()
    else:
        result = locator.get_indices_motifs(1).head()
    print(result)
    print(expected)
    assert (result == expected).all().all()
