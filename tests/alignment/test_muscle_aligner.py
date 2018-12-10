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


"""Test file for testing immuno_probs.alignment.muscle_aligner file."""


from Bio.Align import MultipleSeqAlignment
import pytest

from immuno_probs.alignment.muscle_aligner import MuscleAligner


@pytest.mark.parametrize('infile, cmd, expected', [
    ('tests/test_data/IGH_mus_musculus/ref_genomes/genomicJs.fasta', 'muscle',
     MultipleSeqAlignment),
    pytest.param('tests/test_data/IGH_mus_musculus/ref_genomes/genomicJs.fasta',
                 'fake_command', MultipleSeqAlignment, marks=pytest.mark.xfail)
])
def test_muscle_aligner(infile, cmd, expected):
    """Test if fasta file can be aligned by MUSCLE commandline tool.

    Parameters
    ----------
    infile : string
        A file path to a FASTA file containining the genomic data to align.
    cmd : string
        The MUSCLE terminal command executable location/name.
    expected : MultipleSeqAlignment
        The expected output type MultipleSeqAlignment.

    Raises
    -------
    AssertionError
        If the performed test failed.

    """
    aligner = MuscleAligner(infile=infile, cmd=cmd)
    alignment = aligner.get_muscle_alignment()
    assert isinstance(alignment, expected)
