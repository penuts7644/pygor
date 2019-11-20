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


"""Test file for testing immuno_probs.model.igor_loader file."""


import pytest

from immuno_probs.model.igor_loader import IgorLoader


@pytest.mark.parametrize(
    'infiles, expected',
    [
        (
            [
                'tests/data/human_t_alpha/model_params.txt',
                'tests/data/human_t_alpha/model_marginals.txt',
                'tests/data/human_t_alpha/V_gene_CDR3_anchors.csv',
                'tests/data/human_t_alpha/J_gene_CDR3_anchors.csv'
            ],
            IgorLoader
        )
    ]
)
def test_igor_loader(infiles, expected):
    """Test if fasta file can be aligned by MUSCLE commandline tool.

    Parameters
    ----------
    infiles : list
        A list of file paths to an IGoR model and CDR3 anchor files.
    expected : IgorLoader
        The expected output type IgorLoader.

    Raises
    -------
    AssertionError
        If the performed test failed.

    """
    model = IgorLoader(model_type='alpha', model_params=infiles[0], model_marginals=infiles[1])
    model.set_anchor(gene='V', file=infiles[2])
    model.set_anchor(gene='J', file=infiles[3])
    model.initialize_model()
    assert isinstance(model, expected)
