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


"""Test file for testing immuno_probs.cdr3.olga_container file."""


import pandas
import pytest

from immuno_probs.cdr3.olga_container import OlgaContainer
from immuno_probs.model.igor_loader import IgorLoader


@pytest.mark.parametrize('option, expected', [
    ('generate', [['seq_index', int], ['nt_sequence', str], ['aa_sequence', str],
                  ['gene_choice_v', int], ['gene_choice_j', int]]),
    ('evaluate', pandas.DataFrame(
        [[0, 'TGTGCCAGTAGTATAACAACCCAGGGCTTGTACGAGCAGTACTTC', 0],
         [1, 'TGTGCAGGAATAAACTTTGGAAATGAGAAATTAACCTTT', 6.022455403460228e-08],
         [2, 'TGTGCATTGAACAGAGATGACAAGATCATCTTT', 3.8690138672702246e-07]],
        columns=['seq_index', 'nt_sequence', 'nt_pgen_estimate',])
    )
])
def test_olga_container(option, expected):
    """Test if the container class can generate and evaluate the CDR3's.

    Parameters
    ----------
    option : string
        A option to run in the container ('generate' of 'evaluate').
    expected : pandas.DataFrame
        The expected output pandas.Dataframe with correct columns and values.

    Raises
    -------
    AssertionError
        If the performed test failed.

    """
    params = 'tests/test_data/TRA_homo_sapiens/model_params.txt'
    marginals = 'tests/test_data/TRA_homo_sapiens/model_marginals.txt'
    v_anchors = 'tests/test_data/TRA_homo_sapiens/V_gene_CDR3_anchors.csv'
    j_anchors = 'tests/test_data/TRA_homo_sapiens/J_gene_CDR3_anchors.csv'
    model = IgorLoader(model_params=params, model_marginals=marginals,
                       v_anchors=v_anchors, j_anchors=j_anchors)
    olga_container = OlgaContainer(igor_model=model)
    if option == 'generate':
        result = olga_container.generate(num_seqs=1)
        for i in expected:
            assert isinstance(i[1](result.iloc[0][i[0]]), i[1])
    elif option == 'evaluate':
        pgen_seqs = expected.drop(['nt_pgen_estimate'], axis=1)
        result = olga_container.evaluate(seqs=pgen_seqs)
        for index, row in result.iterrows():
            assert (row['nt_pgen_estimate'] - expected['nt_pgen_estimate'][index]) < 0.0000001
