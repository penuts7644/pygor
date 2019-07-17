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


"""Test file for testing immuno_probs.cdr3.olga_container file."""


import pandas
import pytest
import numpy

from immuno_probs.cdr3.olga_container import OlgaContainer
from immuno_probs.model.igor_loader import IgorLoader


@pytest.mark.parametrize('option, expected', [
    ('generate', [['nt_sequence', str], ['aa_sequence', str],
                  ['v_gene_choice', str], ['j_gene_choice', str]]),
    ('evaluate', pandas.DataFrame(
        [['TGTGCCAGTAGTATAACAACCCAGGGCTTGTACGAGCAGTACTTC', numpy.nan, 0, numpy.nan],
         ['TGTGCAGGAATAAACTTTGGAAATGAGAAATTAACCTTT', numpy.nan, 6.022455403460228e-08, numpy.nan],
         ['TGTGCATTGAACAGAGATGACAAGATCATCTTT', numpy.nan, 3.8690138672702246e-07, numpy.nan]],
        columns=['nt_sequence', 'aa_sequence', 'nt_pgen_estimate', 'aa_pgen_estimate'])
    )
])
def test_olga_container(option, expected):
    """Test if the container class can generate and evaluate the CDR3's.

    Parameters
    ----------
    option : str
        A option to run in the container ('generate' of 'evaluate').
    expected : pandas.DataFrame
        The expected output pandas.Dataframe with correct columns and values.

    Raises
    -------
    AssertionError
        If the performed test failed.

    """
    params = 'tests/data/human_t_alpha/model_params.txt'
    marginals = 'tests/data/human_t_alpha/model_marginals.txt'
    v_anchors = 'tests/data/human_t_alpha/V_gene_CDR3_anchors.csv'
    j_anchors = 'tests/data/human_t_alpha/J_gene_CDR3_anchors.csv'
    model = IgorLoader(model_type='alpha', model_params=params, model_marginals=marginals)
    model.set_anchor(gene='V', file=v_anchors)
    model.set_anchor(gene='J', file=j_anchors)
    model.initialize_model()
    olga_container = OlgaContainer(
        igor_model=model,
        nt_col='nt_sequence',
        nt_p_col='nt_pgen_estimate',
        aa_col='aa_sequence',
        aa_p_col='aa_pgen_estimate',
        v_gene_choice_col='v_gene_choice',
        j_gene_choice_col='j_gene_choice')
    if option == 'generate':
        result = olga_container.generate(num_seqs=1)
        for i in expected:
            assert isinstance(i[1](result.iloc[0][i[0]]), i[1])
    elif option == 'evaluate':
        pgen_seqs = expected.drop(['nt_pgen_estimate', 'aa_pgen_estimate'], axis=1)
        result = olga_container.evaluate(seqs=pgen_seqs, num_threads=1)
        for index, row in result.iterrows():
            assert (row['nt_pgen_estimate'] - expected['nt_pgen_estimate'][index]) < 0.0000001
