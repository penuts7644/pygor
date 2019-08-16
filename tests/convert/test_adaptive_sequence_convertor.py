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


"""Test file for testing immuno_probs.convert.adaptive_sequence_convertor file."""


import pandas
import pytest

from immuno_probs.convert.adaptive_sequence_convertor import AdaptiveSequenceConvertor
from immuno_probs.util.io import read_fasta_as_dataframe, read_separated_to_dataframe


def _process_gene_df(file, nt_col, resolved_col):
    """Read in and process a given reference gene file."""
    gene_df = read_fasta_as_dataframe(
        file=file, col=nt_col, header='info')
    gene_df[nt_col] = gene_df[nt_col].apply(lambda x: x.replace('.', ''))
    gene_df[resolved_col] = pandas.DataFrame(
        list(gene_df['info'].apply(lambda x: x.split('|')[1])))
    gene_df.drop('info', axis=1, inplace=True)
    return gene_df


@pytest.mark.parametrize('value, use_allele, default_allele, expected', [
    ('TRBV5', True, '01', [r'5\*01$']),
    ('TRBV5*09', False, '01', [r'5\*01$']),
    ('TRBV5*09', True, '01', [r'5\*09$']),
    ('TRBV5-02*09', True, '01', [r'5\-2\*09$']),
    ('TRBV5-01*09', True, '01', [r'5\-1\*09$|5\*09$']),
    ('TRBV5-01/4-04*09', True, '01', [r'5\-1\*09$|5\*09$', r'4\-4\*09$']),
])
def test_build_resolved_pattern(value, use_allele, default_allele, expected):
    """Test if correct resolved gene regex patterns are returned.

    Parameters
    ----------
    value : str
        The string value to split and format.
    use_allele : bool
        If True, the allele information from the input genes is used instead
        of the 'default_allele' value.
    default_allele : str
        A default allele value to use when spliting gene choices, and
        'use_allele' option is False.
    expected : pandas.DataFrame
        The expected output pandas.Dataframe with correct columns and values.

    Raises
    -------
    AssertionError
        If the performed test failed.

    """
    asc = AdaptiveSequenceConvertor()
    patterns = asc.build_resolved_pattern(value, use_allele, default_allele)
    for pat in expected:
        assert pat in patterns


@pytest.mark.parametrize('full, partial, expected', [
    ('GATGCTGAAATCACCCAGAGCCCAAGACACAAGATCACAGAGACAGGAAGGCAGGTGACC', 'GAT', 'GAT'),
    ('GATGCTGAAATCACCCAGAGCCCAAGACACAAGATCACAGAGACAGGAAGGCAGGTGACC', 'AGGCA', 'AGGCA'),
])
def test_find_longest_substring(full, partial, expected):
    """Test if correct substring is returned.

    Parameters
    ----------
    full : str
        A full length sequence string.
    partial : str
        A partial length sequence string to compare against the full length
        sequence.
    expected : pandas.DataFrame
        The expected output pandas.Dataframe with correct columns and values.

    Raises
    -------
    AssertionError
        If the performed test failed.

    """
    asc = AdaptiveSequenceConvertor()
    substring = asc.find_longest_substring(full, partial)
    assert substring == expected


@pytest.mark.parametrize('seqs, v_genes, j_genes, expected', [
    ('tests/data/human_t_beta/10_sequence_samples.tsv',
     'tests/data/human_t_beta/ref_genomes/TRBV.fasta',
     'tests/data/human_t_beta/ref_genomes/TRBJ.fasta',
     [6, 4, 2, 6]),
])
def test_convert(seqs, v_genes, j_genes, expected):
    """Test if converted data is returned.

    Parameters
    ----------
    seqs : str
        A filepath to a file containing sequences.
    v_genes : str
        A filepath to a file containing V gene sequences.
    j_genes : str
        A filepath to a file containing J gene sequences.
    expected : pandas.DataFrame
        The expected output pandas.Dataframe with correct columns and values.

    Raises
    -------
    AssertionError
        If the performed test failed.

    """
    seqs = read_separated_to_dataframe(
        file='tests/data/human_t_beta/10_sequence_samples.tsv',
        separator='\t')
    v_genes = _process_gene_df(
        file='tests/data/human_t_beta/ref_genomes/TRBV.fasta',
        nt_col='nt_sequence',
        resolved_col='v_resolved')
    j_genes = _process_gene_df(
        file='tests/data/human_t_beta/ref_genomes/TRBJ.fasta',
        nt_col='nt_sequence',
        resolved_col='j_resolved')
    asc = AdaptiveSequenceConvertor()
    cdr3_df, full_prod_df, full_unprod_df, full_df = asc.convert(
        num_threads=1,
        seqs=seqs,
        ref_v_genes=v_genes,
        ref_j_genes=j_genes,
        row_id_col='row_id',
        nt_col='nt_sequence',
        aa_col='aa_sequence',
        frame_type_col='frame_type',
        cdr3_length_col='cdr3_length',
        v_resolved_col='v_resolved',
        v_gene_choice_col='v_gene_choice',
        j_resolved_col='j_resolved',
        j_gene_choice_col='j_gene_choice',
        use_allele=True,
        default_allele='01')
    assert len(cdr3_df) == expected[0]
    assert len(full_prod_df) == expected[1]
    assert len(full_unprod_df) == expected[2]
    assert len(full_df) == expected[3]
