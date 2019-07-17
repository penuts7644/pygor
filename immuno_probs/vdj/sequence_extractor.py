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


"""Contains SequenceExtractor class for extracting sequences from given DataFrame."""


import re

import pandas
import numpy

from immuno_probs.util.processing import multiprocess_array


class SequenceExtractor(object):
    """Extracts the full length (VDJ for productive, unproductive and the total)
    and CDR3 sequences from a given sequence data file.

    Parameters
    ----------
    ref_v_genes : pandas.DataFrame
        A dataframe containing the reference V gene sequences from IMGT as well
        as V family and V gene names.
    ref_j_genes : pandas.DataFrame
        A dataframe containing the reference J gene sequences from IMGT as well
        as J family and J gene names.
    row_id_col : str
        The name of the column containing the row identiefiers.
    nt_col : str
        The name of the nucleotide sequence column to use.
    aa_col : str
        The name of the aminoacid sequence column to use.
    frame_type_col : str
        The name of the column containing the frame type of the sequence.
    cdr3_length_col : str
        The name of the column specifying the length of the CDR3 sequence.
    v_resolved_col : str
        The name of the column containing the resolved V gene name.
    v_gene_choice_col : str
        The name of the V gene choice column to use.
    j_resolved_col : str
        The name of the column containing the resolved J gene name.
    j_gene_choice_col : str
        The name of the J gene choice column to use.

    Methods
    -------
    extract(num_threads, seqs, use_allele=True, default_allele=None)
        Extract and reconstruct sequences and formats them to ImmunoProbs
        compatible.

    """
    def __init__(self, ref_v_genes, ref_j_genes, row_id_col, nt_col, aa_col,
                 frame_type_col, cdr3_length_col, v_resolved_col,
                 v_gene_choice_col, j_resolved_col, j_gene_choice_col):
        super(SequenceExtractor, self).__init__()
        self.ref_v_genes = ref_v_genes
        self.ref_j_genes = ref_j_genes
        self.col_names = {
            'ROW_ID_COL': row_id_col,
            'NT_COL': nt_col,
            'AA_COL': aa_col,
            'FRAME_TYPE_COL': frame_type_col,
            'CDR3_LENGTH_COL': cdr3_length_col,
            'V_RESOLVED_COL': v_resolved_col,
            'V_GENE_CHOICE_COL': v_gene_choice_col,
            'J_RESOLVED_COL': j_resolved_col,
            'J_GENE_CHOICE_COL': j_gene_choice_col,
        }

    @staticmethod
    def _build_resolved_pattern(value, use_allele, default_allele):
        """Private function to split the resolved gene value in an IMGT formated
        regex pattern.

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

        Returns
        -------
        list
            Containing all the 'family-gene*allele' regex pattern combinations
            for the given value string in IMGT format.

        Notes
        -----
        This function assumes the following: If 'use_allele' is True, the allele
        from the input value is used, if this can't be found, 'default_allele'
        is used for that value anyway. It function always requires a family to
        be found. If there is no gene, only the allele is inserted. If there is
        a gene, family, gene and allele are combined. In addition if that gene
        is number 1 and extra pattern is attached with only the family and
        allele.

        """
        # Setup output list, match correct group and split from allele.
        resolved_list = []
        match = re.match(r'^[a-zA-Z]+[V|D|J](.*)$', value)
        match = match.group(1).split('*')
        allele = default_allele

        # When using allele from file update variable or raise error.
        if len(match) == 2 and use_allele:
            allele = match[1]

        # Get possible gene combinations, loop over each and return final list.
        match = match[0].split('/')
        for gene in match:
            tmp_gene = gene.split('-')

            # Only add when there is gene information.
            if len(tmp_gene) == 2:
                tmp_gene_re = [
                    r'\*'.join([
                        r'\-'.join([
                            tmp_gene[0].lstrip('0'),
                            tmp_gene[1].lstrip('0')
                        ]),
                        allele
                    ])]

                # If gene info is '01', also add version without gene name.
                if tmp_gene[1] == '01':
                    tmp_gene_re.append(r'\*'.join([
                        tmp_gene[0].lstrip('0'),
                        allele
                    ]))
                resolved_list.append(r'$|'.join(tmp_gene_re) + '$')
            else:
                resolved_list.append(r'\*'.join([
                    tmp_gene[0].lstrip('0'),
                    allele
                ]) + '$')
        return resolved_list

    @staticmethod
    def _find_longest_substring(full, partial):
        """Finds the longest overlap between a full length sequences and a
        partial length one.

        Parameters
        ----------
        full : str
            A full length sequence string.
        partial : str
            A partial length sequence string to compare to the full length one.

        Returns
        -------
        str
            The longest substring from the compared input strings.

        """
        t = [[0] * (1 + len(partial)) for i in range(1 + len(full))]
        l, xl = 0, 0
        for x in range(1, 1 + len(full)):
            for y in range(1, 1 + len(partial)):
                if full[x - 1] == partial[y - 1]:
                    t[x][y] = t[x - 1][y - 1] + 1
                    if t[x][y] > l:
                        l = t[x][y]
                        xl = x
                else:
                    t[x][y] = 0
        return full[xl - l: xl]

    def _reassemble_data(self, args):
        """Reassembles the dataframe data by striping the CDR3's and creating the
        full length VDJ sequences using recombination.

        Parameters
        ----------
        args : list
            A collection of arguments containing the dataframe to process for
            the thread and default allele value.

        Returns
        -------
        list
            Four pandas dataframes containing the reassembled data, full length
            productive VDJ sequences, full length unproductive VDJ sequences and
            one with the total full length VDJ sequences.

        """
        # Setup the initial dataframe.
        ary, kwargs = args
        use_allele = kwargs['use_allele']
        default_allele = kwargs['default_allele']
        reassembled_df = pandas.DataFrame(columns=[
            self.col_names['ROW_ID_COL'], self.col_names['NT_COL'], self.col_names['AA_COL'],
            self.col_names['V_GENE_CHOICE_COL'], self.col_names['J_GENE_CHOICE_COL']
        ])
        full_length_prod_df = pandas.DataFrame(
            columns=[self.col_names['ROW_ID_COL'], self.col_names['NT_COL']])
        full_length_unprod_df = pandas.DataFrame(
            columns=[self.col_names['ROW_ID_COL'], self.col_names['NT_COL']])
        full_length_df = pandas.DataFrame(
            columns=[self.col_names['ROW_ID_COL'], self.col_names['NT_COL']])

        # Iterate over the rows with index value.
        for i, row in ary.iterrows():

            # Pre-process the resolved V genes.
            v_gene_choices = []
            imgt_v_gene = pandas.DataFrame()
            if isinstance(row[self.col_names['V_RESOLVED_COL']], str):
                v_resolved = self._build_resolved_pattern(
                    value=row[self.col_names['V_RESOLVED_COL']],
                    use_allele=use_allele, default_allele=default_allele)
                for resolved in v_resolved:
                    imgt_v_gene = self.ref_v_genes.loc[self.ref_v_genes[
                        self.col_names['V_RESOLVED_COL']].str.contains(resolved), :]

                    # Assemble the output dataframe gene choices for the V genes.
                    if not imgt_v_gene.empty:
                        v_gene_choices.append(
                            imgt_v_gene[self.col_names['V_RESOLVED_COL']].values[0])
            if v_gene_choices:
                v_gene_choices = '|'.join(v_gene_choices)
            else:
                v_gene_choices = numpy.nan

            # Pre-process the resolved J genes.
            j_gene_choices = []
            imgt_j_gene = pandas.DataFrame()
            if isinstance(row[self.col_names['J_RESOLVED_COL']], str):
                j_resolved = self._build_resolved_pattern(
                    value=row[self.col_names['J_RESOLVED_COL']],
                    use_allele=use_allele, default_allele=default_allele)
                for resolved in j_resolved:
                    imgt_j_gene = self.ref_j_genes.loc[self.ref_j_genes[
                        self.col_names['J_RESOLVED_COL']].str.contains(resolved), :]

                    # Assemble the output dataframe gene choices for the J genes.
                    if not imgt_j_gene.empty:
                        j_gene_choices.append(
                            imgt_j_gene[self.col_names['J_RESOLVED_COL']].values[0])
            if j_gene_choices:
                j_gene_choices = '|'.join(j_gene_choices)
            else:
                j_gene_choices = numpy.nan

            # Create the trimmed NT sequence (removing primers).
            trimmed_nt_seq = row[self.col_names['NT_COL']][
                (81 - int(row[self.col_names['CDR3_LENGTH_COL']])): 81]

            # Add data row of reassembled data to the dataframe.
            reassembled_df = reassembled_df.append({
                self.col_names['ROW_ID_COL']: i,
                self.col_names['NT_COL']: trimmed_nt_seq,
                self.col_names['AA_COL']: row[self.col_names['AA_COL']],
                self.col_names['V_GENE_CHOICE_COL']: v_gene_choices,
                self.col_names['J_GENE_CHOICE_COL']: j_gene_choices,
            }, ignore_index=True)

            # Create the VDJ full length sequence
            if (not imgt_v_gene.empty and not imgt_j_gene.empty):
                vd_segment = self._find_longest_substring(
                    imgt_v_gene[self.col_names['NT_COL']].values[0], trimmed_nt_seq)
                dj_segment = self._find_longest_substring(
                    imgt_j_gene[self.col_names['NT_COL']].values[0], trimmed_nt_seq)
                split_v = imgt_v_gene[self.col_names['NT_COL']].values[0].rsplit(vd_segment, 1)
                split_j = imgt_j_gene[self.col_names['NT_COL']].values[0].split(dj_segment, 1)
                if (len(split_v[1]) >= len(split_v[0])) or (len(split_j[0]) >= len(split_j[1])):
                    continue
                vdj_sequence = split_v[0] + trimmed_nt_seq + split_j[1]
            else:
                continue

            # Add data row of full length data to the dataframe for productive and
            # unproductive sequences.
            if row[self.col_names['FRAME_TYPE_COL']].lower() == 'in':
                full_length_prod_df = full_length_prod_df.append({
                    self.col_names['ROW_ID_COL']: i,
                    self.col_names['NT_COL']: vdj_sequence
                }, ignore_index=True)
            elif (row[self.col_names['FRAME_TYPE_COL']].lower() == 'out'
                  or row[self.col_names['FRAME_TYPE_COL']].lower() == 'stop'):
                full_length_unprod_df = full_length_unprod_df.append({
                    self.col_names['ROW_ID_COL']: i,
                    self.col_names['NT_COL']: vdj_sequence
                }, ignore_index=True)
            full_length_df = full_length_df.append({
                self.col_names['ROW_ID_COL']: i,
                self.col_names['NT_COL']: vdj_sequence
            }, ignore_index=True)
        return reassembled_df, full_length_prod_df, full_length_unprod_df, full_length_df

    def extract(self, num_threads, seqs, use_allele=True, default_allele=None):
        """Restore and extract the full length VDJ and CDR3 sequences from the
        given dataframe.

        The function needs to reassemble the full length VDJ sequences with the
        given reference V and J gene sequences first.

        Parameters
        ----------
        num_threads : int
            The number of threads to use when processing the data.
        seqs : pandas.DataFrame
            Dataframe containing the sequences that need to be extracted.
        use_allele : bool, optional
            If True, the allele information from the input genes is used instead
            of the 'default_allele' value (default: True).
        default_allele : str, optional
            A default allele value to use when spliting gene choices, and
            'use_allele' option is False (default: None).

        Returns
        -------
        list
            Four pandas dataframes containing the reassembled data, full length
            productive VDJ sequences, full length unproductive VDJ sequences and
            one with the total full length VDJ sequences.

        """
        # Set and perform the multiprocessing task.
        results = multiprocess_array(ary=seqs,
                                     func=self._reassemble_data,
                                     num_workers=num_threads,
                                     use_allele=use_allele,
                                     default_allele=default_allele)
        return results
