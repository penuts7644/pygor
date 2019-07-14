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


import pandas
import numpy

from immuno_probs.util.processing import multiprocess_array


class SequenceExtractor(object):
    """Extracts the full length (VDJ for productive, unproductive and combined)
    and CDR3 sequences from a given data file.

    Parameters
    ----------
    ref_v_genes : pandas.DataFrame
        A dataframe containing the reference V gene sequences from IMGT as well
        as V family and V gene names.
    ref_j_genes : pandas.DataFrame
        A dataframe containing the reference J gene sequences from IMGT as well
        as J family and J gene names.
    nt_col : str
        The name of the nucleotide sequence column to use.
    aa_col : str
        The name of the aminoacid sequence column to use.
    frame_type_col : str
        The name of the column containing the frame type of the sequence.
    cdr3_length_col : str
        The name of the column specifying the length of the CDR3 sequence.
    v_family_col : str
        The name of the column containing the V family name.
    v_gene_col : str
        The name of the column containing the V gene name.
    v_gene_choice_col : str
        The name of the V gene choice column to use.
    j_family_col : str
        The name of the column containing the J family name.
    j_gene_col : str
        The name of the column containing the J gene name.
    j_gene_choice_col : str
        The name of the J gene choice column to use.

    Methods
    -------
    extract(num_threads, seqs)
        Returns a four pandas.DataFrames (CDR3, procuctive, unproductive and
        all) containing the extracted sequences.

    """
    def __init__(self, ref_v_genes, ref_j_genes, i_col, nt_col, aa_col,
                 frame_type_col, cdr3_length_col, v_family_col, v_gene_col,
                 v_gene_choice_col, j_family_col, j_gene_col, j_gene_choice_col):
        super(SequenceExtractor, self).__init__()
        self.ref_v_genes = ref_v_genes
        self.ref_j_genes = ref_j_genes
        self.col_names = {
            'I_COL': i_col,
            'NT_COL': nt_col,
            'AA_COL': aa_col,
            'FRAME_TYPE_COL': frame_type_col,
            'CDR3_LENGTH_COL': cdr3_length_col,
            'V_FAMILY_COL': v_family_col,
            'V_GENE_COL': v_gene_col,
            'V_GENE_CHOICE_COL': v_gene_choice_col,
            'J_FAMILY_COL': j_family_col,
            'J_GENE_COL': j_gene_col,
            'J_GENE_CHOICE_COL': j_gene_choice_col,
        }

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
            the thread.

        Returns
        -------
        DataFrames
            Four pandas datframes containing the reassembled data, full length
            productive VDJ sequences, full length unproductive VDJ sequences and
            one with all full length VDJ sequences.

        Notes
        -----
        The gene name is separated into three sections: family, gene and allele.
        By default, family and gene are required while for the allele '*01' is
        used as default.

        """
        # Setup the initial dataframe.
        ary, _ = args
        reassembled_df = pandas.DataFrame(columns=[
            self.col_names['I_COL'], self.col_names['NT_COL'], self.col_names['AA_COL'],
            self.col_names['V_GENE_CHOICE_COL'], self.col_names['J_GENE_CHOICE_COL']
        ])
        full_length_prod_df = pandas.DataFrame(
            columns=[self.col_names['I_COL'], self.col_names['NT_COL']])
        full_length_unprod_df = pandas.DataFrame(
            columns=[self.col_names['I_COL'], self.col_names['NT_COL']])
        full_length_df = pandas.DataFrame(
            columns=[self.col_names['I_COL'], self.col_names['NT_COL']])

        # Iterate over the rows with index value.
        for i, row in ary.iterrows():

            gene_choice_v = numpy.nan
            imgt_v_gene = pandas.DataFrame()
            if (isinstance(row[self.col_names['V_FAMILY_COL']], str)
                    and isinstance(row[self.col_names['V_GENE_COL']], str)):

                # Pre-process the V gene.
                v_family = row[self.col_names['V_FAMILY_COL']] \
                    .replace('TCR', 'TR').replace('V0', 'V')
                v_gene = row[self.col_names['V_GENE_COL']] \
                    .split('-')[1].split('/')[0].lstrip('0')
                imgt_v_gene = self.ref_v_genes[
                    (self.ref_v_genes[self.col_names['V_FAMILY_COL']] == v_family)
                    & ((self.ref_v_genes[self.col_names['V_GENE_COL']] == v_gene)
                       | ((self.ref_v_genes[self.col_names['V_GENE_COL']].isna())
                          & (v_gene == '1')))].head(1)

                # Assemble the output dataframe gene choices for the V gene.
                if not imgt_v_gene.empty:
                    if (not isinstance(imgt_v_gene[self.col_names['V_GENE_COL']].values[0], str)
                            and v_gene == '1'):
                        gene_choice_v = \
                            imgt_v_gene[self.col_names['V_FAMILY_COL']].values[0] + '*01'
                    else:
                        gene_choice_v = \
                            imgt_v_gene[self.col_names['V_FAMILY_COL']].values[0] + '-' \
                            + imgt_v_gene[self.col_names['V_GENE_COL']].values[0] + '*01'

            gene_choice_j = numpy.nan
            imgt_j_gene = pandas.DataFrame()
            if (isinstance(row[self.col_names['J_FAMILY_COL']], str)
                    and isinstance(row[self.col_names['J_GENE_COL']], str)):

                # Pre-process the J gene.
                j_family = row[self.col_names['J_FAMILY_COL']] \
                    .replace('TCR', 'TR').replace('J0', 'J')
                j_gene = row[self.col_names['J_GENE_COL']] \
                    .split('-')[1].split('/')[0].lstrip('0')
                imgt_j_gene = self.ref_j_genes[
                    (self.ref_j_genes[self.col_names['J_FAMILY_COL']] == j_family)
                    & ((self.ref_j_genes[self.col_names['J_GENE_COL']] == j_gene)
                       | ((self.ref_j_genes[self.col_names['J_GENE_COL']].isna())
                          & (j_gene == '1')))].head(1)

                # Assemble the output dataframe gene choices for the J gene.
                if not imgt_j_gene.empty:
                    if (not isinstance(imgt_j_gene[self.col_names['J_GENE_COL']].values[0], str)
                            and j_gene == '1'):
                        gene_choice_j = \
                            imgt_j_gene[self.col_names['J_FAMILY_COL']].values[0] + '*01'
                    else:
                        gene_choice_j = \
                            imgt_j_gene[self.col_names['J_FAMILY_COL']].values[0] + '-' \
                            + imgt_j_gene[self.col_names['J_GENE_COL']].values[0] + '*01'

            # Create the trimmed NT sequence (removing primers).
            trimmed_nt_seq = row[self.col_names['NT_COL']][
                (81 - int(row[self.col_names['CDR3_LENGTH_COL']])): 81]

            # Add data row of reassembled data to the dataframe.
            reassembled_df = reassembled_df.append({
                self.col_names['I_COL']: i,
                self.col_names['NT_COL']: trimmed_nt_seq,
                self.col_names['AA_COL']: row[self.col_names['AA_COL']],
                self.col_names['V_GENE_CHOICE_COL']: gene_choice_v,
                self.col_names['J_GENE_CHOICE_COL']: gene_choice_j,
            }, ignore_index=True)

            # Create the VDJ full length sequence
            if (not imgt_v_gene.empty
                    and not imgt_j_gene.empty):
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
                    self.col_names['I_COL']: i,
                    self.col_names['NT_COL']: vdj_sequence
                }, ignore_index=True)
            elif (row[self.col_names['FRAME_TYPE_COL']].lower() == 'out'
                  or row[self.col_names['FRAME_TYPE_COL']].lower() == 'stop'):
                full_length_unprod_df = full_length_unprod_df.append({
                    self.col_names['I_COL']: i,
                    self.col_names['NT_COL']: vdj_sequence
                }, ignore_index=True)
            full_length_df = full_length_df.append({
                self.col_names['I_COL']: i,
                self.col_names['NT_COL']: vdj_sequence
            }, ignore_index=True)
        return reassembled_df, full_length_prod_df, full_length_unprod_df, full_length_df

    def extract(self, num_threads, data_df):
        """Restore and extract the full length and CDR3 sequences from the given
        dataframe.

        The function needs to reassemble the full length VDJ sequences with the
        given reference V and J gene sequences first.

        Parameters
        ----------
        num_threads : int
            The number of threads to use when processing the data.
        data_df : pandas.DataFrame
            Dataframe containing the sequences that need to be extracted.

        Returns
        -------
        list
            Containing four pandas.DataFrames: reassembled CDR3 sequences and
            full length sequences (VDJ for productive, unproductive and all).

        """
        # Set and perform the multiprocessing task.
        results = multiprocess_array(ary=data_df,
                                     func=self._reassemble_data,
                                     num_workers=num_threads)
        return results
