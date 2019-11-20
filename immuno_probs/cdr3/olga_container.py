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


"""Contains OlgaContainer class for generating and evaluating CDR3 sequences."""


import olga.sequence_generation as olga_seq_gen
import olga.generation_probability as olga_pgen
import pandas
import numpy

from immuno_probs.util.conversion import nucleotides_to_aminoacids
from immuno_probs.util.processing import multiprocess_array


class OlgaContainer(object):
    """Generates and/or evaluates CDR3 sequences using given an IGoR model.

    Parameters
    ----------
    igor_model : immuno_probs.model.igor_loader.IgorLoader
        IgorLoader object containing the loaded IGoR VJ or VDJ model.
    nt_col : str
        The name of the nucleotide sequence column to use.
    nt_p_col : str
        The name of the nucleotide Pgen column to use.
    aa_col : str
        The name of the aminoacid sequence column to use.
    aa_p_col : str
        The name of the aminoacid Pgen column to use.
    v_gene_choice_col : str
        The name of the V gene choice column to use.
    j_gene_choice_col : str
        The name of the J gene choice column to use.


    Methods
    -------
    generate(num_seqs)
        Returns pandas.DataFrame with nucleotide and aminoacid CDR3 sequences.
    evaluate(seq, num_threads, use_allele=True, default_allele=None)
        Returns the generation probability value for the given sequences.

    """
    def __init__(self, igor_model, nt_col, nt_p_col, aa_col, aa_p_col, v_gene_choice_col, j_gene_choice_col):
        super(OlgaContainer, self).__init__()
        self.igor_model = igor_model
        self.col_names = {
            'NT_COL': nt_col,
            'NT_P_COL': nt_p_col,
            'AA_COL': aa_col,
            'AA_P_COL': aa_p_col,
            'V_GENE_CHOICE_COL': v_gene_choice_col,
            'J_GENE_CHOICE_COL': j_gene_choice_col,
        }

    def generate(self, num_seqs):
        """Generate a given number of CDR3 sequences through OLGA.

        Parameters
        ----------
        num_seqs : int
            An integer specifying the number of sequences to generate.

        Returns
        -------
        pandas.DataFrame
            Containing columns with sequence index, nucleotide CDR3 sequence, amino acid CDR3 sequence, the index of the
            chosen V gene and the index of the chosen J gene.

        Raises
        ------
        TypeError
            When the model type does not equal 'VDJ' or 'VJ'.

        """
        # Create the dataframe and set the generation objects.
        generated_seqs = pandas.DataFrame(
            columns=[self.col_names['NT_COL'], self.col_names['AA_COL'],
                     self.col_names['V_GENE_CHOICE_COL'],
                     self.col_names['J_GENE_CHOICE_COL']])
        seq_gen_model = None
        if self.igor_model.get_type() == "VDJ":
            seq_gen_model = olga_seq_gen.SequenceGenerationVDJ(
                self.igor_model.get_generative_model(),
                self.igor_model.get_genomic_data())
        elif self.igor_model.get_type() == "VJ":
            seq_gen_model = olga_seq_gen.SequenceGenerationVJ(
                self.igor_model.get_generative_model(),
                self.igor_model.get_genomic_data())
        else:
            raise TypeError("OLGA could not create a SequenceGeneration object since model is not of type 'VDJ' or 'VJ'")

        # Generate the sequences, add them to the dataframe and return.
        for _ in range(num_seqs):
            generated_seq = seq_gen_model.gen_rnd_prod_CDR3()
            generated_seqs = generated_seqs.append({
                self.col_names['NT_COL']: generated_seq[0],
                self.col_names['AA_COL']: generated_seq[1],
                self.col_names['V_GENE_CHOICE_COL']: self.igor_model.get_genomic_data().genV[generated_seq[2]][0],
                self.col_names['J_GENE_CHOICE_COL']: self.igor_model.get_genomic_data().genJ[generated_seq[3]][0]
            }, ignore_index=True)
        return generated_seqs

    @staticmethod
    def _locate_genes(genes, ref_genes, use_allele, default_allele):
        """Locates all the given gene values in the reference gene list.

        If a gene family value is specified instead of the whole gene (family + gene identifier), all possible genes within
        that family are located.

        Parameters
        ----------
        genes : list
            Containing gene string values that need to be located.
        ref_genes : list
            Containing reference gene string values.
        use_allele : bool
            If True, the allele information from the input genes is used instead of the 'default_allele' value.
        default_allele : str
            A default allele value to use when spliting gene choices, and 'use_allele' option is False.

        Returns
        -------
        list
            A list with the genes that where located in the reference genes list. If no genes where found, an empty list
            is returned.

        """
        # For each given gene, split up the name into family, gene and allele.
        located_genes = set()
        allele = default_allele
        for name in genes:
            name = name.split('*')
            name[0] = name[0].split('-')
            family, gene, allele = [None] * 3
            if len(name[0]) == 2:
                family, gene = name[0][0], name[0][1]
            else:
                family = name[0][0]
            if len(name) == 2 and use_allele:
                allele = name[1]

            # Collect the subsection of the genes using the reference genes.
            if family and not gene:
                if allele:
                    located_genes.update(
                        [i for i in ref_genes
                         if family in i and '*' + allele in i])
                else:
                    located_genes.update([i for i in ref_genes if family in i])
            elif family and gene:
                if allele:
                    located_genes.update(
                        [i for i in ref_genes
                         if family + '-' + gene in i and '*' + allele in i])
                else:
                    located_genes.update(
                        [i for i in ref_genes if family + '-' + gene in i])
        return list(located_genes)

    def _evaluate(self, args):
        """Private function for evaluating a given nucleotide CDR3 sequence by using OLGA.

        Parameters
        ----------
        args : list
            The arguments from the 'multiprocess_array' function. Consists of an pandas.DataFrame and additional kwargs like
            a GenerationProbability object, the column name containing the nucleotide sequences and value to use as allele information.

        Returns
        -------
        pandas.DataFrame
            Containing columns sequence index number, the generation probability of nucleotide sequence if given and the
            generation probability of aminoacid sequence if given.

        """
        # Set the arguments and pandas.DataFrame.
        ary, kwargs = args
        model = kwargs["model"]
        use_allele = kwargs["use_allele"]
        default_allele = kwargs["default_allele"]
        ref_genes_v = [i[0] for i in self.igor_model.get_genomic_data().genV]
        ref_genes_j = [i[0] for i in self.igor_model.get_genomic_data().genJ]
        pgen_seqs = pandas.DataFrame(
            index=ary.index.tolist(),
            columns=[self.col_names['NT_P_COL'], self.col_names['AA_P_COL']])

        for i, row in ary.iterrows():

            # Evaluate the sequences with V/J gene columns.
            if ((self.col_names['V_GENE_CHOICE_COL'] in ary.columns
                 and isinstance(row[self.col_names['V_GENE_CHOICE_COL']], str))
                    and (self.col_names['J_GENE_CHOICE_COL'] in ary.columns
                         and isinstance(row[self.col_names['J_GENE_CHOICE_COL']], str))):

                # Create all V/J gene combinations for pgen calculation.
                located_v = self._locate_genes(
                    genes=row[self.col_names['V_GENE_CHOICE_COL']].split('|'),
                    ref_genes=ref_genes_v, use_allele=use_allele,
                    default_allele=default_allele)
                located_j = self._locate_genes(
                    genes=row[self.col_names['J_GENE_CHOICE_COL']].split('|'),
                    ref_genes=ref_genes_j, use_allele=use_allele,
                    default_allele=default_allele)
                permutations = [(v, j) for v in located_v for j in located_j]

                # For the nucleotide sequence if exists.
                if (self.col_names['NT_COL'] in ary.columns
                        and isinstance(row[self.col_names['NT_COL']], str)):
                    sum_pgen = 0
                    for v_gene, j_gene in permutations:
                        sum_pgen += model.compute_nt_CDR3_pgen(
                            row[self.col_names['NT_COL']], v_gene, j_gene)
                    pgen_seqs.loc[i, :][self.col_names['NT_P_COL']] = sum_pgen

                # For the amino acid sequence if exists.
                if (self.col_names['AA_COL'] in ary.columns
                        and isinstance(row[self.col_names['AA_COL']], str)):
                    sum_pgen = 0
                    for v_gene, j_gene in permutations:
                        sum_pgen += model.compute_aa_CDR3_pgen(
                            row[self.col_names['AA_COL']], v_gene, j_gene)
                    pgen_seqs.loc[i, :][self.col_names['AA_P_COL']] = sum_pgen

            # If no V/J gene choice column, use less complicated method.
            else:

                # For the nucleotide sequence if exists.
                if (self.col_names['NT_COL'] in ary.columns
                        and isinstance(row[self.col_names['NT_COL']], str)):
                    pgen_seqs.loc[i, :][self.col_names['NT_P_COL']] = \
                        model.compute_nt_CDR3_pgen(row[self.col_names['NT_COL']])

                # For the amino acid sequence if exists.
                if (self.col_names['AA_COL'] in ary.columns
                        and isinstance(row[self.col_names['AA_COL']], str)):
                    pgen_seqs.loc[i, :][self.col_names['AA_P_COL']] = \
                        model.compute_aa_CDR3_pgen(row[self.col_names['AA_COL']])
        return pgen_seqs

    def evaluate(self, seqs, num_threads, use_allele=True, default_allele=None):
        """Evaluate a given nucleotide CDR3 sequences using OLGA.

        This function also checks if the given input sequence file contains the gene index columns for the V and J gene.
        If so, then the V and J gene masks in these columns are used to increase calculation accuracy of the generation
        probabality values.

        Parameters
        ----------
        seqs : pandas.DataFrame
            A pandas dataframe object containing a column with nucleotide CDR3 sequences and/or amino acid sequences.
        num_threads : int
            The number of threads to use when processing the sequences.
        use_allele : bool, optional
            If True, the allele information from the input genes is used instead of the 'default_allele' value (default: True).
        default_allele : str, optional
            A default allele value to use when spliting gene choices, and 'use_allele' option is False (default: None).

        Returns
        -------
        pandas.DataFrame
            Containing the sequence index number, the generation probability of nucleotide sequence if given and the
            generation probability of aminoacid sequence if given.

        Raises
        ------
        TypeError
            When the model type does not equal 'VDJ' or 'VJ'.

        """
        # Set the evaluation objects.
        pgen_model = None
        if self.igor_model.get_type() == "VDJ":
            pgen_model = olga_pgen.GenerationProbabilityVDJ(
                self.igor_model.get_generative_model(),
                self.igor_model.get_genomic_data())
        elif self.igor_model.get_type() == "VJ":
            pgen_model = olga_pgen.GenerationProbabilityVJ(
                self.igor_model.get_generative_model(),
                self.igor_model.get_genomic_data())
        else:
            raise TypeError("OLGA could not create a GenerationProbability object since model is not of type 'VDJ' or 'VJ'")

        # Insert amino acid sequence column if not existent.
        if (self.col_names['NT_COL'] in seqs.columns
                and not self.col_names['AA_COL'] in seqs.columns):
            seqs.insert(seqs.columns.get_loc(self.col_names['NT_COL']) + 1,
                        self.col_names['AA_COL'], numpy.nan)
            seqs[self.col_names['AA_COL']] = seqs[self.col_names['NT_COL']] \
                .apply(nucleotides_to_aminoacids)

        # Use multiprocessing to evaluate the sequences in chunks and return.
        result = multiprocess_array(
            ary=seqs,
            func=self._evaluate,
            num_workers=num_threads,
            model=pgen_model,
            use_allele=use_allele,
            default_allele=default_allele)
        result = pandas.concat(result, axis=0, copy=False)
        return result
