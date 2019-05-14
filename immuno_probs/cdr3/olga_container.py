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


"""OlgaContainer class for generating and evaluating CDR3 sequences with OLGA."""


import olga.sequence_generation as olga_seq_gen
import olga.generation_probability as olga_pgen
import pandas
import numpy

from immuno_probs.util.constant import get_config_data
from immuno_probs.util.conversion import nucleotides_to_aminoacids
from immuno_probs.util.exception import OlgaException
from immuno_probs.util.processing import multiprocess_array


class OlgaContainer(object):
    """Class generating and evaluating CDR3 sequences using an IGoR model.

    Parameters
    ----------
    igor_model : IgorLoader object
        IgorLoader object containing the loaded IGoR VJ or VDJ model.

    Methods
    -------
    generate(num_seqs)
        Returns pandas.DataFrame with NT and AA CDR3 sequences.
    evaluate(seq)
        Returns the generation probability value for the given sequence.

    """
    def __init__(self, igor_model):
        super(OlgaContainer, self).__init__()
        self.igor_model = igor_model

    def generate(self, num_seqs):
        """Generate a given number of CDR3 sequences through OLGA.

        Parameters
        ----------
        num_seqs : int
            An integer specifying the number of sequences to generate.

        Returns
        -------
        pandas.DataFrame
            Containing columns with sequence index, nucleotide CDR3 sequence,
            ammino acid CDR3 sequence, the index of the chosen V gene and the
            index of the chosen J gene.

        """
        # Create the dataframe and set the generation objects.
        generated_seqs = pandas.DataFrame(
            columns=[get_config_data('NT_COL'), get_config_data('AA_COL'),
                     get_config_data('V_GENE_COL'), get_config_data('J_GENE_COL')])
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
            raise OlgaException("OLGA could not create a SequenceGeneration object")

        # Generate the sequences, add them to the dataframe and return.
        for _ in range(num_seqs):
            generated_seq = seq_gen_model.gen_rnd_prod_CDR3()
            generated_seqs = generated_seqs.append({
                get_config_data('NT_COL'): generated_seq[0],
                get_config_data('AA_COL'): generated_seq[1],
                get_config_data('V_GENE_COL'): self.igor_model.get_genomic_data() \
                    .genV[generated_seq[2]][0],
                get_config_data('J_GENE_COL'): self.igor_model.get_genomic_data() \
                    .genJ[generated_seq[3]][0]
            }, ignore_index=True)
        return generated_seqs

    @staticmethod
    def _evaluate(args):
        """Evaluate a given nucleotide CDR3 sequence through OLGA.

        Parameters
        ----------
        args : list
            The arguments from the multiprocess_array function. Consists of an
            pandas.DataFrame and additional kwargs like the
            GenerationProbability object, the column name containing the
            nucleotide sequences and a boolean for usingf V/J masks.

        Returns
        -------
        pandas.DataFrame
            Containing columns sequence index number, the generation probability
            of nucleotide sequence if given and the generation probability of
            aminoacid sequence if given.

        """
        # Set the arguments and pandas.DataFrame.
        ary, kwargs = args
        model = kwargs["model"]
        allele = kwargs["allele"]
        pgen_seqs = pandas.DataFrame(
            index=ary.index.tolist(),
            columns=[get_config_data('NT_P_COL'), get_config_data('AA_P_COL')])

        for i, row in ary.iterrows():

            # Evaluate the sequences with V/J gene columns.
            if ((get_config_data('V_GENE_COL') in ary.columns
                 and isinstance(row[get_config_data('V_GENE_COL')], str))
                    and (get_config_data('J_GENE_COL') in ary.columns
                         and isinstance(row[get_config_data('J_GENE_COL')], str))):
                permutations = [
                    (v_gene, j_gene)
                    for v_gene in row[get_config_data('V_GENE_COL')].split('|')
                    for j_gene in row[get_config_data('J_GENE_COL')].split('|')
                ]

                # For the nucleotide sequence if exists.
                if (get_config_data('NT_COL') in ary.columns
                        and isinstance(row[get_config_data('NT_COL')], str)):
                    sum_pgen = 0
                    for v_gene, j_gene in permutations:
                        sum_pgen += model.compute_nt_CDR3_pgen(
                            row[get_config_data('NT_COL')], v_gene, j_gene)
                    pgen_seqs.loc[i, :][get_config_data('NT_P_COL')] = sum_pgen

                # For the amino acid sequence if exists.
                if (get_config_data('AA_COL') in ary.columns
                        and isinstance(row[get_config_data('AA_COL')], str)):
                    sum_pgen = 0
                    for v_gene, j_gene in permutations:
                        sum_pgen += model.compute_aa_CDR3_pgen(
                            row[get_config_data('AA_COL')], v_gene, j_gene)
                    pgen_seqs.loc[i, :][get_config_data('AA_P_COL')] = sum_pgen

            # If no V/J gene choice column, use less complicated method.
            else:

                # For the nucleotide sequence if exists.
                if (get_config_data('NT_COL') in ary.columns
                        and isinstance(row[get_config_data('NT_COL')], str)):
                    pgen_seqs.loc[i, :][get_config_data('NT_P_COL')] = \
                        model.compute_nt_CDR3_pgen(row[get_config_data('NT_COL')])

                # For the amino acid sequence if exists.
                if (get_config_data('AA_COL') in ary.columns
                        and isinstance(row[get_config_data('AA_COL')], str)):
                    pgen_seqs.loc[i, :][get_config_data('AA_P_COL')] = \
                        model.compute_aa_CDR3_pgen(row[get_config_data('AA_COL')])
        return pgen_seqs

    def evaluate(self, seqs, allele=None):
        """Evaluate a given nucleotide CDR3 sequence through OLGA.

        Parameters
        ----------
        seqs : pandas.DataFrame
            A pandas dataframe object containing column with nucleotide CDR3
            sequences and/or amino acid sequences.
        allele : str, optional
            An allele value to use when calculating the Pgen for each single
            gene combination. If given, only the given allele is used
            (default: use allele info from gene choice column if available).

        Returns
        -------
        pandas.DataFrame
            Containing columns sequence index number, the generation probability
            of nucleotide sequence if given and the generation probability of
            aminoacid sequence if given.

        Notes
        -----
        This function also checks if the given input sequence file contains the
        gene index columns for the V and J gene. If so, then the V and J gene
        masks in these columns are used to incease accuracy of the
        generation probabality.

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
            raise OlgaException("OLGA could not create a GenerationProbability object")

        # Insert amino acid sequence column if not existent.
        if (get_config_data('NT_COL') in seqs.columns
                and not get_config_data('AA_COL') in seqs.columns):
            seqs.insert(seqs.columns.get_loc(get_config_data('NT_COL')) + 1,
                        get_config_data('AA_COL'), numpy.nan)
            seqs[get_config_data('AA_COL')] = seqs[get_config_data('NT_COL')] \
                .apply(nucleotides_to_aminoacids)

        # Use multiprocessing to evaluate the sequences in chunks and return.
        result = multiprocess_array(ary=seqs,
                                    func=self._evaluate,
                                    num_workers=get_config_data('NUM_THREADS'),
                                    model=pgen_model,
                                    allele=allele)
        result = pandas.concat(result, axis=0, copy=False)
        return result


def main():
    """Function to be called when file executed via terminal."""
    print(__doc__)


if __name__ == "__main__":
    main()
