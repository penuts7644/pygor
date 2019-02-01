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


"""OlgaContainer class for generating and evaluating CDR3 sequences with OLGA."""


import olga.sequence_generation as olga_seq_gen
import olga.generation_probability as olga_pgen
import pandas

from immuno_probs.util.constant import get_num_threads
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
            Containing columns with nucleotide CDR3 sequence - 'nt_sequence',
            ammino acid CDR3 sequence - 'aa_sequence', the index of the chosen
            V gene - 'gene_choice_v' and the index of the chosen J
            gene - 'gene_choice_j'.

        """
        # Create the dataframe and set the generation objects.
        generated_seqs = pandas.DataFrame(columns=['nt_sequence', 'aa_sequence',
                                                   'gene_choice_v', 'gene_choice_j'])
        seq_gen_model = None
        if self.igor_model.is_vdj():
            seq_gen_model = olga_seq_gen.SequenceGenerationVDJ(
                self.igor_model.get_generative_model(),
                self.igor_model.get_genomic_data())
        elif self.igor_model.is_vj():
            seq_gen_model = olga_seq_gen.SequenceGenerationVJ(
                self.igor_model.get_generative_model(),
                self.igor_model.get_genomic_data())
        else:
            raise OlgaException("OLGA could not create a SequenceGeneration object")

        # Generate the sequences, add them to the dataframe and return.
        for _ in range(num_seqs):
            generated_seq = seq_gen_model.gen_rnd_prod_CDR3()
            generated_seqs = generated_seqs.append({
                'nt_sequence': generated_seq[0],
                'aa_sequence': generated_seq[1],
                'gene_choice_v': generated_seq[2],
                'gene_choice_j': generated_seq[3],
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
            GenerationProbability object and the name of the column to generate
            probabilities for.

        Returns
        -------
        pandas.DataFrame
            Containing columns with nucleotide CDR3 sequence - <given by function>
            and the generation probability of the sequences - 'pgen_estimate'.

        """
        # Set the arguments and pandas.DataFrame.
        ary, kwargs = args
        model = kwargs["model"]
        column_name = kwargs["column"]
        pgen_seqs = pandas.DataFrame(columns=['nt_sequence', 'pgen_estimate'])

        # Evaluate the sequences, add them to the dataframe and return.
        for _, row in ary.iterrows():
            seq_pgen = model.compute_nt_CDR3_pgen(row[column_name])
            pgen_seqs = pgen_seqs.append({
                column_name: row[column_name],
                'pgen_estimate': seq_pgen,
            }, ignore_index=True)
        return pgen_seqs

    def evaluate(self, seqs):
        """Evaluate a given nucleotide CDR3 sequence through OLGA.

        Parameters
        ----------
        seqs : pandas.DataFrame
            A pandas dataframe object containing column with nucleotide CDR3
            sequences - 'nt_sequence'.

        Returns
        -------
        pandas.DataFrame
            Containing columns with nucleotide CDR3 sequence - 'nt_sequence'
            and the generation probability of the sequence - 'pgen_estimate'.

        """
        # Set the evaluation objects.
        pgen_model = None
        if self.igor_model.is_vdj():
            pgen_model = olga_pgen.GenerationProbabilityVDJ(
                self.igor_model.get_generative_model(),
                self.igor_model.get_genomic_data())
        elif self.igor_model.is_vj():
            pgen_model = olga_pgen.GenerationProbabilityVJ(
                self.igor_model.get_generative_model(),
                self.igor_model.get_genomic_data())
        else:
            raise OlgaException("OLGA could not create a GenerationProbability object")

        # Use multiprocessing to evaluate the sequences in chunks and return.
        result = multiprocess_array(ary=seqs,
                                    func=self._evaluate,
                                    num_workers=get_num_threads(),
                                    model=pgen_model,
                                    column='nt_sequence')
        result = pandas.concat(result, axis=0).reset_index(drop=True)
        result.drop_duplicates(inplace=True)
        return result


def main():
    """Function to be called when file executed via terminal."""
    print(__doc__)


if __name__ == "__main__":
    main()