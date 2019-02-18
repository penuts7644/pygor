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


"""AnchorLocator class for locating CDR3 anchor indices of given sequences."""


import pandas
import numpy

from immuno_probs.util.constant import get_num_threads
from immuno_probs.util.exception import GeneIdentifierException, IndexNotFoundException
from immuno_probs.util.processing import multiprocess_array


class AnchorLocator(object):
    """Class for locating CDR3 anchors within the given nucleotide alignment.

    Parameters
    ----------
    alignment : Bio.AlignIO
        An biopython MUSCLE alignement object from alignment.MuscleAligner.
    gene : string
        A gene identifier, either V or J, specifying the alignement's origin.

    Methods
    -------
    get_indices_motifs(*motifs)
        Returns the indices dictionarys for each of the given motif string.

    Notes
    -----
        This class followes OLGA's CSV formatting standard for the outputed
        CDR3 anchor files. The reference genome alignment's gene name should
        be located in the second position and the function on the fourth
        position in the sequence header (separated by '|').

    """
    def __init__(self, alignment, gene):
        super(AnchorLocator, self).__init__()
        self.alignment = alignment
        self.gene = self._set_gene(gene)
        self.default_motifs = {"V": ["TGT", "TGC"],
                               "J": ["TGG", "TTT", "TTC"]}

    @staticmethod
    def _set_gene(gene):
        """Private setter function for setting the gene identifier.

        Parameters
        ----------
        gene : string
            A gene identifier, either V or J, specifying the alignment's origin.

        Returns
        -------
        str
            The gene character if passing the validation tests.

        Raises
        ------
        GeneIdentifierException
            When gene character is not 'V' or 'J'.

        """
        gene = gene.upper()
        if gene not in ["V", "J"]:
            raise GeneIdentifierException(
                "Gene identifier needs can be either 'V' or 'J'", gene)
        return gene

    @staticmethod
    def _find_conserved_motif_indices(args):
        """Find the most conserved motif region within the MUSCLE alignment.

        This function finds conserved motif regions using the provided V or J
        gene multi-alignment. This is done for each given motif in the input
        list.

        Parameters
        ----------
        args : list
            The arguments from the multiprocess_array function. Consists of an
            list and additional kwargs like the Bio.AlignIO alignment object.

        Returns
        -------
        pandas.DataFrame
            Containing start index values for each sequence identifier in the
            alignment. Each motif has its own row in the dataframe.

        Raises
        ------
        IndexNotFoundException
            When the sequence header is not splitable by '|', if the header
            index 1 (gene name) or if the header index 3 (function) can be
            found.

        """
        # Set the arguments and pandas.DataFrame.
        ary, kwargs = args
        alignment = kwargs["alignment"]
        seq_motif_indices = pandas.DataFrame(columns=['gene', 'anchor_index',
                                                      'function', 'motif'])

        # For each of the motifs in the input array.
        for motif in ary:

            # Loop over alignment (codon len) and collect occurences of motif.
            motif_index_occurances = []
            for i in range(0, alignment.get_alignment_length() - len(motif)):
                motif_counts = numpy.zeros(len(alignment))
                alignment_codon = alignment[:, i:i + len(motif)]

                # For the motif alignment, count motif occurences and add to
                # the counts.
                for seq_record, j in zip(alignment_codon,
                                         range(0, len(alignment_codon))):
                    motif_counts[j] = (seq_record.seq == motif)

                # Calculate average of occurences (between 0 and 1) and add to
                # start index.
                motif_index_occurances.append(
                    float(sum(motif_counts)) / len(alignment_codon))

            # Collect index with highest value attached.
            max_index = numpy.argmax(motif_index_occurances)
            for seq_record in alignment:

                # Only process sequences that contain the motif at the conserved
                # index location.
                if seq_record.seq[max_index:max_index + len(motif)] == motif:
                    try:
                        description_list = seq_record.description.split('|')
                        start_index = len(str(seq_record.seq[0:max_index])
                                          .replace('-', ''))
                        seq_motif_indices = seq_motif_indices.append({
                            'gene': description_list[1],
                            'anchor_index': start_index,
                            'function': description_list[3],
                            'motif': motif,
                        }, ignore_index=True)
                    except IndexError:
                        raise IndexNotFoundException(
                            "FASTA header needs to be separated by '|', " \
                            "needs to have gene name on index 1 and function " \
                            "on index 3", seq_record.description)
        return seq_motif_indices

    def get_indices_motifs(self, *motifs):
        """Collect the conserved indices in the multi-alignment for each motif.

        Parameters
        ----------
        *motif : strings
            One or multiple motif strings to process. If none given, default V/J
            gene motifs are located.

        Returns
        -------
        pandas.DataFrame
            Containing columns with sequence identifiers - 'gene', start index
            values for the anchors - 'anchor_index', sequence function -
            'function' and motifs - 'motif'.

        Notes
        -----
            This function uses the given MUSCLE alignment and gene identifier.
            It locates the most common 'V' (Cysteine - TGT and TGC) or 'J'
            (Tryptophan - TGG, Phenylalanine - TTT and TTC) index that covers
            all sequences in the multi-alignment (Or custom motifs).
            This function uses the NUM_THREADS variable for multiprocessing.

        """
        # Set the motifs arrays and perform the multiprocessing task.
        if not motifs:
            motifs = self.default_motifs[self.gene]
        result = multiprocess_array(ary=motifs,
                                    func=self._find_conserved_motif_indices,
                                    num_workers=get_num_threads(),
                                    alignment=self.alignment)
        result = pandas.concat(result, axis=0, ignore_index=True, copy=False)
        result.drop_duplicates(inplace=True)
        result.reset_index(inplace=True, drop=True)
        return result


def main():
    """Function to be called when file executed via terminal."""
    print(__doc__)


if __name__ == "__main__":
    main()
