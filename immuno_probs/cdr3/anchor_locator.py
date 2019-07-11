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


"""Contains AnchorLocator class for locating CDR3 anchors of given sequences."""


import pandas
import numpy

from immuno_probs.util.processing import multiprocess_array


class AnchorLocator(object):
    """Locates all possibile CDR3 start and end positions within the given
    nucleotide alignment.

    This class followes OLGA's formatting standard for the outputed CDR3 anchor
    files. The gene names in the reference genomic sequence alignment should be
    located in the second position and the function of the gene on the fourth
    position in the sequence header (separated by '|'). The header needs to
    follow IMGT standards.

    Parameters
    ----------
    alignment : Bio.AlignIO
        An biopython MUSCLE alignement object created with the
        alignment.MuscleAligner class.
    gene : str
        A gene identifier, either V or J, specifying the alignment's origin
        gene.

    Methods
    -------
    get_indices_motifs(num_threads, *motifs)
        Returns a pandas.DataFrame containing CDR3 anchors.

    """
    def __init__(self, alignment, gene):
        super(AnchorLocator, self).__init__()
        self.alignment = alignment
        self.gene = self._set_gene(gene)

    @staticmethod
    def _set_gene(gene):
        """Private function for setting the gene identifier value in the class.

        Parameters
        ----------
        gene : str
            A gene identifier, either V or J, specifying the alignment's origin
            gene.

        Returns
        -------
        str
            The input gene character value after passing the validation test.

        Raises
        ------
        ValueError
            When the given gene character does not equal 'V' or 'J'.

        """
        gene = gene.upper()
        if gene not in ["V", "J"]:
            raise ValueError(
                "Gene identifier should be either 'V' or 'J'", gene)
        return gene

    @staticmethod
    def _find_conserved_motif_indices(args):
        """Find the most conserved motif region within the MUSCLE alignment.

        The regions are located for each given motif using the provided V or
        J gene sequence alignment.

        Parameters
        ----------
        args : list
            The arguments from the 'multiprocess_array' function. Consists of an
            list and additional kwargs with the Bio.AlignIO alignment object.

        Returns
        -------
        pandas.DataFrame
            Containing start index values for each sequence identifier in the
            alignment. Each motif has its own row in the dataframe.

        """
        # Set the arguments and pandas.DataFrame.
        ary, kwargs = args
        alignment = kwargs["alignment"]
        seq_motif_indices = pandas.DataFrame(
            columns=['name', 'anchor_index', 'motif'])

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
                    start_index = len(
                        str(seq_record.seq[0:max_index]).replace('-', ''))
                    seq_motif_indices = seq_motif_indices.append({
                        'name': seq_record.description,
                        'anchor_index': start_index,
                        'motif': motif,
                    }, ignore_index=True)
        return seq_motif_indices

    def get_indices_motifs(self, num_threads, *motifs):
        """Collects and returns the CDR3 anchors for each motif from the
        sequence alignment.

        The function locates the most common V (Cysteine - TGT and TGC by
        default) or J (Tryptophan - TGG, Phenylalanine - TTC and TTT by
        default) gene index that best covers all sequences in the alignment.

        Parameters
        ----------
        num_threads : int
            The number of threads to use when processing the motif list.
        *motifs : str
            Various motif strings (combined into a list) to process. At least
            motif needs to be specified.

        Returns
        -------
        pandas.DataFrame
            Containing columns with sequence identifier, start index values for
            the anchors, function of the gene and motifs.

        Raises
        ------
        ValueError
            When no motifs have been specified in the function call.

        Notes
        -----
            This function uses the given MUSCLE alignment and gene identifier
            from the class constructor.

        """
        # Set the motifs arrays and perform the multiprocessing task.
        if not motifs:
            raise ValueError(
                'No motifs are given in function call, one is required')
        result = multiprocess_array(ary=motifs,
                                    func=self._find_conserved_motif_indices,
                                    num_workers=num_threads,
                                    alignment=self.alignment)
        result = pandas.concat(result, axis=0, ignore_index=True, copy=False)
        result.drop_duplicates(inplace=True)
        result.reset_index(inplace=True, drop=True)
        return result
