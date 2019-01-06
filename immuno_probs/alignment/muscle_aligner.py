# ImmunoProbs Python package able to calculate the generation probability of
# V(D)J and CDR3 sequences. Copyright (C) 2018 Wout van Helvoirt

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


"""MuscleAligner class for performing MUSCLE alignments."""


from immuno_probs.util.exception import AlignerException

from Bio.Align.Applications import MuscleCommandline
from StringIO import StringIO
from Bio import AlignIO
from Bio.Application import ApplicationError


class MuscleAligner(object):
    """Class to perform alignments with MUSCLE via biopython's commandline tool.

    Parameters
    ----------
    infile : string
        A file path to a FASTA file containining the genomic data to align.
    **kwargs
        Dict of optional arguments for MuscleCommandline biopython class.

    Methods
    -------
    get_muscle_alignment()
        Returns the generated MUSCLE alignment object.

    """
    def __init__(self, infile, **kwargs):
        super(MuscleAligner, self).__init__()
        self.fasta = infile
        self.kwargs = kwargs
        self.alignment = self._align_fasta()

    def get_muscle_alignment(self):
        """Getter function for collecting the MUSCLE alignment.

        Returns
        -------
        Bio.AlignIO
            A biopython alignment object from the fasta file.

        """
        return self.alignment

    def _align_fasta(self):
        """Uses MUSCLE via commandline to create a multi-alignment from fasta.

        Raises
        ------
        AlignerException
            When the Muscle commandline program returns an error.

        Notes
        -----
            This function uses the given fasta file for creating an alignment.

        """
        try:
            muscle_cline = MuscleCommandline(input=self.fasta, **self.kwargs)
            stdout, _ = muscle_cline()
            return AlignIO.read(StringIO(stdout), "fasta")
        except ApplicationError as err:
            raise AlignerException(err.stderr)


def main():
    """Function to be called when file executed via terminal."""
    print(__doc__)


if __name__ == "__main__":
    main()
