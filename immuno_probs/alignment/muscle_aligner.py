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


"""Contains MuscleAligner class to perform MUSCLE alignments."""


from Bio.Align.Applications import MuscleCommandline
from StringIO import StringIO
from Bio import AlignIO
from Bio.Application import ApplicationError

from immuno_probs.util.exception import AlignerException


class MuscleAligner(object):
    """Performs MUSCLE alignments via biopython's commandline tool.

    Parameters
    ----------
    infile : str
        A file path to a FASTA formatted file containining the genomic sequence
        data that is to be aligned against eachother.
    **kwargs
        Optional arguments used for the MuscleCommandline biopython class. Have
        a look at biopython's documenation for more information on the input
        parameters.

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
        """Collects and returns the computed MUSCLE alignment.

        Returns
        -------
        Bio.AlignIO
            A biopython alignment object containing the alignment made from the
            given input FASTA file.

        """
        return self.alignment

    def _align_fasta(self):
        """Executed MUSCLE via commandline to create a multi-alignment from the
        input FASTA file.

        Raises
        ------
        AlignerException
            When the MUSCLE commandline program returns an error.

        Notes
        -----
            This function uses the FASTA file set in the class constructor for
            creating the alignment.

        """
        try:
            muscle_cline = MuscleCommandline(input=self.fasta, **self.kwargs)
            stdout, _ = muscle_cline()
            return AlignIO.read(StringIO(stdout), "fasta")
        except ApplicationError as err:
            raise AlignerException(err.stderr)
