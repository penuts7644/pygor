# Pygor is part of the IGoR (Inference and Generation of Repertoires) software.
# Pygor Python package can be used to post process files generated by IGoR.
# Copyright (C) 2018 Quentin Marcou & Wout van Helvoirt

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


from Bio.Align.Applications import MuscleCommandline
from StringIO import StringIO
from Bio import AlignIO


class MuscleAligner:
    """Class to perform alignments with MUSCLE via biopython's commandline tool.

    Parameters
    ----------
    fasta_file : string
        A file path to a fasta file containining the genomic data to align.

    """
    def __init__(self, fasta_file):
        super(MuscleAligner, self).__init__()
        self.fasta_file = fasta_file
        self.alignment = None

    def get_muscle_alignment(self):
        """Getter function for collecting the MUSCLE alignment.

        Returns
        -------
        Bio.AlignIO
            A biopython alignment object from the fasta file.

        """
        return self.alignment

    def align(self):
        """Uses MUSCLE via commandline to create a multi-alignment from fasta.

        Notes
        -----
            This function uses the given fasta file for creating an alignment.

        """
        muscle_cline = MuscleCommandline(input=self.fasta_file)
        stdout, stderr = muscle_cline()
        self.alignment = AlignIO.read(StringIO(stdout), "fasta")


def main():
    """Function to be called when file executed via terminal."""
    print(__doc__)


if __name__ == "__main__":
    main()
