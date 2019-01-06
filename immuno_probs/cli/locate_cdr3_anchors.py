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


"""Commandline tool for creating files containing CDR3 anchor indices."""


import os

from immuno_probs.alignment.muscle_aligner import MuscleAligner
from immuno_probs.cdr3.anchor_locator import AnchorLocator
from immuno_probs.util.cli import dynamic_cli_options
from immuno_probs.util.constant import get_working_dir, get_separator
from immuno_probs.util.io import write_dataframe_to_csv


class LocateCdr3Anchors(object):
    """Commandline tool for creating compatible CDR3 anchor CSV files.

    Parameters
    ----------
    subparsers : ArgumentParser
        A subparser object for appending the tool's parser and options.

    Methods
    -------
    run(args)
        Uses the given Namespace commandline arguments to locate the CDR3
        anchors and write them to a file.

    """
    def __init__(self, subparsers):
        super(LocateCdr3Anchors, self).__init__()
        self.subparsers = subparsers
        self._add_options()

    def _add_options(self):
        """Function for adding the parser/options to the input ArgumentParser.

        Notes
        -----
            Uses the class's subparser object for appending the tool's parser
            and options.

        """
        # Create the description and options for the parser.
        description = "This tool creates an alignment from the given " \
            "reference genome FASTA file and seaches the given alignment for " \
            "conserved motif regions. The located regions are written out to " \
            "a CSV file. Note: the FASTA needs to conform to IGMT annotation."
        parser_options = {
            '-ref': {
                'metavar': ('<gene>', '<fasta>'),
                'type': 'str',
                'action': 'append',
                'nargs': 2,
                'required': 'True',
                'help': 'A gene (V or J) followed by a reference genome ' \
                        'FASTA file (IMGT).'
            },
            '--motifs': {
                'type': 'str',
                'nargs': '*',
                'help': "The motifs to look for. If none specified, the " \
                        "default 'V' (Cystein - TGT and TGC) or 'J' " \
                        "(Tryptophan - TGG, Phenylalanine - TTT and TTC)."
            }
        }

        # Add the options to the parser and return the updated parser.
        parser_tool = self.subparsers.add_parser(
            'locate-cdr3-anchors', help=description, description=description)
        parser_tool = dynamic_cli_options(parser=parser_tool,
                                          options=parser_options)


    @staticmethod
    def run(args):
        """Function to execute the commandline tool.

        Parameters
        ----------
        args : Namespace
            Object containing our parsed commandline arguments.

        """
        # Create the directory for the output files.
        directory = os.path.join(get_working_dir(), 'cdr3_anchors')
        if not os.path.isdir(directory):
            os.makedirs(os.path.join(get_working_dir(), 'cdr3_anchors'))

        # Create the alignment and locate the motifs.
        for gene in args.ref:
            aligner = MuscleAligner(infile=gene[1])
            locator = AnchorLocator(alignment=aligner.get_muscle_alignment(),
                                    gene=gene[0])

            if args.motifs is not None:
                anchors_df = locator.get_indices_motifs(args.motifs)
            else:
                anchors_df = locator.get_indices_motifs()

            # Write the pandas dataframe to a CSV file.
            directory, filename = write_dataframe_to_csv(
                dataframe=anchors_df,
                filename='{}_gene_CDR3_anchors'.format(gene[0]),
                directory=directory,
                separator=get_separator())

            print("Written '{}' file to '{}' directory."
                  .format(filename, directory))


def main():
    """Function to be called when file executed via terminal."""
    print(__doc__)


if __name__ == "__main__":
    main()
