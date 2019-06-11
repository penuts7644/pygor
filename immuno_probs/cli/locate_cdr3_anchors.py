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


"""Commandline tool for creating files containing CDR3 anchor indices."""


import os
import sys

import numpy

from immuno_probs.alignment.muscle_aligner import MuscleAligner
from immuno_probs.cdr3.anchor_locator import AnchorLocator
from immuno_probs.util.cli import dynamic_cli_options
from immuno_probs.util.constant import get_config_data
from immuno_probs.util.exception import AlignerException, GeneIdentifierException
from immuno_probs.util.io import copy_to_dir, preprocess_reference_file, \
write_dataframe_to_separated


class LocateCdr3Anchors(object):
    """Commandline tool for creating compatible CDR3 anchor separated files.

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
        description = "Create an alignment for the given reference genome " \
            "FASTA file and seach the given alignment for conserved motif " \
            "regions."
        parser_options = {
            '-ref': {
                'metavar': ('<gene>', '<fasta>'),
                'type': 'str',
                'action': 'append',
                'nargs': 2,
                'required': 'True',
                'help': "A gene (V or J) followed by a reference genome " \
                        "FASTA file. Note: the FASTA reference genome files " \
                        "needs to conform to IGMT annotation (separated by " \
                        "'|' character)."
            },
            '-motif': {
                'type': 'str.upper',
                'action': 'append',
                'help': "The motifs to look for. If none specified, the " \
                        "default 'V' (Cystein - TGT and TGC) or 'J' " \
                        "(Tryptophan - TGG, Phenylalanine - TTC and TTT)."
            }
        }

        # Add the options to the parser and return the updated parser.
        parser_tool = self.subparsers.add_parser(
            'locate-cdr3-anchors', help=description, description=description)
        parser_tool = dynamic_cli_options(parser=parser_tool,
                                          options=parser_options)

    @staticmethod
    def run(args, output_dir):
        """Function to execute the commandline tool.

        Parameters
        ----------
        args : Namespace
            Object containing our parsed commandline arguments.
        output_dir : str
            A directory path for writing output files to.

        """
        # Get the working directory.
        working_dir = get_config_data('WORKING_DIR')

        # Create the alignment and locate the motifs.
        for gene in args.ref:
            sys.stdout.write('Locating CDR3 anchors for {}...'.format(gene[0]))
            filename = preprocess_reference_file(
                os.path.join(working_dir, 'genomic_templates'),
                copy_to_dir(working_dir, gene[1], 'fasta'),
            )
            try:
                aligner = MuscleAligner(infile=filename)
                locator = AnchorLocator(alignment=aligner.get_muscle_alignment(),
                                        gene=gene[0])
            except (AlignerException, GeneIdentifierException) as err:
                sys.stdout.write('\033[91merror\033[0m\n')
                sys.stderr.write(str(err) + '\n')
                return

            if args.motif is not None:
                anchors_df = locator.get_indices_motifs(*args.motif)
            else:
                anchors_df = locator.get_indices_motifs()

            # Modify the dataframe to make it OLGA compliant.
            anchors_df.insert(2, 'function', numpy.nan)
            anchors_df.rename(columns={'name': 'gene'}, inplace=True)
            try:
                anchors_df['gene'], anchors_df['function'] = zip(*anchors_df['gene'].apply(
                    lambda value: (value.split('|')[1], value.split('|')[3])))
            except (IndexError, ValueError):
                sys.stdout.write('\033[91merror\033[0m\n')
                sys.stderr.write("FASTA header needs to be separated by '|', " \
                    "needs to have gene name on index position 1 and function " \
                    "on index position 3: '{}'\n".format(anchors_df['gene']))
                return

            # Apply some filtering to the anchor dataframe.
            anchors_df.drop_duplicates(subset=['gene'], keep='first', inplace=True)
            anchors_df.reset_index(inplace=True, drop=True)

            # Write the pandas dataframe to a separated file with prefix.
            output_prefix = get_config_data('OUT_NAME')
            if not output_prefix:
                output_prefix = 'gene_CDR3_anchors'
            _, filename = write_dataframe_to_separated(
                dataframe=anchors_df,
                filename='{}_{}'.format(gene[0], output_prefix),
                directory=output_dir,
                separator=get_config_data('SEPARATOR'))
            sys.stdout.write("(written '{}' for {} gene)...".format(filename, gene[0]))
            sys.stdout.write('\033[92msuccess\033[0m\n')


def main():
    """Function to be called when file executed via terminal."""
    print(__doc__)


if __name__ == "__main__":
    main()
