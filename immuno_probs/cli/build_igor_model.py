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


"""Commandline tool for creating a custom IGoR V(D)J model."""

import os
import sys
from shutil import copy2

from immuno_probs.model.igor_interface import IgorInterface
from immuno_probs.util.cli import dynamic_cli_options
from immuno_probs.util.constant import get_num_threads, get_working_dir, get_separator
from immuno_probs.util.io import preprocess_input_file, preprocess_reference_file


class BuildIgorModel(object):
    """Commandline tool for creating custom IGoR V(D)J models.

    Parameters
    ----------
    subparsers : ArgumentParser
        A subparser object for appending the tool's parser and options.

    Methods
    -------
    run(args)
        Uses the given Namespace commandline arguments to execute IGoR
        for creating a custom model.

    """
    def __init__(self, subparsers):
        super(BuildIgorModel, self).__init__()
        self.subparsers = subparsers
        self._add_options()

    def _add_options(self):
        """Function for adding the parser and options to the given ArgumentParser.

        Notes
        -----
            Uses the class's subparser object for appending the tool's parser
            and options.

        """
        # Create the description and options for the parser.
        description = "Create a VJ or VDJ model by executing IGoR commandline " \
            "tool via a python subprocess and an initial model parameters."
        parser_options = {
            '-seqs': {
                'metavar': '<fasta/csv>',
                'required': 'True',
                'type': 'str',
                'help': "An input FASTA or CSV file with sequences for " \
                        "training the model. Note: file needs to end on " \
                        "'.fasta' or '.csv'. CSV files need to conform to " \
                        "IGoR standards, 'seq_index' and 'nt_sequence' column."
            },
            '-ref': {
                'metavar': ('<gene>', '<fasta>'),
                'type': 'str',
                'action': 'append',
                'nargs': 2,
                'required': 'True',
                'help': "A gene (V, D or J) followed by a reference genome " \
                        "FASTA file. Note: the FASTA reference genome files " \
                        "needs to conform to IGMT annotation (separated by " \
                        "'|' character)."
            },
            '-init-model': {
                'metavar': '<parameters>',
                'required': 'True',
                'type': 'str',
                'help': "An initial IGoR model parameters txt file."
            },
            '-n-iter': {
                'type': 'int',
                'nargs': '?',
                'default': 1,
                'help': 'The number of inference iterations to perform when ' \
                    'creating the model. (default: %(default)s)'
            }
        }

        # Add the options to the parser and return the updated parser.
        parser_tool = self.subparsers.add_parser(
            'build-igor-model', help=description, description=description)
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
        # Add general igor commands.
        command_list = []
        working_dir = get_working_dir()
        command_list.append(['set_wd', working_dir])
        command_list.append(['threads', str(get_num_threads())])

        # Add sequence and file paths commands.
        ref_list = ['set_genomic']
        for i in args.ref:
            filename = preprocess_reference_file(
                os.path.join(working_dir, 'genomic_templates'), i[1], 1)
            ref_list.append([i[0], filename])
        command_list.append(ref_list)
        command_list.append(['set_custom_model', str(args.init_model)])

        # Add the sequence command after pre-processing of the input file.
        if args.seqs.lower().endswith('.csv'):
            input_seqs = preprocess_input_file(
                os.path.join(working_dir, 'input'), str(args.seqs),
                get_separator(), ';', [0, 1])
            command_list.append(['read_seqs', input_seqs])
        elif args.seqs.lower().endswith('.fasta'):
            command_list.append(['read_seqs', str(args.seqs)])

        # Add alignment commands.
        command_list.append(['align', ['all']])

        # Add inference commands.
        command_list.append(['infer', ['N_iter', str(args.n_iter)]])

        igor_cline = IgorInterface(args=command_list)
        code, _ = igor_cline.call()

        if code != 0:
            print("An error occurred during execution of IGoR " \
                  "command (exit code {})".format(code))
            sys.exit()

        # Write output files to output directory.
        copy2(os.path.join(working_dir, 'inference', 'final_marginals.txt'),
              output_dir)
        print("Written '{}' file to '{}' directory.".format(
            'final_marginals.txt', output_dir))

        copy2(os.path.join(working_dir, 'inference', 'final_parms.txt'),
              output_dir)
        print("Written '{}' file to '{}' directory.".format(
            'final_parms.txt', output_dir))


def main():
    """Function to be called when file executed via terminal."""
    print(__doc__)


if __name__ == "__main__":
    main()
