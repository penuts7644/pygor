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

from Bio import SeqIO

from immuno_probs.model.igor_interface import IgorInterface
from immuno_probs.util.cli import dynamic_cli_options
from immuno_probs.util.constant import get_num_threads, get_working_dir


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
        self.ref_dir = None
        self._add_options()

    def _add_options(self):
        """Function for adding the parser and options to the given ArgumentParser.

        Notes
        -----
            Uses the class's subparser object for appending the tool's parser
            and options.

        """
        # Create the description and options for the parser.
        description = "This tool creates a V(D)J model by executing IGoR " \
            "via a python subprocess. Note: the FASTA reference genome files " \
            "needs to conform to IGMT annotation."
        parser_options = {
            '-seqs': {
                'metavar': '<fasta>',
                'required': 'True',
                'type': 'str',
                'help': 'An input FASTA file with sequences for training ' \
                        'the model.'
            },
            '-ref': {
                'metavar': ('<gene>', '<fasta>'),
                'type': 'str',
                'action': 'append',
                'nargs': 2,
                'required': 'True',
                'help': 'A gene (V, D or J) followed by a reference genome ' \
                        'FASTA file (IMGT).'
            },
            '-model': {
                'metavar': '<parameters>',
                'required': 'True',
                'type': 'str',
                'help': "An initial IGoR model parameters txt file."
            },
            '--n-iter': {
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

    def _format_imgt_reference_fasta(self, fasta):
        """Function for formatting the IMGT reference genome files for IGoR.

        Parameters
        ----------
        fasta : str
            A FASTA file path for a reference genomic template file.

        Returns
        -------
        str
            A string file path to the new reference genomic template file.

        """
        # Create the reference genomic template directory.
        if not os.path.isdir(self.ref_dir):
            os.makedirs(self.ref_dir)

        # Open the fasta file, update the fasta header and write out.
        records = list(SeqIO.parse(str(fasta), "fasta"))
        for rec in records:
            rec.id = rec.description.split('|')[1].split('*')[0]
            rec.description = ""
        updated_path = os.path.join(self.ref_dir, os.path.basename(str(fasta)))
        SeqIO.write(records, str(updated_path), "fasta")
        return updated_path

    def run(self, args):
        """Function to execute the commandline tool.

        Parameters
        ----------
        args : Namespace
            Object containing our parsed commandline arguments.

        """
        # Add general igor commands.
        command_list = []
        directory = get_working_dir()
        command_list.append(['set_wd', str(directory)])
        command_list.append(['threads', str(get_num_threads())])
        self.ref_dir = os.path.join(str(directory), 'genomic_templates')

        # Add sequence and file paths commands.
        if args.ref:
            ref_list = ['set_genomic']
            for i in args.ref:
                filename = self._format_imgt_reference_fasta(i[1])
                ref_list.append([str(i[0]), str(filename)])
            command_list.append(ref_list)
        if args.model:
            command_list.append(['set_custom_model', str(args.model)])
        if args.seqs:
            command_list.append(['read_seqs', str(args.seqs)])

        # Add alignment commands.
        command_list.append(['align', ['all']])

        # Add inference commands.
        if args.n_iter:
            command_list.append(['infer', ['N_iter', str(args.n_iter)]])

        igor_cline = IgorInterface(args=command_list)
        code, _ = igor_cline.call()

        if code != 0:
            print("An error occurred during execution of IGoR " \
                  "command (exit code {})".format(code))
            sys.exit()


def main():
    """Function to be called when file executed via terminal."""
    print(__doc__)


if __name__ == "__main__":
    main()
