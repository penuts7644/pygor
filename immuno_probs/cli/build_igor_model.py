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


"""Commandline tool for creating a custom IGoR V(D)J model."""

import os

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

    @staticmethod
    def run(args):
        """Function to execute the commandline tool.

        Parameters
        ----------
        args : Namespace
            Object containing our parsed commandline arguments.

        """
        # -------------------------------
        # NOTE: IGoR has uses the whole header of the reference genome FASTA
        # files as gene name. The IMGT versions of these files includes more
        # information. We have to rewrite the input FASTA files to make it more
        # consistent. This is a quick fix.
        ref_directory = None
        # End ---------------------------

        # Add general igor commands.
        command_list = []
        if args.set_wd:
            command_list.append(['set_wd', str(args.set_wd)])
            # Can be removed ----------------
            ref_directory = os.path.join(str(args.set_wd), 'genomic_templates')
            # End ---------------------------
        else:
            command_list.append(['set_wd', str(get_working_dir())])
            # Can be removed ----------------
            ref_directory = os.path.join(str(get_working_dir()), 'genomic_templates')
            # End ---------------------------
        if args.threads:
            command_list.append(['threads', str(args.threads)])
        else:
            command_list.append(['threads', str(get_num_threads())])

        # Can be removed ----------------
        if not os.path.isdir(ref_directory):
            os.makedirs(ref_directory)
        # End ---------------------------

        # Add sequence and file paths commands.
        if args.ref:
            ref_list = ['set_genomic']
            for i in args.ref:
                # Can be removed ----------------
                records = list(SeqIO.parse(str(i[1]), "fasta"))
                for rec in records:
                    rec.id = rec.description.split('|')[1]
                    rec.description = ""
                name = os.path.join(ref_directory, os.path.basename(str(i[1])))
                SeqIO.write(records, str(name), "fasta")
                # End ---------------------------
                ref_list.append([str(i[0]), str(name)])
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


def main():
    """Function to be called when file executed via terminal."""
    print(__doc__)


if __name__ == "__main__":
    main()
