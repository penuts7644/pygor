# ImmunoProbs Python package uses a simplified manner for calculating the
# generation probability of V(D)J and CDR3 sequences.
# Copyright (C) 2018 Wout van Helvoirt

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

from immuno_probs.model.igor_interface import IgorInterface
from immuno_probs.util.cli import dynamic_cli_options
from immuno_probs.util.constant import get_num_threads


class CreateIgorModel(object):
    """Commandline tool for creating custom IGoR V(D)J models.

    Parameters
    ----------
    subparsers : ArgumentParser
        A subparser object for appending the tool's parser and options.

    Methods
    -------
    run(args)
        Uses the given Namespace commandline arguments to locate the run IGoR
        for creating a custom model. As of now, still requires an initial model.

    """
    def __init__(self, subparsers):
        super(CreateIgorModel, self).__init__()
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
            "via a python subprocess."
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
                        'FASTA file.'
            },
            '-model-params': {
                'metavar': '<txt>',
                'required': 'True',
                'type': 'str',
                'help': "An initial IGoR model's parameters file."
            },
            '--set-wd': {
                'type': 'str',
                'nargs': '?',
                'help': 'An optional location for creating the IGoR files. ' \
                        'By default, uses the current directory for ' \
                        'written files.'
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
            'create-igor-model', help=description, description=description)
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
        # Add general igor commands.
        command_list = []
        if args.set_wd:
            command_list.append(['set_wd', str(args.set_wd)])
        else:
            command_list.append(['set_wd', str(os.getcwd())])
        if args.threads:
            command_list.append(['threads', str(args.threads)])
        else:
            command_list.append(['threads', str(get_num_threads())])

        # Add sequence and file paths commands.
        if args.ref:
            ref_list = ['set_genomic']
            for i in args.ref:
                ref_list.append([str(i[0]), str(i[1])])
            command_list.append(ref_list)
        if args.model_params:
            command_list.append(['set_custom_model', str(args.model_params)])
        if args.seqs:
            command_list.append(['read_seqs', str(args.seqs)])

        # Add alignment commands.
        command_list.append(['align', ['all']])

        # Add inference commands.
        if args.n_iter:
            command_list.append(['infer', ['N_iter', str(args.n_iter)]])

        igor_cline = IgorInterface(args=command_list)
        code, stdout, stderr, _ = igor_cline.call()

        if code != 0:
            print("An error occureud during for execution of IGoR command: \n")
            print("stderr:\n{}".format(stderr))
            print("stdout:\n{}".format(stdout))


def main():
    """Function to be called when file executed via terminal."""
    print(__doc__)


if __name__ == "__main__":
    main()
