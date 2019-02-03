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


"""Commandline tool for generating V(D)J sequences from and IGoR model."""


import os
import sys

from immuno_probs.cdr3.olga_container import OlgaContainer
from immuno_probs.model.igor_interface import IgorInterface
from immuno_probs.model.igor_loader import IgorLoader
from immuno_probs.util.cli import dynamic_cli_options
from immuno_probs.util.constant import get_num_threads, get_working_dir, get_separator
from immuno_probs.util.io import read_csv_to_dataframe, write_dataframe_to_csv


class GenerateSeqs(object):
    """Commandline tool for generating sequences from and IGoR model.

    Parameters
    ----------
    subparsers : ArgumentParser
        A subparser object for appending the tool's parser and options.

    Methods
    -------
    run(args)
        Uses the given Namespace commandline arguments for generating sequences.

    """
    def __init__(self, subparsers):
        super(GenerateSeqs, self).__init__()
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
        description = "This tool generates V(D)J sequences given an IGoR " \
            "model by executing IGoR via a python subprocess and gnerates " \
            "CDR3 sequences by using the OLGA python package."
        parser_options = {
            '-model': {
                'metavar': ('<parameters>', '<marginals>'),
                'type': 'str',
                'nargs': 2,
                'required': 'True',
                'help': 'A IGoR parameters txt file followed by an IGoR ' \
                        'marginals txt file.'
            },
            '-type': {
                'type': 'str',
                'choices': ['CDR3', 'VDJ'],
                'required': 'True',
                'help': 'The type of sequences to generate. (select one: ' \
                        '%(choices)s)'
            },
            '--generate': {
                'type': 'int',
                'nargs': '?',
                'default': 1,
                'help': 'The number of sequences to generate. (default: ' \
                        '%(default)s)'
            },
            '--anchors': {
                'metavar': ('<v_gene>', '<j_gene>'),
                'type': 'str',
                'nargs': 2,
                'required': '-type=CDR3' in sys.argv or
                            ('-type' in sys.argv and 'CDR3' in sys.argv),
                'help': 'The V and J gene CDR3 anchor files. (required ' \
                        'for -type=CDR3)'
            },
        }

        # Add the options to the parser and return the updated parser.
        parser_tool = self.subparsers.add_parser(
            'generate-seqs', help=description, description=description)
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
        # If the given type of sequences generation is VDJ, use IGoR.
        if args.type == 'VDJ':

            # Add general igor commands.
            command_list = []
            directory = get_working_dir()
            command_list.append(['set_wd', str(directory)])
            command_list.append(['threads', str(get_num_threads())])

            # Add the model command.
            if args.model:
                command_list.append(['set_custom_model', str(args.model[0]),
                                     str(args.model[1])])

            # Add generate command.
            if args.generate:
                command_list.append(['generate', str(args.generate)])

            # Execute IGoR through command line and catch error code.
            igor_cline = IgorInterface(args=command_list)
            code, _ = igor_cline.call()
            if code != 0:
                print("An error occurred during execution of IGoR " \
                      "command (exit code {})".format(code))

            # Merge IGoR generated sequence outputs and remove old files.
            sequence_df = read_csv_to_dataframe(
                filename=os.path.join(
                    directory, 'generated/generated_seqs_werr.csv'),
                separator=';')
            realizations_df = read_csv_to_dataframe(
                filename=os.path.join(
                    directory, 'generated/generated_realizations_werr.csv'),
                separator=';')
            gen_df = sequence_df.merge(realizations_df, on='seq_index')
            directory, filename = write_dataframe_to_csv(
                dataframe=gen_df,
                filename='generated/generated_VDJ_seqs',
                directory=directory,
                separator=';')
            os.remove(os.path.join(
                directory, 'generated/generated_seqs_werr.csv'))
            os.remove(os.path.join(
                directory, 'generated/generated_realizations_werr.csv'))

        # If the given type of sequences generation is CDR3, use OLGA.
        elif args.type == 'CDR3':

            # Create the directory for the output files.
            directory = os.path.join(get_working_dir(), 'generated')
            if not os.path.isdir(directory):
                os.makedirs(os.path.join(get_working_dir(), 'generated'))

            # Load the model, create the sequence generator and generate the sequences.
            model = IgorLoader(
                model_params=args.model[0], model_marginals=args.model[1],
                v_anchors=args.anchors[0], j_anchors=args.anchors[1])
            seq_generator = OlgaContainer(igor_model=model)
            cdr3_seqs_df = seq_generator.generate(num_seqs=args.generate)

            # Write the pandas dataframe to a CSV file.
            directory, filename = write_dataframe_to_csv(
                dataframe=cdr3_seqs_df,
                filename='generated_CDR3_seqs',
                directory=directory,
                separator=get_separator())

            print("Written '{}' file to '{}' directory."
                  .format(filename, directory))


def main():
    """Function to be called when file executed via terminal."""
    print(__doc__)


if __name__ == "__main__":
    main()
