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


"""Commandline tool for evaluating V(D)J sequences using an IGoR model."""


import os
import sys

from immuno_probs.cdr3.olga_container import OlgaContainer
from immuno_probs.model.igor_interface import IgorInterface
from immuno_probs.model.igor_loader import IgorLoader
from immuno_probs.util.cli import dynamic_cli_options
from immuno_probs.util.constant import get_num_threads, get_working_dir, get_separator
from immuno_probs.util.io import read_csv_to_dataframe, write_dataframe_to_csv


class EvaluateSeqs(object):
    """Commandline tool for evaluating sequences using an IGoR model.

    Parameters
    ----------
    subparsers : ArgumentParser
        A subparser object for appending the tool's parser and options.

    Methods
    -------
    run(args)
        Uses the given Namespace commandline arguments for evaluating sequences.

    """
    def __init__(self, subparsers):
        super(EvaluateSeqs, self).__init__()
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
        description = "This tool evaluates V(D)J sequences through IGoR " \
            "commandline python subprocess and CDR3 sequences through OLGA " \
            "python package."
        parser_options = {
            '-seqs': {
                'metavar': '<csv>',
                'required': 'True',
                'type': 'str',
                'help': 'An input CSV file with sequences for evaluation. ' \
                        'Note: uses IGoR file formatting.'
            },
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
            '--anchors': {
                'metavar': ('<v_gene>', '<j_gene>'),
                'type': 'str',
                'nargs': 2,
                'required': '-type=CDR3' in sys.argv or
                            ('-type' in sys.argv and 'CDR3' in sys.argv),
                'help': 'The V and J gene CDR3 anchor files. (required ' \
                        'for -type=CDR3)'
            },
            '--igor-realizations': {
                'metavar': '<csv>',
                'type': 'str',
                'help': 'An optional IGoR realizations file that ' \
                        'corresponds to the given input sequences. ' \
                        '(used when -type=VDJ)'
            },
        }

        # Add the options to the parser and return the updated parser.
        parser_tool = self.subparsers.add_parser(
            'evaluate-seqs', help=description, description=description)
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
        # If the given type of sequences evaluation is VDJ, use IGoR.
        if args.type == 'VDJ':

            # Add general igor commands.
            command_list = []
            directory = get_working_dir()
            command_list.append(['set_wd', str(directory)])
            command_list.append(['threads', str(get_num_threads())])

            # Add the model and sequence commands.
            if args.model:
                command_list.append(['set_custom_model', str(args.model[0]),
                                     str(args.model[1])])
            if args.seqs:
                command_list.append(['read_seqs', str(args.seqs)])

            # Add evaluation commands.
            command_list.append(['evaluate'])
            command_list.append(['output', ['Pgen']])

            igor_cline = IgorInterface(args=command_list)
            code, _ = igor_cline.call()
            if code != 0:
                print("An error occurred during execution of IGoR " \
                      "command (exit code {})".format(code))
                sys.exit()

            # Read in all data frame files, merge them and write out new file.
            sequence_df = read_csv_to_dataframe(
                filename=args.seqs,
                separator=get_separator())
            pgen_df = read_csv_to_dataframe(
                filename=os.path.join(directory, 'output/Pgen_counts.csv'),
                separator=';')
            if args.igor_realizations:
                realizations_df = read_csv_to_dataframe(
                    filename=args.igor_realizations,
                    separator=get_separator())
                gen_df = sequence_df.merge(
                    realizations_df.merge(pgen_df, on='seq_index'),
                    on='seq_index')
            else:
                gen_df = sequence_df.merge(pgen_df, on='seq_index')
            directory, filename = write_dataframe_to_csv(
                dataframe=gen_df,
                filename='output/VDJ_seqs_pgen_estimate',
                directory=directory,
                separator=get_separator())
            os.remove(os.path.join(directory, 'output/Pgen_counts.csv'))

        # If the given type of sequences evaluation is CDR3, use OLGA.
        elif args.type == 'CDR3':

            # Create the directory for the output files.
            directory = os.path.join(get_working_dir(), 'output')
            if not os.path.isdir(directory):
                os.makedirs(os.path.join(get_working_dir(), 'output'))

            # Load the model, create the sequence evaluator and evaluate the sequences.
            model = IgorLoader(
                model_params=args.model[0], model_marginals=args.model[1],
                v_anchors=args.anchors[0], j_anchors=args.anchors[1])
            seq_evaluator = OlgaContainer(igor_model=model)
            sequence_df = read_csv_to_dataframe(filename=args.seqs,
                                                separator=get_separator())
            cdr3_pgen_df = seq_evaluator.evaluate(seqs=sequence_df)

            # Merge IGoR generated sequence outputs and remove old file.
            cdr3_pgen_df_merged = sequence_df.merge(cdr3_pgen_df, on='seq_index')

            # Write the pandas dataframe to a CSV file.
            directory, filename = write_dataframe_to_csv(
                dataframe=cdr3_pgen_df_merged,
                filename='CDR3_seqs_pgen_estimate',
                directory=directory,
                separator=get_separator())

            print("Written '{}' file to '{}' directory."
                  .format(filename, directory))


def main():
    """Function to be called when file executed via terminal."""
    print(__doc__)


if __name__ == "__main__":
    main()
