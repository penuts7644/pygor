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


"""Commandline tool for evaluating V(D)J sequences using an IGoR model."""


import os
import sys

import pandas

from immuno_probs.cdr3.olga_container import OlgaContainer
from immuno_probs.model.default_models import get_default_model_file_paths
from immuno_probs.model.igor_interface import IgorInterface
from immuno_probs.model.igor_loader import IgorLoader
from immuno_probs.util.cli import dynamic_cli_options
from immuno_probs.util.constant import get_num_threads, get_working_dir, get_separator
from immuno_probs.util.io import read_csv_to_dataframe, read_fasta_as_dataframe, write_dataframe_to_csv, preprocess_input_file, preprocess_reference_file


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
        description = "Evaluate VJ or VDJ sequences given a custom IGoR " \
            "model (or build-in) through IGoR commandline tool via python " \
            "subprocess. Or evaluate CDR3 sequences through OLGA."
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
            '-model': {
                'type': 'str',
                'choices': ['tutorial-model', 'human-t-alpha', 'human-t-beta',
                            'human-b-heavy', 'mouse-t-beta'],
                'required': '-custom-model' not in sys.argv,
                'help': "Specify a pre-installed model for evaluation. " \
                        "(required if --custom-model not specified) " \
                        "(select one: %(choices)s)."
            },
            '-ref': {
                'metavar': ('<gene>', '<fasta>'),
                'type': 'str',
                'action': 'append',
                'nargs': 2,
                'required': ('-cdr3' not in sys.argv and '-custom-model' in sys.argv),
                'help': "A gene (V, D or J) followed by a reference genome " \
                        "FASTA file. Note: the FASTA reference genome files " \
                        "needs to conform to IGMT annotation (separated by " \
                        "'|' character). (required for -custom_model without " \
                        "-cdr3)"
            },
            '-type': {
                'type': 'str',
                'choices': ['VDJ', 'VJ'],
                'required': ('-custom-model' in sys.argv),
                'help': 'The type of model to create. (select one: ' \
                        '%(choices)s) (required for -custom_model).'
            },
            '-custom-model': {
                'metavar': ('<parameters>', '<marginals>'),
                'type': 'str',
                'nargs': 2,
                'help': 'A IGoR parameters txt file followed by an IGoR ' \
                        'marginals txt file.'
            },
            '-anchor': {
                'metavar': ('<gene>', '<csv>'),
                'type': 'str',
                'action': 'append',
                'nargs': 2,
                'required': ('-cdr3' in sys.argv and '-custom-model' in sys.argv),
                'help': 'A gene (V or J) followed by a CDR3 anchor CSV file. ' \
                        'Note: need to contain gene in the firts column, ' \
                        'anchor index in the second and gene function in the ' \
                        'third (required for -cdr3 and -custom_model).'
            },
            '-cdr3': {
                'action': 'store_true',
                'help': 'If specified, CDR3 sequences should be evaluated, ' \
                        'else expecting V(D)J input sequences.'
            },
        }

        # Add the options to the parser and return the updated parser.
        parser_tool = self.subparsers.add_parser(
            'evaluate-seqs', help=description, description=description)
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
        # If the given type of sequences evaluation is VDJ, use IGoR.
        if not args.cdr3:

            # Add general IGoR commands.
            command_list = []
            working_dir = get_working_dir()
            command_list.append(['set_wd', working_dir])
            command_list.append(['threads', str(get_num_threads())])

            # Add the model (build-in or custom) command depending on given.
            if args.model:
                files = get_default_model_file_paths(name=args.model)
                command_list.append(['set_custom_model', files['parameters'],
                                     files['marginals']])
                ref_list = ['set_genomic']
                for gene, filename in files['reference'].items():
                    ref_list.append([gene, filename])
                command_list.append(ref_list)
                if args.model == 'tutorial-model':
                    args.seqs = files['seqs']
            elif args.custom_model:
                command_list.append(['set_custom_model', str(args.custom_model[0]),
                                     str(args.custom_model[1])])
                ref_list = ['set_genomic']
                for i in args.ref:
                    filename = preprocess_reference_file(
                        os.path.join(working_dir, 'genomic_templates'), i[1], 1)
                    ref_list.append([i[0], filename])
                command_list.append(ref_list)

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

            # Add evaluation commands.
            command_list.append(['evaluate'])
            command_list.append(['output', ['Pgen']])

            igor_cline = IgorInterface(args=command_list)
            code, _ = igor_cline.call()
            if code != 0:
                print("An error occurred during execution of IGoR " \
                      "command (exit code {})".format(code))
                return

            # Read in all data frame files.
            sequence_df = read_csv_to_dataframe(
                file=args.seqs,
                separator=get_separator())
            vdj_pgen_df = read_csv_to_dataframe(
                file=os.path.join(working_dir, 'output', 'Pgen_counts.csv'),
                separator=';')

            # Merge IGoR generated sequence output dataframes.
            vdj_pgen_df = sequence_df.merge(vdj_pgen_df, on='seq_index')

            # Write the pandas dataframe to a CSV file.
            directory, filename = write_dataframe_to_csv(
                dataframe=vdj_pgen_df,
                filename='VDJ_seqs_pgen_estimate',
                directory=output_dir,
                separator=get_separator())
            print("Written '{}' file to '{}' directory.".format(
                filename, directory))

        # If the given type of sequences evaluation is CDR3, use OLGA.
        elif args.cdr3:

            # Create the directory for the output files.
            working_dir = os.path.join(get_working_dir(), 'output')
            if not os.path.isdir(working_dir):
                os.makedirs(os.path.join(get_working_dir(), 'output'))

            # Load the model and create the sequence evaluator.
            if args.model:
                files = get_default_model_file_paths(name=args.model)
                model = IgorLoader(model_type=files['type'],
                                   model_params=files['parameters'],
                                   model_marginals=files['marginals'])
                model.set_anchor(gene='V', file=files['v_anchors'])
                model.set_anchor(gene='J', file=files['j_anchors'])
                if args.model == 'tutorial-model':
                    args.seqs = files['cdr3']
            elif args.custom_model:
                model = IgorLoader(model_type=args.type,
                                   model_params=args.custom_model[0],
                                   model_marginals=args.custom_model[1])
                for gene in args.anchor:
                    anchor_file = preprocess_input_file(
                        os.path.join(working_dir, 'cdr3_anchors'), str(gene[1]),
                        get_separator(), ',')
                    model.set_anchor(gene=gene[0], file=anchor_file)
            model.initialize_model()
            seq_evaluator = OlgaContainer(igor_model=model)

            # Based on input file type, load in file.
            if args.seqs.lower().endswith('.csv'):
                sequence_df = read_csv_to_dataframe(file=args.seqs,
                                                    separator=get_separator())
            elif args.seqs.lower().endswith('.fasta'):
                sequence_df = read_fasta_as_dataframe(file=args.seqs)

            # Evaluate the sequences.
            cdr3_pgen_df = seq_evaluator.evaluate(seqs=sequence_df)

            # Merge IGoR generated sequence output dataframes.
            cdr3_pgen_df = sequence_df.merge(cdr3_pgen_df, on='seq_index')

            # Write the pandas dataframe to a CSV file.
            directory, filename = write_dataframe_to_csv(
                dataframe=cdr3_pgen_df,
                filename='CDR3_seqs_pgen_estimate',
                directory=output_dir,
                separator=get_separator())
            print("Written '{}' file to '{}' directory.".format(
                filename, directory))


def main():
    """Function to be called when file executed via terminal."""
    print(__doc__)


if __name__ == "__main__":
    main()
