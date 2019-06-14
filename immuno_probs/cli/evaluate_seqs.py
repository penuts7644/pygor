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
import numpy

from immuno_probs.cdr3.olga_container import OlgaContainer
from immuno_probs.model.default_models import get_default_model_file_paths
from immuno_probs.model.igor_interface import IgorInterface
from immuno_probs.model.igor_loader import IgorLoader
from immuno_probs.util.cli import dynamic_cli_options, make_colored
from immuno_probs.util.conversion import nucleotides_to_aminoacids
from immuno_probs.util.constant import get_config_data
from immuno_probs.util.exception import ModelLoaderException, \
GeneIdentifierException, OlgaException
from immuno_probs.util.io import read_separated_to_dataframe, \
read_fasta_as_dataframe, write_dataframe_to_separated, preprocess_separated_file, \
preprocess_reference_file, is_fasta, is_separated, copy_to_dir


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
        description = "Evaluate VDJ or VJ sequences given a custom IGoR " \
            "model (or build-in) through IGoR's commandline tool via python " \
            "subprocess. Or evaluate CDR3 sequences through OLGA."
        parser_options = {
            '-seqs': {
                'metavar': '<fasta/separated>',
                'required': 'True',
                'type': 'str',
                'help': "An input FASTA or separated data file with " \
                        "sequences for training the model."
            },
            '-model': {
                'type': 'str.lower',
                'choices': ['tutorial-model', 'human-t-alpha', 'human-t-beta',
                            'human-b-heavy', 'mouse-t-beta'],
                'required': '-custom-model' not in sys.argv,
                'help': "Specify a pre-installed model for evaluation. " \
                        "(required if -custom-model NOT specified) " \
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
                        "'|' character). (required for -custom-model without " \
                        "-cdr3)"
            },
            '-type': {
                'type': 'str.lower',
                'choices': ['alpha', 'beta', 'light', 'heavy'],
                'required': ('-custom-model' in sys.argv),
                'help': 'The type of model to create. (select one: ' \
                        '%(choices)s) (required for -custom-model).'
            },
            '-custom-model': {
                'metavar': ('<parameters>', '<marginals>'),
                'type': 'str',
                'nargs': 2,
                'help': 'A IGoR parameters file followed by an IGoR ' \
                        'marginals file.'
            },
            '-anchor': {
                'metavar': ('<gene>', '<separated>'),
                'type': 'str',
                'action': 'append',
                'nargs': 2,
                'required': ('-cdr3' in sys.argv and '-custom-model' in sys.argv),
                'help': 'A gene (V or J) followed by a CDR3 anchor separated ' \
                        'data file. Note: need to contain gene in the first ' \
                        'column, anchor index in the second and gene function ' \
                        'in the third (required for -cdr3 and -custom-model).'
            },
            '-cdr3': {
                'action': 'store_true',
                'help': 'If specified, CDR3 sequences should be evaluated, ' \
                        'else expecting V(D)J input sequences.'
            },
            '-use-cdr3-allele': {
                'action': 'store_true',
                'help': "If specified in combination with the '-cdr3' flag, " \
                        "the allele information from the gene choice fields " \
                        "is used to calculate the generation probability " \
                        "(default: allele '*01' is used for each gene)."
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
            working_dir = get_config_data('WORKING_DIR')
            command_list.append(['set_wd', working_dir])
            command_list.append(['threads', str(get_config_data('NUM_THREADS'))])

            # Add the model (build-in or custom) command depending on given.
            sys.stdout.write('Processing genomic reference templates...')
            try:
                if args.model:
                    files = get_default_model_file_paths(name=args.model)
                    model_type = files['type']
                    command_list.append([
                        'set_custom_model',
                        files['parameters'],
                        files['marginals']
                    ])
                    ref_list = ['set_genomic']
                    for gene, filename in files['reference'].items():
                        ref_list.append([gene, filename])
                    command_list.append(ref_list)
                    if args.model == 'tutorial-model':
                        args.seqs = files['seqs']
                elif args.custom_model:
                    model_type = args.type
                    command_list.append([
                        'set_custom_model',
                        copy_to_dir(working_dir, str(args.custom_model[0]), 'txt'),
                        copy_to_dir(working_dir, str(args.custom_model[1]), 'txt'),
                    ])
                    ref_list = ['set_genomic']
                    for i in args.ref:
                        filename = preprocess_reference_file(
                            os.path.join(working_dir, 'genomic_templates'),
                            copy_to_dir(working_dir, i[1], 'fasta'),
                            1
                        )
                        ref_list.append([i[0], filename])
                    command_list.append(ref_list)
                sys.stdout.write(make_colored('success\n', 'green'))
            except IOError as err:
                sys.stdout.write(make_colored('error\n', 'red'))
                sys.stderr.write(make_colored(str(err) + '\n', 'bg-red'))
                return

            # Add the sequence command after pre-processing of the input file.
            sys.stdout.write('Pre-processing input sequence file...')
            try:
                if is_fasta(args.seqs):
                    sys.stdout.write('(FASTA input file extension detected)...')
                    command_list.append([
                        'read_seqs',
                        copy_to_dir(working_dir, str(args.seqs), 'fasta')
                    ])
                elif is_separated(args.seqs, get_config_data('SEPARATOR')):
                    sys.stdout.write('(separated input file type detected)...')
                    try:
                        input_seqs = preprocess_separated_file(
                            os.path.join(working_dir, 'input'),
                            copy_to_dir(working_dir, str(args.seqs), 'csv'),
                            get_config_data('SEPARATOR'),
                            ';',
                            get_config_data('I_COL'),
                            [get_config_data('NT_COL')]
                        )
                        command_list.append(['read_seqs', input_seqs])
                    except KeyError as err:
                        sys.stdout.write(make_colored('error\n', 'red'))
                        sys.stderr.write(make_colored(
                            "Given input sequence file does not have a '{}' " \
                            "column\n".format(get_config_data('NT_COL')), 'bg-red'))
                        return
                else:
                    sys.stdout.write(make_colored('error\n', 'red'))
                    sys.stderr.write(make_colored(
                        'Given input sequence file could not be detected as ' \
                        'FASTA file or separated data type\n', 'bg-red'))
                    return
                sys.stdout.write(make_colored('success\n', 'green'))
            except IOError as err:
                sys.stdout.write(make_colored('error\n', 'red'))
                sys.stderr.write(make_colored(str(err) + '\n', 'bg-red'))
                return

            # Add alignment commands.
            command_list.append(['align', ['all']])

            # Add evaluation commands.
            command_list.append(['evaluate'])
            command_list.append(['output', ['Pgen']])

            # Execute IGoR through command line and catch error code.
            sys.stdout.write('Executing IGoR...')
            igor_cline = IgorInterface(args=command_list)
            exit_code, _, stderr, _ = igor_cline.call()
            if exit_code != 0:
                sys.stdout.write(make_colored('error\n', 'red'))
                sys.stderr.write(make_colored(
                    "An error occurred during execution of IGoR command (exit " \
                    "code {}):\n{}\n".format(exit_code, stderr), 'bg-red'))
                return
            sys.stdout.write(make_colored('success\n', 'green'))

            # Read in all data frame files, based on input file type.
            sys.stdout.write('Processing generation probabilities...')
            try:
                if is_fasta(args.seqs):
                    seqs_df = read_fasta_as_dataframe(file=args.seqs,
                                                      col=get_config_data('NT_COL'))
                elif is_separated(args.seqs, get_config_data('SEPARATOR')):
                    seqs_df = read_separated_to_dataframe(
                        file=args.seqs, separator=get_config_data('SEPARATOR'),
                        index_col=get_config_data('I_COL'))
                full_pgen_df = read_separated_to_dataframe(
                    file=os.path.join(working_dir, 'output', 'Pgen_counts.csv'),
                    separator=';',
                    index_col=get_config_data('I_COL'),
                    cols=['Pgen_estimate'])
                full_pgen_df.rename(columns={'Pgen_estimate': get_config_data('NT_P_COL')},
                                    inplace=True)
                full_pgen_df.loc[:, get_config_data('AA_P_COL')] = numpy.nan
            except IOError as err:
                sys.stdout.write(make_colored('error\n', 'red'))
                sys.stderr.write(make_colored(str(err) + '\n', 'bg-red'))
                return

            # Insert amino acid sequence column if not existent.
            if (get_config_data('NT_COL') in seqs_df.columns
                    and not get_config_data('AA_COL') in seqs_df.columns):
                seqs_df.insert(
                    seqs_df.columns.get_loc(get_config_data('NT_COL')) + 1,
                    get_config_data('AA_COL'), numpy.nan)
                seqs_df[get_config_data('AA_COL')] = seqs_df[get_config_data('NT_COL')] \
                    .apply(nucleotides_to_aminoacids)

            # Merge IGoR generated sequence output dataframes.
            full_pgen_df = seqs_df.merge(full_pgen_df, left_index=True, right_index=True)
            sys.stdout.write(make_colored('success\n', 'green'))

            # Write the pandas dataframe to a separated file.
            sys.stdout.write('Writting file...')
            try:
                output_filename = get_config_data('OUT_NAME')
                if not output_filename:
                    output_filename = 'pgen_estimate_{}'.format(model_type)
                _, filename = write_dataframe_to_separated(
                    dataframe=full_pgen_df,
                    filename=output_filename,
                    directory=output_dir,
                    separator=get_config_data('SEPARATOR'),
                    index_name=get_config_data('I_COL'))
                sys.stdout.write("(written '{}')...".format(filename))
                sys.stdout.write(make_colored('success\n', 'green'))
            except IOError as err:
                sys.stdout.write(make_colored('error\n', 'red'))
                sys.stderr.write(make_colored(str(err) + '\n', 'bg-red'))
                return

        # If the given type of sequences evaluation is CDR3, use OLGA.
        elif args.cdr3:

            # Create the directory for the output files.
            working_dir = os.path.join(get_config_data('WORKING_DIR'), 'output')
            if not os.path.isdir(working_dir):
                os.makedirs(os.path.join(get_config_data('WORKING_DIR'), 'output'))

            # Load the model and create the sequence evaluator.
            sys.stdout.write('Loading model...')
            try:
                if args.model:
                    files = get_default_model_file_paths(name=args.model)
                    model_type = files['type']
                    model = IgorLoader(model_type=model_type,
                                       model_params=files['parameters'],
                                       model_marginals=files['marginals'])
                    args.anchor = [['V', files['v_anchors']],
                                   ['J', files['j_anchors']]]
                    if args.model == 'tutorial-model':
                        args.seqs = files['cdr3']
                    separator = '\t'
                elif args.custom_model:
                    model_type = args.type
                    model = IgorLoader(model_type=model_type,
                                       model_params=args.custom_model[0],
                                       model_marginals=args.custom_model[1])
                    separator = get_config_data('SEPARATOR')
                for gene in args.anchor:
                    anchor_file = preprocess_separated_file(
                        os.path.join(working_dir, 'cdr3_anchors'),
                        str(gene[1]),
                        separator,
                        ','
                    )
                    model.set_anchor(gene=gene[0], file=anchor_file)
                model.initialize_model()
                sys.stdout.write(make_colored('success\n', 'green'))
            except (ModelLoaderException, GeneIdentifierException, IOError) as err:
                sys.stdout.write(make_colored('error\n', 'red'))
                sys.stderr.write(make_colored(str(err) + '\n', 'bg-red'))
                return

            # Based on input file type, load in input file.
            sys.stdout.write('Pre-processing input sequence file...')
            try:
                if is_fasta(args.seqs):
                    sys.stdout.write('(FASTA input file extension detected)...')
                    seqs_df = read_fasta_as_dataframe(file=args.seqs,
                                                      col=get_config_data('NT_COL'))
                elif is_separated(args.seqs, get_config_data('SEPARATOR')):
                    sys.stdout.write('(separated input file type detected)...')
                    seqs_df = read_separated_to_dataframe(
                        file=args.seqs, separator=get_config_data('SEPARATOR'),
                        index_col=get_config_data('I_COL'))
                else:
                    sys.stdout.write(make_colored('error\n', 'red'))
                    sys.stderr.write(make_colored(
                        'Given input sequence file could not be detected as ' \
                        'FASTA file or separated data type\n', 'bg-red'))
                    return
                sys.stdout.write(make_colored('success\n', 'green'))
            except (IOError) as err:
                sys.stdout.write(make_colored('error\n', 'red'))
                sys.stderr.write(make_colored(str(err) + '\n', 'bg-red'))
                return

            # Evaluate the sequences.
            sys.stdout.write('Evaluating sequences...')
            try:
                seq_evaluator = OlgaContainer(igor_model=model)
                if args.use_cdr3_allele:
                    cdr3_pgen_df = seq_evaluator.evaluate(seqs=seqs_df)
                else:
                    cdr3_pgen_df = seq_evaluator.evaluate(
                        seqs=seqs_df, default_allele=get_config_data('ALLELE'))

                # Merge IGoR generated sequence output dataframes.
                cdr3_pgen_df = seqs_df.merge(cdr3_pgen_df, left_index=True, right_index=True)
                sys.stdout.write(make_colored('success\n', 'green'))
            except (OlgaException, IOError) as err:
                sys.stdout.write(make_colored('error\n', 'red'))
                sys.stderr.write(make_colored(str(err) + '\n', 'bg-red'))
                return

            # Write the pandas dataframe to a separated file.
            sys.stdout.write('Writting file...')
            try:
                output_filename = get_config_data('OUT_NAME')
                if not output_filename:
                    output_filename = 'pgen_estimate_{}_CDR3'.format(model_type)
                _, filename = write_dataframe_to_separated(
                    dataframe=cdr3_pgen_df,
                    filename=output_filename,
                    directory=output_dir,
                    separator=get_config_data('SEPARATOR'),
                    index_name=get_config_data('I_COL'))
                sys.stdout.write("(written '{}')...".format(filename))
                sys.stdout.write(make_colored('success\n', 'green'))
            except IOError as err:
                sys.stdout.write(make_colored('error\n', 'red'))
                sys.stderr.write(make_colored(str(err) + '\n', 'bg-red'))
                return


def main():
    """Function to be called when file executed via terminal."""
    print(__doc__)


if __name__ == "__main__":
    main()
