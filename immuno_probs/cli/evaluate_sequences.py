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


import logging
import os
import sys

import pandas
import numpy

from immuno_probs.cdr3.olga_container import OlgaContainer
from immuno_probs.model.default_models import get_default_model_file_paths
from immuno_probs.model.igor_interface import IgorInterface
from immuno_probs.model.igor_loader import IgorLoader
from immuno_probs.util.cli import dynamic_cli_options
from immuno_probs.util.conversion import nucleotides_to_aminoacids
from immuno_probs.util.constant import get_config_data
from immuno_probs.util.io import read_separated_to_dataframe, \
read_fasta_as_dataframe, write_dataframe_to_separated, preprocess_separated_file, \
preprocess_reference_file, is_fasta, is_separated, copy_to_dir


class EvaluateSequences(object):
    """Commandline tool for evaluating sequences using an IGoR model.

    Parameters
    ----------
    subparsers : argparse.ArgumentParser
        A subparser object for appending the tool's parser and options.

    Methods
    -------
    run(args)
        Uses the given Namespace commandline arguments for evaluating sequences.

    """
    def __init__(self, subparsers):
        super(EvaluateSequences, self).__init__()
        self.logger = logging.getLogger(__name__)
        self.subparsers = subparsers
        self._add_options()

    def _add_options(self):
        """Function for adding the parser and options to the given ArgumentParser.

        Notes
        -----
            Uses the class constructor's subparser object for appending the
            tool's parser and options.

        """
        # Create the description and options for the parser.
        description = "Evaluate VDJ or VJ sequences given a custom IGoR " \
            "model (or build-in) through IGoR's commandline tool via python " \
            "subprocess. Or evaluate CDR3 sequences with the model by using " \
            "OLGA."
        parser_options = {
            '-seqs': {
                'metavar': '<fasta/separated>',
                'required': 'True',
                'type': 'str',
                'help': "An input FASTA or separated data file with " \
                        "sequences to evaluate."
            },
            '-model': {
                'type': 'str.lower',
                'choices': get_default_model_file_paths(),
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
                'help': 'If specified (True), CDR3 sequences should be ' \
                        'evaluated, otherwise V(D)J sequences (default: ' \
                        '{}).'.format(
                            get_config_data('EVALUATE', 'EVAL_CDR3', 'bool'))
            },
            '-use-allele': {
                'action': 'store_true',
                'help': "If specified (True), in combination with the '-cdr3' " \
                        " flag, the allele information from the gene choice " \
                        "fields is used to calculate the generation " \
                        "probability (default: {}).".format(
                            get_config_data('EVALUATE', 'USE_ALLELE', 'bool'))
            },
        }

        # Add the options to the parser and return the updated parser.
        parser_tool = self.subparsers.add_parser(
            'evaluate', help=description, description=description)
        parser_tool = dynamic_cli_options(parser=parser_tool,
                                          options=parser_options)

    def run(self, args, output_dir):
        """Function to execute the commandline tool.

        Parameters
        ----------
        args : Namespace
            Object containing our parsed commandline arguments.
        output_dir : str
            A directory path for writing output files to.

        """
        eval_cdr3 = get_config_data('EVALUATE', 'EVAL_CDR3', 'bool')
        if args.cdr3:
            eval_cdr3 = args.cdr3

        # If the given type of sequences evaluation is VDJ, use IGoR.
        if not eval_cdr3:

            # Add general IGoR commands.
            self.logger.info('Setting up initial IGoR command (1/4)')
            command_list = []
            working_dir = get_config_data('COMMON', 'WORKING_DIR')
            command_list.append(['set_wd', working_dir])
            command_list.append(['threads', str(get_config_data('COMMON', 'NUM_THREADS', 'int'))])

            # Add the model (build-in or custom) command depending on given.
            self.logger.info('Processing genomic reference templates (2/4)')
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
            except IOError as err:
                self.logger.error(str(err))
                return

            # Add the sequence command after pre-processing of the input file.
            self.logger.info('Pre-processing input sequence file (3/4)')
            try:
                if is_fasta(args.seqs):
                    self.logger.info('FASTA input file extension detected')
                    command_list.append([
                        'read_seqs',
                        copy_to_dir(working_dir, str(args.seqs), 'fasta')
                    ])
                elif is_separated(args.seqs, get_config_data('COMMON', 'SEPARATOR')):
                    self.logger.info('Separated input file type detected')
                    input_seqs = preprocess_separated_file(
                        os.path.join(working_dir, 'input'),
                        copy_to_dir(working_dir, str(args.seqs), 'csv'),
                        get_config_data('COMMON', 'SEPARATOR'),
                        ';',
                        get_config_data('COMMON', 'I_COL'),
                        [get_config_data('COMMON', 'NT_COL')]
                    )
                    command_list.append(['read_seqs', input_seqs])
                else:
                    self.logger.error(
                        'Given input sequence file could not be detected as ' \
                        'FASTA file or separated data type')
                    return
            except (IOError, KeyError, ValueError) as err:
                self.logger.error(str(err))
                return

            # Add alignment and evealuation commands.
            self.logger.info('Adding additional variables to IGoR command (4/4)')
            command_list.append(['align', ['all']])
            command_list.append(['evaluate'])
            command_list.append(['output', ['Pgen']])

            # Execute IGoR through command line and catch error code.
            self.logger.info('Executing IGoR (this might take a while)')
            try:
                igor_cline = IgorInterface(command=command_list)
                exit_code, _, stderr, _ = igor_cline.call()
                if exit_code != 0:
                    self.logger.error(
                        "An error occurred during execution of IGoR command " \
                        "(exit code %s):\n%s", exit_code, stderr)
                    return
            except OSError as err:
                self.logger.error(str(err))
                return

            # Read in all data frame files, based on input file type.
            self.logger.info('Processing generation probabilities')
            try:
                if is_fasta(args.seqs):
                    seqs_df = read_fasta_as_dataframe(
                        file=args.seqs,
                        col=get_config_data('COMMON', 'NT_COL'))
                elif is_separated(args.seqs, get_config_data('COMMON', 'SEPARATOR')):
                    seqs_df = read_separated_to_dataframe(
                        file=args.seqs,
                        separator=get_config_data('COMMON', 'SEPARATOR'),
                        index_col=get_config_data('COMMON', 'I_COL'))
                full_pgen_df = read_separated_to_dataframe(
                    file=os.path.join(working_dir, 'output', 'Pgen_counts.csv'),
                    separator=';',
                    index_col='seq_index',
                    cols=['Pgen_estimate'])
                full_pgen_df.index.names = [get_config_data('COMMON', 'I_COL')]
                full_pgen_df.rename(
                    columns={'Pgen_estimate': get_config_data('COMMON', 'NT_P_COL')},
                    inplace=True)
                full_pgen_df.loc[:, get_config_data('COMMON', 'AA_P_COL')] = numpy.nan
            except (IOError, KeyError, ValueError) as err:
                self.logger.error(str(err))
                return

            # Insert amino acid sequence column if not existent.
            self.logger.info('Formatting output dataframe')
            if (get_config_data('COMMON', 'NT_COL') in seqs_df.columns
                    and not get_config_data('COMMON', 'AA_COL') in seqs_df.columns):
                seqs_df.insert(
                    seqs_df.columns.get_loc(get_config_data('COMMON', 'NT_COL')) + 1,
                    get_config_data('COMMON', 'AA_COL'), numpy.nan)
                seqs_df[get_config_data('COMMON', 'AA_COL')] = \
                    seqs_df[get_config_data('COMMON', 'NT_COL')].apply(nucleotides_to_aminoacids)

            # Merge IGoR generated sequence output dataframes.
            full_pgen_df = seqs_df.merge(full_pgen_df, left_index=True, right_index=True)

            # Write the pandas dataframe to a separated file.
            try:
                self.logger.info('Writing evaluated data to file system')
                output_filename = get_config_data('COMMON', 'OUT_NAME')
                if not output_filename:
                    output_filename = 'pgen_estimate_{}'.format(model_type)
                _, filename = write_dataframe_to_separated(
                    dataframe=full_pgen_df,
                    filename=output_filename,
                    directory=output_dir,
                    separator=get_config_data('COMMON', 'SEPARATOR'),
                    index_name=get_config_data('COMMON', 'I_COL'))
                self.logger.info("Written '%s'", filename)
            except IOError as err:
                self.logger.error(str(err))
                return

        # If the given type of sequences evaluation is CDR3, use OLGA.
        elif eval_cdr3:

            # Create the directory for the output files.
            working_dir = os.path.join(get_config_data('COMMON', 'WORKING_DIR'), 'output')
            if not os.path.isdir(working_dir):
                os.makedirs(os.path.join(get_config_data('COMMON', 'WORKING_DIR'), 'output'))

            # Load the model and create the sequence evaluator.
            self.logger.info('Loading the IGoR model files')
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
                    separator = get_config_data('COMMON', 'SEPARATOR')
                for gene in args.anchor:
                    anchor_file = preprocess_separated_file(
                        os.path.join(working_dir, 'cdr3_anchors'),
                        str(gene[1]),
                        separator,
                        ','
                    )
                    model.set_anchor(gene=gene[0], file=anchor_file)
                model.initialize_model()
            except (TypeError, OSError, IOError, KeyError, ValueError) as err:
                self.logger.error(str(err))
                return

            # Based on input file type, load in input file.
            self.logger.info('Pre-processing input sequence file')
            try:
                if is_fasta(args.seqs):
                    self.logger.info('FASTA input file extension detected')
                    seqs_df = read_fasta_as_dataframe(
                        file=args.seqs,
                        col=get_config_data('COMMON', 'NT_COL'))
                elif is_separated(args.seqs, get_config_data('COMMON', 'SEPARATOR')):
                    self.logger.info('Separated input file type detected')
                    seqs_df = read_separated_to_dataframe(
                        file=args.seqs,
                        separator=get_config_data('COMMON', 'SEPARATOR'),
                        index_col=get_config_data('COMMON', 'I_COL'))
                else:
                    self.logger.error(
                        'Given input sequence file could not be detected as ' \
                        'FASTA file or separated data type')
                    return
            except (IOError, KeyError, ValueError) as err:
                self.logger.error(str(err))
                return

            # Evaluate the sequences.
            self.logger.info('Evaluating sequences')
            try:
                use_allele = get_config_data('EVALUATE', 'USE_ALLELE', 'bool')
                if args.use_allele:
                    use_allele = args.use_allele
                seq_evaluator = OlgaContainer(
                    igor_model=model,
                    nt_col=get_config_data('COMMON', 'NT_COL'),
                    nt_p_col=get_config_data('COMMON', 'NT_P_COL'),
                    aa_col=get_config_data('COMMON', 'AA_COL'),
                    aa_p_col=get_config_data('COMMON', 'AA_P_COL'),
                    v_gene_choice_col=get_config_data('COMMON', 'V_GENE_CHOICE_COL'),
                    j_gene_choice_col=get_config_data('COMMON', 'J_GENE_CHOICE_COL'))
                cdr3_pgen_df = seq_evaluator.evaluate(
                    seqs=seqs_df,
                    num_threads=get_config_data('COMMON', 'NUM_THREADS', 'int'),
                    use_allele=use_allele,
                    default_allele=get_config_data('EVALUATE', 'DEFAULT_ALLELE'))

                # Merge IGoR generated sequence output dataframes.
                cdr3_pgen_df = seqs_df.merge(cdr3_pgen_df, left_index=True, right_index=True)
            except (TypeError, IOError) as err:
                self.logger.error(str(err))
                return

            # Write the pandas dataframe to a separated file.
            try:
                self.logger.info('Writing evaluated data to file system')
                output_filename = get_config_data('COMMON', 'OUT_NAME')
                if not output_filename:
                    output_filename = 'pgen_estimate_{}_CDR3'.format(model_type)
                _, filename = write_dataframe_to_separated(
                    dataframe=cdr3_pgen_df,
                    filename=output_filename,
                    directory=output_dir,
                    separator=get_config_data('COMMON', 'SEPARATOR'),
                    index_name=get_config_data('COMMON', 'I_COL'))
                self.logger.info("Written '%s'", filename)
            except IOError as err:
                self.logger.error(str(err))
                return


def main():
    """Function to be called when file executed via terminal."""
    print(__doc__)


if __name__ == "__main__":
    main()
