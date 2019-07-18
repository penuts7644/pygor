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


"""Commandline tool for generating V(D)J sequences from and IGoR model."""


import os
import sys

import pandas

from immuno_probs.cdr3.olga_container import OlgaContainer
from immuno_probs.model.default_models import get_default_model_file_paths
from immuno_probs.model.igor_interface import IgorInterface
from immuno_probs.model.igor_loader import IgorLoader
from immuno_probs.util.cli import dynamic_cli_options, make_colored
from immuno_probs.util.conversion import nucleotides_to_aminoacids
from immuno_probs.util.constant import get_config_data
from immuno_probs.util.io import read_separated_to_dataframe, \
write_dataframe_to_separated, preprocess_separated_file, copy_to_dir


class GenerateSequences(object):
    """Commandline tool for generating sequences from and IGoR model.

    Parameters
    ----------
    subparsers : argparse.ArgumentParser
        A subparser object for appending the tool's parser and options.

    Methods
    -------
    run(args)
        Uses the given Namespace commandline arguments for generating sequences.

    """
    def __init__(self, subparsers):
        super(GenerateSequences, self).__init__()
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
        description = "Generate VDJ or VJ sequences given a custom IGoR " \
            "model (or build-in) by executing IGoR's commandline tool via " \
            "python subprocess. Or generate CDR3 sequences from the model by " \
            "using OLGA."
        parser_options = {
            '-model': {
                'type': 'str.lower',
                'choices': get_default_model_file_paths(),
                'required': '-custom-model' not in sys.argv,
                'help': "Specify a pre-installed model for generation. " \
                        "(required if -custom-model NOT specified) " \
                        "(select one: %(choices)s)."
            },
            '-type': {
                'type': 'str.lower',
                'choices': ['alpha', 'beta', 'light', 'heavy'],
                'required': ('-custom-model' in sys.argv),
                'help': 'The type of model to create. (select one: ' \
                        '%(choices)s) (required for -custom-model).'
            },
            '-anchor': {
                'metavar': ('<gene>', '<separated>'),
                'type': 'str',
                'action': 'append',
                'nargs': 2,
                'required': ('-cdr3' in sys.argv and '-custom-model' in sys.argv),
                'help': 'A gene (V or J) followed by a CDR3 anchor separated '
                        'data file. Note: need to contain gene in the first ' \
                        'column, anchor index in the second and gene function ' \
                        'in the third (required for -cdr3 and -custom-model).'
            },
            '-custom-model': {
                'metavar': ('<parameters>', '<marginals>'),
                'type': 'str',
                'nargs': 2,
                'help': 'A IGoR parameters file followed by an IGoR ' \
                        'marginals file.'
            },
            '-n-gen': {
                'type': 'int',
                'nargs': '?',
                'default': 1,
                'help': 'The number of sequences to generate. (default: ' \
                        '%(default)s)'
            },
            '-cdr3': {
                'action': 'store_true',
                'help': 'If specified, CDR3 sequences are generated, else ' \
                        'V(D)J full length sequences.'
            },
        }

        # Add the options to the parser and return the updated parser.
        parser_tool = self.subparsers.add_parser(
            'generate', help=description, description=description)
        parser_tool = dynamic_cli_options(parser=parser_tool,
                                          options=parser_options)

    @staticmethod
    def _process_realizations(data, model, v_gene_choice_col,
                              d_gene_choice_col, j_gene_choice_col):
        """Function for processing an IGoR realization dataframe with indices.

        Parameters
        ----------
        data : pandas.DataFrame
            A pandas dataframe object with the IGoR realization data.
        model : immuno_probs.model.igor_loader.IgorLoader
            Object containing the IGoR model.
        v_gene_choice_col : str
            The name of the V gene choice column to use.
        d_gene_choice_col : str
            The name of the D gene choice column to use.
        j_gene_choice_col : str
            The name of the J gene choice column to use.

        Returns
        -------
        pandas.DataFrame
            A pandas dataframe object with sequence index, the V(D)J gene
            choice columns containing the names of the selected genes.

        """
        # If the suplied model is VDJ, locate important columns and update index values.
        if model.get_type() == "VDJ":
            real_df = pandas.concat([data.filter(regex=("GeneChoice_V_gene_.*")),
                                     data.filter(regex=("GeneChoice_J_gene_.*")),
                                     data.filter(regex=("GeneChoice_D_gene_.*"))],
                                    axis=1, sort=False)
            real_df.columns = [v_gene_choice_col, j_gene_choice_col, d_gene_choice_col]
            real_df[v_gene_choice_col], \
            real_df[j_gene_choice_col], \
            real_df[d_gene_choice_col] = zip(
                *real_df.apply(lambda row: (
                    model.get_genomic_data().genV[int(row[v_gene_choice_col].strip('()'))][0],
                    model.get_genomic_data().genJ[int(row[j_gene_choice_col].strip('()'))][0],
                    model.get_genomic_data().genD[int(row[d_gene_choice_col].strip('()'))][0]
                ), axis=1))

        # Or do the same if the model is VJ.
        elif model.get_type() == "VJ":
            real_df = pandas.concat([data.filter(regex=("GeneChoice_V_gene_.*")),
                                     data.filter(regex=("GeneChoice_J_gene_.*"))],
                                    axis=1, sort=False)
            real_df.columns = [v_gene_choice_col, j_gene_choice_col]
            real_df[v_gene_choice_col], \
            real_df[j_gene_choice_col] = zip(
                *real_df.apply(lambda row: (
                    model.get_genomic_data().genV[int(row[v_gene_choice_col].strip('()'))][0],
                    model.get_genomic_data().genJ[int(row[j_gene_choice_col].strip('()'))][0]
                ), axis=1))
        return real_df

    def run(self, args, output_dir):
        """Function to execute the commandline tool.

        Parameters
        ----------
        args : Namespace
            Object containing our parsed commandline arguments.
        output_dir : str
            A directory path for writing output files to.

        """
        # Setup the necessary column names.
        col_names = {
            'I_COL': get_config_data('I_COL'),
            'NT_COL': get_config_data('NT_COL'),
            'NT_P_COL': get_config_data('NT_P_COL'),
            'AA_COL': get_config_data('AA_COL'),
            'AA_P_COL': get_config_data('AA_P_COL'),
            'V_GENE_CHOICE_COL': get_config_data('V_GENE_CHOICE_COL'),
            'D_GENE_CHOICE_COL': get_config_data('D_GENE_CHOICE_COL'),
            'J_GENE_CHOICE_COL': get_config_data('J_GENE_CHOICE_COL'),
        }

        # If the given type of sequences generation is not CDR3, use IGoR.
        if not args.cdr3:

            # Add general igor commands.
            command_list = []
            working_dir = get_config_data('WORKING_DIR')
            command_list.append(['set_wd', working_dir])
            command_list.append(['threads', str(get_config_data('NUM_THREADS'))])

            # Add the model (build-in or custom) command.
            sys.stdout.write('Processing IGoR model files...')
            try:
                if args.model:
                    files = get_default_model_file_paths(name=args.model)
                    command_list.append([
                        'set_custom_model',
                        files['parameters'],
                        files['marginals']
                    ])
                elif args.custom_model:
                    command_list.append([
                        'set_custom_model',
                        copy_to_dir(working_dir, str(args.custom_model[0]), 'txt'),
                        copy_to_dir(working_dir, str(args.custom_model[1]), 'txt')
                    ])
                sys.stdout.write(make_colored('success\n', 'green'))
            except IOError as err:
                sys.stdout.write(make_colored('error\n', 'red'))
                sys.stderr.write(make_colored(str(err) + '\n', 'bg-red'))
                return

            # Add generate command.
            command_list.append(['generate', str(args.n_gen), ['noerr']])

            # Execute IGoR through command line and catch error code.
            sys.stdout.write('Executing IGoR...')
            try:
                igor_cline = IgorInterface(command=command_list)
                exit_code, _, stderr, _ = igor_cline.call()
                if exit_code != 0:
                    sys.stdout.write(make_colored('error\n', 'red'))
                    sys.stderr.write(make_colored(
                        "An error occurred during execution of IGoR command (exit " \
                        "code {}):\n{}\n".format(exit_code, stderr), 'bg-red'))
                    return
                sys.stdout.write(make_colored('success\n', 'green'))
            except OSError as err:
                sys.stdout.write(make_colored('error\n', 'red'))
                sys.stderr.write(make_colored(str(err), 'bg-red'))
                return

            # Merge the generated output files together (translated).
            sys.stdout.write('Processing sequence realizations...')
            try:
                sequence_df = read_separated_to_dataframe(
                    file=os.path.join(working_dir, 'generated', 'generated_seqs_noerr.csv'),
                    separator=';',
                    index_col='seq_index',
                    cols=['nt_sequence'])
                sequence_df.index.names = [col_names['I_COL']]
                sequence_df.columns = [col_names['NT_COL']]
                sequence_df[col_names['AA_COL']] = sequence_df[col_names['NT_COL']] \
                    .apply(nucleotides_to_aminoacids)
                realizations_df = read_separated_to_dataframe(
                    file=os.path.join(working_dir, 'generated', 'generated_realizations_noerr.csv'),
                    separator=';',
                    index_col='seq_index')
                realizations_df.index.names = [col_names['I_COL']]
                if args.model:
                    files = get_default_model_file_paths(name=args.model)
                    model_type = files['type']
                    model = IgorLoader(model_type=model_type,
                                       model_params=files['parameters'],
                                       model_marginals=files['marginals'])
                elif args.custom_model:
                    model_type = args.type
                    model = IgorLoader(model_type=model_type,
                                       model_params=args.custom_model[0],
                                       model_marginals=args.custom_model[1])
                realizations_df = self._process_realizations(
                    data=realizations_df, model=model,
                    v_gene_choice_col=col_names['V_GENE_CHOICE_COL'],
                    d_gene_choice_col=col_names['D_GENE_CHOICE_COL'],
                    j_gene_choice_col=col_names['J_GENE_CHOICE_COL'])
                full_seqs_df = sequence_df.merge(realizations_df, left_index=True, right_index=True)
                sys.stdout.write(make_colored('success\n', 'green'))
            except (IOError, KeyError, ValueError) as err:
                sys.stdout.write(make_colored('error\n', 'red'))
                sys.stderr.write(make_colored(str(err) + '\n', 'bg-red'))
                return

            # Write the pandas dataframe to a separated file.
            sys.stdout.write('Writting file...')
            try:
                output_filename = get_config_data('OUT_NAME')
                if not output_filename:
                    output_filename = 'generated_seqs_{}'.format(model_type)
                _, filename = write_dataframe_to_separated(
                    dataframe=full_seqs_df,
                    filename=output_filename,
                    directory=output_dir,
                    separator=get_config_data('SEPARATOR'),
                    index_name=col_names['I_COL'])
                sys.stdout.write("(written '{}')...".format(filename))
                sys.stdout.write(make_colored('success\n', 'green'))
            except IOError as err:
                sys.stdout.write(make_colored('error\n', 'red'))
                sys.stderr.write(make_colored(str(err) + '\n', 'bg-red'))
                return

        # If the given type of sequences generation is CDR3, use OLGA.
        elif args.cdr3:

            # Get the working directory.
            working_dir = get_config_data('WORKING_DIR')

            # Load the model, create the sequence generator and generate the sequences.
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
            except (TypeError, OSError, IOError, KeyError, ValueError) as err:
                sys.stdout.write(make_colored('error\n', 'red'))
                sys.stderr.write(make_colored(str(err) + '\n', 'bg-red'))
                return

            # Setup the sequence generator and generate sequences.
            sys.stdout.write('Generating sequences...')
            try:
                seq_generator = OlgaContainer(
                    igor_model=model,
                    nt_col=col_names['NT_COL'],
                    nt_p_col=col_names['NT_P_COL'],
                    aa_col=col_names['AA_COL'],
                    aa_p_col=col_names['AA_P_COL'],
                    v_gene_choice_col=col_names['V_GENE_CHOICE_COL'],
                    j_gene_choice_col=col_names['J_GENE_CHOICE_COL'])
                cdr3_seqs_df = seq_generator.generate(num_seqs=args.n_gen)
                sys.stdout.write(make_colored('success\n', 'green'))
            except (TypeError, IOError) as err:
                sys.stdout.write(make_colored('error\n', 'red'))
                sys.stderr.write(make_colored(str(err) + '\n', 'bg-red'))
                return

            # Write the pandas dataframe to a separated file with.
            sys.stdout.write('Writting file...')
            try:
                output_filename = get_config_data('OUT_NAME')
                if not output_filename:
                    output_filename = 'generated_seqs_{}_CDR3'.format(model_type)
                _, filename = write_dataframe_to_separated(
                    dataframe=cdr3_seqs_df,
                    filename=output_filename,
                    directory=output_dir,
                    separator=get_config_data('SEPARATOR'),
                    index_name=col_names['I_COL'])
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
