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

import pandas

from immuno_probs.cdr3.olga_container import OlgaContainer
from immuno_probs.data.default_models import get_default_model_file_paths
from immuno_probs.model.igor_interface import IgorInterface
from immuno_probs.model.igor_loader import IgorLoader
from immuno_probs.util.cli import dynamic_cli_options
from immuno_probs.util.constant import get_num_threads, get_working_dir, get_separator
from immuno_probs.util.io import read_csv_to_dataframe, write_dataframe_to_csv, preprocess_input_file


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
        description = "Generate VJ or VDJ sequences given a custom IGoR " \
            "model (or build-in) by executing IGoR commandline tool via " \
            "python subprocess. Or generate CDR3 sequences by using the OLGA."
        parser_options = {
            '-model': {
                'type': 'str',
                'choices': ['human-t-alpha', 'human-t-beta', 'human-b-heavy',
                            'mouse-t-beta'],
                'required': '-custom-model' not in sys.argv,
                'help': "Specify a pre-installed model for generation. " \
                        "(required if --custom-model not specified) " \
                        "(select one: %(choices)s)."
            },
            '-type': {
                'type': 'str',
                'choices': ['CDR3', 'VDJ'],
                'required': 'True',
                'help': 'The type of sequences to generate. (select one: ' \
                        '%(choices)s)'
            },
            '-anchors': {
                'metavar': ('<v_gene>', '<j_gene>'),
                'type': 'str',
                'nargs': 2,
                'required': ('-type=CDR3' in sys.argv or
                             ('-type' in sys.argv and 'CDR3' in sys.argv)
                             and '-custom-model' in sys.argv),
                'help': 'The V and J gene CDR3 anchor files. Note: need to ' \
                        'contain gene in the firts column, anchor index in ' \
                        'the second and gene function in the third (required ' \
                        'for -type=CDR3 and -custom_model).'
            },
            '-custom-model': {
                'metavar': ('<parameters>', '<marginals>'),
                'type': 'str',
                'nargs': 2,
                'help': 'A IGoR parameters txt file followed by an IGoR ' \
                        'marginals txt file.'
            },
            '-generate': {
                'type': 'int',
                'nargs': '?',
                'default': 1,
                'help': 'The number of sequences to generate. (default: ' \
                        '%(default)s)'
            },
        }

        # Add the options to the parser and return the updated parser.
        parser_tool = self.subparsers.add_parser(
            'generate-seqs', help=description, description=description)
        parser_tool = dynamic_cli_options(parser=parser_tool,
                                          options=parser_options)

    @staticmethod
    def _process_realizations(data, model):
        """Function for processing an IGoR realization dataframe with indices.

        Parameters
        ----------
        data : pandas.DataFrame
            A pandas dataframe object with the IGoR realization data.
        model : IgorLoader
            Object containing the IGoR model.

        Returns
        -------
        pandas.DataFrame
            A pandas dataframe object with 'seq_index', 'gene_choice_v',
            'gene_choice_j' and optionally 'gene_choice_d' columns containing
            the names of the selected genes.

        """
        # If the suplied model is VDJ, locate important columns and update index values.
        if model.is_vdj():
            real_df = pandas.concat([data[['seq_index']],
                                     data.filter(regex=("GeneChoice_V_gene_.*")),
                                     data.filter(regex=("GeneChoice_J_gene_.*")),
                                     data.filter(regex=("GeneChoice_D_gene_.*"))],
                                    ignore_index=True, axis=1, sort=False)
            real_df.columns = ['seq_index', 'gene_choice_v', 'gene_choice_j', 'gene_choice_d']
            v_gene_names = [V[0] for V in model.get_genomic_data().genV]
            j_gene_names = [J[0] for J in model.get_genomic_data().genJ]
            d_gene_names = [J[0] for J in model.get_genomic_data().genD]
            for i, row in real_df.iterrows():
                real_df.ix[i, 'gene_choice_v'] = v_gene_names[int(row['gene_choice_v'].strip('()'))]
                real_df.ix[i, 'gene_choice_j'] = j_gene_names[int(row['gene_choice_j'].strip('()'))]
                real_df.ix[i, 'gene_choice_d'] = d_gene_names[int(row['gene_choice_d'].strip('()'))]

        # Or do the same if the model is VJ.
        elif model.is_vj():
            real_df = pandas.concat([data[['seq_index']],
                                     data.filter(regex=("GeneChoice_V_gene_.*")),
                                     data.filter(regex=("GeneChoice_J_gene_.*"))],
                                    ignore_index=True, axis=1, sort=False)
            real_df.columns = ['seq_index', 'gene_choice_v', 'gene_choice_j']
            v_gene_names = [V[0] for V in model.get_genomic_data().genV]
            j_gene_names = [J[0] for J in model.get_genomic_data().genJ]
            for i, row in real_df.iterrows():
                real_df.ix[i, 'gene_choice_v'] = v_gene_names[int(row['gene_choice_v'].strip('()'))]
                real_df.ix[i, 'gene_choice_j'] = j_gene_names[int(row['gene_choice_j'].strip('()'))]
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
        # If the given type of sequences generation is VDJ, use IGoR.
        if args.type == 'VDJ':

            # Add general igor commands.
            command_list = []
            working_dir = get_working_dir()
            command_list.append(['set_wd', working_dir])
            command_list.append(['threads', str(get_num_threads())])

            # Add the model (build-in or custom) command.
            if args.model:
                files = get_default_model_file_paths(model_name=args.model)
                command_list.append(['set_custom_model', files['parameters'],
                                     files['marginals']])
            elif args.custom_model:
                command_list.append(['set_custom_model', str(args.custom_model[0]),
                                     str(args.custom_model[1])])

            # Add generate command.
            command_list.append(['generate', str(args.generate), ['noerr']])

            # Execute IGoR through command line and catch error code.
            igor_cline = IgorInterface(args=command_list)
            code, _ = igor_cline.call()
            if code != 0:
                print("An error occurred during execution of IGoR " \
                      "command (exit code {})".format(code))
                return

            # Merge the generated output files together (translated).
            sequence_df = read_csv_to_dataframe(
                file=os.path.join(working_dir, 'generated', 'generated_seqs_noerr.csv'),
                separator=';')
            realizations_df = read_csv_to_dataframe(
                file=os.path.join(working_dir, 'generated', 'generated_realizations_noerr.csv'),
                separator=';')
            if args.model:
                files = get_default_model_file_paths(model_name=args.model)
                model = IgorLoader(model_params=files['parameters'],
                                   model_marginals=files['marginals'])
            elif args.custom_model:
                model = IgorLoader(model_params=args.custom_model[0],
                                   model_marginals=args.custom_model[1])
            realizations_df = self._process_realizations(data=realizations_df,
                                                         model=model)
            vdj_seqs_df = sequence_df.merge(realizations_df, on='seq_index')

            # Write the pandas dataframe to a CSV file.
            directory, filename = write_dataframe_to_csv(
                dataframe=vdj_seqs_df,
                filename='generated_VDJ_seqs',
                directory=output_dir,
                separator=get_separator())
            print("Written '{}' file to '{}' directory.".format(
                filename, directory))

        # If the given type of sequences generation is CDR3, use OLGA.
        elif args.type == 'CDR3':

            # Get the working directory.
            working_dir = get_working_dir()

            # Load the model, create the sequence generator and generate the sequences.
            if args.model:
                files = get_default_model_file_paths(model_name=args.model)
                model = IgorLoader(model_params=files['parameters'],
                                   model_marginals=files['marginals'])
                model.load_anchors(model_params=files['parameters'],
                                   v_anchors=files['v_anchors'],
                                   j_anchors=files['j_anchors'])
            elif args.custom_model:
                model = IgorLoader(model_params=args.custom_model[0],
                                   model_marginals=args.custom_model[1])
                v_anchors = preprocess_input_file(
                    os.path.join(working_dir, 'cdr3_anchors'), str(args.anchors[0]),
                    get_separator(), ',')
                j_anchors = preprocess_input_file(
                    os.path.join(working_dir, 'cdr3_anchors'), str(args.anchors[1]),
                    get_separator(), ',')
                model.load_anchors(model_params=args.custom_model[0],
                                   v_anchors=v_anchors,
                                   j_anchors=j_anchors)
            seq_generator = OlgaContainer(igor_model=model)
            cdr3_seqs_df = seq_generator.generate(num_seqs=args.generate)

            # Write the pandas dataframe to a CSV file.
            directory, filename = write_dataframe_to_csv(
                dataframe=cdr3_seqs_df,
                filename='generated_CDR3_seqs',
                directory=output_dir,
                separator=get_separator())
            print("Written '{}' file to '{}' directory.".format(
                filename, directory))


def main():
    """Function to be called when file executed via terminal."""
    print(__doc__)


if __name__ == "__main__":
    main()
