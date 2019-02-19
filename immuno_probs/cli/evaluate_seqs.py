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
from shutil import copy2

from Bio import SeqIO
import pandas

from immuno_probs.cdr3.olga_container import OlgaContainer
from immuno_probs.data.default_models import get_default_model_file_paths
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
        description = "Evaluate VJ or VDJ sequences given a custom IGoR " \
            "model (or build-in) through IGoR commandline tool via python " \
            "subprocess. Or evaluate CDR3 sequences through OLGA."
        parser_options = {
            '-seqs': {
                'metavar': '<csv>',
                'required': 'True',
                'type': 'str',
                'help': 'An input CSV file with sequences for evaluation. ' \
                        'Note: uses IGoR generated file formatting.'
            },
            '-model': {
                'type': 'str',
                'choices': ['human-t-alpha', 'human-t-beta', 'human-b-heavy',
                            'mouse-t-beta'],
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
                'required': ('-type=VDJ' in sys.argv or
                             ('-type' in sys.argv and 'VDJ' in sys.argv)
                             and '-custom-model' in sys.argv),
                'help': 'A gene (V, D or J) followed by a reference genome ' \
                        'FASTA file. Note: the FASTA reference genome files ' \
                        'needs to conform to IGMT annotation. (required ' \
                        'for -type=VDJ with -custom_model)'
            },
            '-type': {
                'type': 'str',
                'choices': ['CDR3', 'VDJ'],
                'required': 'True',
                'help': 'The type of sequences to generate. (select one: ' \
                        '%(choices)s)'
            },
            '-custom-model': {
                'metavar': ('<parameters>', '<marginals>'),
                'type': 'str',
                'nargs': 2,
                'help': 'A IGoR parameters txt file followed by an IGoR ' \
                        'marginals txt file.'
            },
            '-anchors': {
                'metavar': ('<v_gene>', '<j_gene>'),
                'type': 'str',
                'nargs': 2,
                'required': ('-type=CDR3' in sys.argv or
                             ('-type' in sys.argv and 'CDR3' in sys.argv)
                             and '-custom-model' in sys.argv),
                'help': 'The V and J gene CDR3 anchor files. (required ' \
                        'for -type=CDR3 with -custom_model)'
            },
        }

        # Add the options to the parser and return the updated parser.
        parser_tool = self.subparsers.add_parser(
            'evaluate-seqs', help=description, description=description)
        parser_tool = dynamic_cli_options(parser=parser_tool,
                                          options=parser_options)

    @staticmethod
    def _format_imgt_reference_fasta(directory, fasta):
        """Function for formatting the IMGT reference genome files for IGoR.

        Parameters
        ----------
        directory : str
            A directory path to write the files to.
        fasta : str
            A FASTA file path for a reference genomic template file.

        Returns
        -------
        str
            A string file path to the new reference genomic template file.

        """
        # Create the output directory.
        if not os.path.isdir(directory):
            os.makedirs(directory)

        # Open the fasta file, update the fasta header and write out.
        records = list(SeqIO.parse(str(fasta), "fasta"))
        for rec in records:
            rec.id = rec.description.split('|')[1]
            rec.description = ""
        updated_path = os.path.join(directory, os.path.basename(str(fasta)))
        SeqIO.write(records, updated_path, "fasta")
        return updated_path

    @staticmethod
    def _preprocess_input_seqs(directory, csv):
        """Function for formatting the input sequence file for IGoR.

        Parameters
        ----------
        directory : str
            A directory path to write the files to.
        csv : str
            A CSV formatted file path to precess for IGoR.

        Returns
        -------
        str
            A string file path to the sequence input file.

        """
        # Create the output directory.
        if not os.path.isdir(directory):
            os.makedirs(directory)

        # Open the sequence input file and update the columns.
        sequence_df = read_csv_to_dataframe(
            filename=csv,
            separator=get_separator())
        sequence_df = sequence_df.iloc[:, 0:1]

        # Write the new pandas dataframe to a CSV file.
        directory, filename = write_dataframe_to_csv(
            dataframe=sequence_df,
            filename=os.path.basename(str(csv)),
            directory=directory,
            separator=';')
        return os.path.join(directory, filename)

    def run(self, args, output_dir):
        """Function to execute the commandline tool.

        Parameters
        ----------
        args : Namespace
            Object containing our parsed commandline arguments.
        output_dir : str
            A directory path for writing output files to.

        """
        # If the given type of sequences evaluation is VDJ, use IGoR.
        if args.type == 'VDJ':

            # Add general igor commands.
            command_list = []
            working_dir = get_working_dir()
            command_list.append(['set_wd', working_dir])
            command_list.append(['threads', str(get_num_threads())])

            # Add the model (build-in or custom) command depending on given.
            if args.model:
                files = get_default_model_file_paths(model_name=args.model)
                command_list.append(['set_custom_model', files['parameters'],
                                     files['marginals']])
                ref_list = ['set_genomic']
                for gene, filename in files['reference'].items():
                    ref_list.append([gene, filename])
                command_list.append(ref_list)
            elif args.custom_model:
                command_list.append(['set_custom_model', str(args.custom_model[0]),
                                     str(args.custom_model[1])])
                ref_list = ['set_genomic']
                for i in args.ref:
                    filename = self._format_imgt_reference_fasta(
                        os.path.join(working_dir, 'genomic_templates'), i[1])
                    ref_list.append([i[0], filename])
                command_list.append(ref_list)

            # Add the sequence command after pre-processing of the input file.
            input_seqs = self._preprocess_input_seqs(
                os.path.join(working_dir, 'input'), str(args.seqs))
            command_list.append(['read_seqs', input_seqs])

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
                sys.exit()

            # Read in all data frame files.
            sequence_df = read_csv_to_dataframe(
                filename=args.seqs,
                separator=get_separator())
            vdj_pgen_df = read_csv_to_dataframe(
                filename=os.path.join(working_dir, 'output', 'Pgen_counts.csv'),
                separator=';')

            # Merge IGoR generated sequence output dataframes.
            vdj_pgen_df = sequence_df.merge(vdj_pgen_df, on='seq_index')

            # Write the pandas dataframe to a CSV file.
            directory, filename = write_dataframe_to_csv(
                dataframe=vdj_pgen_df,
                filename=os.path.join('output', 'VDJ_seqs_pgen_estimate'),
                directory=working_dir,
                separator=get_separator())

            # Write output file to output directory.
            copy2(os.path.join(directory, filename), output_dir)
            print("Written '{}' file to '{}' directory.".format(
                filename, output_dir))

        # If the given type of sequences evaluation is CDR3, use OLGA.
        elif args.type == 'CDR3':

            # Create the directory for the output files.
            working_dir = os.path.join(get_working_dir(), 'output')
            if not os.path.isdir(working_dir):
                os.makedirs(os.path.join(get_working_dir(), 'output'))

            # Load the model, create the sequence evaluator and evaluate the sequences.
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
                model.load_anchors(model_params=args.custom_model[0],
                                   v_anchors=args.anchors[0],
                                   j_anchors=args.anchors[1])
            seq_evaluator = OlgaContainer(igor_model=model)
            sequence_df = read_csv_to_dataframe(filename=args.seqs,
                                                separator=get_separator())
            cdr3_pgen_df = seq_evaluator.evaluate(seqs=sequence_df)

            # Merge IGoR generated sequence output dataframes.
            cdr3_pgen_df = sequence_df.merge(cdr3_pgen_df, on='seq_index')

            # Write the pandas dataframe to a CSV file.
            directory, filename = write_dataframe_to_csv(
                dataframe=cdr3_pgen_df,
                filename='CDR3_seqs_pgen_estimate',
                directory=working_dir,
                separator=get_separator())

            # Write output file to output directory.
            copy2(os.path.join(directory, filename), output_dir)
            print("Written '{}' file to '{}' directory.".format(
                filename, output_dir))


def main():
    """Function to be called when file executed via terminal."""
    print(__doc__)


if __name__ == "__main__":
    main()
