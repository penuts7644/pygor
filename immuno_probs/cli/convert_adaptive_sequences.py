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


"""Commandline tool for converting VDJ sequences (and CDR3) from adaptive format."""


import logging
import os

import pandas
import numpy

from immuno_probs.convert.adaptive_sequence_convertor import AdaptiveSequenceConvertor
from immuno_probs.util.cli import dynamic_cli_options
from immuno_probs.util.constant import get_config_data
from immuno_probs.util.io import copy_to_dir, preprocess_reference_file, \
write_dataframe_to_separated, read_fasta_as_dataframe, read_separated_to_dataframe


class ConvertAdaptiveSequences(object):
    """Commandline tool for converting full and CDR3 sequences from a given
    adaptive input sequence file.

    Parameters
    ----------
    subparsers : argparse.ArgumentParser
        A subparser object for appending the tool's parser and options.

    Methods
    -------
    run(args)
        Uses the given Namespace commandline arguments to convert the full
        length (VDJ for productive, unproductive and the total) and CDR3
        sequences from a given adaptive input sequence file.

    """
    def __init__(self, subparsers):
        super(ConvertAdaptiveSequences, self).__init__()
        self.logger = logging.getLogger(__name__)
        self.subparsers = subparsers
        self._add_options()

    def _add_options(self):
        """Function for adding the parser/options to the input ArgumentParser.

        Notes
        -----
            Uses the class constructor's subparser object for appending the
            tool's parser and options.

        """
        # Create the description and options for the parser.
        description = "Converts the full length (VDJ for productive, " \
            "unproductive and the total) and CDR3 sequences from a given " \
            "adaptive input sequence file. The VDJ sequences can be used to " \
            "build a new IGoR model and the CDR3 sequences can be evaluated."
        parser_options = {
            '-seqs': {
                'metavar': '<separated>',
                'required': 'True',
                'type': 'str',
                'help': "An input separated data file with sequences to " \
                        "convert using the defined column names."
            },
            '-ref': {
                'metavar': ('<gene>', '<fasta>'),
                'type': 'str',
                'action': 'append',
                'nargs': 2,
                'required': 'True',
                'help': "A gene (V or J) followed by a reference genome " \
                        "FASTA file. Note: the FASTA reference genome files " \
                        "needs to conform to IGMT annotation (separated by " \
                        "'|' character)."
            },
            '-n-random': {
                'type': 'int',
                'nargs': '?',
                'help': "Number of random sequences (subset) to convert " \
                        "from the given file (default: {}).".format(
                            get_config_data('CONVERT', 'NUM_RANDOM', 'int'))
            },
            '-use-allele': {
                'action': 'store_true',
                'help': "If specified (True), the allele information from " \
                        "the resolved gene fields are used to when " \
                        "reconstructing the gene choices (default: {}).".format(
                            get_config_data('CONVERT', 'USE_ALLELE', 'bool'))
            },
        }

        # Add the options to the parser and return the updated parser.
        parser_tool = self.subparsers.add_parser(
            'convert', help=description, description=description)
        parser_tool = dynamic_cli_options(parser=parser_tool,
                                          options=parser_options)

    @staticmethod
    def _process_gene_df(filename, nt_col, resolved_col):
        """Private function for processing the given gene FASTA file.

        Parameters
        ----------
        filename : str
            File path FASTA of the gene to process.
        nt_col : str
            The column name containing the sequences.
        resolved_col : str
            The column name containing the resolved values.

        """
        # Read in the seperated file as dataframe.
        gene_df = read_fasta_as_dataframe(
            file=filename, col=nt_col, header='info')

        # Modify the dataframe to have gene resolved column, remove the old
        # info column and return.
        gene_df[resolved_col] = pandas.DataFrame(
            list(gene_df['info'].apply(lambda x: x.split('|')[1])))
        gene_df.drop('info', axis=1, inplace=True)
        return gene_df

    def run(self, args, output_dir):
        """Function to execute the commandline tool.

        Parameters
        ----------
        args : Namespace
            Object containing our parsed commandline arguments.
        output_dir : str
            A directory path for writing output files to.

        """
        # Get the working directory.
        working_dir = get_config_data('COMMON', 'WORKING_DIR')

        # Collect and read in the corresponding reference genomic templates.
        self.logger.info('Processing genomic reference templates')
        try:
            for gene in args.ref:
                filename = preprocess_reference_file(
                    os.path.join(working_dir, 'genomic_templates'),
                    copy_to_dir(working_dir, gene[1], 'fasta'),
                )
                if gene[0] == 'V':
                    v_gene_df = self._process_gene_df(
                        filename=filename,
                        nt_col=get_config_data('COMMON', 'NT_COL'),
                        resolved_col=get_config_data('COMMON', 'V_RESOLVED_COL'))
                if gene[0] == 'J':
                    j_gene_df = self._process_gene_df(
                        filename=filename,
                        nt_col=get_config_data('COMMON', 'NT_COL'),
                        resolved_col=get_config_data('COMMON', 'J_RESOLVED_COL'))
        except (IOError, KeyError, ValueError) as err:
            self.logger.error(str(err))
            return

        # Read in the sequence data.
        self.logger.info('Pre-processing input sequence file')
        try:
            seqs_df = read_separated_to_dataframe(
                file=args.seqs,
                separator=get_config_data('COMMON', 'SEPARATOR'),
                cols=[get_config_data('COMMON', 'NT_COL'),
                      get_config_data('COMMON', 'AA_COL'),
                      get_config_data('COMMON', 'FRAME_TYPE_COL'),
                      get_config_data('COMMON', 'CDR3_LENGTH_COL'),
                      get_config_data('COMMON', 'V_RESOLVED_COL'),
                      get_config_data('COMMON', 'J_RESOLVED_COL')])

            # Take a random subsample of sequences in the file.
            n_random = get_config_data('CONVERT', 'NUM_RANDOM', 'int')
            if args.n_random:
                n_random = args.n_random
            if n_random != 0:
                if len(seqs_df) < n_random:
                    self.logger.warn(
                        'Number of random sequences is higher then number of ' \
                        'rows in file, all rows are used')
                    return
        except (IOError, KeyError, ValueError) as err:
            self.logger.error(str(err))
            return

        # Setup the data convertor class and convert data.
        self.logger.info('Converting adaptive file format')
        try:
            use_allele = get_config_data('CONVERT', 'USE_ALLELE', 'bool')
            if args.use_allele:
                use_allele = args.use_allele
            asc = AdaptiveSequenceConvertor()
            cdr3_df, full_prod_df, full_unprod_df, full_df = asc.convert(
                num_threads=get_config_data('COMMON', 'NUM_THREADS', 'int'),
                seqs=seqs_df,
                ref_v_genes=v_gene_df,
                ref_j_genes=j_gene_df,
                row_id_col=get_config_data('COMMON', 'ROW_ID_COL'),
                nt_col=get_config_data('COMMON', 'NT_COL'),
                aa_col=get_config_data('COMMON', 'AA_COL'),
                frame_type_col=get_config_data('COMMON', 'FRAME_TYPE_COL'),
                cdr3_length_col=get_config_data('COMMON', 'CDR3_LENGTH_COL'),
                v_resolved_col=get_config_data('COMMON', 'V_RESOLVED_COL'),
                v_gene_choice_col=get_config_data('COMMON', 'V_GENE_CHOICE_COL'),
                j_resolved_col=get_config_data('COMMON', 'J_RESOLVED_COL'),
                j_gene_choice_col=get_config_data('COMMON', 'J_GENE_CHOICE_COL'),
                use_allele=use_allele,
                default_allele=get_config_data('CONVERT', 'DEFAULT_ALLELE'),
                n_random=n_random)
            cdr3_df.insert(0, get_config_data('COMMON', 'FILE_NAME_ID_COL'),
                           os.path.splitext(os.path.basename(args.seqs))[0])
            full_prod_df.insert(0, get_config_data('COMMON', 'FILE_NAME_ID_COL'),
                                os.path.splitext(os.path.basename(args.seqs))[0])
            full_unprod_df.insert(0, get_config_data('COMMON', 'FILE_NAME_ID_COL'),
                                  os.path.splitext(os.path.basename(args.seqs))[0])
            full_df.insert(0, get_config_data('COMMON', 'FILE_NAME_ID_COL'),
                           os.path.splitext(os.path.basename(args.seqs))[0])
        except KeyError as err:
            self.logger.error(str(err))
            return

        # Copy the output files to the output directory with prefix.
        try:
            self.logger.info('Writing converted files to file system')
            output_prefix = get_config_data('COMMON', 'OUT_NAME')
            if not output_prefix:
                output_prefix = 'converted'
            _, filename_1 = write_dataframe_to_separated(
                dataframe=cdr3_df,
                filename='{}_CDR3'.format(output_prefix),
                directory=output_dir,
                separator=get_config_data('COMMON', 'SEPARATOR'),
                index_name=get_config_data('COMMON', 'I_COL'))
            self.logger.info("Written '%s'", filename_1)
            _, filename_2 = write_dataframe_to_separated(
                dataframe=full_prod_df,
                filename='{}_full_length_productive'.format(output_prefix),
                directory=output_dir,
                separator=get_config_data('COMMON', 'SEPARATOR'),
                index_name=get_config_data('COMMON', 'I_COL'))
            self.logger.info("Written '%s'", filename_2)
            _, filename_3 = write_dataframe_to_separated(
                dataframe=full_unprod_df,
                filename='{}_full_length_unproductive'.format(output_prefix),
                directory=output_dir,
                separator=get_config_data('COMMON', 'SEPARATOR'),
                index_name=get_config_data('COMMON', 'I_COL'))
            self.logger.info("Written '%s'", filename_3)
            _, filename_4 = write_dataframe_to_separated(
                dataframe=full_df,
                filename='{}_full_length'.format(output_prefix),
                directory=output_dir,
                separator=get_config_data('COMMON', 'SEPARATOR'),
                index_name=get_config_data('COMMON', 'I_COL'))
            self.logger.info("Written '%s'", filename_4)
        except IOError as err:
            self.logger.error(str(err))
            return


def main():
    """Function to be called when file executed via terminal."""
    print(__doc__)


if __name__ == "__main__":
    main()
