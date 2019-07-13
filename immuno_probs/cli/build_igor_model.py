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


"""Commandline tool for creating a custom IGoR V(D)J model."""

import os
from shutil import copy2
import sys

from immuno_probs.model.default_models import get_default_model_file_paths
from immuno_probs.model.igor_interface import IgorInterface
from immuno_probs.util.cli import dynamic_cli_options, make_colored
from immuno_probs.util.constant import get_config_data
from immuno_probs.util.io import preprocess_separated_file, \
preprocess_reference_file, is_fasta, is_separated, copy_to_dir


class BuildIgorModel(object):
    """Commandline tool for creating custom IGoR V(D)J models.

    Parameters
    ----------
    subparsers : argparse.ArgumentParser
        A subparser object for appending the tool's parser and options.

    Methods
    -------
    run(args)
        Uses the given Namespace commandline arguments to execute IGoR
        for creating a custom model.

    """
    def __init__(self, subparsers):
        super(BuildIgorModel, self).__init__()
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
        description = "Create a VDJ or VJ model by executing IGoR's " \
            "commandline tool via a python subprocess using default model " \
            "parameters."
        parser_options = {
            '-seqs': {
                'metavar': '<fasta/separated>',
                'required': 'True',
                'type': 'str',
                'help': "An input FASTA or separated data file with sequences for " \
                        "training the model."
            },
            '-ref': {
                'metavar': ('<gene>', '<fasta>'),
                'type': 'str',
                'action': 'append',
                'nargs': 2,
                'required': 'True',
                'help': "A gene (V, D or J) followed by a reference genome " \
                        "FASTA file. Note: the FASTA reference genome files " \
                        "needs to conform to IGMT annotation (separated by " \
                        "'|' character)."
            },
            '-type': {
                'type': 'str.lower',
                'choices': ['alpha', 'beta', 'light', 'heavy'],
                'required': 'True',
                'help': 'The type of model to create. (select one: ' \
                        '%(choices)s).'
            },
            '-n-iter': {
                'type': 'int',
                'nargs': '?',
                'default': 1,
                'help': 'The number of inference iterations to perform when ' \
                        'creating the model. (default: %(default)s)'
            }
        }

        # Add the options to the parser and return the updated parser.
        parser_tool = self.subparsers.add_parser(
            'build', help=description, description=description)
        parser_tool = dynamic_cli_options(parser=parser_tool,
                                          options=parser_options)

    @staticmethod
    def _copy_file_to_output(file, filename, directory):
        """Copies a txt file to the given directory.

        If the file already exists, a number will be appended to the filename.
        The given output directory is created recursively if it does not exist.

        Parameters
        ----------
        file : str
            A string path to the file to copy over to the directory
        filename : str
            Base filename for writting the file, excluding the extension.
        directory : str
            A directory path location to write the file to.

        Returns
        -------
        tuple
            Containing the output directory and the name of the file that has
            been written to disk.

        """
        # Check if the filename is unique, modify name if necessary.
        file_count = 1
        updated_filename = filename

        # Keep modifying the filename until it doesn't exist.
        while os.path.isfile(os.path.join(directory, updated_filename + '.txt')):
            updated_filename = str(filename) + '_' + str(file_count)
            file_count += 1

        # Copy input file to new location and return info.
        copy2(file, os.path.join(directory, updated_filename + '.txt'))
        return (directory, updated_filename + '.txt')

    def run(self, args, output_dir):
        """Function to execute the commandline tool.

        Parameters
        ----------
        args : Namespace
            Object containing our parsed commandline arguments.
        output_dir : str
            A directory path for writing output files to.

        """
        # Add general igor commands.
        command_list = []
        working_dir = get_config_data('WORKING_DIR')
        command_list.append(['set_wd', working_dir])
        command_list.append(['threads', str(get_config_data('NUM_THREADS'))])

        # Add sequence and file paths commands.
        sys.stdout.write('Processing genomic reference templates...')
        try:
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

        # Set the initial model parameters using a build-in model.
        sys.stdout.write('Setting initial model parameters...')
        if args.type in ['beta', 'heavy']:
            command_list.append([
                'set_custom_model',
                get_default_model_file_paths(name='human-t-beta')['parameters']
            ])
        elif args.type in ['alpha', 'light']:
            command_list.append([
                'set_custom_model',
                get_default_model_file_paths(name='human-t-alpha')['parameters']
            ])
        sys.stdout.write(make_colored('success\n', 'green'))

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
                except (KeyError, ValueError) as err:
                    sys.stdout.write(make_colored('error\n', 'red'))
                    sys.stderr.write(make_colored(
                        "Given input sequence file does not have a '{}' column\n" \
                        .format(get_config_data('NT_COL')), 'bg-red'))
                    return
            else:
                sys.stdout.write(make_colored('error\n', 'red'))
                sys.stderr.write(make_colored(
                    'Given input sequence file could not be detected as FASTA ' \
                    'file or separated data type\n', 'bg-red'))
                return
            sys.stdout.write(make_colored('success\n', 'green'))
        except (IOError, KeyError) as err:
            sys.stdout.write(make_colored('error\n', 'red'))
            sys.stderr.write(make_colored(str(err) + '\n', 'bg-red'))
            return

        # Add alignment commands.
        command_list.append(['align', ['all']])

        # Add inference commands.
        command_list.append(['infer', ['N_iter', str(args.n_iter)]])

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

        # Copy the output files to the output directory with prefix.
        sys.stdout.write('Writting files...')
        try:
            output_prefix = get_config_data('OUT_NAME')
            if not output_prefix:
                output_prefix = 'model'
            _, filename_1 = self._copy_file_to_output(
                file=os.path.join(working_dir, 'inference', 'final_marginals.txt'),
                filename='{}_marginals'.format(output_prefix),
                directory=output_dir)
            _, filename_2 = self._copy_file_to_output(
                file=os.path.join(working_dir, 'inference', 'final_parms.txt'),
                filename='{}_params'.format(output_prefix),
                directory=output_dir)
            sys.stdout.write("(written '{}' and '{}')...".format(filename_1, filename_2))
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
