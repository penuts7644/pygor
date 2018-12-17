# ImmunoProbs Python package uses a simplified manner for calculating the
# generation probability of V(D)J and CDR3 sequences.
# Copyright (C) 2018 Wout van Helvoirt

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


from immuno_probs.util.cli import dynamic_cli_options
from immuno_probs.util.io import write_dataframe_to_csv


class CreateIgorModel(object):
    """Commandline tool for creating custom IGoR V(D)J models.

    Parameters
    ----------
    subparsers : ArgumentParser
        A subparser object for appending the tool's parser and options.

    Methods
    -------
    run(args)
        Uses the given Namespace commandline arguments to locate the run IGoR
        for creating a custom model. As of now, still requires an initial model.

    """
    def __init__(self, subparsers):
        super(CreateIgorModel, self).__init__()
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
        description = "This tool creates a V(D)J model by executing IGoR via a " \
            "python subprocess. The resulting files are written to the " \
            "'immuno_probs' directory."
        parser_options = {
            'seqs': {
                'metavar': 'S',
                'type': 'str',
                'help': 'An input FASTA file with sequences for training the model.'
            },
            'v-gene': {
                'metavar': 'V',
                'type': 'str',
                'help': 'Reference genome FASTA for the V gene.'
            },
            'j-gene': {
                'metavar': 'J',
                'type': 'str',
                'help': 'Reference genome FASTA for the J gene.'
            },
            'init-model': {
                'metavar': 'M',
                'type': 'str',
                'help': "An initial IGoR model's parameters file."
            },
            '--d-gene': {
                'type': 'str',
                'nargs': '?',
                'help': 'Optional reference genome FASTA for the D gene if ' \
                    'available.'
            },
            '--set-wd': {
                'type': 'str',
                'nargs': '?',
                'help': 'An optional location for creating the IGoR files. By ' \
                    'default, uses the current directory for written files.'
            },
            '--n-iter': {
                'type': 'int',
                'nargs': '?',
                'default': '1',
                'help': 'The number of inference iterations to perform when ' \
                    'creating the model. (default: %(default)s)'
            }
        }

        # Add the options to the parser and return the updated parser.
        parser_tool = self.subparsers.add_parser('CreateIgorModel',
                                                 help=description, description=description)
        parser_tool = dynamic_cli_options(parser=parser_tool, options=parser_options)


    @staticmethod
    def run(args):
        """Function to execute the commandline tool.

        Parameters
        ----------
        args : Namespace
            Object containing our parsed commandline arguments.

        """
        # TODO


def main():
    """Function to be called when file executed via terminal."""
    print(__doc__)


if __name__ == "__main__":
    main()
