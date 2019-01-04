# ImmunoProbs Python package able to calculate the generation probability of
# V(D)J and CDR3 sequences. Copyright (C) 2018 Wout van Helvoirt

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


from immuno_probs.model.igor_interface import IgorInterface
from immuno_probs.util.cli import dynamic_cli_options
from immuno_probs.util.constant import get_num_threads, get_working_dir


class EvaluateVdjSeqs(object):
    """Commandline tool for evaluating V(D)J sequences using an IGoR model.

    Parameters
    ----------
    subparsers : ArgumentParser
        A subparser object for appending the tool's parser and options.

    Methods
    -------
    run(args)
        Uses the given Namespace commandline arguments to execute IGoR
        for evaluating V(D)J sequences.

    """
    def __init__(self, subparsers):
        super(EvaluateVdjSeqs, self).__init__()
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
            "commandline python subprocess."
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
            }
        }

        # Add the options to the parser and return the updated parser.
        parser_tool = self.subparsers.add_parser(
            'evaluate-vdj-seqs', help=description, description=description)
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
        # Add general igor commands.
        command_list = []
        if args.set_wd:
            command_list.append(['set_wd', str(args.set_wd)])
        else:
            command_list.append(['set_wd', str(get_working_dir())])
        if args.threads:
            command_list.append(['threads', str(args.threads)])
        else:
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
        code, stdout, stderr, _ = igor_cline.call()

        if code != 0:
            print("An error occurred during execution of IGoR command: \n")
            print("stderr:\n{}".format(stderr))
            print("stdout:\n{}".format(stdout))


def main():
    """Function to be called when file executed via terminal."""
    print(__doc__)


if __name__ == "__main__":
    main()
