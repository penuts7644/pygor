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


"""IgorInterface class for interfacing with IGoR's commandline tool."""


import shlex
import subprocess

from immuno_probs.util.exception import SubprocessException


class IgorInterface(object):
    """Interface class for executing commandline processes for IGoR.

    Parameters
    ----------
    args : list
        A list with strings and nested lists that will be build into a
        subprocess command.

    Methods
    -------
    call()
        Call the IGoR program and return stdout and stderr messages.
    get_command()
        Returns the created commandline subprocess string.
    set_command(args)
        Set a new commandline string using a nested list.

    """
    def __init__(self, args):
        super(IgorInterface, self).__init__()
        self.command = self._subprocess_builder(options=args)

    def _subprocess_builder(self, options, level=0):
        """Creates a subprocess command string from ordered input list.

        Parameters
        ----------
        options : list
            A Python nested ordered list with each value being a command/option.
            White spaces will seperate the options. Note: the depth of the
            nested lists will determine the number '-' characters to add in
            front of the first element for each list.
        level : int, optional
            The initial start depth level for indication the number of '-'
            characters to add to the first item in the lists. (default: 0)

        Returns
        -------
        str
            The formatted commandline subprocess as string.

        """
        # Create the commandline string from the options in the list.
        command_str = ""
        for index, val in enumerate(options):
            if isinstance(val, list):
                val = self._subprocess_builder(val, level=level+1)
            if isinstance(val, str):
                if index == 0:
                    val = (level * '-') + val
                command_str += ' ' + val
        return command_str.strip(' ')

    def call(self):
        """Calls IGoR's commandline tool via subprocess.

        Returns
        -------
        tuple
            A tuple containing the exit code, standard out, standard error and
            the executed command.

        Raises
        ------
        SubprocessException
            When the subprocess program execution returns an error.

        Notes
        -----
            This function uses the generated commandline string.

        """
        # Execute the commandline process and return the results.
        updated_command = 'igor ' + self.command
        try:
            process = subprocess.Popen(shlex.split(updated_command))
            (stdout, stderr) = process.communicate()
            return (process.returncode, stdout, stderr, updated_command)
        except OSError as err:
            raise SubprocessException(err)

    def get_command(self):
        """Getter function for collecting the IGoR command.

        Returns
        -------
        str
            A command string for executing IGoR's commandline tool.

        """
        return self.command

    def set_command(self, args):
        """Setter function for setting a IGoR command.

        Parameters
        ----------
        args : list
            A list with strings and nested lists that will be build into a
            subprocess command.

        """
        self.command = self._subprocess_builder(options=args)


def main():
    """Function to be called when file executed via terminal."""
    print(__doc__)


if __name__ == "__main__":
    main()
