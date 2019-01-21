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


"""Test file for testing immuno_probs.model.igor_interface file."""


import pytest

from immuno_probs.model.igor_interface import IgorInterface


@pytest.mark.parametrize('options, expected', [
    (['cmd', '1', '2'], 'cmd 1 2'),
    (['cmd', ['sub', '1', '2'], '3'], 'cmd -sub 1 2 3'),
    (['cmd', ['sub', ['sub-sub']], ['sub', '2', ['sub-sub', '3']]],
     'cmd -sub --sub-sub -sub 2 --sub-sub 3'),
])
def test_igor_interface(options, expected):
    """Test if the IgorInterface class creates the subprocess command properly.

    Parameters
    ----------
    options : list
        A Python nested ordered list with each value being a command/option.
        White spaces will seperate the options. Note: the depth of the
        nested lists will determine the number '-' characters to add in
        front of the first element for each list.
    expected : str
        The correctly formatted commandline string.

    Raises
    -------
    AssertionError
        If the performed test failed.

    """
    igor_cline = IgorInterface(options)
    command = igor_cline.get_command()
    assert command == expected
