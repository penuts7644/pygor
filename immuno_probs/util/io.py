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


"""Contains I/O related functions used in immuno_probs."""


import os
from shutil import copy2

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas


def create_directory_path(directory):
    """Updates and creates given directory path by adding a number at the end.

    Parameters
    ----------
    directory : string
        A directory path location to create recursively.

    Returns
    -------
    str
        The (updated) output directory location path.

    """
    # Check if the directory name is unique, modify name if necessary.
    dir_count = 1
    updated_directory = directory

    # Keep modifying the name until it doesn't exist.
    while os.path.isdir(os.path.join(directory, updated_directory)):
        updated_directory = str(directory) + '_' + str(dir_count)
        dir_count += 1

    # Finally create directory's recursively if not exists.
    if not os.path.isdir(updated_directory):
        os.makedirs(updated_directory)

    return updated_directory


def is_fasta(file):
    """Checks if the input file is valid fasta.

    Parameters
    ----------
    file : string
        Location of the FASTA file to be tested.

    """
    with open(file, "r") as fasta_file:
        return any(SeqIO.parse(fasta_file, "fasta"))


def is_csv(file, separator):
    """Checks if the input file is valid csv.

    Parameters
    ----------
    file : string
        Location of the CSV file to be tested.
    separator : string
        A separator character used for separating the fields in the file.

    """
    dataframe = pandas.read_csv(file, sep=separator, comment='#', header=0, nrows=100)
    return not dataframe.empty


def read_fasta_as_dataframe(file, index_col, seq_col):
    """Creates a pandas.DataFrame from the FASTA file.

    The dataframe contains header name and sequence columns containing the
    corresponding FASTA data.

    Parameters
    ----------
    file : string
        Location of the FASTA file to be read in.
    index_col : string
        The name of the index column.
    seq_col : string
        The name of the FASTA sequence column.

    """
    # Create a dataframe and read in the fasta file.
    fasta_df = pandas.DataFrame(columns=[index_col, seq_col])
    with open(file, 'r') as fasta_file:
        fasta_count = 0
        for _, sequence in SimpleFastaParser(fasta_file):
            fasta_df = fasta_df.append({
                index_col: fasta_count,
                seq_col: sequence.upper()
            }, ignore_index=True)
            fasta_count += 1
    return fasta_df


def read_csv_to_dataframe(file, separator):
    """Read in a CSV file as pandas.DataFrame.

    Parameters
    ----------
    file : string
        File path to be read in as dataframe.
    separator : string
        A separator character used for separating the fields in the file.

    Notes
    -----
        This function uses the global SEPARATOR variable to set the separator
        string for the input CSV file. Comments ('#') in the file are skipped.

    """
    dataframe = pandas.read_csv(file, sep=separator, comment='#', header=0)
    return dataframe


def write_dataframe_to_csv(dataframe, filename, directory, separator):
    """Writes a pandas.DataFrame to a CSV formatted file.

    If the file already exists, a number will be appended to the filename.
    The given output directory is created recursively if it does not exist.
    The column names in the dataframe is used as first line in the csv file.

    Parameters
    ----------
    dataframe : pandas.DataFrame
        The dataframe to be written to the CSV file.
    filename : string
        Base filename for writting the file, excluding the extension.
    directory : string
        A directory path location to create recursively.
    separator : string
        A separator character used for separating the fields in the CSV file.

    Returns
    -------
    tuple
        Containing the output directory and the name of the file that has been
        written to disk.

    """
    # Check if the filename is unique, modify name if necessary.
    file_count = 1
    updated_filename = filename

    # Keep modifying the filename until it doesn't exist.
    while os.path.isfile(os.path.join(directory, updated_filename + '.csv')):
        updated_filename = str(filename) + '_' + str(file_count)
        file_count += 1

    # Write dataframe contents to csv file and return info.
    pandas.DataFrame.to_csv(
        dataframe, path_or_buf=os.path.join(directory, updated_filename + '.csv'),
        sep=separator, index=False, na_rep='na')
    return (directory, updated_filename + '.csv')


def preprocess_csv_file(directory, file, in_sep, out_sep, index_col=None, cols=None):
    """Function for formatting the input sequence file for IGoR.

    Parameters
    ----------
    directory : str
        A directory path to write the file to.
    file : str
        A CSV formatted file path to precess for IGoR.
    in_sep : str
        The input file seperator.
    out_sep : str
        The wanted output file seperator.
    index_col : str, optional
        The name of the index column to use. If specified and given column is
        not found in the dataframe, the index values are generated. (default: No
        index column)
    cols : list, optional
        Containing column names to keep in the output file. The order will
        change the output file column formatting (default: includes all
        columns in the output file).

    Returns
    -------
    str
        A string file path to the newly created file.

    Notes
    -----
        Returns the input file path if no changes will be applied to the file.
        This means, the input seperator and output seperator are equal and the
        columns attribute has not been specified.

    Raises
    ------
    KeyError
        If a column is not found in the input data file.

    """
    # If the seperators are the same and no columns are given, return the input.
    if out_sep == in_sep and cols is None:
        return file

    # Create the output directory.
    if not os.path.isdir(directory):
        os.makedirs(directory)

    # Open the sequence input file.
    sequence_df = read_csv_to_dataframe(
        file=file,
        separator=in_sep)

    # Create new dataframe and add index if specified.
    updated_df = pandas.DataFrame()
    if index_col is not None:
        if index_col in sequence_df.columns:
            updated_df[index_col] = sequence_df.loc[:, index_col]
        else:
            updated_df[index_col] = range(0, sequence_df.size + 1)

    # Add the remainder columns if given.
    if cols is not None:
        for col in cols:
            updated_df[col] = sequence_df.loc[:, [col]]

    # Write the new pandas dataframe to a CSV file.
    directory, filename = write_dataframe_to_csv(
        dataframe=updated_df,
        filename=os.path.basename(str(file)),
        directory=directory,
        separator=out_sep)
    return os.path.join(directory, filename)


def preprocess_reference_file(directory, file, index):
    """Function for formatting the IMGT reference genome files for IGoR.

    Parameters
    ----------
    directory : str
        A directory path to write the file to.
    file : str
        A FASTA file path for a reference genomic template file.
    index : int
        Index of the header line to keep after splitting on '|'.

    Returns
    -------
    str
        A string file path to the new reference FASTA file.

    """
    # Create the output directory.
    if not os.path.isdir(directory):
        os.makedirs(directory)

    # Open the fasta file and update the fasta header.
    records = list(SeqIO.parse(file, "fasta"))
    for rec in records:
        rec.id = rec.description.split('|')[index]
        rec.description = ""

    # Write out the modified file.
    updated_path = os.path.join(directory, os.path.basename(file))
    SeqIO.write(records, updated_path, "fasta")
    return updated_path


def copy_to_dir(directory, file, extension):
    """Function for copying file to directory and modifying the extension.

    Parameters
    ----------
    directory : str
        A directory path to write the file to.
    file : str
        A FASTA file path for a reference genomic template file.
    extension : str
        An extension name, if the same, file will NOT be coppied.

    Returns
    -------
    str
        A string file path to the new file location and name.

    """
    # Check if file name extension if different, or return.
    filename, file_extension = os.path.splitext(os.path.basename(file))
    if file_extension == str('.' + extension):
        return file

    # Copy file to given directory if necessary.
    output_file = os.path.join(directory, filename + '.' + extension)
    copy2(file, output_file)
    return output_file
