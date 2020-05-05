#!/usr/bin/env python3
import os
import warnings
from zipfile import ZipFile
# This file contains functions for common errors that can be checked (e.g. correct input file names, etc.)


def assembly_species_error(key_dir, assembly_name, species):
    """
    Checks that the assembly name and species combination are compatible
    :param key_dir: variant files directory
    :param assembly_name: name of the assembly
    :param species: species name
    :return: (bool) whether the combination exists
    """
    variant_path = os.path.join(key_dir, species)
    variant_assembly_dirs = os.listdir(variant_path)
    assembly_found = False
    for directory in variant_assembly_dirs:
        if assembly_name == directory:
            assembly_found = True

        else:
            pass
    return assembly_found


def input_file_name_error(input_dir, file_list):
    files_not_found = []
    for file in file_list:
        file_path = os.path.join(input_dir, file)
        if os.path.exists(file_path):
            pass
        else:
            files_not_found.append(file)
    return files_not_found


def non_txt_fmt_input_files(input_files_from_dir):
    warn_files = []
    for file in input_files_from_dir:
        if file.endswith(".txt"):
            pass
        else:
            warn_files.append(file)
    return


def non_csv_fmt_conversion_files(key_dir, assembly_name, species):
    """
    Checks that the variant files are csv and have acceptable compression
    :param key_dir: variant files directory
    :param assembly_name: assembly name
    :param species: species
    :return: List of files to exclude based on file extension
    """
    species_path = os.path.join(key_dir, species)
    variant_path = os.path.join(species_path, assembly_name)
    variant_files = os.listdir(variant_path)
    # remove files that do not contain conversion
    variant_exclude_list = []
    for file in variant_files:
        file_split = file.split(".")
        if "conversion" == file_split[2]:
            if assembly_name == file_split[1]:
                if file.endswith(('.csv', '.csv.gz', '.csv.zip', '.csv.bz2', '.csv.xz')):
                    # special zipfile check - should only contain 1 archive
                    if file.endswith('.zip'):
                        full_file_name = os.path.join(key_dir, file)
                        with ZipFile(full_file_name, mode='r', allowZip64=True) as zipped_file:
                            zip_list = ZipFile.namelist(zipped_file)
                            if len(zip_list) != 1:
                                warnings.warn("Zipped archive " + file +
                                              " contains more than one file - skipping", stacklevel=4)
                                variant_exclude_list.append(file)
                            else:
                                pass
                    else:
                        pass
                else:
                    # Check for any tar archives
                    if file.endswith(('.csv.tar.gz', '.csv.tar.bz2')):
                        warnings.warn("Only the following compression types are supported: gzip, bz2, and zip. "
                                      "Unpack tar archive as tar files are not supported.", stacklevel=4)
                    variant_exclude_list.append(file)
            else:
                pass
        else:
            pass
    return variant_exclude_list
