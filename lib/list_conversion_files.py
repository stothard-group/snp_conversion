#!/usr/bin/env python3
import os

# This module finds the files in the variant file directory and outputs the species and assembly options


def list_conversion_files(var_dir):
    """
    Reads variant file dir and outputs species name and assembly options
    :param var_dir: variant file directory
    :return: dict containing species and assembly values
    """
    if os.path.exists(var_dir):
        pass
    else:
        exit("Directory " + var_dir + " does not appear to exist. Please check path")
    species_names = os.listdir(var_dir)
    conversion_dict = {}
    print("Available species and assembly names:")
    for species in species_names:
        species_dirs = os.path.join(var_dir, species)
        assembly_dirs = os.listdir(species_dirs)
        conversion_dict.update({species: assembly_dirs})
    for key in conversion_dict:
        print(""" - {}""".format(key))
        for values in conversion_dict[key]:
            print("""    - {}""".format(values))
    return conversion_dict

