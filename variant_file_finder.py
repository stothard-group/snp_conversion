#!/usr/bin/python
import pandas as pd
import os
import file_parsing as fp
import warnings
import bz2
import gzip
import zipfile
import io


def parse_vars(file_handle):
    skip_lines = 0
    for line in file_handle:
        line = line.rstrip()
        if line.startswith("#"):
            skip_lines = skip_lines + 1
        else:
            break
    return skip_lines


def get_var_df(var_dir, var_file, assembly, species, alt_bool):
    """
    Gets the var file that we are using and returns a df with this info
    :param var_dir: variant file dir
    :param var_file: list of variant files to search
    :param assembly: assembly name
    :param species: species
    :param alt_bool: (bool) whether we are using the regular or alt marker names
    :return: dataframe containing variant file info
    """
    # Read in variant file
    var_species_path = os.path.join(var_dir, species)
    var_assembly_path = os.path.join(var_species_path, assembly)
    filepath = os.path.join(var_assembly_path, var_file)
    header_count = fp.uncompressing(filepath)
    whole_var_df = pd.read_csv(filepath, header=0, skiprows=header_count, compression='infer')
    # Read in position file
    position_name = var_file.replace("conversion", "position")
    position_filepath = os.path.join(var_assembly_path, position_name)
    pos_header_count = fp.uncompressing(position_filepath)
    position_df = pd.read_csv(position_filepath, header=0, skiprows=pos_header_count, compression='infer')
    # Make output dataframe look like the old dataframe
    # Use A and B dataframes instead of ref and alt
    df_A = whole_var_df[whole_var_df["AB"] == "A"].copy()
    df_B = whole_var_df[whole_var_df["AB"] == "B"].copy()
    merged_df_A_B = pd.merge(df_A, df_B, on=["marker_name", "alt_marker_name"])
    if alt_bool is True:
        merged_df_A_B.drop(columns=["marker_name"], inplace=True)
        merged_df_A_B.rename(columns={"alt_marker_name": "marker_name"}, inplace=True)
    else:
        pass
    merged_df_A_B.rename(columns={"marker_name": "Name", "DESIGN_x": "DESIGN_A", "DESIGN_y": "DESIGN_B",
                                   "TOP_x": "TOP_A", "TOP_y": "TOP_B", "PLUS_x": "PLUS_A", "PLUS_y": "PLUS_B",
                                   "FORWARD_x": "FORWARD_A", "FORWARD_y": "FORWARD_B"}, inplace=True)
    if alt_bool is True:
        merged_df_A_B_dropped = merged_df_A_B.drop(columns=["AB_x", "AB_y", "VCF_x", "VCF_y"])
    else:
        merged_df_A_B_dropped = merged_df_A_B.drop(columns=["alt_marker_name", "AB_x", "AB_y", "VCF_x", "VCF_y"])
    merged_df_A_B_dropped["SNP"] = ''
    merged_df_A_B_dropped["BLAST_strand"] = ''
    merged_df_A_B_dropped["Reference_allele_forward_strand"] = ''
    # Add positional info
    if alt_bool is True:
        position_marker_info = position_df[["alt_marker_name", "chromosome", "position"]].copy()
        position_marker_info.rename(columns={"alt_marker_name": "marker_name"}, inplace=True)
    else:
        position_marker_info = position_df[["marker_name", "chromosome", "position"]].copy()
    position_marker_info.rename(columns={"marker_name": "Name", "chromosome": "BLAST_chromosome", "position": "BLAST_position"}, inplace=True)
    merged_positional = pd.merge(merged_df_A_B_dropped, position_marker_info, on="Name")
    merged_positional = merged_positional[["Name", "SNP", "BLAST_chromosome", "BLAST_position", "BLAST_strand",
                                           "Reference_allele_forward_strand", "DESIGN_A", "DESIGN_B", "FORWARD_A",
                                           "FORWARD_B", "PLUS_A", "PLUS_B", "TOP_A", "TOP_B"]]
    merged_positional.fillna(".", inplace=True)
    return merged_positional


def get_vf_header(var_f2_handle):
    """
    Extracts header information from var file. The first 3 lines of each header should be
        #SPECIES=
        #REF=
        #PANEL=
    Although not in this particular order.
    :param var_f2_handle: open var file handle
    :return: dict containing header information for the variant file
    """
    head = [next(var_f2_handle) for x in range(10)]
    header_dict = {}
    for each_line in head:
        if each_line.startswith("#SPECIES="):
            splitline = each_line.split("=")
            header_dict.update({"SPECIES": splitline[1]})
        elif each_line.startswith("#REF="):
            splitline = each_line.split("=")
            header_dict.update({'REF': splitline[1]})
        elif each_line.startswith("#PANEL="):
            splitline = each_line.split("=")
            header_dict.update({'PANEL': splitline[1]})
        else:
            pass
    return header_dict


def var_match(var_files, var_dir, file_df, input_file, converting_file, assembly, species):
    """
    Finds the best matching variant file for the snp panel dataset
    :param var_files: list of variant files in the variant file directory
    :param var_dir: variant file directory name
    :param file_df: dataframe from the snp panel file
    :param input_file: input file name
    :param converting_file: (bool) whether we are converting a file
    :param assembly: assembly name
    :param species: species name
    :return: list containing the single best match (list format is historic)
    """
    matching_files = {}
    matches_list = []
    verbose_log = []
    reg_otherwise_best_matches = {}
    alt_otherwise_best_matches = {}
    header_dict_per_varfile = {}
    missing_sample_dict = {}
    alt_missing_sample_dict = {}
    input_basename = os.path.splitext(input_file)
    # Subset variant files based on assembly name
    subset_var_files = []
    for var_f1 in var_files:
        split_name = var_f1.split(".")
        if split_name[2] == "conversion": # get only the conversion files here
            if split_name[1] == assembly:
                subset_var_files.append(var_f1)
            else:
                pass
    if not subset_var_files:
        exit("No variant files contain the assembly name " + assembly)
    alt_used = False
    # Look only in subset files
    for var_f2 in subset_var_files:
        var_species_dir = os.path.join(var_dir, species)
        var_assembly_dir = os.path.join(var_species_dir, assembly)
        var_f2_path = os.path.join(var_assembly_dir, var_f2)
        # get the header line count
        header_count = fp.uncompressing(var_f2_path)
        # uncompress the file
        filetype = fp.get_uncompressed_file(var_f2_path)
        if filetype == "no match" or filetype == "native_uncompressed":
            filetype = None
        else:
            pass
        # read file into dataframe
        test_var_df = pd.read_csv(var_f2_path, header=0, skiprows=header_count, compression=filetype)
        # read important header information for each file (SPECIES=, PANEL=, REF=)
        var_path_uncompressed = fp.get_uncompressed_file(var_f2_path)
        if var_path_uncompressed == "no match" or var_path_uncompressed == "native_uncompressed":
            # we want to actually open the file here
            with open(var_f2_path, "r") as uncompressed_handle:
                var_file_header_dict = get_vf_header(uncompressed_handle)
        elif var_path_uncompressed == "gzip":
            with gzip.open(var_f2_path, "rt") as gzipped_handle:
                var_file_header_dict = get_vf_header(gzipped_handle)
        elif var_path_uncompressed == 'bz2':
            with bz2.open(var_f2_path, "rt") as bzipped_handle:
                var_file_header_dict =  get_vf_header(bzipped_handle)
        elif var_path_uncompressed == "zip":
            with zipfile.ZipFile(var_f2_path, mode='r', allowZip64=True) as zipped_file:
                zip_list = zipfile.ZipFile.namelist(zipped_file)
                if len(zip_list) != 1:
                    warnings.warn("Zipped archive " + var_f2_path + " contains more than one file.")
                else:
                    with zipped_file.open(zip_list[0], mode='r') as zipped_handle:
                        items_file = io.TextIOWrapper(zipped_handle, encoding='UTF-8', newline='')
                        var_file_header_dict = get_vf_header(items_file)
        else:
            pass
        header_dict_per_varfile.update({var_f2: var_file_header_dict})
        # Get marker names (2 lines per marker)
        marker_names = test_var_df['marker_name'].to_list()
        marker_names_reduced = list(dict.fromkeys(marker_names))
        var_len = len(marker_names_reduced)
        # Get alt marker names
        alt_marker_names = test_var_df['alt_marker_name'].to_list()
        alt_marker_names_reduced = list(dict.fromkeys(alt_marker_names))
        alt_var_len = len(alt_marker_names_reduced)

        # Get SNP panel marker names
        if 'Name' in file_df:
            input_names = file_df['Name'].to_list()
        elif 'SNP Name' in file_df:
            input_names = file_df['SNP Name'].to_list()
        else:
            input_names = []
        # Check if all normal panel markers are found in the var file
        if set(input_names).issubset(set(marker_names_reduced)):
            matching_files.update({var_f2: var_len})
        # If not all normal panel markers are found in var file, find the number of markers that are missing from it
        else:
            missing_samples = set(input_names) - set(marker_names_reduced)
            missing_sample_dict.update({var_f2: missing_samples})
            reg_otherwise_best_matches.update({var_f2: len(missing_samples)})
        # Do the same thing for alt markers
        if set(input_names).issubset(set(alt_marker_names_reduced)):
            matching_files.update({var_f2: alt_var_len})
            alt_used = True
        else:
            missing_alt = set(input_names) - set(alt_marker_names_reduced)
            alt_missing_sample_dict.update({var_f2: missing_alt})
            alt_otherwise_best_matches.update({var_f2: len(missing_alt)})
    # Deal with non-matching files, and print all problem markers to a new file
    if alt_used is False:
        otherwise_best_matches = reg_otherwise_best_matches
    else:
        otherwise_best_matches = alt_otherwise_best_matches
    if not matching_files:
        best_match = min(otherwise_best_matches, key=otherwise_best_matches.get)
        best_match_dict = header_dict_per_varfile[best_match]
        best_match_panel = best_match_dict["PANEL"]
        bad_var_file = input_basename[0] + "_problem_variants.txt"
        # Write problem variants to output file
        with open(bad_var_file, "w") as var_outfile:
            var_missing = missing_sample_dict[best_match]
            var_outfile.write("\n".join([missing for missing in list(var_missing)]))
        var_outfile.close()
        message = "No variant file contains all user-specified variants. The closest panel is " + best_match_panel + ", missing " + str(otherwise_best_matches[best_match]) + \
                  " variant(s). Problem variants have been written to the file " + bad_var_file + "\n"
        exit(message)
    # return the name of the matching variant file
    else:
        keys = list(matching_files.keys())
        if len(keys) == 1:
            matches_list.append(keys[0])
        else:
            match = min(matching_files, key=matching_files.get)
            matches_list.append(match)
            smallest_match_dict = header_dict_per_varfile[match]
            smallest_match_panel = smallest_match_dict["PANEL"]
            if converting_file is not True:
                message = "Multiple variant files match " + input_file + ". The smallest panel was chosen: " + \
                          smallest_match_panel
                print(message)
                verbose_log.append(message)
            else:
                pass
    return matches_list, verbose_log, alt_used
