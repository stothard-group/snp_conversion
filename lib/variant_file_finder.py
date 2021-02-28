#!/usr/bin/env python3
import pandas as pd
import os
from lib.file_parsing import uncompressing, get_uncompressed_file
import warnings
import bz2
import gzip
import zipfile
import io
from functools import reduce
import numpy as np


def parse_vars(file_handle):
    skip_lines = 0
    for line in file_handle:
        line = line.rstrip()
        if line.startswith("#"):
            skip_lines = skip_lines + 1
        else:
            break
    return skip_lines


def get_var_df(var_dir, var_file_list, assembly, species, reg_alt_bool_dict):
    """
    Gets the var file that we are using and returns a df with this info
    :param var_dir: variant file dir
    :param var_file_list: list of variant files to search (this is truly now a list)
    :param assembly: assembly name
    :param species: species
    :param reg_alt_bool_dict: dict with "reg" or "alt" for each var df name
    :return: dataframe containing variant file info
    """
    # Read in variant file
    fully_merged_list = []
    for var_file in var_file_list:
        if reg_alt_bool_dict[var_file] == "reg":
            alt_bool = False
        else:
            alt_bool = True
        var_species_path = os.path.join(var_dir, species)
        var_assembly_path = os.path.join(var_species_path, assembly)
        filepath = os.path.join(var_assembly_path, var_file)
        header_count = uncompressing(filepath)
        whole_var_df = pd.read_csv(
            filepath, header=0, skiprows=header_count, compression="infer"
        )
        # Read in position file
        position_name = var_file.replace("conversion", "position")
        position_filepath = os.path.join(var_assembly_path, position_name)
        pos_header_count = uncompressing(position_filepath)
        position_df = pd.read_csv(
            position_filepath, header=0, skiprows=pos_header_count, compression="infer"
        )

        # Make output dataframe look like the old dataframe
        # Use A and B dataframes instead of ref and alt
        df_A = whole_var_df[whole_var_df["AB"] == "A"].copy()
        df_B = whole_var_df[whole_var_df["AB"] == "B"].copy()
        merged_df_A_B = pd.merge(df_A, df_B, on=["marker_name", "alt_marker_name"])
        if alt_bool is True:
            merged_df_A_B.drop(columns=["marker_name"], inplace=True)
            merged_df_A_B.rename(
                columns={"alt_marker_name": "marker_name"}, inplace=True
            )
        else:
            pass
        merged_df_A_B.rename(
            columns={
                "marker_name": "Name",
                "DESIGN_x": "DESIGN_A",
                "DESIGN_y": "DESIGN_B",
                "TOP_x": "TOP_A",
                "TOP_y": "TOP_B",
                "PLUS_x": "PLUS_A",
                "PLUS_y": "PLUS_B",
                "FORWARD_x": "FORWARD_A",
                "FORWARD_y": "FORWARD_B",
                "VCF_x": "VCF_A",
                "VCF_y": "VCF_B",
                "AB_x": "AB_A",
                "AB_y": "AB_B"
            },
            inplace=True,
        )
        if alt_bool is False:
            merged_df_A_B_dropped = merged_df_A_B.drop(
                columns=["alt_marker_name"]
            )
        else:
            merged_df_A_B_dropped = merged_df_A_B.copy()
        merged_df_A_B_dropped["SNP"] = ""
        merged_df_A_B_dropped["BLAST_strand"] = ""
        merged_df_A_B_dropped["Reference_allele_forward_strand"] = ""
        # Add positional info
        if alt_bool is True:
            position_marker_info = position_df[
                ["alt_marker_name", "chromosome", "position"]
            ].copy()
            position_marker_info.rename(
                columns={"alt_marker_name": "marker_name"}, inplace=True
            )
        else:
            position_marker_info = position_df[
                ["marker_name", "chromosome", "position"]
            ].copy()
        position_marker_info.rename(
            columns={
                "marker_name": "Name",
                "chromosome": "BLAST_chromosome",
                "position": "BLAST_position",
            },
            inplace=True,
        )
        merged_positional = pd.merge(
            merged_df_A_B_dropped, position_marker_info, on="Name"
        )
        merged_positional = merged_positional[
            [
                "Name",
                "SNP",
                "BLAST_chromosome",
                "BLAST_position",
                "BLAST_strand",
                "Reference_allele_forward_strand",
                "DESIGN_A",
                "DESIGN_B",
                "FORWARD_A",
                "FORWARD_B",
                "PLUS_A",
                "PLUS_B",
                "TOP_A",
                "TOP_B",
                "VCF_A",
                "VCF_B",
                "AB_A",
                "AB_B"

            ]
        ]
        merged_positional.fillna(".", inplace=True)
        #print(merged_positional)
        fully_merged_list.append(merged_positional)
    # Will have possibly multiple merged positional files here but all contradictory SNPs should be gone
    ##### DO THINGS HERE TO GET BOTH FILES MERGED
    final_merge = reduce(
        lambda left, right: pd.merge(
            left,
            right,
            on=[
                "Name",
                "SNP",
                "BLAST_chromosome",
                "BLAST_position",
                "BLAST_strand",
                "Reference_allele_forward_strand",
                "DESIGN_A",
                "DESIGN_B",
                "FORWARD_A",
                "FORWARD_B",
                "PLUS_A",
                "PLUS_B",
                "TOP_A",
                "TOP_B",
                "VCF_A",
                "VCF_B",
                "AB_A",
                "AB_B"
            ],
            how="outer",
        ),
        fully_merged_list,
    )
    return final_merge


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
            header_dict.update({"REF": splitline[1]})
        elif each_line.startswith("#PANEL="):
            splitline = each_line.split("=")
            header_dict.update({"PANEL": splitline[1]})
        else:
            pass
    return header_dict


def complete_subset_check(var_dict, input_names, bad_var_file, complete_var_dict):
    """
    Checks whether the user input is found in the entire var df file dataset. Will exit if this is not the case
    Also convert alt names to reg
    :param var_dict: dict containing 2D lists of markers and alt_markers from all variant files
    :param input_names: SNP panel input names
    :param bad_var_file: Name of file to write probes to that are not found in any of the var files
    :param complete_var_dict: dict containing all variant files
    :return: String indicating whether we have all snps in reg, alt, or a both files ("reg", "alt", "both" are options)
    :return: if "both", nested list containing the reg vs alt snps (otherwise an empty list)
    """
    var_probe_reg = []
    var_probe_alt = []
    reg_vs_alt = []
    list_of_dicts = []

    for var_key in var_dict:
        var_probe_reg = var_probe_reg + var_dict[var_key][0]
        var_probe_alt = var_probe_alt + var_dict[var_key][1]
        list_of_dicts.append(complete_var_dict[var_key])
    concat_var_files = pd.concat(list_of_dicts)
    var_probe_both = var_probe_reg + var_probe_alt
    if set(input_names).issubset(set(var_probe_both)):
        if set(input_names).issubset(set(var_probe_reg)):
            reg_alt_both_snps = "reg"
            converted_input_names = input_names
        elif set(input_names).issubset(set(var_probe_alt)):
            reg_alt_both_snps = "alt"
            sub_df = concat_var_files[concat_var_files["alt_marker_name"].isin(input_names)]
            converted_input_names = sub_df["marker_name"].to_list()
        else:
            reg_alt_both_snps = "both"
            reg_snps = [i for i in input_names if i in var_probe_reg]
            converted_input_names = reg_snps.copy()
            alt_snps = [i for i in input_names if i in var_probe_alt]
            sub_df = concat_var_files[concat_var_files["alt_marker_name"].isin(input_names)]
            converted_input_names.extend(sub_df["marker_name"].to_list())
            reg_vs_alt = [reg_snps, alt_snps]


    else:
        reg_al_both_snps = ''
        converted_input_names = []
        missing_snps = list(set(input_names) - set(var_probe_both))
        number_missing = len(missing_snps)
        with open(bad_var_file, "w") as var_outfile:
            var_outfile.write("\n".join(missing_snps))
        var_outfile.close()
        message = (
            str(number_missing)
            + " user-specified SNP probe(s) are not found in any conversion file. Problem variants have been written to the file "
            + bad_var_file
            + "\n"
        )
        exit(message)
    return reg_alt_both_snps, reg_vs_alt, converted_input_names


def check_for_incongruous_matches(matches_list, complete_var_dict, converted_snps):
    """
    Checks for whether there are matches to user probes that are not congruent in multiple selected variant files
    :param matches_list: list of files that user SNPs are found in (selected variant files contributing to the user panel)
    :param complete_var_dict: dict containing all the info for all variant files
    :param converted_snps: user input snps which have been converted to reg format if in alt
    :param reg_vs_alt: 2D list of the reg allele and the alt allele lists
    :return: A concatenated dataframe of all variant files subsetted on user input SNPs that contribute to the user SNP input file
    """
    # Get list of SNPs from the user panel in each var dict
    subset_var_dict = {}
    subset_snp_list = []
    for var_keys in matches_list:
        working_var_df = complete_var_dict[var_keys]
        subset_snps = working_var_df.loc[
            working_var_df["marker_name"].isin(converted_snps)
        ]
        subset_snp_list.append(subset_snps)
        subset_var_dict.update({var_keys: subset_snps})
    merged_var_df = reduce(
        lambda left, right: pd.merge(
            left,
            right,
            on=["marker_name", "VCF", "AB", "TOP", "FORWARD", "DESIGN", "PLUS"],
            how="outer",
        ),
        subset_snp_list,
    )
    concatenated_var_df = pd.concat(subset_snp_list)
    ref_only_df = merged_var_df[merged_var_df["VCF"] == "REF"]
    marker_names_list = ref_only_df["marker_name"].to_list()
    unique, unique_index, unique_counts = np.unique(
        marker_names_list, return_inverse=True, return_counts=True
    )
    count_mask = unique_counts > 1
    marker_duplicates = unique[count_mask]
    if len(marker_duplicates) > 0:
        marker_duplicates_list = "\n".join(marker_duplicates)
        message = (
            "One or more user input SNPs is found in multiple contributing conversion files with conflicting information. "
            "Remove these SNPs from your input file(s):\n" + marker_duplicates_list
        )
        exit(message)
    else:
        pass
    return concatenated_var_df


def best_match_calculation(reg_or_alt, var_df_dict, snp_names, complete_var_dict):
    """
    Performs best match var file calculations. Reg and alt snps are separate
    :param reg_or_alt: (Str) either "reg" or "alt"
    :param var_df_dict: dict containing variant dataframes
    :param snp_names: SNP probes (either reg or alt)
    :return:
    """
    var_missing_count_dict = {}
    var_missing_match = {}
    df_name_to_return = []
    multi_matching_file_dict = {}
    for var in var_df_dict:
        if reg_or_alt == "reg":
            var_names = var_df_dict[var][0]
        else:
            var_names = var_df_dict[var][1]
        missing_list = [x for x in snp_names if x not in var_names]
        missing_count = len(
            missing_list
        )  # if all snps are in the var file, missing count will be 0, otherwise, >0
        match_list = [x for x in snp_names if x in var_names]
        var_missing_count_dict.update({var: missing_count})
        var_missing_match.update({var: [match_list, missing_list]})
        # Check if all normal panel markers are found in the var file
        snp_set = set(snp_names)
        var_set = set(var_names)
        if snp_set.issubset(var_set):
            multi_matching_file_dict.update({var: missing_count})
    # Now we have dicts that have a missing count, and a missing list and match list for all var files
    # We should be able to figure out the best and second, third, etc. best var files for a user snp panel
    best_match = min(
        var_missing_count_dict, key=var_missing_count_dict.get
    )  # Name of the var file that gives the best match
    best_match_missing_val = var_missing_count_dict[best_match]
    df_name_to_return = []
    if best_match_missing_val == 0:
        # EVERYTHING IS AWESOME
        df_name_to_return.append(best_match)
    else:  # The best match does not contain all variants in the user input file
        # :(
        # Subtract the missing values that are matched in best match and iteratively and go until the missing values are 0
        df_name_to_return.append(best_match)
        starting_snp_names = snp_names
        missing_values_found = False
        current_best_match = best_match
        working_var_missing_match = var_missing_match
        while not missing_values_found:
            working_snp_names = list(
                set(starting_snp_names)
                - set(working_var_missing_match[current_best_match][0])
            )
            del working_var_missing_match[current_best_match]
            current_var_missing_count_dict = {}
            new_var_missing_match = {}
            for vars in working_var_missing_match:  # create new dict of matches
                new_missing_list = list(
                    set(working_snp_names) - set(working_var_missing_match[vars][0])
                )
                new_missing_count = len(new_missing_list)
                new_match_list = list(
                    set(snp_names) & set(working_var_missing_match[vars][0])
                )
                current_var_missing_count_dict.update({vars: new_missing_count})
                new_var_missing_match.update({vars: [new_match_list, new_missing_list]})
            current_best_match = min(
                current_var_missing_count_dict, key=current_var_missing_count_dict.get
            )
            current_best_match_missing_val = current_var_missing_count_dict[
                current_best_match
            ]
            df_name_to_return.append(current_best_match)
            if current_best_match_missing_val == 0:
                missing_values_found = True
            else:
                starting_snp_names = working_snp_names
    return df_name_to_return, multi_matching_file_dict


def get_best_match_list(
    reg_alt_both_snps, reg_vs_alt, varfile_dict_reg_alt, input_names, complete_var_dict
):
    """
    Gets a list of the best matching SNP files from user input SNPs, returns a dict of multimatching files and a dict
    for which var dfs are associated with reg or alt values
    :param reg_alt_both_snps: str "reg", "alt", "both"; whether user input snps have regular, alt, or both reg and alt names
    :param reg_vs_alt: if "both", a 2D list containing the SNPs in reg and SNPs in alt; otherwise an empty list
    :param varfile_dict_reg_alt: dict containing reg and alt varfile data
    :param input_names: list of user input snp panel names
    :return:
    """
    user_reg_snps = []
    user_alt_snps = []
    reg_alt_bool_dict = {}
    multi_matching_file_dict_list = []
    if reg_alt_both_snps == "reg":
        df_name_list, multi_matching_file_dict = best_match_calculation(
            "reg", varfile_dict_reg_alt, input_names, complete_var_dict
        )
        multi_matching_file_dict_list.append(multi_matching_file_dict)
        for val in df_name_list:
            reg_alt_bool_dict.update({val: "reg"})
    elif reg_alt_both_snps == "alt":
        df_name_list, multi_matching_file_dict = best_match_calculation(
            "alt", varfile_dict_reg_alt, input_names, complete_var_dict
        )
        multi_matching_file_dict_list.append(multi_matching_file_dict)
        for val in df_name_list:
            reg_alt_bool_dict.update({val: "alt"})
    else:
        df_name_list_reg, multi_matching_file_dict1 = best_match_calculation(
            "reg", varfile_dict_reg_alt, reg_vs_alt[0], complete_var_dict
        )
        multi_matching_file_dict_list.append(multi_matching_file_dict1)
        for val in df_name_list_reg:
            reg_alt_bool_dict.update({val: "reg"})
        df_name_list_alt, multi_matching_file_dict2 = best_match_calculation(
            "alt", varfile_dict_reg_alt, reg_vs_alt[1], complete_var_dict
        )
        multi_matching_file_dict_list.append(multi_matching_file_dict2)
        for val in df_name_list_alt:
            reg_alt_bool_dict.update({val: "alt"})
        df_name_list = df_name_list_reg + df_name_list_alt
    return df_name_list, multi_matching_file_dict_list, reg_alt_bool_dict


def var_match(
    var_files, var_dir, file_df, input_file, converting_file, assembly, species
):
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
        if split_name[2] == "conversion":  # get only the conversion files here
            if split_name[1] == assembly:
                subset_var_files.append(var_f1)
            else:
                pass
    if not subset_var_files:
        exit("No variant files contain the assembly name " + assembly)
    complete_var_dict = {}
    # Look only in subset files
    varfile_dict_reg_alt = {}
    for var_f2 in subset_var_files:
        var_species_dir = os.path.join(var_dir, species)
        var_assembly_dir = os.path.join(var_species_dir, assembly)
        var_f2_path = os.path.join(var_assembly_dir, var_f2)
        # get the header line count
        header_count = uncompressing(var_f2_path)
        # uncompress the file
        filetype = get_uncompressed_file(var_f2_path)
        if filetype == "no match" or filetype == "native_uncompressed":
            filetype = None
        else:
            pass
        # read file into dataframe
        test_var_df = pd.read_csv(
            var_f2_path, header=0, skiprows=header_count, compression=filetype
        )
        # read important header information for each file (SPECIES=, PANEL=, REF=)
        var_path_uncompressed = get_uncompressed_file(var_f2_path)
        if (
            var_path_uncompressed == "no match"
            or var_path_uncompressed == "native_uncompressed"
        ):
            # we want to actually open the file here
            with open(var_f2_path, "r") as uncompressed_handle:
                var_file_header_dict = get_vf_header(uncompressed_handle)
        elif var_path_uncompressed == "gzip":
            with gzip.open(var_f2_path, "rt") as gzipped_handle:
                var_file_header_dict = get_vf_header(gzipped_handle)
        elif var_path_uncompressed == "bz2":
            with bz2.open(var_f2_path, "rt") as bzipped_handle:
                var_file_header_dict = get_vf_header(bzipped_handle)
        elif var_path_uncompressed == "zip":
            with zipfile.ZipFile(var_f2_path, mode="r", allowZip64=True) as zipped_file:
                zip_list = zipfile.ZipFile.namelist(zipped_file)
                if len(zip_list) != 1:
                    warnings.warn(
                        "Zipped archive "
                        + var_f2_path
                        + " contains more than one file."
                    )
                else:
                    with zipped_file.open(zip_list[0], mode="r") as zipped_handle:
                        items_file = io.TextIOWrapper(
                            zipped_handle, encoding="UTF-8", newline=""
                        )
                        var_file_header_dict = get_vf_header(items_file)
        else:
            pass
        header_dict_per_varfile.update({var_f2: var_file_header_dict})
        # Get marker names (2 lines per marker)
        marker_names = test_var_df["marker_name"].to_list()
        marker_names_reduced = list(dict.fromkeys(marker_names))
        # Get alt marker names
        alt_marker_names = test_var_df["alt_marker_name"].to_list()
        alt_marker_names_reduced = list(dict.fromkeys(alt_marker_names))
        complete_var_dict.update({var_f2: test_var_df})
        varfile_dict_reg_alt.update(
            {var_f2: [marker_names_reduced, alt_marker_names_reduced]}
        )
        # Finished creating varfile dict

    # Get SNP panel marker names
    file_columns = list(file_df.columns.values)
    if "Name" in file_columns:
        input_names = file_df["Name"].to_list()
    elif "SNP Name" in file_columns:
        input_names = file_df["SNP Name"].to_list()
    else:
        input_names = []

    # First check if all user input markers are found in the panel files (as a whole); return whether we have "reg", "alt", or "both" markers
    bad_var_file = input_basename[0] + "_problem_variants.txt"
    reg_alt_both_snps, reg_vs_alt, converted_snps = complete_subset_check(
        varfile_dict_reg_alt, input_names, bad_var_file, complete_var_dict
    )

    # Go through list to get the best matches
    df_name_list, multi_matching_file_dict_list, reg_alt_bool_dict = get_best_match_list(
        reg_alt_both_snps,
        reg_vs_alt,
        varfile_dict_reg_alt,
        input_names,
        complete_var_dict,
    )
    if multi_matching_file_dict_list[0]:
        for mmdict in multi_matching_file_dict_list:
            keys = list(mmdict.keys())
            if len(keys) == 1:
                matches_list.append(keys[0])
            else:
                match = min(mmdict, key=mmdict.get)
                matches_list.append(match)
                smallest_match_dict = header_dict_per_varfile[match]
                smallest_match_panel = smallest_match_dict["PANEL"]
                if converting_file is not True:
                    message = (
                        "Multiple variant files match "
                        + input_file
                        + ". The smallest panel was chosen: "
                        + smallest_match_panel
                    )
                    verbose_log.append(message)
                else:
                    pass
    else:
        matches_list = df_name_list
    # Check for incongruous matches
    if len(matches_list) > 1:
        output_dataframe = check_for_incongruous_matches(
            matches_list, complete_var_dict, converted_snps
        )

    return matches_list, verbose_log, reg_alt_bool_dict
