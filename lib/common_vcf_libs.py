#!/usr/bin/env python3
import platform
import os
import pandas as pd
import lib.errors as errors
import warnings
from pandas.errors import ParserError
import time
from lib.file_format_checker import affy_test, TFDP_format_check
from lib.file_parsing import parse_header
import lib.file_conversion as fc
from lib.make_logs import simple_log
from lib.variant_file_finder import var_match, get_var_df
import atexit
import concurrent.futures as cf

# This is a set of common functions to vcf_generator and genotype_concordance
def check_user_platform():
    """
    Determines whether we can perform threading or if we are on a windows machine
    :return: (bool) True if we can thread, otherwise false
    """
    user_platform = platform.system()
    if user_platform == "Linux" or user_platform == "Darwin":
        threading = True
    elif user_platform == "Windows":
        threading = False
    else:
        threading = False
    return threading


def retrieve_user_input_file(panel_path):
    """
    Splits the path and snp file basename from the user input var snp-panel.
    Performs error checking for non text files and files not found.
    :param panel_path: path to snp file or filename
    :return: path to directory containing the snp panel
    :return: snp panel file name
    """
    panel_names = os.path.split(panel_path)
    path_to_file = panel_path
    file_name = panel_names[1]
    # try to open file
    try:
        with open(panel_path) as fi:
            pass
    except FileNotFoundError:
        print("Cannot find file " + panel_path)
    finally:
        fi.close()
        file_array = [panel_path]
    non_txt_file = errors.non_txt_fmt_input_files(file_array)
    if non_txt_file:
        exit("File " + panel_path + " is not a text file.")
    else:
        pass

    return path_to_file, file_name


def retrieve_all_variant_files(variant_dir, assembly, species):
    """
    Gets a list of all possible variant files (checking for an allowed assembly name) and excludes files in the variant
    directory that are not .csv format
    :param variant_dir: directory containing variant files (default: variant_position_files
    :param assembly: user-inputted name of assembly (e.g. UMD3.1)
    :return: list of all acceptable variant files
    """
    variant_species_dir = os.path.join(variant_dir, species)
    variant_assembly_dir = os.path.join(variant_species_dir, assembly)
    variant_files = os.listdir(variant_assembly_dir)
    variant_exclude = errors.non_csv_fmt_conversion_files(
        variant_dir, assembly, species
    )
    if variant_exclude:
        warnings.warn(
            "The following variant files are not .csv files are are excluded from analysis: "
            + ", ".join(variant_exclude),
            stacklevel=3,
        )
        for wrong_var in variant_exclude:
            variant_files.remove(wrong_var)  # remove non-csv variant files
    else:
        pass
    return variant_files


def parse_AFFY_file(affy_file_path, file_name, logfile):
    """
    Reads in known affymetrix file and returns affymetrix header dict and affymetrix data dict
    :param affy_file_path: file path to the user input SNP panel (in affy format)
    :param file_name: name of SNP panel for writing to log file
    :param logfile: log file
    :return: dict containing affy data (genotype, not AB data)
    """
    affy_df = pd.read_csv(affy_file_path, sep="\t", mangle_dupe_cols=True)
    # make dataframe look like the illumina one (remove AB columns)
    header_row = list(set(affy_df.columns))
    header_row.remove("probeset_id")
    affy_header_dict = {}
    for value in header_row:
        new_affy_df = affy_df[["probeset_id", value]].copy()
        ab_list = ["AA", "AB", "BB", "NoCall"]
        test_ab_df = ~new_affy_df.iloc[:, 1:].isin(ab_list)
        not_ab = sum(test_ab_df.sum(axis=1))
        affy_header_dict.update({value: not_ab})
    copy_affy_df = affy_df.copy()  # keep a copy just in case
    for keys in affy_header_dict:
        if affy_header_dict[keys] == 0:
            affy_df.drop(columns=keys, inplace=True)
    affy_df.rename(columns={"probeset_id": "Name"}, inplace=True)
    return affy_df


def long_file_to_array(illumina_df, logfile, file_name):
    """
    Parses a Long format file and creates a dict where each value is a dataframe for each 'file type'
    :param illumina_df: long file dataframe
    :param logfile: logfile name
    :param file_name: name of the snp panel file
    :return: dict containing all file type-specific dataframes in a long file
    """
    log_array = []
    type_list = ["PLUS"]
    df_dict = {}
    for ty in type_list:
        out1, out2 = fc.gen_long_output_col_names(ty)
        # Make sure allele type exists
        if out1 not in illumina_df.columns or out2 not in illumina_df.columns:
            message = (
                "Columns " + out1 + " and " + out2 + " might not be in the input file"
            )
            warnings.warn(message, stacklevel=4)
            log_array.append(message)
        else:
            sub_long_df = illumina_df[["SNP Name", "Sample ID", out1, out2]]
            # print new df by sample and save to dict by "file type"
            sample_list = sub_long_df["Sample ID"]
            unique_samples = list(set(sample_list))
            pos1 = sub_long_df.columns.values[2]
            pos2 = sub_long_df.columns.values[3]
            output_df = sub_long_df[["SNP Name"]].copy()
            output_df = output_df.rename(columns={"SNP Name": "Name"})
            sub_sample_df_list = []
            for sample in unique_samples:
                sub_sub_df = sub_long_df[sub_long_df["Sample ID"] == sample].copy()
                sub_sub_df.rename(columns={"SNP Name": "Name"}, inplace=True)
                sub_sample_df = sub_sub_df[["Name"]].copy()
                sub_sample_df[sample] = sub_sub_df[[pos1, pos2]].apply(
                    lambda x: "".join(x), axis=1
                )
                sub_sample_df_list.append(sub_sample_df)
            output_df_2 = reduce(
                lambda x, y: pd.merge(x, y, on="Name", how="left"), sub_sample_df_list
            )
            df_dict.update({ty: output_df_2})
    if logfile is not None:
        logging = simple_log(log_array, file_name, logfile)
    return df_dict


def update_snp_panel_names(panel_df):
    """
     Adds .1 values to this dataframe as the rest of the script depends on them
    :param panel_df: panel dataframe to update
    :return: updated dataframe
    """
    column_list = list(panel_df)
    column_list.remove("Name")
    update_col_dict = {}
    for colname in column_list:
        new_name = colname + ".1"
        update_col_dict.update({colname: new_name})
    panel_df.rename(columns=update_col_dict, inplace=True)
    return panel_df


def check_input_snp_panel(
    snp_panel,
    file_name,
    file_type,
    variant_files,
    assembly,
    variant_directory,
    can_we_thread,
    n_threads,
    logfile,
    species,
):
    """
    Function uses the file_format_check module to determine whether the SNP panel has a correct input format, or
    determines the input format if the format is "mixed" or unknown.
    :param snp_panel: path to the user input SNP panel
    :param file_name: name of user input SNP panel
    :param file_type: user-specified input file type
    :param variant_files: list of acceptable variant files
    :param assembly: user-inputted assembly to use
    :param variant_directory: dir containing variant files
    :param can_we_thread: (bool) test of whether we can thread
    :param n_threads: number of threads for threading
    :param logfile: log file name to pass to functions
    :param species: species name for variant file finding
    :return: determined filetype, if the file is correct (correct file), name of log file
    """
    log_array = []
    timestr = time.strftime("%H:%M:%S")
    logfile_text = timestr + " ..... Checking the format of the SNP panel"
    print("Checking SNP panel format")
    log_array.append(logfile_text)
    # declare some variables to make things easier later
    panel_dataframe = None
    long_file_df_dict = {}
    # test for affymetrix file
    if file_type == "AFFY" or file_type == "AFFY-PLUS":
        affy_flag = True
    else:
        affy_flag = False
    affy_flag = affy_test(snp_panel, file_name, file_type, affy_flag)
    if affy_flag is True:
        # Parse affymetrix file
        logfile_text = "Parsing AFFY file"
        log_array.append(logfile_text)
        affy_data_dict = parse_AFFY_file(snp_panel, file_name, logfile)
        panel_dataframe = affy_data_dict
    else:  # not an affy file
        # Read in header info of Illumina file
        header_row, header_dict = parse_header(snp_panel)
        # Read in data from Illumina file
        illumina_data_df = None
        try:
            illumina_data_df = pd.read_csv(snp_panel, skiprows=header_row, sep="\t")
        except ParserError:
            exit("First line might not contain column names - check formatting")
        if illumina_data_df is not None:
            if file_type != "LONG":
                logfile_text = "Parsing " + file_type + " file"
                log_array.append(logfile_text)
                illumina_data_df.rename(
                    columns={"Unnamed: 0": "Name"}, inplace=True
                )
                if file_type != "PLUS":
                    converting_file = False
                    var_list, log_text, reg_alt_bool_dict = var_match(
                        variant_files,
                        variant_directory,
                        illumina_data_df,
                        file_name,
                        converting_file,
                        assembly,
                        species,
                    )
                    for text in log_text:
                        log_array.append(text)
                    var_df = get_var_df(
                        variant_directory, var_list, assembly, species, reg_alt_bool_dict
                    )
                    matrix_type = file_type
                    output_type = "PLUS"
                    converted_df, column_names, logfile = fc.split_and_convert(
                        illumina_data_df,
                        var_df,
                        matrix_type,
                        output_type,
                        can_we_thread,
                        n_threads,
                        file_name,
                        logfile,
                    )
                    # Reorder converted dataframe to match the original order of samples and snp names
                    specified_out = "PLUS"
                    reordered_df, logfile = fc.reorder_converted_df(
                        converted_df,
                        column_names,
                        illumina_data_df,
                        specified_out,
                        output_type,
                        file_name,
                        logfile,
                    )
                    reordered_df.rename({"": "Name"}, axis=1, inplace=True)
                    panel_dataframe = reordered_df
                else:
                    illumina_data_df.rename({"": "Name"}, axis=1, inplace=True)
                    panel_dataframe = illumina_data_df

                # add .1 values to this dict as the rest of the script depends on them
            else:
                # make a dict containing the innards of the Long file and pass PLUS df to panel_dataframe, if it exists
                logfile_text = "Parsing Long file"
                log_array.append(logfile_text)
                long_file_df_dict = long_file_to_array(
                    illumina_data_df, logfile, file_name
                )
                if "PLUS" in long_file_df_dict:
                    panel_dataframe = long_file_df_dict["PLUS"]
                else:
                    panel_dataframe = illumina_data_df
    # Find corresponding var file
    converting_file = False
    if panel_dataframe is None:
        exit("something went wrong creating the snp panel dataframe")
    timestr = time.strftime("%H:%M:%S")
    logfile_text = timestr + " ..... Finding the matching variant file "
    log_array.append(logfile_text)
    var_list, log_text, reg_alt_bool_dict = var_match(
        variant_files,
        variant_directory,
        panel_dataframe,
        file_name,
        converting_file,
        assembly,
        species,
    )
    for text in log_text:
        log_array.append(text)
    var_df = get_var_df(variant_directory, var_list, assembly, species, reg_alt_bool_dict)
    #  Check if PLUS format exists, and if not, convert long fwd format to plus (or whatever format exists)
    if affy_flag is False:
        if file_type == "LONG":
            if "PLUS" not in long_file_df_dict:
                # make PLUS dataframe # SHOVE FILE INTO FILE CONVERSION AND CONVERT TO PLUS FORMAT
                # use long_as_matrix
                df_for_conversion, matrix_type, logfile = fc.long_as_matrix(
                    illumina_data_df, file_name, logfile
                )
                out_type_for_long = "PLUS"
                converted_df, column_names, logfile = fc.split_and_convert(
                    df_for_conversion,
                    var_df,
                    matrix_type,
                    out_type_for_long,
                    can_we_thread,
                    n_threads,
                    file_name,
                    logfile,
                )
                reordered_df, logfile = fc.reorder_converted_df(
                    converted_df,
                    column_names,
                    illumina_data_df,
                    out_type_for_long,
                    out_type_for_long,
                    file_name,
                    logfile,
                )
                reordered_df.rename({"": "Name"}, axis=1, inplace=True)
                panel_dataframe = reordered_df
        else:
            pass
    # Convert affymetrix file to AFFY-PLUS
    else:
        if file_type == "AFFY":
            matrix_type = "FWD"
            converted_df, column_names, logfile = fc.split_and_convert(
                panel_dataframe,
                var_df,
                matrix_type,
                "PLUS",
                can_we_thread,
                n_threads,
                file_name,
                logfile,
            )
            reordered_df, logfile = fc.reorder_converted_df(
                converted_df,
                column_names,
                panel_dataframe,
                "FWD",
                "PLUS",
                file_name,
                logfile,
            )
            reordered_df.rename(columns={"": "Name"}, inplace=True)
            panel_dataframe = reordered_df
        else:
            pass

    # Quick test user SNP panel
    #is_mixed = True
    #fmt = "PLUS"
    #logfile_text = "Checking the panel file format relative to the variant file"
    #log_array.append(logfile_text)
    #format_check, format_log_out = TFDP_format_check(
    #    panel_dataframe, var_df, fmt, file_name, is_mixed, logfile
    #)
    #if format_check != 0:
    #    message = (
    #        "Quick format check may have found inaccuracies: run check_format module"
    #    )
    #    log_array.append(message)
    #    atexit.register(simple_log, log_array, file_name, format_log_out)
    #    exit(message)
    #else:
    #    logfile_text = "SNP panel file is correctly formatted"
    #    log_array.append(logfile_text)
    format_check = 0
    format_log_out = logfile
    if affy_flag is False:
        panel_dataframe_updated = update_snp_panel_names(panel_dataframe)
    else:
        panel_dataframe_updated = panel_dataframe
    if logfile is not None:
        logging = simple_log(log_array, file_name, logfile)
    return format_check, format_log_out, panel_dataframe_updated, var_df


def subsampling_panel_and_varframe(panel_dataframe, variant_sub_dataframe, sample):
    positional_df = panel_dataframe[["Name", sample]]

    positional_df = pd.merge(
        left=positional_df, right=variant_sub_dataframe, on="Name", how="left"
    )
    # Remove rows where there is no positional info (BLAST_chromosome and BLAST_position have '.' values)
    pos_df_nullvals_only = positional_df[positional_df.BLAST_chromosome == "."]
    pos_df_nullvals_removed = positional_df[positional_df.BLAST_chromosome != "."]
    pos_df_nullvals_removed2 = pos_df_nullvals_removed[
        pos_df_nullvals_removed.BLAST_position != "."
    ]
    pos_df_list = [pos_df_nullvals_removed2, pos_df_nullvals_only]
    pos_df_dict = {sample: pos_df_list}
    return pos_df_dict


def get_snp_panel_positional_info(
    panel_dataframe, var_dataframe, can_we_thread, n_threads, logfile, file_name
):
    """
    Creates a dataframe for each animal containing BLAST chromosome & position info from the provided var dataframe
    :param panel_dataframe: SNP user panel df obtained from the function check_input_snp_panel
    :param var_dataframe: matching variant file df obtained from the function check_input_snp_panel and
    associated functions var_match and get_var_df
    :param can_we_thread: (bool) is threading allowed on machine
    :param n_threads: number of threads to use
    :param logfile: logfile
    :param file_name: name of snp file for logging purposes
    :return: hash of animal-specific dataframes containing position info
    """
    log_array = []
    timestr = time.strftime("%H:%M:%S")
    logfile_text = (
        timestr + " ..... Getting chromosome and position information for SNPs "
    )
    log_array.append(logfile_text)
    samples = list(panel_dataframe.columns)
    samples.remove("Name")
    variant_sub_dataframe = var_dataframe[
        ["Name", "BLAST_chromosome", "BLAST_position"]
    ]
    positional_info_dict = {}
    results_list = []
    if can_we_thread is True:
        with cf.ProcessPoolExecutor(max_workers=n_threads) as executor:
            pos_df_dict = {
                executor.submit(
                    subsampling_panel_and_varframe,
                    panel_dataframe,
                    variant_sub_dataframe,
                    each_sample,
                ): each_sample
                for each_sample in samples
            }
            for n in cf.as_completed(pos_df_dict):
                data = n.result()
                results_list.append(data)
        for each_result in results_list:
            key_list = list(each_result.keys())
            each_sample = key_list[0]
            positional_info_dict.update({each_sample: each_result[each_sample]})
    else:
        for each_sample in samples:
            pos_df_dict = subsampling_panel_and_varframe(
                panel_dataframe, variant_sub_dataframe, each_sample
            )
            positional_info_dict.update({each_sample: pos_df_dict[each_sample]})

    # for each_sample in samples:
    #    positional_df = panel_dataframe[['Name', each_sample]]
    #
    #        positional_df = pd.merge(left=positional_df, right=variant_sub_dataframe, on='Name', how='left')
    #        # Remove rows where there is no positional info (BLAST_chromosome and BLAST_position have '.' values)
    #        pos_df_nullvals_only = positional_df[positional_df.BLAST_chromosome == '.']
    #        pos_df_nullvals_removed = positional_df[positional_df.BLAST_chromosome != '.']
    #        pos_df_nullvals_removed = pos_df_nullvals_removed[pos_df_nullvals_removed.BLAST_position != '.']
    #        positional_info_dict.update({each_sample: [pos_df_nullvals_removed, pos_df_nullvals_only]})
    if logfile is not None:
        logging = simple_log(log_array, file_name, logfile)
    return positional_info_dict, logfile
