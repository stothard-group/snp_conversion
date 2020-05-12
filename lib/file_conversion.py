#!/usr/bin/env python3

# Converts between formats - This version uses the information from the variant file only, not the manifest files
# Version 2, to try to make the program more modular and deal with missing long data better

import argparse
import os
import sys
import lib.file_format_checker as fc
import lib.file_parsing as fp
import pandas as pd
import warnings
import time
import concurrent.futures as cf
import lib.variant_file_finder as vff
import lib.make_logs as make_logs
import lib.errors as errors
from functools import reduce
from toolz import interleave
from pandas.errors import ParserError
import platform
import numpy as np
import atexit


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


def get_input_snp_files(input_dir, file_list):
    """
    Gets a list of input SNP panel files, also checks if the files are text files
    :param input_dir: directory containing input files
    :param file_list: list of files supplied by user (may be None)
    :return: list of files to convert
    """
    if file_list is None:
        #  get all files in directory
        snp_files_array = os.listdir(input_dir)
    else:
        # get all files in directory using file list
        snp_files_array = file_list.split(",")
        # check that these files actually exist in the input directory
        input_not_exists = errors.input_file_name_error(input_dir, snp_files_array)
        if input_not_exists:
            exit(
                "File(s) "
                + ", ".join(input_not_exists)
                + " are not found in "
                + input_dir
            )
        else:
            pass
    # check to exclude any non .txt files
    non_txt_files = errors.non_txt_fmt_input_files(snp_files_array)
    if non_txt_files:
        warnings.warn(
            "The following files are not .txt files are are excluded from analysis: "
            + ", ".join(non_txt_files),
            stacklevel=3,
        )
        for non_txt in non_txt_files:
            snp_files_array.remove(non_txt)
    else:
        pass
    # check that there are files left in snp_files_array and conversion dir
    if snp_files_array:
        pass
    else:
        exit("No matrix files in " + input_dir)
    return snp_files_array


def generate_outfile_name(file, output_suffix, logfile):
    """
    Get the name of the output file for the converted panel
    :param file: file name
    :param output_suffix: user-supplied suffix
    :param logfile: logfile
    :return: name of output file
    """
    log_array = []
    names = file.rpartition(".")
    basename = names[0]
    outfile_name = basename + "." + output_suffix + ".txt"
    if os.path.exists(outfile_name):
        warning_text = (
            "Output file " + outfile_name + " already exists - deleting output file"
        )
        warnings.warn(warning_text, stacklevel=3)
        os.remove(outfile_name)
        log_array.append(warning_text)
    if logfile is not None:
        logging = make_logs.simple_log(log_array, file, logfile)
    return outfile_name, logfile


def read_input_file_to_df(file_path, file, file_type, logfile):
    """
    Reads input file to dataframe
    :param file_path: path to input file
    :param file: just the input file name
    :param file_type: user-supplied input file type
    :param logfile: logfile
    :return: dataframe containing input panel
    :return: AB dataframe for writing
    :return: number of header rows to skip
    :return: header dict for writing to output file
    """
    # Reading in file
    timestr = time.strftime("%H:%M:%S")
    logfile_text = timestr + " ..... Reading in input file"
    log_array = [logfile_text]
    affy_flag = False
    if file_type == "affymetrix":
        affy_flag = True
    # check for affymetrix file
    affy_flag = fc.affy_test(file_path, file, file_type, affy_flag)
    if file_type == "affymetrix" or affy_flag is True:
        affy_df = pd.read_csv(file_path, sep="\t", mangle_dupe_cols=True)
        # make dataframe look like the illumina one (remove AB columns)
        header_row = list(set(affy_df.columns))
        header_row.remove("probeset_id")
        # quick check to make sure header row is really affy format
        head_check = fc.affy_head_check(header_row)
        if head_check is False:
            message = "Affymetrix file may not be properly formatted: check columns"
            log_array.append(message)
            warnings.warn(message, stacklevel=3)
        else:
            pass
        affy_header_dict = {}
        for value in header_row:
            new_affy_df = affy_df[["probeset_id", value]].copy()
            ab_list = ["AA", "AB", "BB", "NoCall"]
            test_ab_df = ~new_affy_df.iloc[:, 1:].isin(ab_list)
            not_ab = sum(test_ab_df.sum(axis=1))
            affy_header_dict.update({value: not_ab})
        affy_AB_df = affy_df.copy()  # keep a copy just in case
        for keys in affy_header_dict:
            if affy_header_dict[keys] == 0:
                affy_df.drop(columns=keys, inplace=True)
            else:
                affy_AB_df.drop(columns=keys, inplace=True)
        affy_df.rename(columns={"probeset_id": "SNP Name"}, inplace=True)
        df = affy_df.copy()
        header_row_num = int
        header_dict = {}
    else:  # file_type not affymetrix
        header_row_num, header_dict = fp.parse_header(file_path)
        try:
            df = pd.read_csv(file_path, skiprows=header_row_num, sep="\t")
        except ParserError:
            message = "First line might not contain column names - check formatting"
            log_array.append(message)
            atexit.register(make_logs.simple_log, log_array, file, logfile)
            sys.exit()
        affy_AB_df = None
        if file_type != "LONG":
            df.rename(columns={"Unnamed: 0": "SNP Name"}, inplace=True)
    if logfile is not None:
        logging = make_logs.simple_log(log_array, file, logfile)
    return df, affy_AB_df, header_row_num, header_dict, logfile


def format_check_parse_errors(correct_file, determined_ft, file_type, file, logfile):
    """
    Parses the outcome of file_format_check
    :param correct_file: (bool) whether the file is formatted correctly
    :param determined_ft: filetype as determined by FFC
    :param file_type: user input file type
    :param file: user input panel file
    :param logfile: logfile
    :return: (bool) whether we can continue with the conversion
    """
    log_array = []
    if correct_file is True:
        if determined_ft == file_type or file_type == "mixed":
            pass  # everything is cool and we just want to close the file, and split it
        elif determined_ft != file_type and file_type != "mixed":
            message = "The user-input file type does not match the file type determined by the check_format module"
            warnings.warn(message, stacklevel=4)
            log_array.append(message)
        else:
            pass
    else:
        # Do not do conversion if file does not have 100% correct SNPs
        message = "Could not complete conversion of " + file
        log_array.append(message)
        atexit.register(make_logs.simple_log, log_array, file, logfile)
        sys.exit()
    if logfile is not None:
        logging = make_logs.simple_log(log_array, file, logfile)
    return True, logfile


def gen_long_output_col_names(out_type):
    out_sample1 = "Pos1"
    out_sample2 = "Pos2"
    if out_type == "TOP":
        out_sample1 = "Allele1 - Top"
        out_sample2 = "Allele2 - Top"
    elif out_type == "FWD":
        out_sample1 = "Allele1 - Forward"
        out_sample2 = "Allele2 - Forward"
    elif out_type == "AB":
        out_sample1 = "Allele1 - AB"
        out_sample2 = "Allele2 - AB"
    elif out_type == "DESIGN":
        out_sample1 = "Allele1 - Design"
        out_sample2 = "Allele2 - Design"
    elif out_type == "PLUS":
        out_sample1 = "Allele1 - Plus"
        out_sample2 = "Allele2 - Plus"
    else:
        message = (
            "Allele type is not Top, Forward, Design, or Plus. Remove before proceeding"
        )
        exit(message)
    return out_sample1, out_sample2


def write_header(outfile, h_dict):
    """
    Writes the header data to the output file
    :param outfile: outfile name
    :param h_dict: header dict
    :return:
    """
    with open(outfile, "a") as output_write:
        output_write.write("[Header]\n")
        for member in h_dict:
            if member != "Content":
                output_write.write(member + "\t" + h_dict[member] + "\n")
            else:
                output_write.write(member + "\t\t" + h_dict[member] + "\n")
        output_write.write("[Data]\n")
    output_write.close()
    return True


def easy_long_conversion(
    file_df, out_1, out_2, outfile_name, header_dict, file, logfile
):
    """
    Performs an easy long -> matrix conversion, as long as the output type is found within the long file
    :param file_df: dataframe containing snp input
    :param out_1: output type 1 from gen_long_output_col_names
    :param out_2: output type 2 from gen_long_output_col_names
    :param outfile_name: name of output file
    :param header_dict: header dict for writing the new file
    :param file: file name for logging purposes
    :param logfile: logfile
    :return: (bool) The output file has been written (true)
    """
    # Performing easy long file conversion
    timestr = time.strftime("%H:%M:%S")
    logfile_text = timestr + " ..... Performing file conversion"
    log_array = [logfile_text]

    sub_long_df = file_df[["SNP Name", "Sample ID", out_1, out_2]]
    sample_list = sub_long_df["Sample ID"].to_numpy()
    unique_samples = np.unique(sample_list)
    mapping = {}
    sample_position_dict = {}
    for i in unique_samples:
        mapping[i] = np.where(sample_list == i)[0]
    for each in mapping:
        sample_position_dict.update({each: mapping[each][0]})
    sorted_unique_samples = sorted(unique_samples, key=sample_position_dict.get)
    pos1 = sub_long_df.columns.values[2]
    pos2 = sub_long_df.columns.values[3]
    sub_sample_df_list = []
    for sample in sorted_unique_samples:
        sub_sample_df = sub_long_df[["SNP Name"]].copy()
        sub_sample_df[sample] = sub_long_df[sub_long_df["Sample ID"] == sample][
            [pos1, pos2]
        ].apply(lambda x: "".join(x), axis=1)
        sub_sample_df.dropna(inplace=True)
        sub_sample_df_list.append(sub_sample_df)
    output_df = reduce(
        lambda left, right: pd.merge(left, right, on=["SNP Name"], how="outer"),
        sub_sample_df_list,
    )
    output_df.fillna("--", inplace=True)
    output_df.rename(columns={"SNP Name": ""}, inplace=True)
    # Write header
    timestr = time.strftime("%H:%M:%S")
    logfile_text = timestr + " ..... Writing data to output file"
    log_array.append(logfile_text)
    header_file = write_header(outfile_name, header_dict)
    # Write file
    col_list = list(output_df)
    output_df.to_csv(outfile_name, index=None, mode="a", header=col_list, sep="\t")
    if logfile is not None:
        logging = make_logs.simple_log(log_array, file, logfile)
    return True, logfile


def long_as_matrix(file_df, file, logfile):
    """
    Creates a matrix of any input type within the LONG file (preferentially Forward) to be used for conversion to a
    matrix type that for whatever reason is not in the long file
    :param file_df: input file dataframe
    :param file: file name for logging purposes
    :param logfile: logfile
    :return: new matrix dataframe
    :return: matrix input type (FWD, TOP, etc.)
    """
    # Long input file parsing for conversion
    timestr = time.strftime("%H:%M:%S")
    logfile_text = timestr + " ..... Parsing long file for conversion"
    log_array = [logfile_text]

    # Look for FWD format first
    if "Allele1 - Forward" in file_df.columns:
        df_type = "FWD"
    elif "Allele1 - Top" in file_df.columns:
        df_type = "TOP"
    elif "Allele1 - Design" in file_df.columns:
        df_type = "DESIGN"
    elif "Allele1 - Plus" in file_df.columns:
        df_type = "PLUS"
    else:
        message = "Input file column names are not the expected format"
        log_array.append(message)
        atexit.register(make_logs.simple_log, log_array, file, logfile)
        sys.exit()
    in_1, in_2 = gen_long_output_col_names(df_type)
    matrix_df = file_df[["SNP Name", "Sample ID", in_1, in_2]].copy()
    sample_list = matrix_df["Sample ID"].to_numpy()
    unique_samples = np.unique(sample_list)
    new_matrix = matrix_df[["SNP Name"]].copy()
    matrix_df_list = []
    for samples in unique_samples:
        sub_matrix = matrix_df[matrix_df["Sample ID"] == samples].copy()
        new_matrix[samples] = sub_matrix[[in_1, in_2]].apply(
            lambda x: "".join(x), axis=1
        )
        sub_matrix_2 = new_matrix[["SNP Name", samples]].copy()
        sub_matrix_2.dropna(inplace=True)
        # print(sub_matrix_2)
        matrix_df_list.append(sub_matrix_2)
    matrix_out = reduce(
        lambda left, right: pd.merge(left, right, on=["SNP Name"], how="outer"),
        matrix_df_list,
    )
    if logfile is not None:
        logging = make_logs.simple_log(log_array, file, logfile)
    return matrix_out, df_type, logfile


def get_df_colnames(in_type, out_type2):
    """
    Generates column names for the conversion dataframe based on the input and desired output types
    :param in_type: the determined_ft or matrix_type
    :param out_type2: the desired filetype
    :return: the names of the input and output filetypes
    """
    if in_type == "FWD":
        input_a = "FORWARD_A"
        input_b = "FORWARD_B"
    elif in_type == "TOP":
        input_a = "TOP_A"
        input_b = "TOP_B"
    elif in_type == "DESIGN":
        input_a = "DESIGN_A"
        input_b = "DESIGN_B"
    elif in_type == "AB":
        input_a = "NA"
        input_b = "NA"
    elif in_type == "PLUS":
        input_a = "PLUS_A"
        input_b = "PLUS_B"
    else:
        input_a = "NA"
        input_b = "NA"
    if out_type2 == "TOP":
        output_a = "TOP_A"
        output_b = "TOP_B"
    elif out_type2 == "FWD":
        output_a = "FORWARD_A"
        output_b = "FORWARD_B"
    elif out_type2 == "PLUS":
        output_a = "PLUS_A"
        output_b = "PLUS_B"
    elif out_type2 == "DESIGN":
        output_a = "DESIGN_A"
        output_b = "DESIGN_B"
    elif out_type2 == "AB":
        output_a = "NA"
        output_b = "NA"
    elif out_type2 == "LONG":
        output_a = "NA"
        output_b = "NA"
    else:
        output_a = "NA"
        output_b = "NA"
    return input_a, input_b, output_a, output_b


def conversion(dataframe, var_dataframe, in_t, out_t):
    """
    Performs the actual genotype conversion
    :param dataframe: input dataframe
    :param var_dataframe: variant dataframe we are using for conversion
    :param in_t: input type
    :param out_t: output type
    :return: list containing [snp name, sample, output allele 1, output allele 2] for each row in the input df
    """
    input_a, input_b, output_a, output_b = get_df_colnames(in_t, out_t)
    first_col = ["."]
    second_col = ["."]
    tuples = [tup for tup in dataframe.set_index("SNP Name").itertuples()]
    out_first_second = []
    sample = dataframe.columns[1]
    for tup in tuples:
        name = tup[0]
        cA = tup[2]
        cB = tup[3]

        if in_t == "AB" or out_t == "AB":
            # do stuff here to deal with AB
            if in_t == "AB":
                small_df_out_a = var_dataframe.loc[
                    var_dataframe["Name"] == name, [output_a]
                ]
                small_df_out_b = var_dataframe.loc[
                    var_dataframe["Name"] == name, [output_b]
                ]
                if cA != "-" and cB != "-":
                    first_col = small_df_out_a.values[0]
                    second_col = small_df_out_b.values[0]
                else:
                    first_col = "-"
                    second_col = "-"
            elif out_t == "AB":
                small_df_in_a = var_dataframe.loc[
                    var_dataframe["Name"] == name, [input_a]
                ]
                small_df_in_b = var_dataframe.loc[
                    var_dataframe["Name"] == name, [input_b]
                ]
                if cA == cB:
                    # figure out which is A or B
                    if cA == small_df_in_a.values and cB == small_df_in_a.values:
                        first_col = "A"
                        second_col = "A"
                    elif cA == small_df_in_b.values and cB == small_df_in_b.values:
                        first_col = "B"
                        second_col = "B"
                    elif cA == "-" or cB == "-":
                        first_col = "-"
                        second_col = "-"
                else:
                    first_col = "A"
                    second_col = "B"

        else:
            small_df_in_a = var_dataframe.loc[var_dataframe["Name"] == name, input_a]
            small_df_in_b = var_dataframe.loc[var_dataframe["Name"] == name, input_b]
            small_df_out_a = var_dataframe.loc[var_dataframe["Name"] == name, output_a]
            small_df_out_b = var_dataframe.loc[var_dataframe["Name"] == name, output_b]
            if small_df_in_a.empty and small_df_in_b.empty:
                in_a = "-"
                in_b = "-"
            else:
                in_a = small_df_in_a.values[0]
                in_b = small_df_in_b.values[0]
            if small_df_out_a.empty and small_df_out_b.empty:
                out_a = "-"
                out_b = "-"
            else:
                out_a = small_df_out_a.values[0]
                out_b = small_df_out_b.values[0]
            if cA != "-" or cB != "-":
                if cA == in_a or cA == in_b:
                    if cA == in_a:
                        first_col = out_a
                    elif cA == in_b:
                        first_col = out_b
                else:
                    first_col = "-"
                if cB == in_a or cB == in_b:
                    if cB == in_b:
                        second_col = out_b
                    elif cB == in_a:
                        second_col = out_a
                else:
                    second_col = "-"

            else:
                first_col = "-"
                second_col = "-"
        if first_col == "." and second_col == ".":
            first_col = "-"
            second_col = "-"
        out_first_second.append(
            [name, sample, first_col, second_col]
        )  # Get the program to just print out the found values
    return out_first_second


def list_to_concat_df(output_type, results_list):
    """
    Concatenates a list of dataframes into a single dataframe
    :param output_type: output file type
    :param results_list: list of dataframes to concatenate
    :return: a single, concatenated dataframe
    """
    df_list_concat = []
    out_sample1, out_sample2 = gen_long_output_col_names(output_type)
    for li in results_list:
        mini_df = pd.DataFrame(
            li, columns=["SNP Name", "Sample ID", out_sample1, out_sample2]
        )
        df_list_concat.append(mini_df)
    output_dataframe = pd.concat(df_list_concat, axis=0, ignore_index=True)
    return output_dataframe


def split_and_convert(
    df_for_conversion,
    var_df,
    matrix_type,
    output_type,
    can_we_thread,
    n_threads,
    file,
    logfile,
):
    """
    Splits input dataframe into one dataframe per sample, and parallel converts these dataframes
    :param df_for_conversion: dataframe to be converted
    :param var_df: variant dataframe used for conversion
    :param matrix_type: usually input file type, but could be FWD for affymetrix or matrix-converted long file
    :param output_type: desired output type
    :param can_we_thread: (bool) do we have threading capability
    :param n_threads: how many threads to use
    :param file: file name for logging
    :param logfile: logfile
    :return: converted dataframe
    :return: column names of original dataframe
    """
    # Splitting file by sample
    timestr = time.strftime("%H:%M:%S")
    logfile_text = timestr + " ..... Splitting file by sample"
    log_array = [logfile_text]

    col0_names = list(df_for_conversion.columns)
    if "SNP Name" in col0_names:
        col0_names.remove("SNP Name")
    elif "Name" in col0_names:
        col0_names.remove("Name")
        df_for_conversion.rename(columns={"Name": "SNP Name"}, inplace=True)
    in_a_name, in_b_name, out_a_name, out_b_name = get_df_colnames(
        matrix_type, output_type
    )
    split_df_list = []
    # print(col0_names)
    for col0 in col0_names:
        name = df_for_conversion[["SNP Name", col0]].copy()
        name["{}".format(in_a_name)], name["{}".format(in_b_name)] = zip(
            *name[col0].apply(lambda y: list(y))
        )
        split_df_list.append(name)
        # Write DF to file
        # name_temp_file = file_path + "_" + col0
        # name.to_csv(name_temp_file, sep='\t', index=False)
        # long_temp_list.append(name_temp_file)
    results_dict = {}
    results_array = []
    if output_type == "LONG":
        vals = ["TOP", "FWD", "AB", "DESIGN", "PLUS"]
        message = "Converting to LONG file"
        log_array.append(message)
    else:
        vals = [output_type]
        message = "Converting to " + output_type + " file"
        log_array.append(message)
    for out_type in vals:
        results = []
        results_list = []
        if can_we_thread is True:
            with cf.ProcessPoolExecutor(max_workers=n_threads) as executor:
                allele = {
                    executor.submit(conversion, df, var_df, matrix_type, out_type): df
                    for df in split_df_list
                }
                results.append(allele)
                for n in cf.as_completed(allele):
                    data = n.result()
                    results_list.append(data)
        else:
            for df in split_df_list:
                data = conversion(df, var_df, matrix_type, out_type)
                results_list.append(data)
        output_results = list_to_concat_df(out_type, results_list)
        output_results = output_results.sort_values(by=["Sample ID", "SNP Name"])
        results_dict.update({out_type: output_results})
        results_array.append(output_results)
    timestr = time.strftime("%H:%M:%S")
    message = timestr + " ..... Conversion complete"
    log_array.append(message)
    # If there are multiple output types (out_type == LONG), then concatenate these
    for results in results_array:
        results.set_index(["SNP Name", "Sample ID"], inplace=True)
    output_df = pd.concat(
        results_array, sort=False, axis=1, levels=[["SNP Name", "Sample ID"]]
    )
    # Remove any markers that do not have associated data of the new type
    output_df.replace("", np.nan, inplace=True)
    rm_marker_list = output_df[output_df.isnull().any(axis=1)].index.to_list()
    output_df.fillna("-", inplace=True)
    new_m_list = []
    for marker in rm_marker_list:
        new_m_list.append(marker[0])
    rd_rm_marker_list = list(dict.fromkeys(new_m_list))
    if rd_rm_marker_list:
        marker_str = "\n".join(rd_rm_marker_list)
        message = (
            "The following markers do not have information associated with the output type "
            + output_type
            + " in the variant position file and have been converted to a non-informational position:\n"
            + marker_str
        )
        log_array.append(message)
        print(message)
    else:
        pass
    if logfile is not None:
        logging = make_logs.simple_log(log_array, file, logfile)
    return output_df, col0_names, logfile


def reorder_converted_df(
    converted_df, col0_names, file_df, specified_out, output_type, file, logfile
):
    """
    Reorders converted dataframe to match the original file in terms of samples and snps
    :param converted_df: dataframe that has been converted
    :param col0_names: column names of original dataframe
    :param file_df: original file dataframe
    :param specified_out: user-specified output format (may be AFFY-PLUS)
    :param output_type: type of output file
    :param file: input file name for logging purposes
    :param logfile: logfile
    :return: reordered dataframe
    """
    # Reorder output to match input file
    timestr = time.strftime("%H:%M:%S")
    logfile_text = (
        timestr + " ..... Reorder to match SNP and sample order of input file"
    )
    log_array = [logfile_text]

    new_df_list = []
    converted_df.reset_index(inplace=True)
    if output_type == "LONG":
        for v in col0_names:
            min_df = converted_df.loc[converted_df["Sample ID"] == v]
            new_df_list.append(min_df)
        new_output_df = pd.concat(new_df_list, ignore_index=True)
        new_output_df["GC Score"] = "-"
        new_output_df.rename(columns={"Name": "SNP Name"}, inplace=True)
        final_output_df = new_output_df.copy()
    else:
        not_long_df = file_df[["SNP Name"]].copy()
        converted_df.rename(columns={"SNP Name": ""}, inplace=True)
        for v in col0_names:
            min_df = converted_df.loc[converted_df["Sample ID"] == v]
            new_df_list.append(min_df)
        new_output_df = pd.concat(new_df_list, ignore_index=True)
        # get allele1 and allele2 names
        pos1 = new_output_df.columns.values[2]
        pos2 = new_output_df.columns.values[3]

        new_output_df["output"] = new_output_df[[pos1, pos2]].apply(
            lambda x: "".join(x), axis=1
        )
        col1 = new_output_df.columns[2]
        col2 = new_output_df.columns[3]
        df_list = []
        for col0 in col0_names:
            new_output_df2 = new_output_df[new_output_df["Sample ID"] == col0]
            new_output_df2 = new_output_df2.rename(columns={"output": col0})
            new_output_df2.drop(["Sample ID", col1, col2], axis=1, inplace=True)
            new_output_df2.rename(
                {new_output_df2.columns[0]: "SNP Name"}, axis=1, inplace=True
            )
            new_output_df2.set_index("SNP Name", inplace=True)
            # change all '--' columns in new df to '---' to match affy format
            if specified_out == "AFFY-PLUS":
                new_output_df2 = new_output_df2.replace("--", "---")
            else:
                pass
            df_list.append(new_output_df2)
        pre_final_output = reduce(
            lambda left, right: pd.merge(left, right, on=["SNP Name"], how="outer"),
            df_list,
        )
        pre_final_output.reset_index(inplace=True)
        pre_final_output.rename({"SNP Name": ""}, axis=1, inplace=True)
        final_output_df = pre_final_output.copy()
    if logfile is not None:
        logging = make_logs.simple_log(log_array, file, logfile)
    return final_output_df, logfile


def write_empty_header(outfile_name):
    """
    Writes an empty header if for some reason there is no header in original data file
    :param outfile_name: name of output file
    :return: (bool) True indicating that empty header has been written to output file
    """
    with open(outfile_name, "a") as empty_header_output:
        empty_header_output.write("[Header]\n[Data]\n")
    empty_header_output.close()
    return True


def write_to_outfile(
    reordered_df, specified_out, affy_AB_df, header_dict, file, outfile_name, logfile
):
    """
    Writes the converted dataframe to an output file (can be any matrix or long file format)
    :param reordered_df: dataframe that is converted, reordered, and ready to be written to outfile
    :param specified_out: specified output format by user
    :param affy_AB_df: dict of AB values from affymetrix file (if non-affymetrix input, this var is None)
    :param header_dict: dict containing header information for writing to output file
    :param file: file name for logging purposes
    :param outfile_name: name of the output file to write to
    :param logfile: logfile
    :return: (bool) True indicating that the outfile has been written
    """
    # Writing to output file
    timestr = time.strftime("%H:%M:%S")
    logfile_text = timestr + " ..... Writing converted file to output"
    log_array = [logfile_text]

    if specified_out != "AFFY-PLUS":  # if the file out is illumina
        # Write header
        if header_dict:
            header_file = write_header(outfile_name, header_dict)
        else:
            header_file = write_empty_header(outfile_name)

        # Write dataframes to outfile
        col_list = list(reordered_df)
        reordered_df.to_csv(
            outfile_name, index=None, mode="a", header=col_list, sep="\t"
        )
    else:  # If the file output is an affy format
        # Write header
        col_list = list(reordered_df)
        new_list = ["probeset_id"]
        for sam in col_list:
            new_list.append(sam)
            new_list.append(sam)
        # Write dataframe as affy formatted file
        # use affy_AB_df, change final_output_df 'Name' back to probeset_id
        reordered_df.rename(
            columns={reordered_df.columns[0]: "probeset_id"}, inplace=True
        )
        merged_df = reduce(
            lambda left, right: pd.merge(left, right, on=["probeset_id"], how="outer"),
            [reordered_df, affy_AB_df],
        )[interleave([affy_AB_df, reordered_df])]
        merged_df = merged_df.loc[:, ~merged_df.columns.duplicated()]
        merged_df.columns = [col.replace(".1", "") for col in merged_df.columns]
        merged_df.to_csv(outfile_name, index=None, mode="a", header=True, sep="\t")
    if logfile is not None:
        logging = make_logs.simple_log(log_array, file, logfile)
    return True, logfile


########################################################################################################################


def convert_file(
    input_dir,
    file_list,
    file_type,
    specified_out,
    output_suffix,
    n_threads,
    is_verbose,
    conversion_dir,
    summarize,
    assembly,
    tabular,
    species,
    plink,
):
    """
    Main file conversion module. Converts a list of files between formats.
    :param input_dir: directory containing input files to convert
    :param file_list: list of files to convert in input directory
    :param file_type: input file type
    :param specified_out: output file type (specified by user)
    :param output_suffix: suffix to append to output file(s)
    :param n_threads: number of threads to use for conversion
    :param is_verbose: (bool) write messages to log file
    :param conversion_dir: directory containing variant files to perform conversion with
    :param assembly: assembly name to use
    :param species: species name for variant file
    :return: True (converted file)
    """
    # Error checking and getting the necessary files
    # Error check species and assembly name information
    species_assembly_correlation = errors.assembly_species_error(
        conversion_dir, assembly, species
    )
    if species_assembly_correlation:
        pass
    else:
        exit(
            "The assembly "
            + assembly
            + " is incompatible with species "
            + species
            + ". See README for acceptable "
            "assemblies."
        )
    variant_species_dir = os.path.join(conversion_dir, species)
    variant_assembly_dir = os.path.join(variant_species_dir, assembly)
    variant_files = os.listdir(variant_assembly_dir)

    # Exclude any files without .csv ending (but with correct assembly
    variant_exclude = errors.non_csv_fmt_conversion_files(
        conversion_dir, assembly, species
    )
    if variant_exclude:
        warnings.warn(
            "The following variant files are not .csv files are are excluded from analysis: "
            + ", ".join(variant_exclude),
            stacklevel=3,
        )
        for wrong_var in variant_exclude:
            variant_files.remove(wrong_var)
        if not variant_files:
            exit("No variant conversion files in " + variant_assembly_dir)
    else:
        pass

    # Determine the operating system of the user
    can_we_thread = check_user_platform()

    # Get input files
    snp_files_array = get_input_snp_files(input_dir, file_list)

    # Set some variables
    if specified_out == "AFFY-PLUS":
        output_type = "PLUS"
    else:
        output_type = specified_out

    # Convert each file
    for file in snp_files_array:

        # Get log file name if verbose logging is true
        if is_verbose:
            timestr = time.strftime("%Y%m%d%H%M%S")
            log_suffix = "-" + timestr + ".log"
            log_input = make_logs.get_logname(log_suffix, file)
        else:
            log_input = None

        # Create outfile name
        outfile_name, logfile1 = generate_outfile_name(file, output_suffix, log_input)

        # Read file into dataframe
        file_path = os.path.join(input_dir, file)
        (
            file_df,
            affy_AB_dict,
            header_row,
            header_dict,
            logfile2,
        ) = read_input_file_to_df(file_path, file, file_type, logfile1)

        # Get the correct variant file
        timestr = time.strftime("%H:%M:%S")
        logfile_text_1 = timestr + " ..... Finding correct variant file"
        log_array_1 = [logfile_text_1]
        if logfile2 is not None:
            logging = make_logs.simple_log(log_array_1, file, logfile2)

        converting_file = True
        var_list, logfile2_text, alt_bool = vff.var_match(
            variant_files,
            conversion_dir,
            file_df,
            file,
            converting_file,
            assembly,
            species,
        )
        var_df = vff.get_var_df(
            conversion_dir, var_list[0], assembly, species, alt_bool
        )

        # Check file format
        timestr = time.strftime("%H:%M:%S")
        logfile_text_2 = timestr + " ..... Running file format check"
        log_array_2 = [logfile2_text, logfile_text_2]
        if logfile2 is not None:
            logging = make_logs.simple_log(log_array_2, file, logfile2)

        get_snp_pan = (
            False  ######## CHANGE THIS IF YOU WANT TO PUT SNP PANEL OPTION IN HERE
        )
        internal_sum = False
        format_check_plink = False
        determined_ft, correct_file, logfile3 = fc.file_format_check(
            input_dir,
            file,
            file_type,
            get_snp_pan,
            is_verbose,
            conversion_dir,
            logfile2,
            assembly,
            internal_sum,
            tabular,
            species,
            format_check_plink,
        )
        file_is_correct, logfile4 = format_check_parse_errors(
            correct_file, determined_ft, file_type, file, logfile3
        )

        # Actually do format conversion
        # Deal with long-to-matrix conversion if necessary
        if determined_ft == "LONG":
            # See if we can do an easy long conversion using the data we already have
            out_1, out_2 = gen_long_output_col_names(output_type)
            if out_1 in file_df.columns:
                long_easy_output, logfile5 = easy_long_conversion(
                    file_df, out_1, out_2, outfile_name, header_dict, file, logfile4
                )  # need to test whether the output file is correctly formatted
                df_for_conversion = None
                matrix_type = None
                logfile8 = logfile5
            else:
                # Convert the long file to a matrix-like array for conversion
                df_for_conversion, matrix_type, logfile5 = long_as_matrix(
                    file_df, file, logfile4
                )
        else:
            df_for_conversion = file_df
            if determined_ft == "affymetrix":
                matrix_type = "FWD"
            else:
                matrix_type = determined_ft
        # Split file based on columns into a dict of dataframes and convert
        if df_for_conversion is not None:
            converted_df, column_names, logfile6 = split_and_convert(
                df_for_conversion,
                var_df,
                matrix_type,
                output_type,
                can_we_thread,
                n_threads,
                file,
                logfile4,
            )
            # Reorder converted dataframe to match the original order of samples and snp names
            reordered_df, logfile7 = reorder_converted_df(
                converted_df,
                column_names,
                file_df,
                specified_out,
                output_type,
                file,
                logfile6,
            )
            # Write files reordered_df to output file
            matrix_or_long_hard_output, logfile8 = write_to_outfile(
                reordered_df,
                specified_out,
                affy_AB_dict,
                header_dict,
                file,
                outfile_name,
                logfile7,
            )
        else:
            logfile8 = logfile4
        # Check file format of converted file
        timestr = time.strftime("%H:%M:%S")
        logfile_text_3 = timestr + " ..... Checking for correct conversion"
        log_array_3 = [logfile_text_3]
        if logfile8 is not None:
            logging = make_logs.simple_log(log_array_3, file, logfile8)

        new_input_dir = "."
        get_snp_panel = False
        correct_output_type, correct_output, logfile9 = fc.file_format_check(
            new_input_dir,
            outfile_name,
            specified_out,
            get_snp_panel,
            is_verbose,
            conversion_dir,
            logfile8,
            assembly,
            summarize,
            tabular,
            species,
            plink,
        )
        # Messaging based on whether conversion occurred properly
        log_array = []
        if correct_output is True:
            message = "File " + file + " was converted properly"
            print(message)
            log_array.append(message)
        else:
            message = "The user output file type does not match the file type determined by the file format checker"
            log_array.append(message)
            warnings.warn(message, stacklevel=4)
        # Write to log file if verbose
        if logfile9 is not None:
            logging = make_logs.simple_log(log_array, file, logfile9)
    return True


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Performs file conversion between Illumina matrix and long formats"
    )
    parser.add_argument(
        "--input-dir", type=str, help="directory containing input files"
    )
    parser.add_argument(
        "--file-list",
        type=str,
        help="[optional] comma-separated list of files in the input directory",
    )
    parser.add_argument(
        "--input-format",
        type=str,
        default="mixed",
        choices=["TOP", "FWD", "AB", "PLUS", "DESIGN", "LONG", "mixed"],
        help="Type of file(s) expected: 'TOP', 'FWD', 'AB', 'PLUS', 'DESIGN', 'LONG', or 'mixed'",
    )

    parser.add_argument(
        "--output-format",
        type=str,
        choices=["TOP", "FWD", "AB", "PLUS", "DESIGN", "LONG"],
        help="Type of file(s) to be created: 'TOP', 'FWD', 'AB', 'PLUS', 'DESIGN', 'LONG'"
        "Only one type of output can be specified at a time"
        "PLUS converts the data to the forward strand of the reference genome"
        "LONG refers to the long-format Illumina file",
    )
    parser.add_argument("--output-name", type=str, help="Suffix of output file")
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=2,
        required=False,
        help="[optional] Number of threads to use during conversion (default = 2)",
    )
    parser.add_argument(
        "-v",
        "--verbose-logging",
        action="store_true",
        default=False,
        required=False,
        help="[optional] Write output to both STDOUT and log file",
    )
    parser.add_argument(
        "--key-dir",
        type=str,
        default="variant_position_files",
        help="Directory containing genotype conversion key files (default = variant_position_files)",
    )
    parser.add_argument(
        "--assembly",
        type=str,
        required=True,
        help="Assembly name which must be included in genotype conversion file name",
    )
    parser.add_argument(
        "--species",
        required=True,
        type=str,
        choices=["bos_taurus", "sus_scrofa"],
        help="Organism name",
    )
    parser.add_argument(
        "-s",
        "--summary",
        action="store_true",
        default=False,
        required=False,
        help="Summarize converted SNP file in *_summary.txt file",
    )
    parser.add_argument(
        "--tabular",
        action="store_true",
        default=True,
        required=False,
        help="Output summary file in tabular format (default: False)",
    )
    parser.add_argument(
        "--plink",
        action="store_true",
        default=False,
        required=False,
        help="Creates PLINK flat files (PED and MAP) (default: False)",
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    in_dir = args.input_dir
    f_list = args.file_list
    f_type = args.input_format
    o_type = args.output_format
    o_suffix = args.output_name
    threads = args.threads
    vb = args.verbose_logging
    conv_dir = args.key_dir
    assembly_val = args.assembly

    out = convert_file(
        in_dir,
        f_list,
        f_type,
        o_type,
        o_suffix,
        threads,
        vb,
        conv_dir,
        args.summary,
        assembly_val,
        args.tabular,
        args.species,
        args.plink,
    )
