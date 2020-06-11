#!/usr/bin/env python3
import argparse
import os
import sys
from lib.file_parsing import check_header_data_congruency
from lib.file_parsing import parse_header
import pandas as pd
import warnings
import numpy as np
from lib.variant_file_finder import var_match, get_var_df
import lib.make_logs as make_logs
from lib.write_summary import write_summary
import lib.errors as errors
from pandas.errors import ParserError
import lib.plink_output as make_plink
import time
import atexit

#####################################################################################################################
# This value controls the minimum fraction of correct SNPs required for the program to determine the file input type,
# if the input type is not given (or is 'mixed')
minimum_correct_snp_fraction = 0.95


#####################################################################################################################


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


def get_snp_files(file_list, input_dir, logfile):
    """
    Gets user input files and removes ones that do not end in .txt
    :param file_list: user list of input files (may be None)
    :param input_dir: directory containing input files
    :param logfile: logfile name (may be None)
    :return: list of snp files
    """
    if file_list is None:
        snp_files_array = os.listdir(input_dir)
    else:
        # get all files in directory using file list
        snp_files = file_list
        snp_files_array = snp_files.split(",")
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
    # Check that these are all .txt files and exclude any non-txt files
    non_txt_files = errors.non_txt_fmt_input_files(snp_files_array)
    if non_txt_files:
        message = (
                "The following files are not .txt files are are excluded from analysis: "
                + ", ".join(non_txt_files)
        )
        warnings.warn(message, stacklevel=4)
        for non_txt in non_txt_files:
            snp_files_array.remove(non_txt)
    else:
        pass
    # check that there are files left in snp_files_array
    if snp_files_array:
        pass
    else:
        exit("No matrix files in " + input_dir)
    return snp_files_array, logfile


def get_variant_files(conversion_dir, assembly, species):
    """
    Gets a list of possible variant files using assembly and species information
    :param conversion_dir: variant position files directory
    :param assembly: assembly name
    :param species: species name
    :return: list of variant files
    """

    # Get variant position files
    variant_species_dir = os.path.join(conversion_dir, species)
    variant_assembly_dir = os.path.join(variant_species_dir, assembly)
    variant_files = os.listdir(variant_assembly_dir)
    # Check that the assembly and species combination works
    species_error = errors.assembly_species_error(conversion_dir, assembly, species)
    if not species_error:
        exit("Cannot find assembly name in any variant file in " + conversion_dir)
    # Exclude any files without .csv ending (but with correct assembly)
    variant_exclude = errors.non_csv_fmt_conversion_files(
        conversion_dir, assembly, species
    )
    if variant_exclude:
        for wrong_var in variant_exclude:
            variant_files.remove(wrong_var)
        if not variant_files:
            exit("No variant conversion files in " + variant_assembly_dir)
    else:
        pass
    return variant_files


def affy_test(file_path, file, file_type, affy_flag):
    """
    Determines whether the input file might be AFFY format, especially important if input format is mixed
    :param file_path: path to input snp file
    :param file: input file name
    :param file_type: input file type
    :param affy_flag: (bool) True if AFFY file was specified
    :return: (Bool) whether or not file is affy format
    """
    with open(file_path) as f:
        line = f.readline()
        line = line.rstrip("\n")
        line_array = line.split("\t")
        if "probeset_id" in line:
            if not affy_flag and not file_type == "mixed":
                warnings.warn(
                    "File "
                    + file
                    + " appears to be in AFFY format, but input format is "
                    + file_type,
                    stacklevel=4,
                )
            affy_flag = True
        elif (line == "[Header]" or line == "[Data]") and file_type == "AFFY":
            warnings.warn(
                "File "
                + file
                + " appears to be in Illumina format, but input format is AFFY",
                stacklevel=4,
            )
            affy_flag = False
        elif (len(line_array) >= 3) and line_array[1] == line_array[
            2
        ]:  # This breaks when there is header info with no [Header] designation
            if not affy_flag and not file_type == "mixed":
                warnings.warn(
                    "File "
                    + file
                    + " appears to be in AFFY format, but input format is "
                    + file_type,
                    stacklevel=4,
                )
            affy_flag = True
        else:
            pass
        f.close()
    return affy_flag


def affy_head_check(header):
    """
    Checks the header format to see if AFFY has AB columns
    :param header: header line of AFFY file (first line)
    :return: (Bool) Whether header has AB format columns as well (True)
    """
    simple_list = []
    header_correct = None
    for sample in header:
        if not sample.endswith(".1"):
            simple_list.append(sample)
    for simple in simple_list:
        new_simple = simple + ".1"
        if new_simple in header:
            header_correct = True
        else:
            header_correct = False
    return header_correct


def parse_affy_file(file_path, file, logfile):
    """
    Parses AFFY file input
    :param file_path: path to input file
    :param file: input file name
    :param logfile: logfile
    :return: AFFY dict, affy_dataframe, logfile
    """
    timestr = time.strftime("%H:%M:%S")
    logfile_text = timestr + " ..... Reading in AFFY input file"
    log_array = [logfile_text]
    # Test affy forward values here
    affy_df = pd.read_csv(file_path, sep="\t", mangle_dupe_cols=True)
    # make dataframe look like the illumina one (remove AB columns)
    header_row = list(set(affy_df.columns))
    header_row.remove("probeset_id")
    # quick check to make sure header row is really affy format (2 columns for every sample)
    head_check = affy_head_check(header_row)
    if head_check is False:
        message = "Affymetrix file may not have AB columns"
        warnings.warn(message, stacklevel=4)
        if logfile is not None:
            log_array.append(message)
            logging = make_logs.simple_log(log_array, file, logfile)
    else:
        pass
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
    if logfile is not None:
        logging = make_logs.simple_log(log_array, file, logfile)
    return affy_header_dict, affy_df, logfile


def read_illumina_input(file_path, header_row, file, logfile):
    """
    Reads the Illumina file into a dataframe
    :param file_path: path to input file
    :param header_row: rows to skip for header info
    :param file: file name for logging
    :param logfile: logfile
    :return: Illumina dataframe, (bool) whether file is a long file, logfile
    """
    timestr = time.strftime("%H:%M:%S")
    logfile_text = timestr + " ..... Reading in Illumina input file"
    log_array = [logfile_text]
    with open(file_path, "r") as input_file:
        # Make sure input file will not throw a dataframe error (from an incorrect header structure)
        try:
            df = pd.read_csv(input_file, skiprows=header_row, sep="\t")
        except ParserError:
            message = "First line might not contain column names - check formatting"
            log_array.append(message)
            atexit.register(make_logs.simple_log, log_array, file, logfile)
            sys.exit()
        long_ft = False
        df_columns_long_check = list(df.columns)
        for type_of_file in df_columns_long_check:
            if "Allele1" in type_of_file:
                long_ft = True
        if long_ft is False:
            # Put in a check here to make sure there are the same number of column labels as columns
            if "Unnamed: 0" not in df.columns:
                df.index.names = ["Unnamed: 0"]
                df.reset_index(inplace=True)
            df.rename(columns={"Unnamed: 0": "Name"}, inplace=True)
    input_file.close()
    if logfile is not None:
        logging = make_logs.simple_log(log_array, file, logfile)
    return df, long_ft, logfile


def get_long_df_dict(long_df, header_row, header_dict, numbers_exist, file, logfile):
    """
    Creates a dict containing all dataframes (sample, pos1, pos2) where keys are data types in long file
    :param long_df: dataframe of Long Illumina data
    :param header_row: number of rows to skip in input file
    :param header_dict: dict containing header information
    :param numbers_exist: num_snps and num_samples exist
    :param file: input file name
    :param logfile: logfile
    :return: dict containing all dataframes in long file
    """
    timestr = time.strftime("%H:%M:%S")
    logfile_text = timestr + " ..... Parsing Long file"
    log_array = [logfile_text]
    # make sure header info accurately reflects file
    if header_row != 0 and numbers_exist:
        samples_list = list(set(long_df["Sample ID"]))
        num_samples = len(samples_list)
        num_snps = len(long_df["SNP Name"].unique())
        header_congruency = check_header_data_congruency(
            header_dict, num_samples, num_snps
        )
        if header_congruency:
            append_text = "File " + file
            log_array.append(append_text)
        for each_warning in header_congruency:
            log_array.append(each_warning)
    else:
        pass
    # do stuff here to get the long file into arrays
    type_list = ["TOP", "FWD", "AB", "DESIGN", "PLUS"]
    df_dict = {}
    for ty in type_list:
        out1, out2 = gen_long_output_col_names(ty)
        # Make sure allele type exists
        if out1 not in long_df.columns or out2 not in long_df.columns:
            message = (
                    "Columns " + out1 + " and " + out2 + " might not be in the input file"
            )
            warnings.warn(message, stacklevel=4)
        else:
            sub_long_df = long_df[["SNP Name", "Sample ID", out1, out2]]
            # print new df by sample and save to dict by "file type"
            sample_list = sub_long_df["Sample ID"]
            unique_samples = list(set(sample_list))
            pos1 = sub_long_df.columns.values[2]
            pos2 = sub_long_df.columns.values[3]
            output_df = sub_long_df[["SNP Name"]].copy()
            output_df = output_df.rename(columns={"SNP Name": "Name"})
            for sample in unique_samples:
                sub_sample_df = pd.DataFrame()
                sub_sample_df[sample] = sub_long_df[[pos1, pos2]].apply(
                    lambda x: "".join(x), axis=1
                )
                output_df[sample] = sub_sample_df[sample].apply(lambda v: v)
                df_dict.update({ty: output_df})
    if logfile is not None:
        logging = make_logs.simple_log(log_array, file, logfile)
    return df_dict, logfile


def check_illumina_df_header(
        illumina_df, header_row, header_dict, numbers_exist, file, logfile
):
    """
    Checks Illumina dataframe header
    :param illumina_df: dataframe of matrix Illumina data
    :param header_row: number of rows to skip in input file
    :param header_dict: dict containing header information
    :param numbers_exist: num_snps and num_samples exist
    :param file: input file name
    :param logfile: logfile
    """
    log_array = []
    if header_row != 0 and numbers_exist:
        num_samples = len(illumina_df.columns) - 1
        num_snps = len(illumina_df.index)
        header_congruency = check_header_data_congruency(
            header_dict, num_samples, num_snps
        )
        if header_congruency:
            append_text = "File " + file
            log_array.append(append_text)
        for each_warning in header_congruency:
            log_array.append(each_warning)
    if logfile is not None:
        logging = make_logs.simple_log(log_array, file, logfile)
    illumina_df_checked = True
    return illumina_df_checked, logfile


def TFDP_format_check(dataframe, var_dataframe, fmt, in_file, is_mixed, log):
    """
    Beefy format checking module. This does all the work, comparing input dataframe to variant dataframe
    :param dataframe: dataframe from input file
    :param var_dataframe: variant dataframe
    :param fmt: dataframe format (one of FWD, PLUS, TOP, DESIGN)
    :param in_file: input file name
    :param is_mixed: (Bool) whether the input type was set by the user to "mixed" or unspecified
    :param log: logfile
    :return: (int) a count of the non-matching positions, logfile
    """
    # The A/B that are variable should be equal to TOP_A and TOP_B, respectively, if file_format is TOP
    # The A/B that are variable should be equal to FORWARD_A and FORWARD_B, respectively, if file_format is FWD
    # Create dataframe with ID, TOP_A, and TOP_B
    col_to_keep = []
    if fmt == "TOP":
        a_column = "TOP_A"
        b_column = "TOP_B"
        keep = ["TOP_A", "TOP_B"]
    elif fmt == "FWD":
        a_column = "FORWARD_A"
        b_column = "FORWARD_B"
        keep = ["FORWARD_A", "FORWARD_B"]
    elif fmt == "PLUS":
        a_column = "PLUS_A"
        b_column = "PLUS_B"
        keep = ["PLUS_A", "PLUS_B"]
    elif fmt == "DESIGN":
        a_column = "DESIGN_A"
        b_column = "DESIGN_B"
        keep = ["DESIGN_A", "DESIGN_B"]
    for i in keep:
        col_to_keep.append(i)
    # Keep only the database entries that are in the user input file
    if "SNP Name" in list(dataframe.columns):
        dataframe.rename(columns={"SNP Name": "Name"}, inplace=True)

    namelist = list(dataframe["Name"])
    subset_var_dataframe = var_dataframe.query("Name in @namelist").copy()
    # print(subset_var_dataframe)
    col_to_keep.append("Name")
    var_df_cols = list(var_dataframe)
    col_to_drop = list(set(var_df_cols) - set(col_to_keep))
    subset_var_dataframe.drop(col_to_drop, inplace=True, axis=1)
    # print(subset_var_dataframe)

    # Split each sample into a separate dataframe, and separate [NN] into two columns for comparison
    column_matches = []
    col0_names = list(dataframe.columns)
    col0_names.remove("Name")
    all_non_matches = 0
    log_exists = log
    for col0 in col0_names:
        name = dataframe[["Name", col0]].copy()
        name["{}".format(a_column)], name["{}".format(b_column)] = zip(
            *name[col0].apply(lambda x: list(x))
        )
        # Keep only rows that have differences (not AA/TT/CC/GG/--)
        nn_list = ["AA", "TT", "CC", "GG", "II", "DD", "--", "---"]
        name_differences_df = name[~name[col0].isin(nn_list)].copy()
        name_differences_df.drop(col0, axis=1, inplace=True)
        # Keep only homozygous rows (AA/TT/CC/GG)
        nn_list_2 = ["AA", "TT", "CC", "GG"]
        name_same_df = name[name[col0].isin(nn_list_2)].copy()
        name_same_df.drop(col0, axis=1, inplace=True)
        # Subset var dataframe again, based on having dropped the rows above (for heterozygous only)
        indiv_namelist = list(name_differences_df["Name"])
        indiv_subset_var_dataframe = subset_var_dataframe.query(
            "Name in @indiv_namelist"
        )

        # Compare dataframe with user input (name_differences_df) to database data (indiv_subset_var_dataframe)
        compare_df = pd.merge(
            indiv_subset_var_dataframe,
            name_differences_df,
            on=["Name", a_column, b_column],
            how="left",
            indicator="Exist",
        )
        compare_df["Exist"] = np.where(compare_df.Exist == "both", 0, 1)
        # check if reversing the user input (AB -> BA) matches the information in the variant dataframe, and update the
        # Exist column accordingly
        user_df_reversed = name_differences_df.copy()
        a_column_new = b_column + "_new"
        b_column_new = a_column + "_new"
        user_df_reversed.rename(
            columns={a_column: a_column_new, b_column: b_column_new}, inplace=True
        )
        user_df_reversed.columns = user_df_reversed.columns.str.replace("_new", "")
        reverse_compare_df = pd.merge(
            indiv_subset_var_dataframe,
            user_df_reversed,
            on=["Name", a_column, b_column],
            how="left",
            indicator="Exist",
        )
        reverse_compare_df["Exist"] = np.where(reverse_compare_df.Exist == "both", 0, 1)
        # Merge "Exist" columns, and if either is 0, then super_exist = 0
        super_compare_df = pd.merge(
            compare_df,
            reverse_compare_df,
            on=["Name", a_column, b_column, "Exist"],
            how="left",
            indicator="SuperExist",
        )
        super_compare_df["SuperExist"] = np.where(
            super_compare_df.SuperExist == "both", 1, 0
        )
        super_compare_df.drop(columns="Exist", inplace=True)
        # Write to log if there are inconsistent values
        if (
                not super_compare_df[super_compare_df["SuperExist"] == 1].empty
                and is_mixed is False
        ):
            # check that we still only have allowed values
            write_log = make_logs.inconsistent_values(
                super_compare_df[super_compare_df["SuperExist"] == 1],
                col0,
                name_differences_df,
                in_file,
                log_exists,
            )
            log_exists = write_log
        match = super_compare_df["SuperExist"].sum()
        column_matches.append(match)
        # Check that homozygous alleles are still allowed based on the var dataframe
        name_same_df.set_index("Name", inplace=True)
        svd_reindexed = subset_var_dataframe.set_index("Name")

        # merge svd_reindexed with name_same_df on name column and keep identical names
        # make a master merge that finds if there are any 'wrong' homozygous values not found in the var df at all
        merged_svd_namesame = pd.merge(
            left=name_same_df, right=svd_reindexed, how="left", on="Name"
        )
        merged_svd_namesame["merge_A"] = np.where(
            (merged_svd_namesame.iloc[:, 0] == merged_svd_namesame.iloc[:, 2]),
            merged_svd_namesame.iloc[:, 2],
            0,
        )
        merged_svd_namesame["merge_B"] = np.where(
            (merged_svd_namesame.iloc[:, 0] == merged_svd_namesame.iloc[:, 3]),
            merged_svd_namesame.iloc[:, 3],
            0,
        )
        merged_svd_namesame["master_merge"] = np.where(
            (merged_svd_namesame["merge_A"] == merged_svd_namesame["merge_B"]), 1, 0
        )
        # make a new df that only has rows where master_merge == 1
        bad_merge_df = merged_svd_namesame[merged_svd_namesame.master_merge != 0]
        if bad_merge_df.empty is False:
            for idx, row in bad_merge_df.iterrows():
                column_matches.append(1)
                hom_allele = row[0]
                var_alleles = [row[2], row[3]]
                if is_mixed is False:
                    write_hom_log = make_logs.inconsistent_values_homozygous(
                        in_file, log_exists, hom_allele, var_alleles, idx, col0
                    )
                    log_exists = write_hom_log
                else:
                    pass
                acceptable_alleles = ["A", "T", "G", "C", "I", "D", "-"]
                if hom_allele not in acceptable_alleles:
                    exit(
                        "Unexpected allele "
                        + hom_allele
                        + " found at "
                        + name_same_df[idx]
                    )

    all_non_matches = sum(column_matches)
    return all_non_matches, log_exists


def check_affy_format(
        affy_df,
        var_df,
        specified_file_type,
        file,
        file_type,
        summary_inconsistency_value,
        logfile,
):
    """
    Check format of an AFFY dataframe using TFDR
    :param affy_df: AFFY input dataframe
    :param var_df: variant dataframe
    :param specified_file_type: file type specified by user
    :param file: file name for logging
    :param file_type: same as specified file type, unless AFFY-PLUS in which case this is PLUS
    :param summary_inconsistency_value: number of inconsistencies in dataframe
    :param logfile: logfile
    :return: determined filetype, (bool) whether file is correctly formatted, number of inconsistencies, logfile
    """
    timestr = time.strftime("%H:%M:%S")
    logfile_text = timestr + " ..... Checking Affymetrix file"
    log_array = [logfile_text]
    if file_type == "mixed":
        is_mixed = True
    else:
        is_mixed = False
    if specified_file_type == "AFFY-PLUS":
        affy_fmt = "PLUS"
        test_both = False
    elif specified_file_type == "AFFY":
        affy_fmt = "FWD"
        test_both = False
    else:
        test_both = True
    if is_mixed or test_both:

        affy_fmt = "FWD"
        format_check_fwd, logfile_name = TFDP_format_check(
            affy_df, var_df, affy_fmt, file, is_mixed, logfile
        )
        affy_fmt = "PLUS"
        format_check_plus, logfile_name = TFDP_format_check(
            affy_df, var_df, affy_fmt, file, is_mixed, logfile
        )
        if format_check_fwd == 0:
            summary_inconsistency_value = format_check_fwd
            affy_fmt = "FWD"
            format_check = format_check_fwd
        elif format_check_plus == 0:
            summary_inconsistency_value = format_check_plus
            affy_fmt = "PLUS"
            format_check = format_check_plus
        else:
            summary_inconsistency_value = format_check_fwd
            affy_fmt = "FWD"
            format_check = format_check_fwd

    else:
        format_check, logfile_name = TFDP_format_check(
            affy_df, var_df, affy_fmt, file, is_mixed, logfile
        )
        summary_inconsistency_value = format_check
    if format_check > 0:
        if is_mixed is False:
            out_string = (
                "Affymetrix file is incorrectly formatted - see log for details."
            )
        else:
            out_string = "Affymetrix file is incorrectly formatted - re-run with --input-format AFFY for details"
        warnings.warn(out_string, stacklevel=4)
        filetype = None
        correct_format = False
    else:  # format check is correct
        if affy_fmt == "PLUS":
            out_string = "File " + file + " is correctly formatted in Affymetrix PLUS format"
            print(out_string)
        elif affy_fmt == "FWD":
            out_string = "File " + file + " is correctly formatted in Affymetrix format"
            print(out_string)
        log_array.append(out_string)
        filetype = "AFFY"
        correct_format = True
    if logfile is not None:
        logging = make_logs.simple_log(log_array, file, logfile)
    return filetype, correct_format, summary_inconsistency_value, logfile


def AB_check(dataframe, file, is_mixed, log):
    """
    Checks AB format matrix dataframes
    :param dataframe: dataframe containing AB data
    :param file: input file name
    :param is_mixed: (bool) treat this dataframe as if type is unknown
    :param log: logfile
    :return: inconsistencies count, logfile (ab_warn)
    """
    ab_list = ["AA", "AB", "BB", "--", "NoCall"]
    new_df = ~dataframe.iloc[:, 1:].isin(ab_list)
    format_check = sum(new_df.sum(axis=1))
    ab_warn = None
    if format_check > 0 and is_mixed is False:
        ab_warn = make_logs.ab_warning(dataframe, new_df, file, log)
    return format_check, ab_warn


def long_format_consistency_check(long_df, variant_df):
    """
    Checks for inequivalency in Long file.
    :param long_df: dataframe containing all long data (not the df dict)
    :param variant_df: variant dataframe
    :return: count of consistent values, count of inconsistent values, dict containing the inequivalent lines for logging
    """
    # Consistency check: pick first 2 columns, figure out if they are AA/BB/AB, make sure that the other ones match this
    long_df.set_index("SNP Name", inplace=True)
    # reorder dataframe
    long_df_cols = list(long_df.columns.values)
    new_cols = [
        "Sample ID",
        "Allele1 - Top",
        "Allele2 - Top",
        "Allele1 - Forward",
        "Allele2 - Forward",
        "Allele1 - AB",
        "Allele2 - AB",
        "Allele1 - Design",
        "Allele2 - Design",
        "Allele1 - Plus",
        "Allele2 - Plus",
    ]
    extra_cols = list(set(long_df_cols) - set(new_cols))
    if "GC Score" in long_df_cols:
        extra_cols.remove("GC Score")
    if extra_cols:
        extra_list = ", ".join(extra_cols)
        message = "Columns " + extra_list + " have been excluded from analysis"
        warnings.warn(message, stacklevel=4)
    else:
        pass
    new_cols_excluded = []
    for col_name in new_cols:
        if col_name in long_df.columns:
            new_cols_excluded.append(col_name)
    col_name_dict = {
        "Allele1 - Forward": "FORWARD",
        "Allele1 - Top": "TOP",
        "Allele1 - AB": "AB",
        "Allele1 - Design": "DESIGN",
        "Allele1 - Plus": "PLUS",
    }
    sorted_long_df = long_df[new_cols_excluded]
    new_cols_no_sample_id = new_cols_excluded
    new_cols_no_sample_id.remove("Sample ID")
    new_cols_half = [x for x in new_cols_no_sample_id if "2" not in x]
    consistency_count = 0
    inconsistency_count = 0
    return_inequiv_dict = {}
    array_len = len(new_cols_excluded)
    array_list = list(range(1, array_len))
    for index, row in sorted_long_df.iterrows():
        row_array = [row[i] for i in array_list]
        var_row = variant_df.loc[variant_df["Name"] == index, :]
        long_df_top_A = row_array[0]
        long_df_top_B = row_array[1]
        var_df_top_A = var_row["TOP_A"].values[0]
        var_df_top_B = var_row["TOP_B"].values[0]
        # print(long_df_top_A, long_df_top_B, var_df_top_A, var_df_top_B)
        a_or_b = ""
        top_empty = False
        if all(x == "-" for x in row_array):
            continue
        else:
            if long_df_top_A == long_df_top_B and long_df_top_A != "-":
                if long_df_top_A == var_df_top_A:
                    a_or_b = "A"
                elif long_df_top_A == var_df_top_B:
                    a_or_b = "B"
                else:
                    inconsistency_count = inconsistency_count + 1
                    return_inequiv_dict.update({index: row})
                    continue
            elif long_df_top_A != long_df_top_B:
                a_or_b = "AB"
            elif long_df_top_A == "-" and long_df_top_B == "-":
                top_empty = True
        if top_empty is True:
            warnings.warn(
                "No genotypes for TOP strand at "
                + index
                + ", "
                + row_array[0]
                + "- this could indicate an error in the data"
            )
        if a_or_b == "A":
            a_val = "A"
            b_val = "A"
        elif a_or_b == "B":
            a_val = "B"
            b_val = "B"
        elif a_or_b == "AB":
            a_val = "A"
            b_val = "B"
        else:
            a_val = "-"
            b_val = "-"
        sub_var_row = []

        for col2 in new_cols_half:
            if col2 == "Allele1 - AB":
                var_row_a_val = a_val
                var_row_b_val = b_val
            else:
                var_row_a_val = var_row[
                    "{}_{}".format(col_name_dict[col2], a_val)
                ].values[0]
                var_row_b_val = var_row[
                    "{}_{}".format(col_name_dict[col2], b_val)
                ].values[0]
            sub_var_row.append(var_row_a_val)
            sub_var_row.append(var_row_b_val)
        # If user input row has dashes in any columns, delete the corresponding columns from the variant row
        non_dash_index = []
        for element in range(0, len(row_array)):
            if row_array[element] != "-":
                non_dash_index.append(element)
        # make a new list with only the positions that are not dashes
        new_sub_var_row = []
        new_row_array = []
        for pos in non_dash_index:
            new_sub_var_row.append(sub_var_row[pos])
            new_row_array.append(row_array[pos])
        # Check for consistency between var file list and user input list
        if new_sub_var_row == new_row_array:
            consistency_count = consistency_count + 1
        else:
            inconsistency_count = inconsistency_count + 1
            # make a dict to return of inequivalent lines
            return_inequiv_dict.update({index: row})
    return consistency_count, inconsistency_count, return_inequiv_dict


def check_long_format(long_df_dict, long_df, var_df, file, logfile):
    """
    Checks all dataframes within the long format dict, and checks internal consistency of each row
    :param long_df_dict: dict containing all long data, separated by data type (keys)
    :param long_df: dataframe containing all long data (not long_df_dict)
    :param var_df: variant dataframe
    :param file: file name for logging
    :param logfile: logfile
    :return: determined file type, (bool) whether file is correctly formatted, number of inconsistencies, number of long equivalent rows, number of long inequivalent rows, logfile
    """
    # Set some variables
    is_mixed = False
    correct_format = False
    correct_count = []
    total_inconsistencies = 0
    inequivalent_cols = []
    timestr = time.strftime("%H:%M:%S")
    logfile_text = timestr + " ..... Checking Illumina LONG file"
    log_array = [logfile_text]
    # Run TFDP check on each dataframe within the long df dict
    for types in long_df_dict:
        if types == "AB":
            format_check, logfile = AB_check(
                long_df_dict[types], file, is_mixed, logfile
            )
            if format_check > 0:
                warning = "Columns for type " + types + " have unexpected values"
                print(warning)
                log_array.append(warning)
                correct_count.append(1)
                inequivalent_cols.append(types)
                total_inconsistencies = total_inconsistencies + format_check
            else:
                correct_count.append(0)
        else:
            format_check, logfile_name = TFDP_format_check(
                long_df_dict[types], var_df, types, file, is_mixed, logfile
            )
            if format_check > 0:
                warning = "Columns for type " + types + " have unexpected values"
                warnings.warn(warning, stacklevel=4)
                log_array.append(warning)
                correct_count.append(1)
                inequivalent_cols.append(types)
                total_inconsistencies = total_inconsistencies + format_check
            else:
                correct_count.append(0)
    # check internal consistency of each row
    long_equivalency, long_inequivalency, inequiv_dict = long_format_consistency_check(
        long_df, var_df
    )
    if long_inequivalency != 0:
        # Write inequivalencies to log file
        logfile = make_logs.long_inequivalency(file, logfile, inequiv_dict)
        # Warn inequivalencies
        warning = "One or more genotypes are inequivalent within a row"
        warnings.warn(warning, stacklevel=4)
        log_array.append(warning)
    if sum(correct_count) == 0:
        filetype = "LONG"
        correct_format = True
        if long_inequivalency == 0:
            print("File " + file + " is correctly formatted in Illumina LONG format")
        elif long_inequivalency != 0:
            print("File " + file + " may have incorrect genotype data")
        summary_inconsistency_value = 0
    else:
        warning = "Long format file contains errors"
        print(warning)
        log_array.append(warning)
        summary_inconsistency_value = total_inconsistencies
        filetype = None
    if logfile is not None:
        logging = make_logs.simple_log(log_array, file, logfile)
    return (
        filetype,
        correct_format,
        summary_inconsistency_value,
        long_equivalency,
        long_inequivalency,
        logfile,
    )


def check_illumina_matrix_format(
        illumina_df, var_df, file_type, file, is_mixed, header_row, header_dict, logfile
):
    """
    Puts Illumina matrix df through AB_check and TFDP
    :param illumina_df: illumina matrix dataframe
    :param var_df: variant dataframe
    :param file_type: user-specified input file type; if this is not mixed, the program assumes you know the file type
    :param file: filename for logging
    :param is_mixed: (bool) consider the file mixed for logging purposes
    :param header_row: number of rows to skip for header
    :param header_dict: dict containing illumina header info
    :param logfile: logfile
    :return: determined filetype, (bool) whether file is correctly formatted, number of inconsistencies, logfile
    """
    timestr = time.strftime("%H:%M:%S")
    logfile_text = timestr + " ..... Checking Illumina matrix file"
    log_array = [logfile_text]
    format_check_dict = {}
    filetype_out = None
    correct_format = False
    out_string = ""
    # If we "know" what the file type is to start with
    if file_type != "mixed":
        if file_type == "AB":
            format_check, logfile = AB_check(
                illumina_df, file, is_mixed, logfile
            )  # format check should == 0 (no non-AA/AB/BB/-- values in AB file)
            summary_inconsistency_value = format_check
            format_check_dict.update({"AB": format_check})
            if format_check != 0:
                out_string = "AB file is incorrectly formatted - see log for details"
            else:
                out_string = (
                        "File " + file + " is correctly formatted in Illumina AB format"
                )
                filetype_out = "AB"
                correct_format = True
        elif file_type == "TOP":
            format_check, logfile_name = TFDP_format_check(
                illumina_df, var_df, file_type, file, is_mixed, logfile
            )
            format_check_dict.update({"TOP": format_check})
            summary_inconsistency_value = format_check
            if format_check > 0:
                out_string = "TOP file is incorrectly formatted - see log for details"
            else:
                out_string = (
                        "File " + file + " is correctly formatted in Illumina TOP format"
                )
                filetype_out = "TOP"
                correct_format = True
        elif file_type == "FWD":
            format_check, logfile_name = TFDP_format_check(
                illumina_df, var_df, file_type, file, is_mixed, logfile
            )
            summary_inconsistency_value = format_check
            format_check_dict.update({"FWD": format_check})
            if format_check > 0:
                out_string = "FWD file is incorrectly formatted - see log for details"
            else:
                out_string = (
                        "File "
                        + file
                        + " is correctly formatted in Illumina FORWARD format"
                )
                filetype_out = "FWD"
                correct_format = True
        elif file_type == "DESIGN":
            format_check, logfile_name = TFDP_format_check(
                illumina_df, var_df, file_type, file, is_mixed, logfile
            )
            summary_inconsistency_value = format_check
            format_check_dict.update({"DESIGN": format_check})
            if format_check > 0:
                out_string = (
                    "DESIGN file is incorrectly formatted - see log for details"
                )
            else:
                out_string = (
                        "File " + file + " is correctly formatted in Illumina DESIGN format"
                )
                correct_format = True
                filetype_out = "DESIGN"
        elif file_type == "PLUS":
            format_check, logfile_name = TFDP_format_check(
                illumina_df, var_df, file_type, file, is_mixed, logfile
            )
            summary_inconsistency_value = format_check
            format_check_dict.update({"PLUS": format_check})
            if format_check > 0:
                out_string = "PLUS file is incorrectly formatted - see log for details"
            else:
                out_string = (
                        "File " + file + " is correctly formatted in Illumina PLUS format"
                )
                correct_format = True
                filetype_out = "PLUS"
        else:
            summary_inconsistency_value = None
        output_value = format_check_dict[file_type]
        print(out_string)
        log_array.append(out_string)
    # If file type IS mixed, then we have to try to guess
    else:
        # Warn if there are too few SNPs to make an assignment
        if len(illumina_df.index) < 50:
            out_string = (
                "There may not be enough SNPs to determine matrix format accurately"
            )
            warnings.warn(out_string, stacklevel=4)
            log_array.append(out_string)
        is_mixed = True
        # Test formats
        try_AB_format, logfile = AB_check(illumina_df, file, is_mixed, logfile)
        try_TOP_format, logfile = TFDP_format_check(
            illumina_df, var_df, "TOP", file, is_mixed, logfile
        )
        try_FWD_format, logfile = TFDP_format_check(
            illumina_df, var_df, "FWD", file, is_mixed, logfile
        )
        try_PLUS_format, logfile = TFDP_format_check(
            illumina_df, var_df, "PLUS", file, is_mixed, logfile
        )
        try_DESIGN_format, logfile = TFDP_format_check(
            illumina_df, var_df, "DESIGN", file, is_mixed, logfile
        )
        if try_AB_format == 0:
            message = "File " + file + " is correctly formatted in Illumina AB format"
            print(message)
            log_array.append(message)
            filetype_out = "AB"
            correct_format = True
            summary_inconsistency_value = try_AB_format

        elif try_TOP_format == 0:
            message = "File " + file + " is correctly formatted in Illumina TOP format"
            print(message)
            log_array.append(message)
            filetype_out = "TOP"
            correct_format = True
            summary_inconsistency_value = try_TOP_format
        elif try_FWD_format == 0:
            message = "File " + file + " is correctly formatted in Illumina FWD format"
            print(message)
            log_array.append(message)
            filetype_out = "FWD"
            correct_format = True
            summary_inconsistency_value = try_FWD_format
        elif try_PLUS_format == 0:
            message = "File " + file + " is correctly formatted in Illumina PLUS format"
            print(message)
            log_array.append(message)
            filetype_out = "PLUS"
            correct_format = True
            summary_inconsistency_value = try_PLUS_format
        elif try_DESIGN_format == 0:
            message = (
                    "File " + file + " is correctly formatted in Illumina DESIGN format"
            )
            print(message)
            log_array.append(message)
            filetype_out = "DESIGN"
            correct_format = True
            summary_inconsistency_value = try_DESIGN_format
        # Find best-fitting file format
        else:
            list_of_formats = [
                try_AB_format,
                try_FWD_format,
                try_TOP_format,
                try_PLUS_format,
                try_DESIGN_format,
            ]
            x = min(list_of_formats, key=float)
            # Get the minimum number of SNPs that have to be correct to determine format
            if header_row != 0:
                n_snps = int(header_dict["Num SNPs"])
            else:
                n_snps = len(list(illumina_df.index))
            min_snps = round(minimum_correct_snp_fraction * n_snps)
            if try_AB_format == x and try_AB_format < min_snps:
                message = (
                        "File "
                        + file
                        + " may be in AB format with "
                        + str(try_AB_format)
                        + " inconsistent SNP(s)"
                )
                print(message)
                log_array.append(message)
                filetype_out = "AB"
                correct_format = False
                summary_inconsistency_value = try_AB_format
            elif try_TOP_format == x and try_TOP_format < min_snps:
                message = (
                        "File "
                        + file
                        + " may be in TOP format with "
                        + str(try_TOP_format)
                        + " inconsistent SNP(s)"
                )
                print(message)
                log_array.append(message)
                filetype_out = "TOP"
                correct_format = False
                summary_inconsistency_value = try_TOP_format
            elif try_FWD_format == x and try_FWD_format < min_snps:
                message = (
                        "File "
                        + file
                        + " may be in FWD format with "
                        + str(try_FWD_format)
                        + " inconsistent SNP(s)"
                )
                print(message)
                log_array.append(message)
                filetype_out = "FWD"
                correct_format = False
                summary_inconsistency_value = try_FWD_format
            elif try_PLUS_format == x and try_PLUS_format < min_snps:
                message = (
                        "File "
                        + file
                        + " may be in PLUS format with "
                        + str(try_PLUS_format)
                        + " inconsistent SNP(s)"
                )
                print(message)
                log_array.append(message)
                filetype_out = "PLUS"
                correct_format = False
                summary_inconsistency_value = try_PLUS_format
            elif try_DESIGN_format == x and try_DESIGN_format < min_snps:
                message = (
                        "File "
                        + file
                        + " may be in DESIGN format with "
                        + str(try_PLUS_format)
                        + " inconsistent SNP(s)"
                )
                print(message)
                log_array.append(message)
                filetype_out = "DESIGN"
                correct_format = False
                summary_inconsistency_value = try_DESIGN_format
            else:
                message = (
                        "File type for "
                        + file
                        + " could not be determined: too many SNPs with inconsistent formatting"
                )
                print(message)
                filetype_out = None
                summary_inconsistency_value = None
                log_array.append(message)
    if logfile is not None:
        logging = make_logs.simple_log(log_array, file, logfile)
    return filetype_out, correct_format, summary_inconsistency_value, logfile


def snp_panel(var_ls, f):
    """
    Gets the SNP panel name and returns to user
    :param var_ls: list of variant files (should be 1)
    :param f: input file name
    :return: message text for log file (printed here)
    """
    # Get possible SNP panels
    error_text = []
    if len(var_ls) == 1:
        var = var_ls[0]
        message = "SNPs from file " + f + " are found in variant file " + var
        print(message)
        error_text.append(message)
    elif len(var_ls) > 1:
        var_string = ", ".join(var_ls)
        message = (
                "SNPs from file "
                + f
                + " are found in variant files "
                + var_string
        )
        print(message)
        error_text.append(message)
    else:
        # Program should quit before seeing this message
        message = "No conversion file contains all input SNPs"
        print(message)
        error_text.append(message)
    return error_text


#####################################################################################################################
def file_format_check(
        input_dir,
        file_list,
        specified_file_type,
        get_snp_panel,
        verbose_logging,
        conversion_dir,
        return_log,
        assembly,
        summarize,
        tabular,
        species,
        make_ped_map,
):
    """
    Main file format check function. Determines (or checks) the file format according to the variant position files
    :param input_dir: directory containing matrix files
    :param file_list: list of input files for processing (may be None)
    :param specified_file_type: input file type (default: mixed)
    :param get_snp_panel: (Bool) output snp panel for user
    :param verbose_logging: (Bool) write log to file
    :param conversion_dir: directory containing variant position file (default: variant_position_files)
    :param return_log: logfile name (will have timestamp information), may be supplied by another module
    :param assembly: genome assembly name
    :param summarize: (Bool) for use in conversion - summarize converted SNP file
    :param tabular: Type of summary file to create (tabular, pretty, basic)
    :param species: species name for variant position file
    :param make_ped_map: (Bool) make PLINK PED and MAP files
    :return: determined file type, correct format (Bool), output (empty) of writing log file
    """
    ##### Set some initial variables #####
    correct_format = False
    determined_ft = None
    log_input = None

    # Get all input files
    snp_files_array, logfile = get_snp_files(file_list, input_dir, return_log)

    # Get variant position files
    variant_files = get_variant_files(conversion_dir, assembly, species)

    # Check each file
    for file in snp_files_array:
        # Set some initial variables
        print("Checking file " + file)
        if verbose_logging is True:
            if return_log:
                timestr = time.strftime("%Y%m%d%H%M%S")
                log_suffix = "-" + timestr + ".log"
                log_input = make_logs.get_logname(log_suffix, file)
            else:
                log_input = None
        else:
            log_input = None
        # Deal with AFFY-PLUS as an input type
        if specified_file_type == "AFFY-PLUS":
            file_type = "PLUS"
        else:
            file_type = specified_file_type  # from here, only use file_type var
        # Get affy_flag
        affy_flag = False
        if file_type == "AFFY" or specified_file_type == "AFFY-PLUS":
            affy_flag = True
        file_path = os.path.join(input_dir, file)
        error_log_text = []
        summary_inconsistency_value = int
        correct_format = False  # set correct format as false initially
        long_df_dict = {}
        illumina_matrix_df = False

        # quick test for whether this still might be an AFFY file
        affy_flag = affy_test(file_path, file, file_type, affy_flag)

        # If affy_flag is True, parse AFFY file
        if affy_flag is True:
            affy_header_dict, affy_df, log_input = parse_affy_file(
                file_path, file, log_input
            )
            generic_input_df = affy_df
            long_ft = False
            illumina_df = None
            header_row = 0
            header_dict = None
        # If an Illumina format
        else:
            # Parse header of Illumina file
            header_row, header_dict = parse_header(file_path)
            numbers_exist = False
            if "Num Samples" in header_dict and "Num SNPs" in header_dict:
                numbers_exist = True
            # Read in Illumina file
            illumina_df, long_ft, log_input = read_illumina_input(
                file_path, header_row, file, log_input
            )
            # If long file, create dict containing all dataframes by type
            if file_type == "LONG" or long_ft is True:
                long_df_dict, log_input = get_long_df_dict(
                    illumina_df, header_row, header_dict, numbers_exist, file, log_input
                )
                generic_input_df = next(iter(long_df_dict.values()))
            # If not long file, use dataframe as is
            else:
                illumina_head_check, log_input = check_illumina_df_header(
                    illumina_df, header_row, header_dict, numbers_exist, file, log_input
                )
                generic_input_df = illumina_df
                illumina_matrix_df = True
        # Get corresponding variant file
        converting_file = False
        var_list, mod_verbose_log, reg_alt_bool_dict = var_match(
            variant_files,
            conversion_dir,
            generic_input_df,
            file,
            converting_file,
            assembly,
            species,
        )
        for text in mod_verbose_log:
            error_log_text.append(text)
        if log_input is not None:
            logging = make_logs.simple_log(error_log_text, file, logfile)
        var_df = get_var_df(
            conversion_dir, var_list, assembly, species, reg_alt_bool_dict
        )

        # If Affymetrix file, check in TFDP
        if affy_flag is True:
            (
                determined_ft,
                correct_format,
                summary_inconsistency_value,
                log_input,
            ) = check_affy_format(
                generic_input_df,
                var_df,
                specified_file_type,
                file,
                file_type,
                summary_inconsistency_value,
                log_input,
            )
            # Write summary file
            if summarize is True:
                equiv = 0
                inequiv = 0
                if specified_file_type == "AFFY-PLUS":
                    sum_type = "PLUS"
                else:  # specified type is "AFFY"
                    sum_type = "FWD"
                sum_file = write_summary(
                    generic_input_df,
                    sum_type,
                    file,
                    summary_inconsistency_value,
                    equiv,
                    inequiv,
                    tabular,
                )
            else:
                pass
        else:
            pass
        # If Illumina Long file, run TFDP on all dataframes
        if file_type == "LONG" or long_ft is True:
            (
                determined_ft,
                correct_format,
                summary_inconsistency_value,
                long_equivalency,
                long_inequivalency,
                log_input,
            ) = check_long_format(long_df_dict, illumina_df, var_df, file, log_input)
            if summarize is True:
                sum_file = write_summary(
                    illumina_df,
                    "LONG",
                    file,
                    summary_inconsistency_value,
                    long_equivalency,
                    long_inequivalency,
                    tabular,
                )
        else:
            pass
        # If Illumina matrix format (format can be known or mixed), run AB_check and TFDP on dataset
        if illumina_matrix_df is True:
            is_mixed = False
            (
                determined_ft,
                correct_format,
                summary_inconsistency_value,
                logfile,
            ) = check_illumina_matrix_format(
                illumina_df,
                var_df,
                file_type,
                file,
                is_mixed,
                header_row,
                header_dict,
                log_input,
            )
            # Write summary file
            if summarize is True:
                equiv = 0
                inequiv = 0
                sum_file = write_summary(
                    illumina_df,
                    file_type,
                    file,
                    summary_inconsistency_value,
                    equiv,
                    inequiv,
                    tabular,
                )
            else:
                pass
        else:
            pass
        # Get possible SNP panels if requested
        if get_snp_panel is True:
            snp_message = snp_panel(var_list, file)
            for item in snp_message:
                error_log_text.append(item)
        # Write plink files if required
        if correct_format is True:
            if make_ped_map is True:
                basename = os.path.splitext(file)
                ped_file = make_plink.create_ped_file(
                    basename[0], generic_input_df, error_log_text
                )
                map_file = make_plink.create_map_file(
                    basename[0], generic_input_df, var_df, species, error_log_text
                )

        # Write final log files
        if log_input is not None:
            logging = make_logs.simple_log(error_log_text, file, logfile)
    return determined_ft, correct_format, log_input


if __name__ == "__main__":
    # Get user input with argparse
    parser = argparse.ArgumentParser(description="Checks format of input files")
    parser.add_argument("input", type=str, help="directory containing input files")
    parser.add_argument(
        "--file-list",
        type=str,
        help="[Optional] comma-separated list of files in the input directory",
    )
    parser.add_argument(
        "--input-format",
        type=str,
        choices=["TOP", "FWD", "AB", "DESIGN", "mixed", "LONG", "PLUS", "AFFY"],
        help="Type of file(s) expected: 'TOP', 'FWD', 'AB', 'DESIGN', 'LONG', 'PLUS', or 'mixed' (Illumina) or 'AFFY' (Affymetrix)",
    )

    parser.add_argument(
        "--get-snp-panel",
        action="store_true",
        default=False,
        help="[Optional] Will determine which genotype conversion key files contain all SNPs in the input",
    )
    parser.add_argument(
        "-v",
        "--verbose-logging",
        action="store_true",
        default=False,
        required=False,
        help="[Optional] Write output to both STDOUT and log file",
    )
    parser.add_argument(
        "--conversion",
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
        default=False,
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
    log_file = None
    file_format_check(
        args.input,
        args.file_list,
        args.input_format,
        args.get_snp_panel,
        args.verbose_logging,
        args.conversion,
        log_file,
        args.assembly,
        args.summary,
        args.tabular,
        args.species,
        args.plink,
    )
