#!/usr/bin/python
import argparse
import os
import sys
import file_parsing as fp
import pandas as pd
import warnings
import numpy as np
import variant_file_finder as vff
import make_logs
import file_conversion as fc
import write_summary as summ
import errors
from pandas.errors import ParserError
import plink_output as make_plink


#####################################################################################################################
# This value controls the minimum fraction of correct SNPs required for the program to determine the file input type,
# if the input type is not given (or is 'mixed')
minimum_correct_snp_fraction = 0.95
#####################################################################################################################


def snp_panel(var_ls, f):
    # Get possible SNP panels
    error_text = []
    if len(var_ls) == 1:
        var = var_ls[0]
        message = "SNPs from file " + f + " are found in variant file " + var
        print(message)
        error_text.append(message)
    elif len(var_ls) > 1:
        var_string = ', '.join(var_ls)
        message = "SNPs from file " + f + " are found in variant files " + var_string + \
                  "\n" + var_ls[0] + " was used for analysis"
        print(message)
        error_text.append(message)
    else:
        # Program should quit before seeing this message
        message = "No conversion file contains all input SNPs"
        print(message)
        error_text.append(message)
    return error_text


def affy_head_check(header):
    simple_list = []
    header_correct = None
    for sample in header:
        if not sample.endswith('.1'):
            simple_list.append(sample)
    for simple in simple_list:
        new_simple = simple + ".1"
        if new_simple in header:
            header_correct = True
        else:
            header_correct = False
    return header_correct


def affy_test(file_path, file, file_type, affy_flag):
    with open(file_path) as f:
        line = f.readline()
        line = line.rstrip("\n")
        line_array = line.split("\t")
        if 'probeset_id' in line:
            if not affy_flag and not file_type == 'mixed':
                warnings.warn("File " + file + " appears to be in affymetrix format, but input format is " + file_type,
                              stacklevel=4)
            affy_flag = True
        elif (line == '[Header]' or line == '[Data]') and file_type == 'affymetrix':
            warnings.warn("File " + file + " appears to be in Illumina format, but input format is affymetrix",
                          stacklevel=4)
            affy_flag = False
        elif (len(line_array) >= 3) and line_array[1] == line_array[2]:  # This breaks when there is header info with no [Header] designation
            if not affy_flag and not file_type == 'mixed':
                warnings.warn("File " + file + " appears to be in affymetrix format, but input format is " + file_type,
                              stacklevel=4)
            affy_flag = True
        else:
            pass
        f.close()
    return(affy_flag)


def AB_check(dataframe, file, is_mixed, log):
    ab_list = ['AA', 'AB', 'BB', '--', 'NoCall']
    new_df = ~dataframe.iloc[:, 1:].isin(ab_list)
    result = sum(new_df.sum(axis=1))
    ab_warn = None
    if result > 0 and is_mixed is False:
        ab_warn = make_logs.ab_warning(dataframe, new_df, file, log)
    return result, ab_warn


def TFDP_format_check(dataframe, var_dataframe, fmt, in_file, is_mixed, log):
    # The A/B that are variable should be equal to TOP_A and TOP_B, respectively, if file_format is TOP
    # The A/B that are variable should be equal to FORWARD_A and FORWARD_B, respectively, if file_format is FWD
    # Create dataframe with ID, TOP_A, and TOP_B
    col_to_keep = []
    if fmt == 'TOP':
        a_column = 'TOP_A'
        b_column = 'TOP_B'
        keep = ['TOP_A', 'TOP_B']
    elif fmt == 'FWD':
        a_column = 'FORWARD_A'
        b_column = 'FORWARD_B'
        keep = ['FORWARD_A', 'FORWARD_B']
    elif fmt == 'PLUS':
        a_column = 'PLUS_A'
        b_column = 'PLUS_B'
        keep = ['PLUS_A', 'PLUS_B']
    elif fmt == 'DESIGN':
        a_column = 'DESIGN_A'
        b_column = 'DESIGN_B'
        keep = ['DESIGN_A', 'DESIGN_B']
    for i in keep:
        col_to_keep.append(i)
    # Keep only the database entries that are in the user input file
    if 'SNP Name' in list(dataframe.columns):
        dataframe.rename(columns={'SNP Name': 'Name'}, inplace=True)

    namelist = list(dataframe['Name'])
    subset_var_dataframe = var_dataframe.query('Name in @namelist').copy()
    #print(subset_var_dataframe)
    col_to_keep.append('Name')
    var_df_cols = list(var_dataframe)
    col_to_drop = list(set(var_df_cols) - set(col_to_keep))
    subset_var_dataframe.drop(col_to_drop, inplace=True, axis=1)
    #print(subset_var_dataframe)

    # Split each sample into a separate dataframe, and separate [NN] into two columns for comparison
    column_matches = []
    col0_names = list(dataframe.columns)
    col0_names.remove('Name')
    all_non_matches = 0
    log_exists = log
    for col0 in col0_names:
        name = dataframe[['Name', col0]].copy()
        name['{}'.format(a_column)], name['{}'.format(b_column)] = zip(*name[col0].apply(lambda x: list(x)))
        # Keep only rows that have differences (not AA/TT/CC/GG/--)
        nn_list = ['AA', 'TT', 'CC', 'GG', 'II', 'DD', '--', '---']
        name_differences_df = name[~name[col0].isin(nn_list)].copy()
        name_differences_df.drop(col0, axis=1, inplace=True)
        # Keep only homozygous rows (AA/TT/CC/GG)
        nn_list_2 = ['AA', 'TT', 'CC', 'GG']
        name_same_df = name[name[col0].isin(nn_list_2)].copy()
        name_same_df.drop(col0, axis=1, inplace=True)
        # Subset var dataframe again, based on having dropped the rows above (for heterozygous only)
        indiv_namelist = list(name_differences_df['Name'])
        indiv_subset_var_dataframe = subset_var_dataframe.query('Name in @indiv_namelist')

        # Compare dataframe with user input (name_differences_df) to database data (indiv_subset_var_dataframe)
        compare_df = pd.merge(indiv_subset_var_dataframe, name_differences_df, on=['Name', a_column, b_column],
                              how='left', indicator='Exist')
        compare_df['Exist'] = np.where(compare_df.Exist == 'both', 0, 1)
        # check if reversing the user input (AB -> BA) matches the information in the variant dataframe, and update the
        # Exist column accordingly
        user_df_reversed = name_differences_df.copy()
        a_column_new = b_column + "_new"
        b_column_new = a_column + "_new"
        user_df_reversed.rename(columns={a_column: a_column_new,
                                         b_column: b_column_new}, inplace=True)
        user_df_reversed.columns = user_df_reversed.columns.str.replace("_new", "")
        reverse_compare_df = pd.merge(indiv_subset_var_dataframe, user_df_reversed, on=['Name', a_column, b_column],
                              how='left', indicator='Exist')
        reverse_compare_df['Exist'] = np.where(reverse_compare_df.Exist == 'both', 0, 1)
        # Merge "Exist" columns, and if either is 0, then super_exist = 0
        super_compare_df = pd.merge(compare_df, reverse_compare_df, on=['Name', a_column, b_column, 'Exist'],
                              how='left', indicator='SuperExist')
        super_compare_df['SuperExist'] = np.where(super_compare_df.SuperExist == 'both', 1, 0)
        super_compare_df.drop(columns='Exist', inplace=True)
        # Write to log if there are inconsistent values
        if not super_compare_df[super_compare_df['SuperExist'] == 1].empty and is_mixed is False:
            # check that we still only have allowed values
            write_log = make_logs.inconsistent_values(super_compare_df[super_compare_df['SuperExist'] == 1], col0,
                                                      name_differences_df, in_file, log_exists)
            log_exists = write_log
        match = super_compare_df['SuperExist'].sum()
        column_matches.append(match)
        # Check that homozygous alleles are still allowed based on the var dataframe
        name_same_df.set_index("Name", inplace=True)
        svd_reindexed = subset_var_dataframe.set_index("Name")

        # merge svd_reindexed with name_same_df on name column and keep identical names
        # make a master merge that finds if there are any 'wrong' homozygous values not found in the var df at all
        merged_svd_namesame = pd.merge(left=name_same_df, right=svd_reindexed, how='left', on='Name')
        merged_svd_namesame['merge_A'] = np.where((merged_svd_namesame.iloc[:, 0] == merged_svd_namesame.iloc[:, 2]),
                                                  merged_svd_namesame.iloc[:, 2], 0)
        merged_svd_namesame['merge_B'] = np.where((merged_svd_namesame.iloc[:, 0] == merged_svd_namesame.iloc[:, 3]),
                                                  merged_svd_namesame.iloc[:, 3], 0)
        merged_svd_namesame['master_merge'] = np.where((merged_svd_namesame['merge_A'] ==
                                                        merged_svd_namesame['merge_B']), 1, 0)
        # make a new df that only has rows where master_merge == 1
        bad_merge_df = merged_svd_namesame[merged_svd_namesame.master_merge != 0]
        if bad_merge_df.empty is False:
            for idx, row in bad_merge_df.iterrows():
                column_matches.append(1)
                hom_allele = row[0]
                var_alleles = [row[2], row[3]]
                if is_mixed is False:
                    write_hom_log = make_logs.inconsistent_values_homozygous(in_file, log_exists, hom_allele,
                                                                             var_alleles, idx, col0)
                    log_exists = write_hom_log
                else:
                    pass
                acceptable_alleles = ['A', 'T', 'G', 'C', 'I', 'D', '-']
                if hom_allele not in acceptable_alleles:
                    exit("Unexpected allele " + hom_allele + " found at " + name_same_df[idx])

    all_non_matches = sum(column_matches)
    return all_non_matches, log_exists


def long_format_consistency_check(long_df, variant_df):
    # Consistency check: pick first 2 columns, figure out if they are AA/BB/AB, make sure that the other ones match this
    long_df.set_index("SNP Name", inplace=True)
    # reorder dataframe
    long_df_cols = list(long_df.columns.values)
    new_cols = ['Sample ID', 'Allele1 - Top', 'Allele2 - Top', 'Allele1 - Forward',
                'Allele2 - Forward', 'Allele1 - AB', 'Allele2 - AB', 'Allele1 - Design',
                'Allele2 - Design', 'Allele1 - Plus', 'Allele2 - Plus']
    extra_cols = list(set(long_df_cols) - set(new_cols))
    if 'GC Score' in long_df_cols:
        extra_cols.remove('GC Score')
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
    col_name_dict = {'Allele1 - Forward': 'FORWARD',
                     'Allele1 - Top': 'TOP',
                     'Allele1 - AB': 'AB',
                     'Allele1 - Design': 'DESIGN',
                     'Allele1 - Plus': 'PLUS'}
    sorted_long_df = long_df[new_cols_excluded]
    new_cols_no_sample_id = new_cols_excluded
    new_cols_no_sample_id.remove('Sample ID')
    new_cols_half = [x for x in new_cols_no_sample_id if "2" not in x]
    consistency_count = 0
    inconsistency_count = 0
    return_inequiv_dict = {}
    array_len = len(new_cols_excluded)
    array_list = list(range(1, array_len))
    for index, row in sorted_long_df.iterrows():
        row_array = [row[i] for i in array_list]
        var_row = variant_df.loc[variant_df['Name'] == index, :]
        long_df_top_A = row_array[0]
        long_df_top_B = row_array[1]
        var_df_top_A = var_row['TOP_A'].values[0]
        var_df_top_B = var_row['TOP_B'].values[0]
        #print(long_df_top_A, long_df_top_B, var_df_top_A, var_df_top_B)
        a_or_b = str
        top_empty = False
        if all(x == '-' for x in row_array):
            continue
        else:
            if long_df_top_A == long_df_top_B and long_df_top_A != '-':
                if long_df_top_A == var_df_top_A:
                    a_or_b = 'A'
                elif long_df_top_A == var_df_top_B:
                    a_or_b = 'B'
                else:
                    inconsistency_count = inconsistency_count + 1
                    return_inequiv_dict.update({index: row})
                    continue
            elif long_df_top_A != long_df_top_B:
                a_or_b = 'AB'
            elif long_df_top_A == '-' and long_df_top_B == '-':
                top_empty = True
        if top_empty is True:
            warnings.warn("No genotypes for TOP strand at " + index + ", " + row_array[0] +
                          "- this could indicate an error in the data")
        if a_or_b == 'A':
            a_val = 'A'
            b_val = 'A'
        elif a_or_b == 'B':
            a_val = 'B'
            b_val = 'B'
        elif a_or_b == 'AB':
            a_val = 'A'
            b_val = 'B'
        else:
            a_val = '-'
            b_val = '-'
        sub_var_row = []

        for col2 in new_cols_half:
            if col2 == 'Allele1 - AB':
                var_row_a_val = a_val
                var_row_b_val = b_val
            else:
                var_row_a_val = var_row['{}_{}'.format(col_name_dict[col2], a_val)].values[0]
                var_row_b_val = var_row['{}_{}'.format(col_name_dict[col2], b_val)].values[0]
            sub_var_row.append(var_row_a_val)
            sub_var_row.append(var_row_b_val)
        #sub_var_row = [var_row['{}_{}'.format(col_name_dict[new_cols[0]], a_val)].values[0], var_row['TOP_{}'.format(b_val)].values[0],
        #               var_row['FORWARD_{}'.format(a_val)].values[0], var_row['FORWARD_{}'.format(b_val)].values[0],
        #               a_val, b_val,
        #               var_row['DESIGN_{}'.format(a_val)].values[0], var_row['DESIGN_{}'.format(b_val)].values[0],
        #               var_row['PLUS_{}'.format(a_val)].values[0], var_row['PLUS_{}'.format(b_val)].values[0]]

        # If user input row has dashes in any columns, delete the corresponding columns from the variant row
        non_dash_index = []
        for element in range(0, len(row_array)):
            if row_array[element] != '-':
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

    ###################################################################################################################
    # Get files based on user input


def file_format_check(input_d, file_list, specified_file_type, get_snp_panel, verbose_logging, conversion_dir, return_log,
                      assembly, summarize, tabular, species, make_ped_map):
    input_dir = input_d
    if file_list is None:
        #  get all files in directory (warn and exclude files that are not in .txt format)
        snp_files_array = os.listdir(input_dir)
    else:
        # get all files in directory using file list
        snp_files = file_list
        snp_files_array = snp_files.split(',')
        # check that these files actually exist in the input directory
        input_not_exists = errors.input_file_name_error(input_d, snp_files_array)
        if input_not_exists:
            exit("File(s) " + ", ".join(input_not_exists) + " are not found in " + input_d)
        else:
            pass
    # Check that these are all .txt files
    # check to exclude any non .txt files
    non_txt_files = errors.non_txt_fmt_input_files(snp_files_array)
    if non_txt_files:
        warnings.warn(
            "The following files are not .txt files are are excluded from analysis: " + ', '.join(non_txt_files),
            stacklevel=4)
        for non_txt in non_txt_files:
            snp_files_array.remove(non_txt)
    else:
        pass
    # check that there are files left in snp_files_array
    if snp_files_array:
        pass
    else:
        exit("No matrix files in " + input_d)

    # Deal with Affy-Plus as an input type
    if specified_file_type == 'AFFY-PLUS':
        file_type = 'PLUS'
    else:
        file_type = specified_file_type



    # Get variant position files
    variant_species_dir = os.path.join(conversion_dir, species)
    variant_assembly_dir = os.path.join(variant_species_dir, assembly)
    variant_files = os.listdir(variant_assembly_dir)
    # Check that the assembly and species combination works
    species_error = errors.assembly_species_error(conversion_dir, assembly, species)
    if not species_error:
        exit("Cannot find assembly name in any variant file in " + conversion_dir)
    # Exclude any files without .csv ending (but with correct assembly)
    variant_exclude = errors.non_csv_fmt_conversion_files(conversion_dir, assembly, species)
    if variant_exclude:
        warnings.warn(
            "The following variant files are not .csv files are are excluded from analysis: " + ', '.join(
                variant_exclude),
            stacklevel=3)
        for wrong_var in variant_exclude:
            variant_files.remove(wrong_var)
        if not variant_files:
            exit("No variant conversion files in " + variant_assembly_dir)
    else:
        pass

    # Get affy_flag
    affy_flag = False
    if file_type == 'affymetrix' or specified_file_type == 'AFFY-PLUS':
        affy_flag = True

    # Read in files and parse header info
    correct_format = False
    for file in snp_files_array:
        print("Checking file " + file)
        file_path = os.path.join(input_dir, file)
        error_log_text = []
        summary_inconsistency_value = int
        inequivalent_cols = []
        # quick test for whether this still might be an affymetrix file
        affy_flag = affy_test(file_path, file, file_type, affy_flag)
        if affy_flag is True:
            # Test affy forward values here
            affy_df = pd.read_csv(file_path, sep='\t', mangle_dupe_cols=True)
            # make dataframe look like the illumina one (remove AB columns)
            header_row = list(set(affy_df.columns))
            header_row.remove('probeset_id')
            # quick check to make sure header row is really affy format (2 columns for every sample)
            head_check = affy_head_check(header_row)
            if head_check is False:
                message = "Affymetrix file may not be properly formatted: check columns"
                warnings.warn(message, stacklevel=4)
                error_log_text.append(message)
            else:
                pass
            affy_header_dict = {}
            for value in header_row:
                new_affy_df = affy_df[['probeset_id', value]].copy()
                ab_list = ['AA', 'AB', 'BB', 'NoCall']
                test_ab_df = ~new_affy_df.iloc[:, 1:].isin(ab_list)
                not_ab = sum(test_ab_df.sum(axis=1))
                affy_header_dict.update({value: not_ab})
            copy_affy_df = affy_df.copy()  # keep a copy just in case
            for keys in affy_header_dict:
                if affy_header_dict[keys] == 0:
                    affy_df.drop(columns=keys, inplace=True)

            # Get variant file using variant_file_finder
            ##
            converting_file = False
            affy_df.rename(columns={'probeset_id': 'Name'}, inplace=True)
            var_list, mod_verbose_log, alt_bool = vff.var_match(variant_files, conversion_dir, affy_df, file, converting_file,
                                                      assembly, species)
            for text in mod_verbose_log:
                error_log_text.append(text)
            var_df = vff.get_var_df(conversion_dir, var_list[0], assembly, species, alt_bool)

            # Pass the dataframe to TFDP_format_check
            if file_type == 'mixed':
                is_mixed = True
            else:
                is_mixed = False
            if specified_file_type == 'AFFY-PLUS':
                affy_fmt = 'PLUS'
            elif specified_file_type == 'affymetrix':
                affy_fmt = 'FWD'
            else:
                affy_fmt = 'FWD'  # most likely mixed format (not from the check post-file-conversion)
            format_check, logfile_name = TFDP_format_check(affy_df, var_df, affy_fmt, file, is_mixed, return_log)
            summary_inconsistency_value = format_check
            if format_check > 0:
                if is_mixed is False:
                    out_string = "Affymetrix file is incorrectly formatted - see log for details."
                else:
                    out_string = "Affymetrix file is incorrectly formatted - re-run with --input-format affymetrix for details"
                warnings.warn(out_string, stacklevel=4)
                filetype = None
            else:  # format check is correct
                if head_check is True:
                    out_string = "File " + file + " is correctly formatted"
                    print(out_string)
                    error_log_text.append(out_string)
                else:
                    pass
                filetype = 'affymetrix'
                correct_format = True
                if is_mixed is True:
                    if head_check is True:
                        out_string = "File " + file + " is in Affymetrix format"

                    else:
                        out_string = "File " + file + " appears to be in Affymetrix format, but is missing expected columns"
                    print(out_string)
                    error_log_text.append(out_string)
            # Get possible SNP panels
            if get_snp_panel is True:
                snp_message = snp_panel(var_list, file)
                for item in snp_message:
                    error_log_text.append(item)
            if verbose_logging is True:
                if logfile_name is not None:
                    log_name = logfile_name
                else:
                    if return_log is None:
                        log_name = None
                    else:
                        log_name = return_log
                v_log = make_logs.simple_log(error_log_text, file, log_name)
            else:
                v_log = None
            # Write summary file
            if summarize is True:
                equiv = 0
                inequiv = 0
                sum_file = summ.write_summary(affy_df, 'FWD', file, summary_inconsistency_value, equiv, inequiv, tabular)
            else:
                pass

        else:  # Not affymetrix -> some sort of Illumina format
            file_path = os.path.join(input_dir, file)
            header_row, header_dict = fp.parse_header(file_path)
            if 'Num Samples' in header_dict and 'Num SNPs' in header_dict:
                numbers_exist = True
            else:
                numbers_exist = False
            # Make new dataframes
            with open(file_path, 'r') as input_file:
                # Make sure input file will not throw a dataframe error (from an incorrect header structure)
                try:
                    df = pd.read_csv(input_file, skiprows=header_row, sep="\t")
                except ParserError:
                    exit("First line might not contain column names - check formatting")

                if file_type != 'LONG':
                    # Put in a check here to make sure there are the same number of column labels as columns
                    if 'Unnamed: 0' not in df.columns:
                        df.index.names = ['Unnamed: 0']
                        df.reset_index(inplace=True)
                    df.rename(columns={'Unnamed: 0': 'Name'}, inplace=True)

                #
                # Check for long and if so, create df to get the right var file
                df_columns_long_check = list(df.columns)
                long_ft = False
                for type_of_file in df_columns_long_check:
                    if 'Allele1' in type_of_file:
                        long_ft = True
                if file_type == 'LONG' or long_ft is True:
                    # make sure header info accurately reflects file
                    if header_row != 0 and numbers_exist:
                        samples_list = list(set(df['Sample ID']))
                        num_samples = len(samples_list)
                        num_snps = len(df['SNP Name'].unique())
                        header_congruency = fp.check_header_data_congruency(header_dict, num_samples, num_snps)
                        if header_congruency:
                            append_text = "File " + file
                            error_log_text.append(append_text)
                        for each_warning in header_congruency:
                            error_log_text.append(each_warning)
                    else:
                        pass
                    # do stuff here to get the long file into arrays
                    type_list = ['TOP', 'FWD', 'AB', 'DESIGN', 'PLUS']
                    df_dict = {}
                    for ty in type_list:
                        out1, out2 = fc.gen_long_output_col_names(ty)
                        # Make sure allele type exists
                        if out1 not in df.columns or out2 not in df.columns:
                            message = "Columns " + out1 + " and " + out2 + " might not be in the input file" ##### check long file equivalencies here
                            warnings.warn(message, stacklevel=4)
                        else:
                            sub_long_df = df[['SNP Name', 'Sample ID', out1, out2]]
                            # print new df by sample and save to dict by "file type"
                            sample_list = sub_long_df['Sample ID']
                            unique_samples = list(set(sample_list))
                            pos1 = sub_long_df.columns.values[2]
                            pos2 = sub_long_df.columns.values[3]
                            output_df = sub_long_df[['SNP Name']].copy()
                            output_df = output_df.rename(columns={'SNP Name': 'Name'})
                            for sample in unique_samples:
                                sub_sample_df = pd.DataFrame()
                                sub_sample_df[sample] = sub_long_df[[pos1, pos2]].apply(lambda x: ''.join(x), axis=1)
                                output_df[sample] = sub_sample_df[sample].apply(lambda v: v)
                                df_dict.update({ty: output_df})
                    # get var df
                    converting_file = False
                    df_var = next(iter(df_dict.values()))
                    var_list, mod_verbose_log, alt_bool = vff.var_match(variant_files, conversion_dir, df_var, file, converting_file, assembly, species)
                    for text in mod_verbose_log:
                        error_log_text.append(text)
                    var_df = vff.get_var_df(conversion_dir, var_list[0], assembly, species, alt_bool)
                    # run file format check on each dataframe
                    is_mixed = False
                    correct_count = []
                    ab_log_name = None
                    logfile_name = None
                    total_inconsistencies = 0
                    for types in df_dict:
                        if types == 'AB':
                            if logfile_name is not None:
                                return_log = logfile_name
                            format_check, ab_log_name = AB_check(df_dict[types], file, is_mixed, return_log)
                            if format_check > 0:
                                warning = 'Columns for type ' + types + ' have unexpected values'
                                print(warning)
                                error_log_text.append(warning)
                                correct_count.append(1)
                                inequivalent_cols.append(types)
                                total_inconsistencies = total_inconsistencies + format_check
                            else:
                                correct_count.append(0)
                        else:
                            if ab_log_name is not None:
                                return_log = ab_log_name
                            format_check, logfile_name = TFDP_format_check(df_dict[types], var_df, types, file, is_mixed, return_log)
                            if format_check > 0:
                                warning = 'Columns for type ' + types + ' have unexpected values'
                                warnings.warn(warning, stacklevel=4)
                                error_log_text.append(warning)
                                correct_count.append(1)
                                inequivalent_cols.append(types)
                                total_inconsistencies = total_inconsistencies + format_check
                            else:
                                correct_count.append(0)
                    # check internal consistency of each row
                    long_equivalency, long_inequivalency, inequiv_dict = long_format_consistency_check(df, var_df)
                    if long_inequivalency != 0:
                        # Write inequivalencies to log file
                        if verbose_logging is True:
                            if logfile_name is not None:
                                log_name = logfile_name
                            else:
                                if return_log is None:
                                    log_name = None
                                else:
                                    log_name = return_log
                            v_log = make_logs.long_inequivalency(file, log_name, inequiv_dict)
                            logfile_name = v_log
                        else:
                            v_log = None
                        # Warn inequivalencies
                        warning = "One or more genotypes are inequivalent within a row"
                        warnings.warn(warning, stacklevel=4)
                        error_log_text.append(warning)

                        # Get possible SNP panels
                    if get_snp_panel is True:
                        snp_message = snp_panel(var_list, file)
                        for item in snp_message:
                            error_log_text.append(item)

                    if sum(correct_count) == 0:
                        filetype = 'LONG'
                        correct_format = True
                        if long_inequivalency == 0:
                            if specified_file_type == 'mixed':
                                print("File " + file + " is in long format")
                            else:
                                print("File " + file + " is correctly formatted")
                        elif long_inequivalency != 0:
                            print("File " + file + " may have incorrect genotype data")

                        summary_inconsistency_value = 0
                    else:
                        warning = "Long format file contains errors"
                        print(warning)
                        error_log_text.append(warning)
                        summary_inconsistency_value = total_inconsistencies

                        filetype = None
                    if verbose_logging is True:
                        if ab_log_name is not None and logfile_name is not None:
                            if ab_log_name == logfile_name:
                                log_name = logfile_name
                        else:
                            if ab_log_name is not None:
                                log_name = ab_log_name
                            elif logfile_name is not None:
                                log_name = logfile_name
                            else:
                                if return_log is None:
                                    log_name = None
                                else:
                                    log_name = return_log
                        v_log = make_logs.simple_log(error_log_text, file, log_name)
                    else:
                        v_log = None
                    if summarize is True:
                        sum_file = summ.write_summary(df, 'LONG', file, summary_inconsistency_value, long_equivalency,
                                                      long_inequivalency, tabular)
                    else:
                        pass

                else:  # Not a LONG format file, but some Illumina matrix format
                    # make sure header info accurately reflects file
                    if header_row != 0 and numbers_exist:
                        num_samples = len(df.columns) - 1
                        num_snps = len(df.index)
                        header_congruency = fp.check_header_data_congruency(header_dict, num_samples, num_snps)
                        if header_congruency:
                            append_text = "File " + file
                            error_log_text.append(append_text)
                        for each_warning in header_congruency:
                            error_log_text.append(each_warning)
                    # Get the right var file
                    converting_file = False
                    var_list, mod_verbose_log, alt_bool = vff.var_match(variant_files, conversion_dir, df, file, converting_file, assembly, species)
                    for text in mod_verbose_log:
                        error_log_text.append(text)
                    var_df = vff.get_var_df(conversion_dir, var_list[0], assembly, species, alt_bool)
                    # Check AB file format
                    is_mixed = False
                    ab_log_name = None
                    logfile_name = None
                    if file_type == "AB":
                        format_check, ab_log_name = AB_check(df, file, is_mixed, return_log)  # format check should == 0 (no non-AA/AB/BB/-- values in AB file)
                        summary_inconsistency_value = format_check
                        if format_check != 0:
                            warning = "AB file is incorrectly formatted - see log for details"
                            warnings.warn(warning, stacklevel=4)
                            filetype = None

                        else:
                            out_string = "File " + file + " is correctly formatted"
                            print(out_string)
                            error_log_text.append(out_string)
                            filetype = 'AB'
                            correct_format = True
                    elif file_type == "TOP":
                        format_check, logfile_name = TFDP_format_check(df, var_df, file_type, file, is_mixed, return_log)
                        summary_inconsistency_value = format_check
                        if format_check > 0:
                            warning = "TOP file is incorrectly formatted - see log for details"
                            warnings.warn(warning, stacklevel=4)
                            filetype = None

                        else:
                            out_string = "File " + file + " is correctly formatted"
                            print(out_string)
                            error_log_text.append(out_string)
                            filetype = 'TOP'
                            correct_format = True
                    elif file_type == 'FWD':
                        format_check, logfile_name = TFDP_format_check(df, var_df, file_type, file, is_mixed, return_log)
                        summary_inconsistency_value = format_check
                        if format_check > 0:
                            out_string = "FWD file is incorrectly formatted - see log for details"
                            warnings.warn(out_string, stacklevel=4)
                            filetype = None

                        else:
                            out_string = "File " + file + " is correctly formatted"
                            print(out_string)
                            error_log_text.append(out_string)
                            filetype = 'FWD'
                            correct_format = True
                    elif file_type == 'DESIGN':
                        format_check, logfile_name = TFDP_format_check(df, var_df, file_type, file, is_mixed, return_log)
                        summary_inconsistency_value = format_check
                        if format_check > 0:
                            out_string = "DESIGN file is incorrectly formatted - see log for details"
                            warnings.warn(out_string, stacklevel=4)
                            filetype = None

                        else:
                            out_string = "File " + file + " is correctly formatted"
                            print(out_string)
                            error_log_text.append(out_string)
                            correct_format = True
                            filetype = 'DESIGN'
                    elif file_type == 'PLUS':
                        format_check, logfile_name = TFDP_format_check(df, var_df, file_type, file, is_mixed, return_log)
                        summary_inconsistency_value = format_check
                        if format_check > 0:
                            out_string = "PLUS file is incorrectly formatted - see log for details"
                            warnings.warn(out_string, stacklevel=4)
                            filetype = None

                        else:
                            out_string = "File " + file + " is correctly formatted"
                            print(out_string)
                            error_log_text.append(out_string)
                            correct_format = True
                            filetype = 'PLUS'

                    else:  # mix of format types or unknown format type
                        # Warn if there are too few SNPs to make an assignment
                        if len(df.index) < 50:
                            out_string = "There may not be enough SNPs to determine matrix format accurately"
                            warnings.warn(out_string, stacklevel=4)
                            error_log_text.append(out_string)

                        is_mixed = True
                        # Test formats
                        try_AB_format, AB_log_name = AB_check(df, file, is_mixed, return_log)
                        try_TOP_format, logfile_name = TFDP_format_check(df, var_df, "TOP", file, is_mixed, return_log)
                        try_FWD_format, logfile_name = TFDP_format_check(df, var_df, "FWD", file, is_mixed, return_log)
                        try_PLUS_format, logfile_name = TFDP_format_check(df, var_df, "PLUS", file, is_mixed, return_log)
                        try_DESIGN_format, logfile_name = TFDP_format_check(df, var_df, 'DESIGN', file, is_mixed, return_log)
                        if try_AB_format == 0:
                            message = "File " + file + " is in AB format"
                            print(message)
                            error_log_text.append(message)
                            filetype = 'AB'
                            correct_format = True
                            summary_inconsistency_value = try_AB_format

                        elif try_TOP_format == 0:
                            message = "File " + file + " is in TOP format"
                            print(message)
                            error_log_text.append(message)
                            filetype = 'TOP'
                            correct_format = True
                            summary_inconsistency_value = try_TOP_format
                        elif try_FWD_format == 0:
                            message = "File " + file + " is in FWD format"
                            print(message)
                            error_log_text.append(message)
                            filetype = 'FWD'
                            correct_format = True
                            summary_inconsistency_value = try_FWD_format
                        elif try_PLUS_format == 0:
                            message = "File " + file + " is in PLUS format"
                            print(message)
                            error_log_text.append(message)
                            filetype = 'PLUS'
                            correct_format = True
                            summary_inconsistency_value = try_PLUS_format
                        elif try_DESIGN_format == 0:
                            message = "File " + file + " is in DESIGN format"
                            print(message)
                            error_log_text.append(message)
                            filetype = 'DESIGN'
                            correct_format = True
                            summary_inconsistency_value = try_DESIGN_format

                        else:
                            list_of_formats = [try_AB_format, try_FWD_format, try_TOP_format, try_PLUS_format, try_DESIGN_format]
                            x = min(list_of_formats, key=float)
                            # Get the minimum number of SNPs that have to be correct to determine format
                            if header_row != 0:
                                n_snps = int(header_dict['Num SNPs'])
                            else:
                                n_snps = len(list(df.index))
                            min_snps = round(minimum_correct_snp_fraction * n_snps)
                            if try_AB_format == x and try_AB_format < min_snps:
                                message = "File " + file + " may be in AB format with " + str(try_AB_format) + \
                                          " inconsistent SNP(s)"
                                print(message)
                                error_log_text.append(message)
                                filetype = 'AB'
                                correct_format = False
                                summary_inconsistency_value = try_AB_format
                            elif try_TOP_format == x and try_TOP_format < min_snps:
                                message = "File " + file + " may be in TOP format with " + str(try_TOP_format) + \
                                          " inconsistent SNP(s)"
                                print(message)
                                error_log_text.append(message)
                                filetype = 'TOP'
                                correct_format = False
                                summary_inconsistency_value = try_TOP_format
                            elif try_FWD_format == x and try_FWD_format < min_snps:
                                message = "File " + file + " may be in FWD format with " + str(try_FWD_format) + \
                                          " inconsistent SNP(s)"
                                print(message)
                                error_log_text.append(message)
                                filetype = 'FWD'
                                correct_format = False
                                summary_inconsistency_value = try_FWD_format
                            elif try_PLUS_format == x and try_PLUS_format < min_snps:
                                message = "File " + file + " may be in PLUS format with " + str(try_PLUS_format) + \
                                          " inconsistent SNP(s)"
                                print(message)
                                error_log_text.append(message)
                                filetype = 'PLUS'
                                correct_format = False
                                summary_inconsistency_value = try_PLUS_format
                            elif try_DESIGN_format == x and try_DESIGN_format < min_snps:
                                message = "File " + file + " may be in DESIGN format with " + str(try_PLUS_format) + \
                                          " inconsistent SNP(s)"
                                print(message)
                                error_log_text.append(message)
                                filetype = 'DESIGN'
                                correct_format = False
                                summary_inconsistency_value = try_DESIGN_format
                            else:
                                message = "File type for " + file + \
                                          " could not be determined: too many SNPs with unclear formatting"
                                print(message)
                                file_type = None
                                error_log_text.append(message)
                    # Write PED and MAP files
                    if correct_format is True:
                        if make_ped_map is True:
                            basename = os.path.splitext(file)
                            ped_file = make_plink.create_ped_file(basename[0], df, error_log_text)
                            map_file = make_plink.create_map_file(basename[0], df, var_df, species, error_log_text)
                    # Get possible SNP panels
                    if get_snp_panel is True:
                        snp_message = snp_panel(var_list, file)
                        for item in snp_message:
                            error_log_text.append(item)
                    if verbose_logging is True:
                        if ab_log_name is not None:
                            log_name = ab_log_name
                        elif logfile_name is not None:
                            log_name = logfile_name
                        else:
                            if return_log is None:
                                log_name = None
                            else:
                                log_name = return_log
                        v_log = make_logs.simple_log(error_log_text, file, log_name)
                    else:
                        v_log = None
                    # Write summary file
                    if summarize is True:
                        equiv = 0
                        inequiv = 0
                        sum_file = summ.write_summary(df, file_type, file, summary_inconsistency_value, equiv, inequiv, tabular)
                    else:
                        pass

    return filetype, correct_format, v_log


if __name__ == "__main__":
    # Get user input with argparse
    parser = argparse.ArgumentParser(description="Checks format of input files")
    parser.add_argument(
        "input",
        type=str,
        help="directory containing input files"
    )
    parser.add_argument(
        "--file-list",
        type=str,
        help="[optional] comma-separated list of files in the input directory"
    )
    parser.add_argument(
        "--input-format",
        type=str,
        choices=['TOP', 'FWD', 'AB', 'DESIGN', 'mixed', 'LONG', 'PLUS',  'affymetrix'],
        help="Type of file(s) expected: 'TOP', 'FWD', 'AB', 'DESIGN', 'LONG', 'PLUS', or 'mixed' (Illumina) or 'affymetrix' (Affymetrix)"
    )

    parser.add_argument(
        "--get-snp-panel",
        action="store_true",
        default=False,
        help="[optional] Will determine which genotype conversion key files contain all SNPs in the input"
    )
    parser.add_argument(
        "-v",
        "--verbose-logging",
        action="store_true",
        default=False,
        required=False,
        help="[optional] Write output to both STDOUT and log file"

    )
    parser.add_argument(
        '--key-dir',
        type=str,
        default="variant_position_files",
        help="Directory containing genotype conversion key files (default = variant_position_files)"
    )
    parser.add_argument(
        '--assembly',
        type=str,
        required=True,
        help="Assembly name which must be included in genotype conversion file name"
    )
    parser.add_argument(
        '--species',
        required=True,
        type=str,
        choices=["bos_taurus", "sus_scrofa"],
        help="Organism name"
    )
    parser.add_argument(
        '-s',
        '--summary',
        action="store_true",
        default=False,
        required=False,
        help="Summarize converted SNP file in *_summary.txt file"
    )
    parser.add_argument(
        '--tabular',
        action="store_true",
        default=False,
        required=False,
        help="Output summary file in tabular format (default: False)"
    )
    parser.add_argument(
        '--plink',
        action="store_true",
        default=False,
        required=False,
        help="Creates PLINK flat files (PED and MAP) (default: False)"
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    log_file = None
    file_format_check(args.input, args.file_list, args.input_format, args.get_snp_panel, args.verbose_logging, args.key_dir, log_file, args.assembly, args.summary, args.tabular, args.species, args.plink)



