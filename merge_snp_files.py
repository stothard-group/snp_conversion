#!/usr/bin/python
import os
import sys
import pandas as pd
import file_parsing as fp
import functools
import file_conversion as conversion
import argparse
from functools import reduce
import numpy as np


# This module can merge multiple files of the same format. If merging files of different formats, convert them to the
# same format first.


def update_header_dict(entire_h, df, in_type):
    gsgt_version = []
    processing_date = []
    content = []
    for keys in entire_h:
        each_dict = entire_h[keys]
        # Use dict.get to make sure that the value exists, or replace with an empty string
        gsgt_version.append(each_dict.get('GSGT Version', ''))
        processing_date.append(each_dict.get('Processing Date', ''))
        content.append(each_dict.get('Content', ''))
    gsgt_str = ", ".join(filter(None, list(set(gsgt_version))))  # this line makes a list of the entries in
    # gsgt_version, removes empty values, and then joins remaining values with a comma
    proc_date_str = ", ".join(filter(None, list(set(processing_date))))
    content_str = ", ".join(filter(None, list(set(content))))
    new_header_info = ({'GSGT Version': gsgt_str, 'Processing Date': proc_date_str, 'Content': content_str})
    # Count long vs matrix files differently
    if in_type != "LONG":
        num_snps = len(df.index)
        num_samples = len(df.columns) - 1
    else:
        snps = df['SNP Name']
        num_snps = len(list(set(snps)))
        samples = df['Sample ID']
        num_samples = len(list(set(samples)))
    new_header_info.update({"Num SNPs": str(num_snps), "Total SNPs": str(num_snps), "Num Samples": str(num_samples),
                            "Total Samples": str(num_samples)})
    return new_header_info


def files_to_dfs(snp_array, in_dir):
    f_dict = {}
    entire_head_dict = {}
    for f in snp_array:
        file_path = os.path.join(in_dir, f)
        header_row, header_dict = fp.parse_header(file_path)
        entire_head_dict.update({f: header_dict})
        with open(file_path, 'r') as input_file:
            df = pd.read_csv(input_file, skiprows=header_row, sep="\t")
            df.rename(columns={'Unnamed: 0': 'Name'}, inplace=True)
            f_dict.update({f: df})
        input_file.close()
    return f_dict, entire_head_dict


def merge_matrix_files(input_dir, input_files, input_type, outfile_name):
    # Remove output file if it already exists
    if os.path.exists(outfile_name):
        os.remove(outfile_name)
    log_file = None
    # Get list of files to merge
    if input_files is None:
        #  get all files in directory
        snp_files_array = os.listdir(input_dir)
    else:
        # get all files in directory using file list
        snp_files_array = input_files.split(',')
    file_dict = {}
    file_count = 0

    # IF files are in matrix format
    if input_type != "LONG":
        # Read files into dataframes
        file_dict, entire_header_dict = files_to_dfs(snp_files_array, input_dir)
        # Start merging dataframes
        all_dfs = []
        for key in file_dict:
            all_dfs.append(file_dict[key])
        new_df = functools.reduce(lambda df1, df2: df1.merge(df2, "outer"), all_dfs)
        # Deal with NaNs
        new_df.fillna('--', inplace=True)
        # Update header information and write to new file
        header_out = update_header_dict(entire_header_dict, new_df, input_type)
        header_file = conversion.write_header(outfile_name, header_out)
        # Write dataframes to outfile
        col_list = list(new_df)
        new_df.to_csv(outfile_name, index=None, mode='a', header=col_list, sep='\t')
    # IF files are in LONG format
    else:
        long_file_dict, long_entire_head_dict = files_to_dfs(snp_files_array, input_dir)
        all_dfs = []
        for key in long_file_dict:
            all_dfs.append(long_file_dict[key])
        new_long_df = pd.concat(all_dfs, ignore_index=True)
        header_out = update_header_dict(long_entire_head_dict, new_long_df, input_type)
        header_files = conversion.write_header(outfile_name, header_out)
        # Write new DF to outfile
        col_list = list(new_long_df)
        new_long_df.to_csv(outfile_name, index=None, mode='a', header=col_list, sep='\t')
    return True


def merge_affy_files(input_dir, input_files, input_type, outfile_name):
    # Remove output file if it already exists
    if os.path.exists(outfile_name):
        os.remove(outfile_name)
    log_file = None
    # Get list of files to merge
    if input_files is None:
        #  get all files in directory
        snp_files_array = os.listdir(input_dir)
    else:
        # get all files in directory using file list
        snp_files_array = input_files.split(',')
    df_list = []
    if input_type == 'affymetrix':
        for affy_file in snp_files_array:
            affy_file_path = os.path.join(input_dir, affy_file)
            affy_df = pd.read_csv(affy_file_path, sep='\t', mangle_dupe_cols=True)
            df_list.append(affy_df)
        merged_df = reduce(lambda left, right: pd.merge(left, right, on=['probeset_id'], how='outer'), df_list)
        # Detect issues with duplicated columns
        all_columns = list(merged_df.columns)
        cols_without_dups = []
        for test_col in all_columns:
            if "_x" in test_col or "_y" in test_col:
                pass
            else:
                cols_without_dups.append(test_col)
        merged_final_df = merged_df[cols_without_dups]
        if any('_x' in n for n in all_columns) and any('_y' in m for m in all_columns):
            x_list = []
            y_list = []
            for col_name in merged_df.columns:
                if "_x" in col_name:
                    x_list.append(col_name)
                elif "_y" in col_name:
                    y_list.append(col_name)
                else:
                    pass
            # For each duplicated column, combine them so all new values are added to a single column
            for each_x in x_list:
                xy_df = pd.DataFrame()
                basename = each_x.replace('_x', '')
                matching_y = basename + '_y'
                xy_df['probeset_id'] = merged_df['probeset_id'].copy()
                xy_df[each_x] = merged_df[each_x].copy()
                xy_df[matching_y] = merged_df[matching_y].copy()
                xy_df.reset_index().groupby('index').max()
                # double check that the columns are truly identical
                compare_df = xy_df.merge(xy_df, on=['probeset_id', each_x, matching_y], how='left', indicator='Exist')
                compare_df['Exist'] = np.where(compare_df.Exist == 'both', 0, 1)
                if compare_df[compare_df['Exist'] == 1].empty:
                    pass
                else:
                    exit("There are duplicate columns in the files which cannot be properly merged")
                # add one of them back to the dataframe and delete the x and y versions
                # get index positions of x and y versions
                x_pos, y_pos = [merged_df.columns.get_loc(each_x), merged_df.columns.get_loc(matching_y)]
                xy_df.rename(columns={each_x: basename}, inplace=True)
                xy_df.drop(matching_y, axis=1, inplace=True)
                merged_final_df = merged_final_df.merge(xy_df, on='probeset_id', how='outer')
                final_col_list = list(merged_final_df.columns)
                final_col_list.remove(basename)
                final_col_list.insert(x_pos, basename)
                merged_final_df = merged_final_df.reindex(columns=final_col_list)
        else:
            pass  # no duplicated columns to merge
        # Deal with NaNs
        col_list = []
        for col in merged_final_df.columns:
            if ".1" not in col:
                col_list.append(col)
        merged_final_df.fillna('---', inplace=True)
        final_df = merged_final_df.apply(lambda x: x.replace('---', 'NoCall') if x.name in col_list else x)
        final_df.columns = [col.replace('.1', '') for col in final_df.columns]

        final_df.to_csv(outfile_name, index=None, mode='a', header=True, sep='\t')


if __name__ == "__main__":
    # Get user input with argparse
    parser = argparse.ArgumentParser(description="Merges two files of the same format")
    parser.add_argument(
        "--input-dir",
        type=str,
        default='.',
        help="Directory containing input files. Default is the current working directory"
    )
    parser.add_argument(
        "--input-format",
        type=str,
        choices=['TOP', 'FWD', 'AB', 'DESIGN', 'PLUS', 'LONG', 'affymetrix'],
        required=True,
        help="Type of files: 'TOP', 'FWD', 'AB', 'DESIGN', 'PLUS', 'LONG', or  'affymetrix'. All files to be merged "
             "must be in the same format. "
             "If they are in PLUS format, they must be derived from the same genome (i.e. the same genotype conversion "
             "key)."
             "'affymetrix' format can be either the native FORWARD format or the converted PLUS format."
    )
    parser.add_argument(
        '--output-name',
        type=str,
        help="Name of merged file"
    )
    parser.add_argument(
        "--file-list",
        type=str,
        help="[optional] comma-separated list of files in the input directory to be merged"
    )


    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    log_file = None
    if args.input_format == 'affymetrix':
        output = merge_affy_files(args.input_dir, args.file_list, args.input_format, args.output_name)
    else:
        output = merge_matrix_files(args.input_dir, args.file_list, args.input_format, args.output_name)
