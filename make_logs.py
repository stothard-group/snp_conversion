#!/usr/bin/python
import os
import time
import pandas as pd

# This module contains functions for outputting information to .log files


def get_logname(suffix, input_file):
    if input_file.endswith(".txt"):
        log_name = input_file.replace(".txt", suffix, 1)
    elif input_file.endswith(".csv"):
        log_name = input_file.replace(".csv", suffix, 1)
    else:
        log_name = input_file + suffix
    return log_name


def inconsistent_values_homozygous(in_file, log_exists, user_allele, var_alleles, index, col0):
    # open log file in append mode (if new file, write header information)
    if log_exists is not None:
        log_file_name = log_exists
    else:
        timestr = time.strftime("%Y%m%d%H%M%S")
        log_suffix = "-" + timestr + ".log"
        log_file_name = get_logname(log_suffix, in_file)
    if not os.path.exists(log_file_name):
        with open(log_file_name, 'a+') as log_write_header:
            log_write_header.write("Inconsistent values:\n")
            header_string = "Sample\tName\tUser Input\tConversion Key\n"
            log_write_header.write(header_string)
        log_write_header.close()
    else:
        pass
    expected_out = var_alleles[0] + " or " + var_alleles[1]
    output_list = [col0, index, user_allele, expected_out]
    output_string = "\t".join(output_list)
    with open(log_file_name, 'a+') as log_output:
        log_output.write(output_string + "\n")
    log_output.close()
    return log_file_name


def inconsistent_values(var_df, sample, user_df, in_file, log_exists):
    # open log file in append mode (if new file, write header information)
    if log_exists is not None:
        log_file_name = log_exists
    else:
        timestr = time.strftime("%Y%m%d%H%M%S")
        log_suffix = "-" + timestr + ".log"
        log_file_name = get_logname(log_suffix, in_file)
    if not os.path.exists(log_file_name):
        with open(log_file_name, 'a+') as log_write_header:
            log_write_header.write("Inconsistent values:\n")
            header_string = "Sample\tName\tUser Input\tConversion Key\n"
            log_write_header.write(header_string)
        log_write_header.close()
    else:
        pass
    # write inconsistent values to log in format
    # sample - SNP - user input - correct input (based on variant file)
    # get variant file values
    output_array = var_df.values.tolist().pop(0)
    output_array.pop()
    output_array.insert(0, sample)
    # get user values
    user_values = user_df[user_df['Name'] == output_array[1]]
    user_list = user_values.values.tolist().pop(0)
    acceptable_alleles = ['A', 'T', 'G', 'C', 'I', 'D', '-']
    # check that we have acceptable alleles
    if user_list[1] not in acceptable_alleles:
        if log_exists is None:
            os.remove(log_file_name)
        else:
            with open(log_file_name, 'a+') as log_out:
                log_out.write("Unexpected allele " + user_list[1] + " for " + output_array[1] + " in sample " + output_array[0])
            log_out.close()
        exit("Unexpected allele " + user_list[1] + " for " + output_array[1] + " in sample " + output_array[0])
    elif user_list[2] not in acceptable_alleles:
        if log_exists is None:
            os.remove(log_file_name)
        else:
            with open(log_file_name, 'a+') as log_out:
                log_out.write("Unexpected allele " + user_list[1] + " for " + output_array[1] + " in sample " + output_array[0])
            log_out.close()
        exit("Unexpected allele " + user_list[2] + " for " + output_array[1] + " in sample " + output_array[0])
    else:
        pass
    # find the incorrect value and make list to write to output
    with open(log_file_name, 'a+') as log_output:
        if user_list[1] == output_array[2] or user_list[1] == output_array[3]:
            pass
        else:
            expected_out = output_array[2] + " or " + output_array[3]
            output_list = [output_array[0], output_array[1], user_list[1], expected_out]
            output_string = "\t".join(output_list)
            log_output.write(output_string + "\n")
        if user_list[2] == output_array[2] or user_list[2] == output_array[3]:
            pass
        else:
            expected_out = output_array[2] + " or " + output_array[3]
            output_list = [output_array[0], output_array[1], user_list[2], expected_out]
            output_string = "\t".join(output_list)
            log_output.write(output_string + "\n")
            pass
    log_output.close()
    return log_file_name


def long_inequivalency(in_file, log_exists, inequiv_info):
    # open log file in append mode (if new file, write header information)
    if log_exists is not None:
        log_file_name = log_exists
    else:
        timestr = time.strftime("%Y%m%d%H%M%S")
        log_suffix = "-" + timestr + ".log"
        log_file_name = get_logname(log_suffix, in_file)
    if not os.path.exists(log_file_name):
        with open(log_file_name, 'a+') as log_write_header:
            log_write_header.write("Inequivalent values:\n")
        log_write_header.close()
    else:
        with open(log_file_name, 'a+') as log_write_header:
            log_write_header.write("\nInequivalent values:\n")
        log_write_header.close()
    # Get reformat inequiv dict and write to log file
    inequiv_df = pd.DataFrame()
    for key in inequiv_info:
        inequiv_df = pd.concat([inequiv_df, inequiv_info[key].to_frame().T], sort=False)
    inequiv_df.index.name = 'SNP Name'
    inequiv_df.reset_index(inplace=True)
    with open(log_file_name, 'a+') as log_output:
        inequiv_df.to_csv(log_output, sep='\t', index=False)
        log_output.close()
    return log_file_name


def ab_warning(dataframe, new_df, filename, log_name):
    if log_name is not None:
        log_file_name = log_name
    else:
        timestr = time.strftime("%Y%m%d%H%M%S")
        log_suffix = "-" + timestr + ".log"
        log_file_name = get_logname(log_suffix, filename)
    # check if file exists and write header if it doesn't
    if not os.path.exists(log_file_name):
        with open(log_file_name, 'a+') as log_write_header:
            log_write_header.write("Inconsistent values:\n")
            header_string = "Sample\tName\tUser Input\n"
            log_write_header.write(header_string)
        log_write_header.close()
    else:
        pass
    # Get lines that do not match AB format
    cols = list(dataframe)
    ab_list = ['AA', 'AB', 'BB', '--']
    cols.remove('Name')
    with open(log_file_name, 'a+') as ab_log:
        for col in cols:
            sample_index = dataframe.columns.get_loc(col)
            res = dataframe[new_df[col]]
            row_list = res.values.tolist()
            for row in row_list:
                output_list = [col, row[0]]
                for value in row:
                    if value not in ab_list and value == row[sample_index]:
                        output_str = "\t".join(output_list)
                        output_str = output_str + "\t" + value
                        ab_log.write(output_str + "\n")
    ab_log.close()
    return log_file_name


# Creates a simple log file and writes to it (for verbose_logging)


def simple_log(message, input_file, log_name):
    if log_name is not None:
        log_file_name = log_name
    else:
        timestr = time.strftime("%Y%m%d%H%M%S")
        log_suffix = "-" + timestr + ".log"
        log_file_name = get_logname(log_suffix, input_file)
    message_text = "\n".join(message)
    message_text = message_text + "\n"
    with open(log_file_name, 'a+') as simple:
        simple.write(message_text)
    simple.close()
    return log_file_name


