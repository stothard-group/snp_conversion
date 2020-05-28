#!/usr/bin/env python3

import os
import lib.errors as errors
import atexit
from lib.make_logs import simple_log, get_logname
import pandas as pd
import lib.common_vcf_libs as gt
import allel
import argparse
import sys
import time
import platform
import concurrent.futures as cf
import gzip
from collections import OrderedDict

# This module performs a concordance check of a SNP panel with genotype data in a vcf file


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


def find_missing_samples(line, samples_strip):
    """
    Figures out if there are samples in the user SNP panel that are not in the VCF file and reports these in a list
    :param line: the header line of the VCF file (#CHROM POS ID REF ALT...)
    :param samples_strip: list of samples found in the SNP panel file
    :return: list of samples not found in the VCF file
    """
    line.rstrip("\n")
    line_split = line.split("\t")
    vcf_samples = line_split[9:]
    missing_samples = []
    kept_samples = []
    for sam in samples_strip:
        if sam in vcf_samples:
            kept_samples.append(sam)
        else:
            missing_samples.append(sam)
    return missing_samples, kept_samples


def gt_from_sample_columns(
    vcf_df, genotype_path, can_we_thread, n_threads, logfile, file_name, log_array
):
    """
    Extracts only the genotype information from the vcf dataframe and returns a new dataframe with non-GT fields removed
    :param vcf_df: dataframe containing vcf data produced by function get_genotypes_from_vcf
    :param genotype_path: path to the genotype file
    :param can_we_thread: (bool) is threading allowed on machine
    :param n_threads: number of threads
    :param logfile: logfile
    :param file_name: name of snp panel file for logging puposes
    :param log_array: array of log text if we quit here
    :return: vcf dataframe with non-GT fields removed from each sample, also FORMAT column is removed as it is no
    longer necessary
    """
    print("Parsing Genotype GT info")

    sample_list = (list(vcf_df.columns))[8:]
    read_in_genotype = allel.read_vcf(
        genotype_path, fields=["calldata/GT"], samples=sample_list
    )
    # Reshape array based on length and number of samples
    gen_shape = read_in_genotype["calldata/GT"].shape
    rows = gen_shape[0]
    columns = gen_shape[1] * gen_shape[2]
    if gen_shape[2] != 2:
        message = (
            "More than 2 alleles detected at a locus - organism may not be polyploid"
        )
        print(message)
        log_array.append(message)
        atexit.register(simple_log, log_array, file_name, logfile)
        exit()
    reshaped = read_in_genotype["calldata/GT"].reshape(rows, columns)
    df = pd.DataFrame(reshaped)
    # Relabel df columns
    sample_to_df_dict = {}
    counter = 0
    for each_sample in sample_list:
        sam1 = each_sample + "_A1"
        sam2 = each_sample + "_A2"
        sample_to_df_dict.update({counter: sam1})
        counter = counter + 1
        sample_to_df_dict.update({counter: sam2})
        counter = counter + 1
    df.rename(columns=sample_to_df_dict, inplace=True)
    copied_vcf_df = vcf_df.drop(columns=sample_list)
    # Merge genotype data with original vcf
    genotype_df = pd.merge(
        left=copied_vcf_df, right=df, how="outer", left_index=True, right_index=True
    )
    # Check for NaN after merge, check that genotype_df is the same length as the array value gen_shape[0] (aka rows)
    nan_values = genotype_df.isnull().values.any()
    if nan_values or len(genotype_df.index) != rows:
        message = "Something went wrong in extracting genotype data, potentially with scikit-allel"
        print(message)
        log_array.append(message)
        atexit.register(simple_log, log_array, file_name, logfile)
        exit()
    genotype_df.drop(columns=["FORMAT"], inplace=True)
    # genotype_df.to_csv("outfile_test.txt", sep='\t')
    return genotype_df


def find_indels(ref_allele, alt_alleles):
    """
    Figures out if there are indels using the VCF ref and alt alleles, and makes a dict that has allele: D/I
    :param ref_allele: the vcf file reference allele
    :param alt_alleles: the vcf file alt alleles
    :return: dict containing the alt allele's index value and whether we have an I or a D
    """
    ref = len(ref_allele)
    alt_indel_dict = {}
    for al in alt_alleles:
        alt = len(al)
        if ref > alt:
            alt_indel_dict.update({alt_alleles.index(al): "D"})
        elif alt > ref:
            alt_indel_dict.update({alt_alleles.index(al): "I"})
        else:
            assert ref != alt, (
                "There may be an indel in the vcf where ref and alt alleles have the same number of "
                "bases"
            )
    return alt_indel_dict


def convert_allele(row, sample):
    """
    Performs the row-specific conversion of numeric genotypes to alleles, for each sample
    :param row: one row of the genotype dataframe produced by the function gt_from_sample_columns
    :param sample: sample name that we are working on
    :return: dataframe column for each sample with converted genotype
    """
    # print(row)
    single_vals = [int(i) for i in row[sample].split("/")]
    allele_ref = row["REF"]
    alleles_alt = row["ALT"]
    output_allele = ""
    if "," in alleles_alt:
        alt_list = alleles_alt.split(",")
    else:
        alt_list = [alleles_alt]
    # Get dict of indels
    alt_indels = False
    for alt_val in alt_list:
        if len(alt_val) > 1:
            alt_indels = True
    # Deal with * in alt alleles (indicating missing values due to overlapping deletions) as we go through the list
    for n in range(len(alt_list)):
        if n == "*":
            alt_list[n] = "D"
    if len(allele_ref) > 1 or alt_indels is True:
        indel_dict = find_indels(allele_ref, alt_list)
    else:
        indel_dict = {}
    # Run conversion
    if row[sample] == "-1/-1":
        output_allele = "--"
    elif row[sample] == "0/0":
        output_allele = allele_ref + allele_ref
    elif row[sample] == "0/1":
        if 0 not in indel_dict:
            output_allele = allele_ref + alt_list[0]
        else:
            output_allele = allele_ref + indel_dict[0]
    elif row[sample] == "1/1":
        if 0 not in indel_dict:
            output_allele = alt_list[0] + alt_list[0]
        else:
            output_allele = indel_dict[0] + indel_dict[0]
    elif single_vals[0] >= 2 or single_vals[1] >= 2:
        out1 = None
        out2 = None
        counter = 2
        n_alleles = len(alleles_alt)
        if single_vals[0] == 0:
            out1 = allele_ref
        elif single_vals[0] == 1:
            if 0 not in indel_dict:
                out1 = alt_list[0]
            else:
                out1 = indel_dict[0]
        else:
            while counter <= n_alleles:
                if single_vals[0] == counter:
                    minus_count = counter - 1
                    if minus_count not in indel_dict:
                        out1 = alt_list[counter - 1]
                    else:
                        out1 = indel_dict[counter - 1]
                counter = counter + 1
                # print(counter, n_alleles)
        counter = 2
        while counter <= n_alleles:
            if single_vals[1] == counter:
                minus_count = counter - 1
                if minus_count not in indel_dict:
                    out2 = alt_list[counter - 1]
                else:
                    out2 = indel_dict[counter - 1]
            counter = counter + 1
            # print(counter, n_alleles)
        if out1 is not None and out2 is not None:
            output_allele = out1 + out2
        else:
            exit("error converting multiple output alleles")

    else:
        exit("Fatal error - genotype code is not recognized " + row[sample])
    return output_allele


def numeric_to_allele_conversion(geno_parsed_df):
    """
    Converts genotypes from numeric format (-1,0,1,etc.) to alleles in REF and ALT, returns converted genotype dataframe
    :param geno_parsed_df: dataframe containing subsampled vcf information and individual genotypes from get_genotypes_from_vcf
    function
    :return: genotype dataframe with both alleles and numbers
    """
    sample_list = (list(geno_parsed_df.columns))[7:]
    basenames = []
    for sample in sample_list:
        sample_basename = sample.rsplit("_", 1)[0]
        if sample_basename not in basenames:
            basenames.append(sample_basename)
    for each_base in basenames:
        # Combine the sample dataframes first
        pos1 = each_base + "_A1"
        pos2 = each_base + "_A2"
        geno_parsed_df[each_base] = (
            geno_parsed_df[pos1]
            .astype(str)
            .str.cat(geno_parsed_df[pos2].astype(str), sep="/")
        )
        geno_parsed_df.drop(columns=[pos1, pos2], inplace=True)
    for each_base2 in basenames:
        new_sample_name = each_base2 + "_allele"
        print("Allele conversion for sample " + each_base2)
        geno_parsed_df[new_sample_name] = geno_parsed_df.apply(
            lambda row: convert_allele(row, each_base2), axis=1
        )
    geno_parsed_df_finished = geno_parsed_df.copy()
    return geno_parsed_df_finished


def get_genotypes_from_vcf(
    genotype_path,
    subsampled_gt_df,
    full_gt_df,
    can_we_thread,
    n_threads,
    logfile,
    file_name,
):
    """
    Extracts genotypes from user-input vcf file
    :param genotype_path: path to the genotype file
    :param subsampled_gt_df: subsampled vcf file that has only positions in the user SNP file
    :param full_gt_df: full vcf dataframe
    :param can_we_thread: (bool) is threading allowed on this machine
    :param n_threads: number of threads to use
    :param logfile: logfile
    :return: dataframe containing normal vcf info + genotypes (no format info), where each genotype has its own column.
    This is a subsampled dataframe containing only the positions found in the user SNP file
    :return: Logfile
    """
    # Parse sample columns to remove non-genotype information
    log_array = []
    timestr = time.strftime("%H:%M:%S")
    logfile_text = timestr + " ..... Converting genotype information from VCF file "
    log_array.append(logfile_text)

    gt_parsed_df = gt_from_sample_columns(
        full_gt_df,
        genotype_path,
        can_we_thread,
        n_threads,
        logfile,
        file_name,
        log_array,
    )
    # Remove the rows that are not in the subsampled df before doing allele conversion
    subsampled_positions = subsampled_gt_df[["#CHROM", "POS"]]
    subsampled_gt_parsed_df = pd.merge(
        left=gt_parsed_df, right=subsampled_positions, on=["#CHROM", "POS"], how="inner"
    )
    # Convert numeric genotypes to alleles
    subsampled_converted_df = numeric_to_allele_conversion(subsampled_gt_parsed_df)
    if logfile is not None:
        logging = simple_log(log_array, file_name, logfile)
    return subsampled_converted_df, logfile


def check_filter_values(qual_val, filt_val, vcf_df, log_array, logfile, file_name):
    """
    If the user has requested filtering, this checks that the filter values are appropriate.
    The quality value must be a float or int, and the filt_val must be one or more values from the FILTER column
    of the vcf file
    :param qual_val: the QUAL value to filter on
    :param filt_val: the list containing one or more FILTER values to filter on
    :param vcf_df: dataframe of entire vcf file
    :param log_array: array of log messages to append here if we end up exiting
    :param logfile: logfile
    :param file_name: name of snp panel file for logging purposes
    :return: Bool value of whether the quality values are appropriate
    """
    # Deal with qual values
    if qual_val is not None:
        if isinstance(qual_val, int) or isinstance(qual_val, float):
            pass
        else:
            message = "QUAL filtering value must be a float or int"
            log_array.append(message)
            atexit.register(simple_log, log_array, file_name, logfile)
            exit()

    else:
        pass
    # Deal with filter values
    possible_filt_vals = vcf_df.FILTER.unique()
    if filt_val:
        filt_val_check = set(filt_val).issubset(possible_filt_vals)
        if not filt_val_check:
            message = "One or more values for filtering from FILTER are not valid"
            log_array.append(message)
            atexit.register(simple_log, log_array, file_name, logfile)
            exit()
    else:
        pass
    return True


def handle_vcf_file(
    file_handle, vcf_header_count, samples_stripped, log_array, file_name, logfile
):
    """
    reads a filehandle that is not compressed (or the gzip-parsed version)
    :param file_handle: open genotype file (VCF) filehandle
    :return: missing and kept samples
    """
    missing_samples = []
    kept_samples = []
    for line in file_handle:
        if line.startswith("##"):
            vcf_header_count = vcf_header_count + 1
        elif line.startswith("#CHROM"):
            # Handle missing samples
            missing_samples, kept_samples = find_missing_samples(line, samples_stripped)
            break
        else:
            logfile_text = (
                "VCF file does not contain a line starting with '#CHROM', signalling the beginning of "
                "the data. Exiting..."
            )
            log_array.append(logfile_text)
            atexit.register(simple_log, log_array, file_name, logfile)
            exit()
    return vcf_header_count, missing_samples, kept_samples


def subsample_vcf_by_position(
    genotype_file_input,
    dict_of_panels,
    samples,
    filter_conditions,
    qual_crit,
    filter_crit,
    logfile,
    file_name,
):
    """
    Subsamples the vcf file to only include the positions that are represented by the SNP panels, to increase speed
    :param genotype_file_input: vcf file
    :param dict_of_panels: dict containing one panel per animal, where each df contains the genotype and positional info
    :param filter_conditions: bool for whether we should filter
    :param qual_crit: quality value to filter on, if requested
    :param filter_crit: filter value to filter on, if requested
    :param logfile: logfile
    :param file_name: snp panel file for logging purposes
    :return: subsampled vcf dataframe with only relevant columns
    :return: full vcf dataframe
    :return: list of missing samples
    :return: number of vcf sites NOT represented in the SNP panel
    :return: logfile
    """
    log_array = []
    timestr = time.strftime("%H:%M:%S")
    logfile_text = timestr + " ..... Subsampling VCF file to reduce computation "
    log_array.append(logfile_text)
    sample_panel_list = next(iter(dict_of_panels.values()))
    sample_panel = sample_panel_list[0]
    # change dtype of sample panel BLAST_chromosome
    sample_panel.BLAST_chromosome.astype(object)
    # Read in genotype file to dataframe
    # Strip '.1' from user samples (only if affymetrix)
    samples_stripped = [sam.replace(".1", "") for sam in samples]
    # Get line number of vcf data start
    vcf_header_count = 0
    if genotype_file_input.endswith(".gz"):
        with gzip.open(genotype_file_input, "rt", errors="ignore") as gz_gt_file:
            vcf_header_count, missing_samples, kept_samples = handle_vcf_file(
                gz_gt_file,
                vcf_header_count,
                samples_stripped,
                log_array,
                file_name,
                logfile,
            )
        gz_gt_file.close()
    else:
        with open(genotype_file_input, "r",) as gt_file:
            vcf_header_count, missing_samples, kept_samples = handle_vcf_file(
                gt_file,
                vcf_header_count,
                samples_stripped,
                log_array,
                file_name,
                logfile,
            )
        gt_file.close()
    # Read vcf to dataframe and skip lines until the data
    normal_vcf_cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "FORMAT"]
    columns_to_use = normal_vcf_cols + kept_samples
    vcf_to_df = pd.read_csv(
        genotype_file_input,
        compression="infer",
        skiprows=vcf_header_count,
        sep="\t",
        usecols=columns_to_use,
        dtype={"#CHROM": object},
    )
    # Do a quick check if the user has requested filtering
    if filter_conditions:
        timestr = time.strftime("%H:%M:%S")
        logfile_text = timestr + " ..... Subsampling VCF file to reduce computation "
        log_array.append(logfile_text)
        filter_check = check_filter_values(
            qual_crit, filter_crit, vcf_to_df, log_array, logfile, file_name
        )
    else:
        pass
    vcf_for_merge = vcf_to_df.copy()
    vcf_for_merge.rename(
        columns={"#CHROM": "BLAST_chromosome", "POS": "BLAST_position"}, inplace=True
    )
    # sample_panel['BLAST_position'] = sample_panel['BLAST_position'].astype(int)
    # 'merge' the sample panel with the vcf file on chromosome and blast position information
    merged_vcf_panel = pd.merge(
        left=sample_panel,
        right=vcf_for_merge,
        on=["BLAST_chromosome", "BLAST_position"],
        how="outer",
        indicator=True,
    )
    # Count number of vcf positions that are not in the user snp panel
    vcf_pos_only = merged_vcf_panel[merged_vcf_panel["_merge"] == "right_only"].shape[0]
    merged_vcf_panel_both = merged_vcf_panel[
        merged_vcf_panel["_merge"] == "both"
    ].copy()
    merged_vcf_panel_both.rename(
        columns={"BLAST_chromosome": "#CHROM", "BLAST_position": "POS"}, inplace=True
    )
    subsampled_df = merged_vcf_panel_both[columns_to_use].copy()
    if logfile is not None:
        logging = simple_log(log_array, file_name, logfile)
    return subsampled_df, vcf_to_df, missing_samples, vcf_pos_only, logfile


def get_discordant_name_list(both_list, forward_df, reverse_df):
    """
    This gets a list of Names (probes) that are discordant; i.e. are found in both snp and vcf file based on positional
    info, but are not in the forward or reverse dfs
    :param both_list: list of names (probes) in both the snp panel and vcf file based on positional info
    :param forward_df: panel and vcf file merged on positions and FORWARD genotype
    :param reverse_df: panel and vcf file merged on positions and REVERSE genotype
    :return: List of ids that are discordant
    """
    forward_list = forward_df["Name"].tolist()
    reverse_list = reverse_df["Name"].tolist()
    combined_list = forward_list + sorted(set(reverse_list) - set(forward_list))
    discordant_list = list(set(both_list) - set(combined_list))
    return discordant_list


def quality_filtering(
    vcf_to_filter, panel_to_filter, qual_val, filter_val, log_array, logfile, file_name
):
    """
    Filters the vcf file (and the panel file) based on the quality score and/or the filter value
    :param vcf_to_filter: vcf file that has been subsampled (prior to panel concordancy analysis in compare_genotypes)
    :param panel_to_filter: snp panel (prior to panel concordancy analysis in compare_genotypes)
    :param qual_val: user-input quality value to filter on
    :param filter_val: one or more values in the Filter column; can be any value in the VCF file
    :param log_array: array of logfile statements if we exit here
    :param logfile: logfile
    :param file_name: name of snp panel file for logging purposes
    :return: filtered vcf file
    :return: filtered panel file
    :return: number of positions excluded based on filtering criteria
    """
    # perform vcf filtering
    vcf_filtered_on_qual = vcf_to_filter[vcf_to_filter["QUAL"] > qual_val].copy()
    vcf_filtered_qual_and_filter = vcf_filtered_on_qual[
        vcf_filtered_on_qual["FILTER"].isin(filter_val)
    ].copy()
    # subtract the number of positions in the filtered vcf file from the original
    vcf_filtered_length = len(vcf_filtered_qual_and_filter.index)
    vcf_unfiltered_length = len(vcf_to_filter.index)
    num_filtered_pos = vcf_unfiltered_length - vcf_filtered_length
    # filter the snp panel file
    panel_columns = list(panel_to_filter.columns)
    merge_filtered = pd.merge(
        left=vcf_filtered_qual_and_filter,
        right=panel_to_filter,
        left_on=["#CHROM", "POS"],
        right_on=["BLAST_chromosome", "BLAST_position"],
        how="inner",
    )
    panel_filtered = merge_filtered[panel_columns]
    # check that the vcf filtered df is the same length as the panel filtered df
    if len(panel_filtered.index) != len(vcf_filtered_qual_and_filter.index):
        message = "Quality filtering may not have been performed properly"
        log_array.append(message)
        atexit.register(simple_log, log_array, file_name, logfile)
        exit()
    return vcf_filtered_qual_and_filter, panel_filtered, num_filtered_pos


def comparison(
    animal,
    SNP_panel_dict,
    vcf_converted,
    filter_bool,
    qual_tofilter,
    filt_tofilter,
    log_array,
    logfile,
    file_name,
):
    """
    Function that actually performs the comparison - separated so it can be run on a parallelized per-sample basis
    :param animal: name of the animal/sample we are doing comparison with
    :param SNP_panel_dict: dict containing the user SNP panels
    :param vcf_converted: dataframe of converted, subsampled vcf file (also contains non _allele positions with numeric
    genotypes)
    :param filter_bool: Bool to filter samples
    :param qual_tofilter: quality value to filter on
    :param filt_tofilter: filter value to filter on
    :param log_array: arracy containing log messagees to write to log file if program exits during quality_filtering
    :param logfile: logfile
    :param file_name: name of input file for logging purposes
    :return:
    1. correct homozygous matches
    2. correct (forward) heterozygous matches
    3. reverse heterozygous matches
    4. dashes in SNP file
    5. discordant positions
    6. number of filtered positions, if filtering was done (otherwise 'None')
    """
    animal_name = animal.replace(".1", "")
    working_snp_panel = SNP_panel_dict[animal][0]
    animal_id = animal_name + "_allele"
    # working_snp_panel['BLAST_position'] = working_snp_panel['BLAST_position'].astype(int)
    # add a reversal column here for the SNP panel df
    working_snp_panel["reversal"] = working_snp_panel.loc[:, animal].apply(
        lambda x: x[::-1]
    )
    working_vcf = vcf_converted[
        ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", animal_id]
    ].copy()
    # Do filtering here, so we can compute only filtered positions and provide a count of excluded positions
    if filter_bool is True:
        new_working_vcf, new_working_snp_panel, filtered_pos_int = quality_filtering(
            working_vcf,
            working_snp_panel,
            qual_tofilter,
            filt_tofilter,
            log_array,
            logfile,
            file_name,
        )
        working_vcf.update(new_working_vcf)
        working_snp_panel.update(new_working_snp_panel)
        filtered_pos_val = str(filtered_pos_int)
    else:
        filtered_pos_val = None
    # subset the snp panel
    merged_snp_vcf_panel = pd.merge(
        left=working_vcf,
        right=working_snp_panel,
        left_on=["#CHROM", "POS"],
        right_on=["BLAST_chromosome", "BLAST_position"],
        how="outer",
        indicator=True,
    )
    # Count number of snp positions that are not in the vcf and get list of panel-specific IDs
    panel_only_not_vcf = merged_snp_vcf_panel[
        merged_snp_vcf_panel["_merge"] == "right_only"
    ]
    # Count number of positions that are in both the panel and vcf and get list of both IDs
    # This also subsamples the panel with the correct name/chrom/position information for panel and vcf file
    panel_and_vcf_df = merged_snp_vcf_panel[merged_snp_vcf_panel["_merge"] == "both"]
    panel_and_vcf_list = panel_and_vcf_df["Name"].tolist()
    subset_snp_panel = panel_and_vcf_df[
        ["Name", animal, "BLAST_chromosome", "BLAST_position", "reversal"]
    ]
    # Start doing the merges to actually assess concordance

    merged_panel_and_vcf_forward = pd.merge(
        left=subset_snp_panel,
        right=working_vcf,
        how="inner",
        left_on=["BLAST_chromosome", "BLAST_position", animal],
        right_on=["#CHROM", "POS", animal_id],
    )
    merged_panel_and_vcf_reverse = pd.merge(
        left=subset_snp_panel,
        right=working_vcf,
        how="inner",
        left_on=["BLAST_chromosome", "BLAST_position", "reversal"],
        right_on=["#CHROM", "POS", animal_id],
    )
    discordant_merge_list = get_discordant_name_list(
        panel_and_vcf_list, merged_panel_and_vcf_forward, merged_panel_and_vcf_reverse
    )
    # If we have discordant values, make a discordant df, otherwise, make an empty df
    # The discordant df also has AA/-- and --/AA panel/VCF values
    discordant_df = panel_and_vcf_df[
        panel_and_vcf_df["Name"].isin(discordant_merge_list)
    ]

    # Find any dashes (-- or ---) in the snp panel file, and find corresponding positions in vcf file
    if "--" in subset_snp_panel[animal].values:
        snp_dashes = subset_snp_panel[subset_snp_panel[animal] == "--"]
    elif "---" in subset_snp_panel[animal].values:
        snp_dashes = subset_snp_panel[subset_snp_panel[animal] == "---"]
        snp_dashes.replace({"---": "--"}, inplace=True)
        dashes_path = animal + "_snp_dashes.txt"
        snp_dashes.to_csv(dashes_path, sep="\t")
    else:
        snp_dashes = pd.DataFrame()
    if not snp_dashes.empty:
        dashes_merge_whole = pd.merge(
            left=snp_dashes,
            right=working_vcf,
            how="left",
            left_on=["BLAST_chromosome", "BLAST_position", animal],
            right_on=["#CHROM", "POS", animal_id],
        )
        dashes_merge_whole.dropna(inplace=True)
        dashes_merge = dashes_merge_whole[
            dashes_merge_whole[animal_id].isin(["--"])
        ].copy()
    else:
        dashes_merge = pd.DataFrame()

    # Get the correctly-merged homozygous dataframe
    nn_list = ["AA", "TT", "CC", "GG"]
    correct_homozygous_df = merged_panel_and_vcf_forward[
        merged_panel_and_vcf_forward[animal_id].isin(nn_list)
    ].copy()

    # Get the correctly-merged heterozygous dataframe (remove homozygous positions)
    homozygous_name_list = correct_homozygous_df["Name"].tolist()
    correct_heterozygous_df = merged_panel_and_vcf_forward[
        ~merged_panel_and_vcf_forward["Name"].isin(homozygous_name_list)
    ].copy()

    # Get the reverse-merged heterozygous dataframe (remove homozygous positions)
    reverse_heterozygous_df = merged_panel_and_vcf_reverse[
        ~merged_panel_and_vcf_reverse["Name"].isin(homozygous_name_list)
    ].copy()

    # Pass back dataframes as values in a dict
    # correct_homozygous_df: AA = AA
    # correct_heterozygous_df: AB = AB
    # reverse_heterozygous_df: AB = BA
    # dashes_merge: -- = --
    # discordant_df: AB = AA, --/AA, AA/--
    info_to_return = {
        animal_name: [
            correct_homozygous_df,
            correct_heterozygous_df,
            reverse_heterozygous_df,
            dashes_merge,
            discordant_df,
            filtered_pos_val,
        ]
    }
    return info_to_return


def compare_genotypes(
    vcf_converted,
    SNP_panel_dict,
    missing_samples,
    filter_bool,
    qual_tofilter,
    filt_tofilter,
    can_we_thread,
    n_threads,
    logfile,
    file_name,
):
    """
    Performs the comparison between the user SNP panel and the subsampled VCF dataframe. Written on a per-animal basis
    so that eventually this could be parallelized with futures
    :param vcf_converted: dataframe of converted, subsampled vcf file (also contains non _allele positions with numeric genotypes)
    :param SNP_panel_dict: dict containing user SNP panels
    :param missing_samples: list of samples missing from the vcf
    :param filter_bool: Bool to filter samples
    :param qual_tofilter: quality value to filter on
    :param filt_tofilter: filter value to filter on
    :param can_we_thread: (bool) is threading allowed on this machine
    :param n_threads: number of threads to use
    :param logfile: logfile
    :param file_name: name of snp panel file for logging purposes
    :return: dict of animals where values are
    1. correct homozygous matches
    2. correct (forward) heterozygous matches
    3. reverse heterozygous matches
    4. dashes in SNP file
    5. discordant positions
    6. number of filtered positions, if filtering was done (otherwise 'None')
    :return: logfile
    """
    log_array = []
    timestr = time.strftime("%H:%M:%S")
    logfile_text = timestr + " ..... Beginning genotype comparison "
    log_array.append(logfile_text)
    for missing in missing_samples:
        if missing in vcf_converted:
            del SNP_panel_dict[missing]
    # remove missing from the vcf before iterating through them
    present_animals_list = []
    for animal in SNP_panel_dict:
        animal_id = animal.replace(".1", "_allele")
        if animal_id in vcf_converted.columns:
            present_animals_list.append(animal)
        else:
            pass
    if filter_bool is True:
        message = "Filtering VCF file on user-defined values"
        print(message)
        log_array.append(message)

    return_dict = {}
    # Parallelize here
    if can_we_thread is True:
        results_list = []
        with cf.ProcessPoolExecutor(max_workers=n_threads) as executor:
            animal_df_dict = {
                executor.submit(
                    comparison,
                    animal,
                    SNP_panel_dict,
                    vcf_converted,
                    filter_bool,
                    qual_tofilter,
                    filt_tofilter,
                    log_array,
                    logfile,
                    file_name,
                ): animal
                for animal in present_animals_list
            }
            for n in cf.as_completed(animal_df_dict):
                data = n.result()
                results_list.append(data)
        for each_result in results_list:
            key_list = list(each_result.keys())
            each_sample = key_list[0]
            animal_name = each_sample.replace(".1", "")
            return_dict.update({animal_name: each_result[animal_name]})
    else:  # can_we_thread is false
        for animal in present_animals_list:
            animal_name = animal.replace(".1", "")
            result = comparison(
                animal,
                SNP_panel_dict,
                vcf_converted,
                filter_bool,
                qual_tofilter,
                filt_tofilter,
                log_array,
                logfile,
                file_name,
            )
            return_dict.update({animal_name: result[animal_name]})

    if logfile is not None:
        logging = simple_log(log_array, file_name, logfile)
    return return_dict, logfile


def get_metastats_numbers(
    samp_list_total,
    comp_dict,
    panel_only_samples,
    sample_snp_panel,
    num_vcf_pos,
    num_shared_pos,
):
    """
    Gets all the values needed for writing to the metastatistics output file
    :param samp_list_total: complete list of samples in the snp panel (all values)
    :param comp_dict: dict where keys are samples that are in the panel and vcf
    :param panel_only_samples: list of animals that are only in the snp panel (not vcf)
    :param sample_snp_panel: a random snp panel from the panel dict
    :param num_vcf_pos: number of positions in the vcf file
    :param num_shared_pos: number of shared positions in sample and vcf file
    :return full_sample_count: total sample count in snp panel
    :return snp_in_vcf_count: samples in panel and vcf (count)
    :return snp_in_vcf_percentage: percent of samples in panel and vcf
    :return both_string: string containing all animals in both panel and vcf
    :return snp_only_string: string containing animals only the panel (not vcf)
    :return number_panel_positions: number of positions in the SNP panel (total)
    :return number_pos_in_vcf: number of positions in the vcf file (total)
    :return number_shared_pos: number of positions shared by panel and vcf
    :return percent_shared_positions_panel: shared positions as a percent of the panel
    :return percent_shared_positions_vcf: shared positions as a percent of the vcf
    :return filtered_pos: number of positions that were filtered based on quality or filter scores
    """
    print("Calculating metastatistics")
    full_sample_count = str(len(samp_list_total))  # total number of samples
    snp_in_vcf_count = str(len(comp_dict))  # number of samples in snp panel and vcf
    snp_in_vcf_percentage = str(
        round((len(comp_dict) / len(samp_list_total) * 100), 0)
    )  # percent of samples in snp panel
    # and vcf compared to total samples
    OD_comp_dict = OrderedDict(comp_dict)
    both_string = ", ".join(
        list(OD_comp_dict.keys())
    )  # list - all animals in both panel and vcf
    snp_only_string = ", ".join(
        panel_only_samples
    )  # list - animals only in panel (not vcf)
    number_panel_positions = str(
        len(sample_snp_panel.index)
    )  # total number of SNP panel positions
    number_pos_in_vcf = str(num_vcf_pos)  # number of positions in vcf
    number_shared_pos = str(
        num_shared_pos
    )  # number of positions shared by panel and vcf
    percent_shared_positions_panel = str(
        round((num_shared_pos / len(sample_snp_panel.index) * 100), 0)
    )  # percent of
    # positions shared by panel and vcf, relative to the panel
    int_percent_shared_positions_vcf = round(
        (num_shared_pos / num_vcf_pos) * 100, 0
    )  # percent of positions shared by
    # panel and vcf, relative to the vcf
    if int_percent_shared_positions_vcf < 0.1:
        percent_shared_positions_vcf = "<0.1"
    else:
        percent_shared_positions_vcf = str(int_percent_shared_positions_vcf)
    # number of filtered positions
    first_compdict_entry = next(iter(comp_dict.values()))
    filtered_pos = first_compdict_entry[5]

    return [
        full_sample_count,
        snp_in_vcf_count,
        snp_in_vcf_percentage,
        both_string,
        snp_only_string,
        number_panel_positions,
        number_pos_in_vcf,
        number_shared_pos,
        percent_shared_positions_panel,
        percent_shared_positions_vcf,
        filtered_pos,
    ]


def count_ref_or_alt(hom_df, animal_snp_name):
    """
    Counts how many homozygous SNP panel alleles are the reference base or the alternate base
    :param hom_df: the homozygous match dataframe
    :param animal_snp_name: the SNP panel name of the animal which we split and compare to REF or ALT from vcf
    :return: a list, [number of homozygous reference alleles, number of homozygous alternate alleles]
    """
    hom_df["A1"], hom_df["A2"] = zip(*hom_df[animal_snp_name].apply(lambda x: list(x)))
    hom_df["ref_true"] = (hom_df["A1"] == hom_df["REF"]).astype(int)
    ref_true_count = hom_df["ref_true"].sum()
    hom_df["alt_true"] = hom_df.apply(lambda row: row.A1 in row.ALT, axis=1).astype(int)
    alt_true_count = hom_df["alt_true"].sum()
    ref_alt_list = [ref_true_count, alt_true_count]
    return ref_alt_list


def count_discordant_types(disc_df, animal_snp_name, animal_vcf_name):
    """
    Counts the ways in which discordant alleles exist
    :param disc_df: the dataframe of all informational discordant positions
    :param animal_snp_name: the name of the column containing the SNP panel alleles
    :param animal_vcf_name: the name of the column containing the VCF panel alleles
    :return: two dicts, one for homozygous stats and heterozygous stats
    PANEL vs VCF
    hom:
    hom_ref_reversal_count = AA vs BB | BB vs AA
    hom_to_het_ref_count = AA vs AB | BB vs AB
    hom_to_hom_not_vcf_count = CC vs AA/BB/AB
    het:
    het_to_hom_vcf_count = AB vs AA | AB vs BB
    het_not_vcf_count = AC vs AA/BB/AB | CD vs AA/BB/AB
    """
    # First remove any '--' from the discordant df
    if not disc_df.empty:
        # disc_df.replace('--', np.nan, inplace=True)
        # disc_df.dropna(inplace=True)
        disc_df = disc_df[disc_df[animal_snp_name] != "--"]
        disc_df = disc_df[disc_df[animal_vcf_name] != "--"]

    if not disc_df.empty:  # If there is a discordant dataframe
        disc_df["panel_A1"], disc_df["panel_A2"] = zip(
            *disc_df[animal_snp_name].apply(lambda x: list(x))
        )
        disc_df["vcf_A1"], disc_df["vcf_A2"] = zip(
            *disc_df[animal_vcf_name].apply(lambda x: list(x))
        )
        hom_disc_df = disc_df[disc_df["panel_A1"] == disc_df["panel_A2"]].copy()
        het_disc_df = disc_df[disc_df["panel_A1"] != disc_df["panel_A2"]].copy()
        if hom_disc_df.empty:
            hom_ref_reversal_count = 0
            hom_to_het_ref_count = 0
            hom_to_hom_not_vcf_count = 0
        else:
            # discordant and homozygous panel to homozygous ref/alt vcf (alleles == both ref, or alleles == both alt) (AA vs BB)
            hom_disc_df["hom_ref_reversal"] = (
                (hom_disc_df["vcf_A1"] == hom_disc_df["vcf_A2"])
                & (
                    (hom_disc_df["panel_A1"] == hom_disc_df["REF"])
                    | (hom_disc_df.apply(lambda r: r.panel_A1 in r.ALT, axis=1))
                )
            ).astype(int)
            hom_ref_reversal_count = hom_disc_df["hom_ref_reversal"].sum()
            # discordant and homozygous panel to heterozygous vcf (alleles == REF/ALT) (AA vs AB)
            hom_disc_df["hom_to_het_ref"] = (
                (hom_disc_df["vcf_A1"] != hom_disc_df["vcf_A2"])
                & (
                    (hom_disc_df["panel_A1"] == hom_disc_df["REF"])
                    | (hom_disc_df.apply(lambda r: r.panel_A1 in r.ALT, axis=1))
                )
            ).astype(int)
            hom_to_het_ref_count = hom_disc_df["hom_to_het_ref"].sum()
            # discordant and homozygous to hom not ref or alt (CC vs AA/AB/BB)
            hom_disc_df["hom_to_hom_not_vcf"] = (
                (hom_disc_df["panel_A1"] != hom_disc_df["REF"])
                & (hom_disc_df.apply(lambda r: r.panel_A1 not in r.ALT, axis=1))
            ).astype(int)
            hom_to_hom_not_vcf_count = hom_disc_df["hom_to_hom_not_vcf"].sum()
        if het_disc_df.empty:
            het_to_hom_vcf_count = 0
            het_not_vcf_count = 0
        else:
            # discordant and heterozygous to homozygous (alt or ref) (AB vs AA)
            het_disc_df["het_to_hom_vcf"] = (
                (het_disc_df["vcf_A1"] == het_disc_df["vcf_A2"])
                & (
                    (
                        (het_disc_df["panel_A1"] == het_disc_df["REF"])
                        & (het_disc_df.apply(lambda r: r.panel_A2 in r.ALT, axis=1))
                    )
                    | (
                        (het_disc_df.apply(lambda r: r.panel_A1 in r.ALT, axis=1))
                        & (het_disc_df["panel_A2"] == het_disc_df["REF"])
                    )
                )
            ).astype(int)
            het_to_hom_vcf_count = het_disc_df["het_to_hom_vcf"].sum()
            # discordant and heterozygous not in ref or alt (AC or CD vs ?? doesn't matter)
            het_disc_df["het_not_vcf"] = (
                (
                    (
                        (het_disc_df["panel_A1"] != het_disc_df["REF"])
                        & (het_disc_df.apply(lambda r: r.panel_A1 not in r.ALT, axis=1))
                    )
                    | (
                        (het_disc_df["panel_A2"] != het_disc_df["REF"])
                        & (het_disc_df.apply(lambda r: r.panel_A2 not in r.ALT, axis=1))
                    )
                )
            ).astype(int)
            het_not_vcf_count = het_disc_df["het_not_vcf"].sum()
        hom_dict = {
            "hom_to_hom_vcf": hom_ref_reversal_count,
            "hom_to_het_vcf": hom_to_het_ref_count,
            "hom_to_not_vcf": hom_to_hom_not_vcf_count,
        }
        het_dict = {
            "het_to_hom_vcf": het_to_hom_vcf_count,
            "het_to_not_vcf": het_not_vcf_count,
        }
    else:  # if there is no discordant dataframe, set all discordant counts to 0
        hom_dict = {
            "hom_to_hom_vcf": 0,
            "hom_to_het_vcf": 0,
            "hom_to_not_vcf": 0,
        }
        het_dict = {"het_to_hom_vcf": 0, "het_to_not_vcf": 0}

    return hom_dict, het_dict


def calc_percentages(stats_list):
    """
    Calculates the appropriate percentages for each statistic
    :param stats_list: list of statistics from generate_animal_stats
    :return: new list containing stats with percentages
    """
    # informational_number, noninformational_number, homozygous_ref_informational_number,
    # homozygous_alt_informational_number, heterozygous_informational_number, total_concordant_number,
    # total_true_concordant_number, total_discordant_number, homozygous_ref_concordant_number,
    # homozygous_alt_concordant_number, snp_hom_vcf_hom, snp_hom_vcf_het, snp_hom_no_vcf,
    # heterozygous_concordant_whole_number, heterozygous_concordant_true_number,
    # heterozygous_concordant_reverse_number, snp_het_vcf_hom, snp_het_no_vcf
    informational_number = stats_list[0]
    noninformational_number = stats_list[1]
    total_positions = informational_number + noninformational_number
    informational_percentage = round((informational_number / total_positions) * 100, 2)
    informational_val = (
        str(informational_number) + " (" + str(informational_percentage) + "%)"
    )
    noninformational_percentage = round(
        (noninformational_number / total_positions) * 100, 2
    )
    noninformational_val = (
        str(noninformational_number) + " (" + str(noninformational_percentage) + "%)"
    )

    # Of informational percentage
    homozygous_ref_info_number = stats_list[2]
    homozygous_ref_percentage = round(
        (homozygous_ref_info_number / informational_number) * 100, 2
    )
    homozygous_ref_val = (
        str(homozygous_ref_info_number) + " (" + str(homozygous_ref_percentage) + "%)"
    )
    homozygous_alt_info_number = stats_list[3]
    homozygous_alt_info_percentage = round(
        (homozygous_alt_info_number / informational_number) * 100, 2
    )
    homozygous_alt_val = (
        str(homozygous_alt_info_number)
        + " ("
        + str(homozygous_alt_info_percentage)
        + "%)"
    )
    heterozygous_informational_number = stats_list[4]
    heterozygous_percentage = round(
        (heterozygous_informational_number / informational_number) * 100, 2
    )
    heterozygous_val = (
        str(heterozygous_informational_number)
        + " ("
        + str(heterozygous_percentage)
        + "%)"
    )

    total_concordant = stats_list[5]
    total_concordant_percentage = round(
        (total_concordant / informational_number) * 100, 2
    )
    total_concordant_val = (
        str(total_concordant) + " (" + str(total_concordant_percentage) + "%)"
    )
    total_true_concordant = stats_list[6]
    total_true_concordant_percentage = round(
        (total_true_concordant / informational_number) * 100, 2
    )
    total_true_concordant_val = (
        str(total_true_concordant) + " (" + str(total_true_concordant_percentage) + "%)"
    )
    total_discordant = stats_list[7]
    total_discordant_percentage = round(
        (total_discordant / informational_number) * 100, 2
    )
    total_discordant_val = (
        str(total_discordant) + " (" + str(total_discordant_percentage) + "%)"
    )

    # Of concordant percentage
    homozygous_ref_concordant = stats_list[8]
    hom_ref_concordant_percentage = round(
        (homozygous_ref_concordant / total_concordant) * 100, 2
    )
    hom_ref_conc_val = (
        str(homozygous_ref_concordant)
        + " ("
        + str(hom_ref_concordant_percentage)
        + "%)"
    )
    hom_alt_concordant = stats_list[9]
    hom_alt_concordant_percentage = round(
        (hom_alt_concordant / total_concordant) * 100, 2
    )
    hom_alt_conc_val = (
        str(hom_alt_concordant) + " (" + str(hom_alt_concordant_percentage) + "%)"
    )

    heterozygous_concordant_whole = stats_list[13]
    het_conc_whole_percentage = round(
        (heterozygous_concordant_whole / total_concordant) * 100, 2
    )
    het_conc_whole_val = (
        str(heterozygous_concordant_whole)
        + " ("
        + str(het_conc_whole_percentage)
        + "%)"
    )
    heterozygous_concordant_true = stats_list[14]
    het_conc_true_percentage = round(
        (heterozygous_concordant_true / total_concordant) * 100, 2
    )
    het_conc_true_val = (
        str(heterozygous_concordant_true) + " (" + str(het_conc_true_percentage) + "%)"
    )
    heterozygous_concordant_reverse = stats_list[15]
    het_conc_rev_percentage = round(
        (heterozygous_concordant_reverse / total_concordant) * 100, 2
    )
    het_conc_rev_val = (
        str(heterozygous_concordant_reverse)
        + " ("
        + str(het_conc_rev_percentage)
        + "%)"
    )

    # Of discordant percentage
    snp_hom_vcf_hom = stats_list[10]
    snp_hom_vcf_het = stats_list[11]
    snp_hom_no_vcf = stats_list[12]

    snp_het_vcf_hom = stats_list[16]
    snp_het_no_vcf = stats_list[17]

    if total_discordant == 0:
        snp_hom_vcf_hom_percentage = 0
        snp_hom_vcf_het_percentage = 0
        snp_hom_no_vcf_percentage = 0
        snp_het_vcf_hom_percentage = 0
        snp_het_no_vcf_percentage = 0
        snp_hom_vcf_hom_val = "0 (0.0%)"
        snp_hom_vcf_het_val = "0 (0.0%)"
        snp_hom_no_vcf_val = "0 (0.0%)"
        snp_het_vcf_hom_val = "0 (0.0%)"
        snp_het_no_vcf_val = "0 (0.0%)"
    else:
        snp_hom_vcf_hom_percentage = round(
            (snp_hom_vcf_hom / total_discordant) * 100, 2
        )
        snp_hom_vcf_hom_val = (
            str(snp_hom_vcf_hom) + " (" + str(snp_hom_vcf_hom_percentage) + "%)"
        )
        snp_hom_vcf_het_percentage = round(
            (snp_hom_vcf_het / total_discordant) * 100, 2
        )
        snp_hom_vcf_het_val = (
            str(snp_hom_vcf_het) + " (" + str(snp_hom_vcf_het_percentage) + "%)"
        )
        snp_hom_no_vcf_percentage = round((snp_hom_no_vcf / total_discordant) * 100, 2)
        snp_hom_no_vcf_val = (
            str(snp_hom_no_vcf) + " (" + str(snp_hom_no_vcf_percentage) + "%)"
        )
        snp_het_vcf_hom_percentage = round(
            (snp_het_vcf_hom / total_discordant) * 100, 2
        )
        snp_het_vcf_hom_val = (
            str(snp_het_vcf_hom) + " (" + str(snp_het_vcf_hom_percentage) + "%)"
        )
        snp_het_no_vcf_percentage = round((snp_het_no_vcf / total_discordant) * 100, 2)
        snp_het_no_vcf_val = (
            str(snp_het_no_vcf) + " (" + str(snp_het_no_vcf_percentage) + "%)"
        )

    values_with_percentages = [
        informational_val,
        noninformational_val,
        homozygous_ref_val,
        homozygous_alt_val,
        heterozygous_val,
        total_concordant_val,
        total_true_concordant_val,
        total_discordant_val,
        hom_ref_conc_val,
        hom_alt_conc_val,
        het_conc_whole_val,
        het_conc_true_val,
        het_conc_rev_val,
        snp_hom_vcf_hom_val,
        snp_hom_vcf_het_val,
        snp_hom_no_vcf_val,
        snp_het_vcf_hom_val,
        snp_het_no_vcf_val,
    ]

    percentages = [
        informational_percentage,
        noninformational_percentage,
        homozygous_ref_percentage,
        homozygous_alt_info_percentage,
        heterozygous_percentage,
        total_concordant_percentage,
        total_true_concordant_percentage,
        total_discordant_percentage,
        hom_ref_concordant_percentage,
        hom_alt_concordant_percentage,
        het_conc_whole_percentage,
        het_conc_true_percentage,
        het_conc_rev_percentage,
        snp_hom_vcf_hom_percentage,
        snp_hom_vcf_het_percentage,
        snp_hom_no_vcf_percentage,
        snp_het_vcf_hom_percentage,
        snp_het_no_vcf_percentage,
    ]
    return values_with_percentages, percentages


def generate_animal_stats(
    animal_basename,
    homozygous,
    het_forward,
    het_reverse,
    dashes,
    discordant,
    log_array,
    logfile,
    file_name,
):
    """
    Computes the shared marker, informational position, and concordancy statistics for a single animal.
    :param animal_basename: name of the current animal (without .1 or _allele)
    :param homozygous: df containing homozygous matches
    :param het_forward: df containing heterozygous forward matches
    :param het_reverse: df containing heterozygous reverse matches
    :param dashes: df containing dash matches
    :param discordant: df containing discordant matches (plus --/AA and AA/--)
    :param log_array: array of log statements if we exit here
    :param logfile: logfile
    :param file_name: name of snp panel file for logging purposes
    :return: list of stats for a single animal
    0: informational positions count
    1: non-informational positions count
    2: homozygous reference informational number
    3: homozygous alternate informational number
    4: heterozygous informational number
    5: total concordant (het reversal allowed)
    6: true concordant (no het reversal allowed)
    7: total discordant
    8: homozygous reference concordant
    9: homozygous alternate concordant
    10: panel hom and vcf hom (discordant AA vs BB or BB vs AA)
    11: panel hom and vcf het (discordant AA vs AB or BB vs AB)
    12: panel hom distinct from vcf (discordant CC vs AA/AB/BB)
    13: heterozygous concordant (whole)
    14: heterozygous concordant (true only)
    15: heterozygous concordant (reverse only)
    16: panel het to vcf hom (discordant AB vs AA or AB vs BB)
    17: panel het distinct from vcf (discordant BC/CD vs AA/AB/BB)
    """
    print("Calculating statistics for animal " + animal_basename)
    total_homozygous = len(homozygous.index)
    total_heterozygous_forward = len(het_forward.index)
    total_heterozygous_reverse = len(het_reverse.index)
    total_dashes_in_both = len(dashes.index)
    total_discordant_including_noninfo = len(discordant.index)

    # Create dict containing only informational discordant positions
    snp_animal_name = animal_basename + ".1"
    vcf_animal_name = animal_basename + "_allele"
    disc_noninformational = discordant.loc[
        (discordant[snp_animal_name] == "--") | (discordant[vcf_animal_name] == "--")
    ]
    disc_informational = discordant.loc[
        (discordant[snp_animal_name] != "--") & (discordant[vcf_animal_name] != "--")
    ]

    total_disc_informational = len(disc_informational.index)
    total_disc_noninformational = len(disc_noninformational.index)

    # Calculate informational content
    # hom + hetF + hetR + discordant - disc(AA/--) - disc(--/AA)
    informational_number = (
        total_homozygous
        + total_heterozygous_forward
        + total_heterozygous_reverse
        + total_disc_informational
    )  # 0

    # Calculuate non-informational content
    # dashes + disc(AA/--) + disc(--/AA)
    noninformational_number = total_dashes_in_both + total_disc_noninformational  # 1

    # Of Informational
    # Calculate homozygous reference, alternate, and heterozygous counts
    ref_alt_count_list = count_ref_or_alt(homozygous, snp_animal_name)
    homozygous_ref_informational_number = ref_alt_count_list[0]  # 2
    homozygous_alt_informational_number = ref_alt_count_list[1]  # 3
    heterozygous_informational_number = (
        total_heterozygous_forward + total_heterozygous_reverse
    )  # 4
    if total_homozygous != (
        homozygous_ref_informational_number + homozygous_alt_informational_number
    ):
        message = "Error determining the number of homozygous reference versus alt matching alleles"
        print(message)
        log_array.append(message)
        atexit.register(simple_log, log_array, file_name, logfile)
        exit()

    # Concordance statistics
    total_concordant_number = (
        total_homozygous + total_heterozygous_forward + total_heterozygous_reverse
    )  # 5
    total_true_concordant_number = total_homozygous + total_heterozygous_forward  # 6
    total_het_reversal_number = total_heterozygous_reverse
    total_discordant_number = total_disc_informational  # 7

    # Concordance statistics broken down by homozygous/heterozygous
    homozygous_ref_concordant_number = homozygous_ref_informational_number  # 8
    homozygous_alt_concordant_number = homozygous_alt_informational_number  # 9
    discordant_hom, discordant_het = count_discordant_types(
        disc_informational, snp_animal_name, vcf_animal_name
    )
    snp_hom_vcf_hom = discordant_hom["hom_to_hom_vcf"]  # 10
    snp_hom_vcf_het = discordant_hom["hom_to_het_vcf"]  # 11
    snp_hom_no_vcf = discordant_hom["hom_to_not_vcf"]  # 12
    heterozygous_concordant_whole_number = heterozygous_informational_number  # 13
    heterozygous_concordant_true_number = total_heterozygous_forward  # 14
    heterozygous_concordant_reverse_number = total_heterozygous_reverse  # 15
    snp_het_vcf_hom = discordant_het["het_to_hom_vcf"]  # 16
    snp_het_no_vcf = discordant_het["het_to_not_vcf"]  # 17

    stats_list = [
        informational_number,
        noninformational_number,
        homozygous_ref_informational_number,
        homozygous_alt_informational_number,
        heterozygous_informational_number,
        total_concordant_number,
        total_true_concordant_number,
        total_discordant_number,
        homozygous_ref_concordant_number,
        homozygous_alt_concordant_number,
        snp_hom_vcf_hom,
        snp_hom_vcf_het,
        snp_hom_no_vcf,
        heterozygous_concordant_whole_number,
        heterozygous_concordant_true_number,
        heterozygous_concordant_reverse_number,
        snp_het_vcf_hom,
        snp_het_no_vcf,
    ]
    # Compute percentages
    stats_values_with_percentages, stats_percentages = calc_percentages(stats_list)
    str_stats_list = [str(s) for s in stats_list]
    stats_text = "\t".join(str_stats_list)
    return stats_list, stats_values_with_percentages


def generate_statistics(
    comparison_dict,
    full_samples_list,
    nonvcf_samples,
    vcfonly_positions,
    snp_panel_dict,
    vcf_pos_number,
    panel_vcf_pos_number,
    out_pfx,
    logfile,
    file_name,
):
    """
    Generates the concordance statistics per animal, then per chromosome, filter value, annotation status; for
    homozygous and heterozygous matches
    :param comparison_dict: dict containing an array of dataframes per animal, from compare_genotypes:
    1. correct homozygous matches
    2. correct (forward) heterozygous matches
    3. reverse heterozygous matches
    4. dashes in SNP file
    5. discordant positions
    6. number of filtered positions, if filtering was done (otherwise None), given as a string
    :param full_samples_list: list of all samples in the user input snp panel
    :param nonvcf_samples: list of samples that are only in the SNP panel and not the vcf file
    :param vcfonly_positions: list of positions that are only in the vcf file and not found in the SNP panel
    :param snp_panel_dict: dict containing all the snp panels
    :param vcf_pos_number: the number of positions in the full vcf
    :param panel_vcf_pos_number: number of positions shared by panel and vcf
    :param out_pfx: output prefix for metastats file
    :param logfile: logfile
    :param file_name: snp file name for logging purposes
    :return: dict of statistics per animal
    """
    log_array = []
    timestr = time.strftime("%H:%M:%S")
    logfile_text = timestr + " ..... Computing concordance statistics "
    log_array.append(logfile_text)
    # Create the metastatistics file first
    metastatistics_filepath = out_pfx + "_metastatistics.txt"
    sample_panel_list = next(iter(snp_panel_dict.values()))
    sample_panel = sample_panel_list[0]
    metastats_numbers = get_metastats_numbers(
        full_samples_list,
        comparison_dict,
        nonvcf_samples,
        sample_panel,
        vcf_pos_number,
        panel_vcf_pos_number,
    )
    with open(metastatistics_filepath, "w") as metastats:
        metastats.write(
            """There are {0} animals in SNP panel.
Of these, {1} ({2}%) are present in the genotype file.
Animals in SNP panel and genotype file: {3}
Animals in SNP panel but not genotype file: {4}""".format(
                metastats_numbers[0],
                metastats_numbers[1],
                metastats_numbers[2],
                metastats_numbers[3],
                metastats_numbers[4],
            )
        )
        metastats.write(
            """
There are {0} probes found in the SNP panel and {1} positions found in the genotype VCF file.
{2} positions are found in both the SNP panel and genotype file, based on genomic location.
This shared value represents {3}% of all SNP panel positions and {4}% of all genotype file positions.
""".format(
                metastats_numbers[5],
                metastats_numbers[6],
                metastats_numbers[7],
                metastats_numbers[8],
                metastats_numbers[9],
            )
        )
        if metastats_numbers[10] is not None:
            filtered_num = metastats_numbers[10]
            filtered_num = str(filtered_num)
            metastats.write(
                """
User-directed filtering was applied, which removed {} positions.
""".format(
                    filtered_num
                )
            )
    metastats.close()

    # Compute the statistics for each animal
    stats_per_animal_dict = {}
    for animal in sorted(comparison_dict.keys()):
        homozygous_dict = comparison_dict[animal][0]
        heterozygous_forward = comparison_dict[animal][1]
        heterozygous_reverse = comparison_dict[animal][2]
        matching_dashes = comparison_dict[animal][3]
        discordant_df = comparison_dict[animal][4]
        filtered_count = comparison_dict[animal][5]
        panel_name = animal + ".1"
        unknown_pos_df = snp_panel_dict[panel_name][
            1
        ]  # not sure what to do with this here... add to animal-spc stats?
        (
            animal_specific_stats,
            animal_specific_stats_with_percentages,
        ) = generate_animal_stats(
            animal,
            homozygous_dict,
            heterozygous_forward,
            heterozygous_reverse,
            matching_dashes,
            discordant_df,
            log_array,
            logfile,
            file_name,
        )
        animal_specific_stats.insert(0, panel_vcf_pos_number)
        animal_specific_stats_with_percentages.insert(0, panel_vcf_pos_number)
        stats_per_animal_dict.update({animal: animal_specific_stats_with_percentages})

    # Write animal stats to file (tab formatted)
    # make an output dataframe with the correct rows
    output_index = [
        "Markers shared",
        "Informational",
        "Non-informational",
        "Homozygous reference",
        "Homozygous alternate",
        "Heterozygous",
        "Total concordant",
        "True concordant",
        "Total discordant",
        "Homozygous reference concordant",
        "Homozygous alternate concordant",
        "Discordant: panel hom. & vcf hom.",
        "Discordant: panel hom. & vcf het.",
        "Discordant: panel hom. & distinct from vcf",
        "Heterozygous concordant (total)",
        "Heterozygous concordant (true only)",
        "Heterozygous concordant (gt reversed only)",
        "Discordant: panel het. & vcf hom.",
        "Discordant: panel het. & distinct from vcf",
    ]
    output_stats_df = pd.DataFrame(index=output_index)
    for each_animal in stats_per_animal_dict:
        output_stats_df[each_animal] = stats_per_animal_dict[each_animal]
    if logfile is not None:
        logging = simple_log(log_array, file_name, logfile)
    return output_stats_df, logfile


def extract_discordant_positions(comparison_dict, logfile, file_name):
    """
    Writes discordant positions to a separate file, if requested with the --extract-discordant flag
    :param comparison_dict: dict containing an array of dataframes per animal, from compare_genotypes:
    1. correct homozygous matches
    2. correct (forward) heterozygous matches
    3. reverse heterozygous matches
    4. dashes in SNP file
    5. discordant positions
    6. number of filtered positions, if filtering was done (otherwise None), given as a string
    :param logfile: logfile
    :param file_name: name of snp panel file for logging purposes
    :return: True value, logfile
    """
    log_array = []
    timestr = time.strftime("%H:%M:%S")
    logfile_text = timestr + " ..... Writing discordant genotypes to file "
    log_array.append(logfile_text)
    discordant_array = []
    for animal in comparison_dict:
        discordant_df = comparison_dict[animal][4].copy()
        animal_allele = animal + "_allele"
        animal_p1 = animal + ".1"
        discordant_df_nodash = discordant_df[
            discordant_df[animal_allele] != "--"
        ].copy()
        discordant_df_nodash.rename(
            columns={
                animal_allele: "VCF genotype",
                animal_p1: "Panel genotype",
                "#CHROM": "CHROM",
            },
            inplace=True,
        )
        discordant_df_nodash["Sample"] = animal
        discordant_array.append(discordant_df_nodash)
    concat_discordant = pd.concat(discordant_array)
    # concat_discordant['POS'].round()
    reduce_discordant = concat_discordant[
        ["CHROM", "POS", "Name", "Sample", "VCF genotype", "Panel genotype"]
    ].copy()
    sorted_discordant = reduce_discordant.sort_values(
        by=["CHROM", "POS", "Name", "Sample"], ascending=True
    )
    disc_output_name = file_name.replace(".txt", "_discordant.txt")
    sorted_discordant.to_csv(disc_output_name, index=False, sep="\t")
    if logfile is not None:
        logging = simple_log(log_array, file_name, logfile)
    return True, logfile


def split_output_text(stats_df):
    """
    Splits the text in the statistics df so that the explanatory text can be inserted easily for the tabular and pretty
    output files
    :param stats_df: statistics dataframe
    :return: list of dataframes containing the separated statistics
    """
    markers_shared = stats_df.loc[["Markers shared"], :]
    shared_breakdown = stats_df.loc[["Informational", "Non-informational"], :]
    informational_breakdown = stats_df.loc[
        ["Homozygous reference", "Homozygous alternate", "Heterozygous"], :
    ]
    concordant_total = stats_df.loc[
        ["Total concordant", "True concordant", "Total discordant"], :
    ]
    concordant_rest = stats_df.loc[
        [
            "Homozygous reference concordant",
            "Homozygous alternate concordant",
            "Discordant: panel hom. & vcf hom.",
            "Discordant: panel hom. & vcf het.",
            "Discordant: panel hom. & distinct from vcf",
            "Heterozygous concordant (total)",
            "Heterozygous concordant (true only)",
            "Heterozygous concordant (gt reversed only)",
            "Discordant: panel het. & vcf hom.",
            "Discordant: panel het. & distinct from vcf",
        ],
        :,
    ]
    split_frames = [
        markers_shared,
        shared_breakdown,
        informational_breakdown,
        concordant_total,
        concordant_rest,
    ]
    return split_frames


def stat_to_file(statistics_df, file_names, output_t, out_pfx, logfile, file_name):
    """
    Writes the statistics data to a file, in either tab-separated or "pretty" format
    :param statistics_df: dataframe containing all animal statistics
    :param file_names: list of [snp file basename, vcf file basename] for metastats file
    :param output_t: output type, either 'basic', 'tabular' or 'pretty'
    :param out_pfx: name of the output file
    :param logfile: logfile
    :param file_name: snp file name for logging purposes
    :return: True
    """
    log_array = []
    timestr = time.strftime("%H:%M:%S")
    logfile_text = timestr + " ..... Writing statistics to file "
    log_array.append(logfile_text)
    out_filename = out_pfx + "_statistics.txt"
    # Make a basic output option (just write df to file)
    if output_t == "basic":
        statistics_df.to_csv(out_filename, sep="\t")
    # Make a tabular format but with more information than basic
    elif output_t == "tabular":
        with open(out_filename, "w") as stat_out:
            stat_out.write(
                """Statistics for SNP panel {0} and genotype VCF file {1}
Only animals with markers in both the SNP panel file and VCF file are shown. 
See the associated metastatistics file for more information.\n
""".format(
                    file_names[0], file_names[1]
                )
            )
        stat_out.close()
        output_text_split = split_output_text(statistics_df)
        markers_shared = output_text_split[0]
        markers_shared.to_csv(out_filename, sep="\t", mode="a", header=True)
        with open(out_filename, "a") as stat_out:
            stat_out.write("\nOf the shared markers:\n")
        stat_out.close()
        shared_breakdown = output_text_split[1]
        shared_breakdown.to_csv(out_filename, sep="\t", mode="a", header=True)
        with open(out_filename, "a") as stat_out:
            stat_out.write("\nOf markers with genotype information:\n")
        stat_out.close()
        informational_breakdown = output_text_split[2]
        informational_breakdown.to_csv(out_filename, sep="\t", mode="a", header=True)
        with open(out_filename, "a") as stat_out:
            stat_out.write("\nConcordancy analysis:\n")
        stat_out.close()
        concordant_total = output_text_split[3]
        concordant_total.to_csv(out_filename, sep="\t", mode="a", header=True)
        with open(out_filename, "a") as stat_out:
            stat_out.write("\n")
        stat_out.close()
        concordant_rest = output_text_split[4]
        concordant_rest.to_csv(out_filename, sep="\t", mode="a", header=True)
    # Make a pretty format that is not tabular
    elif output_t == "pretty":
        with open(out_filename, "w") as stat_out:
            stat_out.write(
                """Statistics for SNP panel {0} and genotype VCF file {1}
Only animals with markers in both the SNP panel file and VCF file are shown. 
See the associated metastatistics file for more information.\n\nMarkers shared:\n
""".format(
                    file_names[0], file_names[1]
                )
            )
            output_text_split = split_output_text(statistics_df)
            marker_string = output_text_split[0].to_string(justify="center")
            stat_out.write(marker_string)
            stat_out.write("\n\nOf the shared markers:\n")
            breakdown_string = output_text_split[1].to_string(justify="center")
            stat_out.write(breakdown_string)
            stat_out.write("\n\nOf markers with genotype information:\n")
            informational_string = output_text_split[2].to_string(justify="center")
            stat_out.write(informational_string)
            stat_out.write("\n\nConcordancy analysis:\n")
            concordant_total_string = output_text_split[3].to_string(justify="center")
            stat_out.write(concordant_total_string)
            stat_out.write("\n")
            concordant_rest_string = output_text_split[4].to_string(justify="center")
            stat_out.write(concordant_rest_string)
        stat_out.close()

    else:
        message = "Unrecognized file output type"
        log_array.append(message)
    # Add percentage explanation
    with open(out_filename, "a") as additional_stat_out:
        additional_stat_out.write(
            """

Percentages are calculated as follows:
Informational, Non-informational: percentage of Markers shared
Homozygous reference, Homozygous alternate, Heterozygous: percentage of Informational
Total concordant, True concordant, Total discordant: percentage of Informational
Homozygous and heterozygous concordant: percentage of Total concordant
Discordant subtypes: percentage of Total discordant"""
        )
    additional_stat_out.close()
    if logfile is not None:
        logging = simple_log(log_array, file_name, logfile)
    return True


# Old options for testing
# input_name = "G_CCGP_downsampled.txt"
# input_name = "G_CCGP_incorrect_homozygous2.txt"
# input_name = "G_CCGP_very_downsampled.txt"
# input_name = "G_CCGP_wrongrow1.txt"
# variant_dir_path = "variant_position_files"
# assembly_input = "ARS-UCD1.2"
# log_input = None
# file_type_input = 'affymetrix'
# genotype_input = os.path.join('genotyping', 'SNPs_reduced_7M.recode.vcf')
# output_type = 'pretty'
# output_name = 'test_statistics.txt'
# filter_conditions = True
# make sure these get default values if they are not set
# qual_crit = 0
# filter_crit = ['PASS']

###################################################################################3
def concordance_analysis(
    variant_dir_path,
    assembly_input,
    snp_path,
    file_type_input,
    verbose_logging,
    genotype_input,
    filter_bool,
    qual_crit,
    filter_crit,
    output_name,
    output_type,
    extract_discordant,
    n_threads,
    species,
):
    # Get user panel files and variant files, and check that these exist
    variant_file_list = gt.retrieve_all_variant_files(
        variant_dir_path, assembly_input, species
    )
    if not variant_file_list:
        exit("No variant conversion files in " + variant_dir_path)

    # Get SNP panel basename (file)
    panel_path, snp_file = gt.retrieve_user_input_file(snp_path)

    # Get log file name if verbose logging is true
    if verbose_logging:
        timestr = time.strftime("%Y%m%d%H%M%S")
        log_suffix = "-" + timestr + ".log"
        log_input = get_logname(log_suffix, snp_file)
    else:
        log_input = None

    # Check for threading
    can_we_thread = gt.check_user_platform()

    # Perform concordance analysis
    # Check inputs
    # Check file format of input file
    format_check_val, logfile_out1, panel_df, variant_file_df = gt.check_input_snp_panel(
        panel_path,
        snp_file,
        file_type_input,
        variant_file_list,
        assembly_input,
        variant_dir_path,
        can_we_thread,
        n_threads,
        log_input,
        species,
    )
    # Format check has passed
    # Create a dict of per-animal dataframes that have BLAST chromosome & position info from var file
    position_dict, logfile_out2 = gt.get_snp_panel_positional_info(
        panel_df, variant_file_df, can_we_thread, n_threads, logfile_out1, snp_file
    )
    # First subsample vcf file using the panel-specific position info
    samples_list = list(panel_df.columns)
    samples_list.remove("Name")
    if file_type_input != "AFFY":
        update_col_dict = {}
        for colname in samples_list:
            new_name = colname + ".1"
            update_col_dict.update({colname: new_name})
        panel_df.rename(columns=update_col_dict, inplace=True)
    (
        subsampled_vcf,
        full_vcf,
        samples_not_in_vcf,
        vcf_only_positions,
        logfile_out3,
    ) = subsample_vcf_by_position(
        genotype_input,
        position_dict,
        samples_list,
        filter_bool,
        qual_crit,
        filter_crit,
        logfile_out2,
        snp_file,
    )
    # Use scikit_allel to read vcf file and get genotype info (for subsampled vcf file)
    converted_dataframe, logfile_out4 = get_genotypes_from_vcf(
        genotype_input,
        subsampled_vcf,
        full_vcf,
        can_we_thread,
        n_threads,
        logfile_out3,
        snp_file,
    )

    # Compare the user SNP panel values to the VCF values on a per-animal basis
    compared_dict, logfile_out5 = compare_genotypes(
        converted_dataframe,
        position_dict,
        samples_not_in_vcf,
        filter_bool,
        qual_crit,
        filter_crit,
        can_we_thread,
        n_threads,
        logfile_out4,
        snp_file,
    )

    # Perform statistics on the merged dataframes (homozygous, heterozygous-forward, heterozyous-reversed, non-matching)
    full_vcf_num_pos = len(full_vcf.index)
    shared_pos_number = len(subsampled_vcf.index)
    statistics_df, logfile_out6 = generate_statistics(
        compared_dict,
        samples_list,
        samples_not_in_vcf,
        vcf_only_positions,
        position_dict,
        full_vcf_num_pos,
        shared_pos_number,
        output_name,
        logfile_out5,
        snp_file,
    )

    # Write discordant positions to file, if specified
    if extract_discordant:
        write_disc, logfile_out7 = extract_discordant_positions(
            compared_dict, logfile_out6, snp_file
        )
    else:
        logfile_out7 = logfile_out6
    # Get SNP panel and VCF file basenames for stats
    snp_basename = os.path.basename(snp_path)
    vcf_basename = os.path.basename(genotype_input)
    input_names = [snp_basename, vcf_basename]
    # Write statistics to file
    concordance_out = stat_to_file(
        statistics_df, input_names, output_type, output_name, logfile_out7, snp_file
    )
    print("Statistics complete")
    return concordance_out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Determines concordance between SNP panel genotypes and VCF file "
        "genotypes"
    )
    parser.add_argument(
        "--snp-panel",
        type=str,
        required=True,
        help="Name or path to the SNP panel file. Can be Illumina or Affymetrix format",
    )
    parser.add_argument(
        "--panel-type",
        type=str,
        required=True,
        choices=["TOP", "FWD", "AB", "PLUS", "DESIGN", "LONG", "AFFY"],
        help="Type of panel file: 'TOP', 'FWD', 'AB', 'PLUS', 'DESIGN', 'LONG', 'AFFY'",
    )
    parser.add_argument(
        "--vcf-file",
        type=str,
        required=True,
        help="Name or path to the VCF file containing genotype information",
    )
    parser.add_argument(
        "--conversion",
        type=str,
        required=False,
        default="variant_position_files",
        help="Directory containing genotype conversion key files (default directory: variant_position_files)",
    )
    parser.add_argument(
        "--assembly",
        type=str,
        required=True,
        help="Assembly name (see README for all available choices)",
    )
    parser.add_argument(
        "--species",
        required=True,
        type=str,
        help="Species name (see README for all available choices)",
    )
    parser.add_argument(
        "--filter-vcf",
        action="store_true",
        default=False,
        required=False,
        help="[Optional] Use this flag to specify that the VCF file should be filtered on QUAL and/or FILTER values. "
        "Must be used in conjunction with --qual and/or --filter parameters",
    )
    parser.add_argument(
        "--qual",
        type=float,
        required=False,
        default=None,
        help="[Optional] Only perform concordance analysis using variants with QUAL scores higher than this value. "
        "Must be used in conjunction with the --filter-vcf flag",
    )
    parser.add_argument(
        "--filter",
        type=str,
        required=False,
        default=None,
        help="[Optional] Only perform concordance analysis using variants with these value(s) in the FILTER field. "
        "Can be a single value or a comma-separated list. Must be used in conjunction with --filter-vcf flag",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=False,
        default="concordance",
        help="[Optional] Output prefix to append to '_metastatistics.txt' and '_statistics.txt' files (default: "
        "'concordance')",
    )
    parser.add_argument(
        "--output-type",
        type=str,
        required=False,
        choices=["basic", "tabular", "pretty"],
        default="tabular",
        help="[Optional] Type of output for statistics file: basic (tsv, stats only), tabular (tsv with extra info), "
        "pretty (nice formatting) (default: tabular)",
    )
    parser.add_argument(
        "-v",
        "--verbose-logging",
        action="store_true",
        required=False,
        default=False,
        help="[Optional] Write program steps and errors to a log file",
    )
    parser.add_argument(
        "--extract-discordant",
        action="store_true",
        required=False,
        default=False,
        help="[Optional] Write discordant positions to a file called [snp-panel]_discordant.txt",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=2,
        required=False,
        help="[Optional] Number of threads to use during conversion (default = 2)",
    )
if __name__ == "__main__":

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    # Error check input values
    species_assembly_correlation = errors.assembly_species_error(
        args.conversion, args.assembly, args.species
    )
    if not species_assembly_correlation:
        exit(
            "The assembly "
            + args.assembly
            + " is incompatible with species "
            + args.species
            + ". See README for acceptable "
            "assemblies."
        )

    if args.filter_vcf:
        filter_or_qual_set = False
        if args.qual is not None:
            filter_or_qual_set = True
        if args.filter is not None:
            filter_or_qual_set = True
        if not filter_or_qual_set:
            exit(
                "The --filter-vcf flag was invoked, but no FILTER or QUAL filter values were provided"
            )

    try:
        vcf_test = allel.read_vcf(args.vcf)
    except RuntimeError as rt_inst:
        print("VCF file could not be opened: " + str(rt_inst))
        exit()

    # Parse filter value(s)
    if args.filter is not None:
        filter_str = args.filter
        filter_vals = filter_str.split(",")
    else:
        filter_vals = []

    output = concordance_analysis(
        args.conversion,
        args.assembly,
        args.snp_panel,
        args.panel_type,
        args.verbose_logging,
        args.vcf,
        args.filter_vcf,
        args.qual,
        filter_vals,
        args.output,
        args.output_type,
        args.extract_discordant,
        args.threads,
        args.species,
    )
