#!/usr/bin/env python3
import pandas as pd
import sys
import lib.common_vcf_libs as gt
import time
import datetime
from lib.make_logs import simple_log, get_logname
from lib.variant_file_finder import var_match
import os
from lib.file_parsing import uncompressing
import numpy as np
from functools import reduce
import argparse
import gzip


# This module creates a vcf file from a SNP panel file


def get_ref_alt_values(
    position_dict,
    variant_file_list,
    logfile,
    snp_file,
    conversion_dir,
    assembly,
    species,
):
    """
    Creates the start of an output dataframe with chromosome and position information, and ref and alt data from
    variant position files
    :param position_dict: hash of animal-specific dataframes containing position info
    :param variant_file_list: list of variant files
    :param logfile: logfile
    :param snp_file: file name for logging
    :param conversion_dir: dir containing variant_position_files
    :param assembly: assembly name
    :param species: species name
    :return: dataframe with position and ALT/REF allele info
    """
    generic_position_df_list = next(iter(position_dict.values()))
    generic_position_df = generic_position_df_list[0]
    # Read in var file ref and alt data
    converting_file = True
    var_match_list_out, mod_verbose_log, reg_alt_bool_dict = var_match(
        variant_file_list,
        conversion_dir,
        generic_position_df,
        snp_file,
        converting_file,
        assembly,
        species,
    )
    concat_df_list = []
    for var_file in var_match_list_out:
        if reg_alt_bool_dict[var_file] == "reg":
            alt_bool = False
        else:
            alt_bool = True
        var_species_path = os.path.join(conversion_dir, species)
        var_assembly_path = os.path.join(var_species_path, assembly)
        filepath = os.path.join(var_assembly_path, var_file)
        header_count = uncompressing(filepath)
        whole_var_df = pd.read_csv(
            filepath, header=0, skiprows=header_count, compression="infer"
        )

        df_ref = whole_var_df[whole_var_df["VCF"] == "REF"].copy()
        df_alt = whole_var_df[whole_var_df["VCF"] == "ALT"].copy()
        merged_df_ref_alt = pd.merge(df_ref, df_alt, on=["marker_name", "alt_marker_name"])
        if alt_bool is True:
            ref_alt_df = merged_df_ref_alt[["alt_marker_name", "PLUS_x", "PLUS_y"]].copy()
            add_ref_alt_df = pd.merge(
                generic_position_df,
                ref_alt_df,
                left_on="Name",
                right_on="alt_marker_name",
                how="outer",
            )
            merged_df_ref_alt.drop(columns=["marker_name"], inplace=True)
            merged_df_ref_alt.rename(
                columns={"alt_marker_name": "marker_name"}, inplace=True
            )
        else:
            ref_alt_df = merged_df_ref_alt[["marker_name", "PLUS_x", "PLUS_y"]].copy()
            add_ref_alt_df = pd.merge(
                generic_position_df,
                ref_alt_df,
                left_on="Name",
                right_on="marker_name",
                how="outer",
            )
        combined_df = add_ref_alt_df[
            ["Name", "BLAST_chromosome", "BLAST_position", "PLUS_x", "PLUS_y"]
        ].copy()
        combined_df.rename(
            columns={
                "Name": "ID",
                "BLAST_chromosome": "#CHROM",
                "BLAST_position": "POS",
                "PLUS_x": "REF",
                "PLUS_y": "ALT",
            },
            inplace=True,
        )
        arranged_cols = ["#CHROM", "POS", "ID", "REF", "ALT"]
        rearranged_df = combined_df[arranged_cols].copy()
        rearranged_df.dropna(inplace=True)
        concat_df_list.append(rearranged_df)
    merged_df_final = reduce(
        lambda left, right: pd.merge(left, right, how="outer", on=["#CHROM", "POS", "ID", 'REF', 'ALT']), concat_df_list
    )
    return merged_df_final, logfile


def create_vcf_genotypes(ref_alt_df, position_dict, discard_snp, filename, logfile):
    """
    Converts bases in each animal dataframe to a numeric genotype format: 0/0 hom ref, 0/1 het, 1/1 hom alt.
    All genotypes are unphased. Output: dataframe with initial vcf info, FORMAT column, and one col per animal with
    only genotype information
    :param ref_alt_df: initial vcf dataframe containing position, ID, ref and alt values
    :param position_dict: dict containing animal dataframes, one df per animal
    :param discard_snp: (bool) discard snp marker names as ID columns (becomes '.')
    :param filename: name of snp file for logging
    :param logfile: logfile
    :return: vcf dataframe (CHROM, POS, ID, REF, ALT, FORMAT, animals...)
    """
    log_array = []
    message = "Converting genotype data to VCF style"
    log_array.append(message)
    vcf_df_list = []
    for animal in position_dict:
        animal_df = position_dict[animal][0]

        # split sample col and assess each value
        panelname_1 = animal + "_A1"
        panelname_2 = animal + "_A2"
        # Add ref and alt code values
        # 0/0 = 10
        # 1/1 = 11
        # 0/1 = 12
        # ./. = 0 (NA and indel)
        animal_refalt = pd.merge(
            left=ref_alt_df,
            right=animal_df,
            how="outer",
            left_on=["#CHROM", "POS", "ID"],
            right_on=["BLAST_chromosome", "BLAST_position", "Name"],
        )
        animal_refalt[panelname_1], animal_refalt[panelname_2] = zip(
            *animal_refalt[animal].apply(lambda x: list(x))
        )
        animal_refalt["GT_REF"] = np.where(
            (animal_refalt["REF"] == animal_refalt[panelname_1])
            & (animal_refalt["REF"] == animal_refalt[panelname_2]),
            10,
            np.nan,
        )
        animal_refalt["GT_ALT"] = np.where(
            (animal_refalt["ALT"] == animal_refalt[panelname_1])
            & (animal_refalt["ALT"] == animal_refalt[panelname_2]),
            11,
            np.nan,
        )
        animal_refalt["GT_HET"] = np.where(
            (
                (animal_refalt["REF"] == animal_refalt[panelname_1])
                & (animal_refalt["ALT"] == animal_refalt[panelname_2])
            )
            | (
                (animal_refalt["REF"] == animal_refalt[panelname_2])
                & (animal_refalt["ALT"] == animal_refalt[panelname_1])
            ),
            12,
            np.nan,
        )
        animal_refalt["GT_NA"] = np.where(
            (animal_refalt[panelname_1] == "-") & (animal_refalt[panelname_2] == "-"),
            0,
            np.nan,
        )
        animal_refalt["GT_INDEL"] = np.where(
            ((animal_refalt[panelname_1] == "I") & (animal_refalt[panelname_2] == "I"))
            | (
                (animal_refalt[panelname_1] == "D")
                & (animal_refalt[panelname_2] == "D")
            ),
            0,
            np.nan,
        )
        animal_refalt.fillna(0, inplace=True)
        animal_refalt["GT_NUM"] = (
            animal_refalt["GT_REF"]
            + animal_refalt["GT_ALT"]
            + animal_refalt["GT_HET"]
            + animal_refalt["GT_NA"]
            + animal_refalt["GT_INDEL"]
        )
        # make df with rows where GT_num are >= 13
        refalt_fail = animal_refalt[animal_refalt["GT_NUM"] >= 13]
        if refalt_fail.empty:
            pass
        else:
            print(
                "Something went wrong during genotype conversion. Printing problem values: "
            )
            print(refalt_fail)
            exit()
        animal_name = animal.replace(".1", "")
        animal_refalt[animal_name] = animal_refalt["GT_NUM"].map(
            {10: "0/0", 11: "1/1", 12: "0/1", 0: "./."}
        )
        animal_refalt.drop(
            [
                "Name",
                animal,
                "BLAST_chromosome",
                "BLAST_position",
                panelname_1,
                panelname_2,
                "GT_REF",
                "GT_ALT",
                "GT_HET",
                "GT_NA",
                "GT_INDEL",
                "GT_NUM",
            ],
            axis=1,
            inplace=True,
        )
        vcf_df_list.append(animal_refalt)
    vcf_df_working = reduce(
        lambda left, right: pd.merge(
            left, right, on=["#CHROM", "POS", "ID", "REF", "ALT"], how="outer"
        ),
        vcf_df_list,
    )
    vcf_df_working.insert(loc=5, column="QUAL", value=".")
    vcf_df_working.insert(loc=6, column="FILTER", value=".")
    vcf_df_working.insert(loc=7, column="INFO", value=".")
    vcf_df_working.insert(loc=8, column="FORMAT", value="GT")
    # Replace ID values if keep_snp is false
    if discard_snp:
        vcf_df_working.ID = "."
    else:
        pass
    vcf_out = vcf_df_working.astype({"POS": "int"})

    ## TODO: sort values in CHROM column
    if logfile:
        logfile = simple_log(log_array, filename, logfile)
    return vcf_out, logfile


def write_vcf_file(
    vcf_df, output_name, assembly, species, compress, commandline, filename, logfile
):
    """
    Write the vcf dataframe to an output file - function also adds preamble
    :param vcf_df: vcf dataframe
    :param output_name: output name for writing (may be empty)
    :param assembly: assembly name
    :param species: species name
    :param compress: compress the file with gzip
    :param commandline: command line call for writing to vcf file metadata
    :param filename: snp filename for writing
    :param logfile: logfile
    :return: file written to output
    """
    log_array = []
    message = "Writing VCF file"
    log_array.append(message)
    if not output_name:
        split_name = os.path.splitext(filename)
        output = split_name[0] + ".vcf"
    else:
        output = output_name
    if compress:
        output = output + ".gz"
        compression_type = "gzip"
    else:
        pass
        compression_type = None
    # Delete file if it already exists
    if os.path.exists(output):
        os.remove(output)
    # Get date info
    currentDT = datetime.datetime.now()
    date_info = currentDT.strftime("%Y%m%d")
    # Get command values
    commandline_str = " ".join(commandline)
    comm = commandline_str.replace("./", "")
    header_string = """##fileformat=VCFv4.2
##filedate={}
##source=vcf_generator.py
##reference={}
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##vcf_generator_commandline.vcf_generator=<{}>
""".format(
        date_info, assembly, comm
    )
    if compress:
        with gzip.open(output, "ab") as output_fh:
            output_fh.write(header_string.encode())
        output_fh.close()
    else:
        with open(output, "a") as output_fh:
            output_fh.write(header_string)
        output_fh.close()
    # Write DF to output
    vcf_df.to_csv(output, sep="\t", index=False, mode="a", compression=compression_type)
    if logfile:
        logfile = simple_log(log_array, filename, logfile)
    return output, logfile


def vcf_generator(
    snp_path,
    conversion_dir,
    species,
    assembly,
    panel_type,
    output_name,
    n_threads,
    verbose_logging,
    compress,
    discard_snp,
):
    """
    Uses data from a SNP panel file to generate a VCF file
    :param snp_path: path to (or name of) SNP panel file
    :param conversion_dir: directory containing variant files (variant_position_files)
    :param species: species name
    :param assembly: assembly name
    :param panel_type: Type of panel file: 'TOP', 'FWD', 'AB', 'PLUS', 'DESIGN', 'LONG', 'affymetrix'
    :param output_name: full name of output file (default: input_basename.vcf)
    :param n_threads: number of threads to use, default: 2
    :param verbose_logging: (bool) whether to create log file
    :param compress: (bool) compress the final vcf file with gzip
    :param discard_snp: (bool) discard the snp maker names in the ID column
    :return: vcf file
    """
    # Get some initial variables
    # Check for threading
    can_we_thread = gt.check_user_platform()
    # Get SNP panel basename (file)
    panel_path, snp_file = gt.retrieve_user_input_file(snp_path)
    # Set logging info
    if verbose_logging:
        timestr = time.strftime("%Y%m%d%H%M%S")
        log_suffix = "-" + timestr + ".log"
        log_input = get_logname(log_suffix, snp_file)
    else:
        log_input = None

    # Get user panel files and variant files, and check that these exist
    variant_file_list = gt.retrieve_all_variant_files(conversion_dir, assembly, species)
    if not variant_file_list:
        exit("No variant conversion files in " + conversion_dir)

    # Check file format of input file - this also converts the panel_df to plus
    (
        format_check_val,
        logfile_out1,
        panel_df,
        variant_file_df,
    ) = gt.check_input_snp_panel(
        panel_path,
        snp_file,
        panel_type,
        variant_file_list,
        assembly,
        conversion_dir,
        can_we_thread,
        n_threads,
        log_input,
        species,
    )
    # Create a dict of per-animal dataframes that have BLAST chromosome & position info from var file
    position_dict, logfile_out2 = gt.get_snp_panel_positional_info(
        panel_df, variant_file_df, can_we_thread, n_threads, logfile_out1, snp_file
    )

    # Get ref and alt values from variant file dataframe (need to read in original variant files)
    ref_alt_df, logfile_out3 = get_ref_alt_values(
        position_dict,
        variant_file_list,
        logfile_out2,
        snp_file,
        conversion_dir,
        assembly,
        species,
    )

    # Convert to vcf genotypes
    vcf_df, logfile_out4 = create_vcf_genotypes(
        ref_alt_df, position_dict, discard_snp, snp_file, logfile_out3
    )

    # Write to vcf file
    commandline = sys.argv
    vcf_to_file, logfile_out5 = write_vcf_file(
        vcf_df,
        output_name,
        assembly,
        species,
        compress,
        commandline,
        snp_file,
        logfile_out4,
    )


if __name__ == "__main__":
    # Get user input with argparse
    parser = argparse.ArgumentParser(
        description="Converts SNP panel file to VCF format file"
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
        choices=["TOP", "FWD", "AB", "DESIGN", "LONG", "PLUS", "AFFY"],
        help="Type of file(s) expected: 'TOP', 'FWD', 'AB', 'DESIGN', 'LONG', 'PLUS', or 'AFFY' (Affymetrix)",
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
        "--output",
        required=False,
        help="Full output file name to store VCF results (default: [input basename].vcf",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=2,
        required=False,
        help="[Optional] Number of threads to use if conversion to PLUS format is required (default = 2)",
    )
    parser.add_argument(
        "-c",
        "--compress",
        action="store_true",
        default=False,
        help="Compress VCF output file with gzip (not bgzip)",
    )
    parser.add_argument(
        "--discard-snp",
        action="store_true",
        default=False,
        help="Set ID field values to '.' (default: ID field contains SNP marker names)",
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    vcf_generator(
        args.snp_panel,
        args.conversion,
        args.species,
        args.assembly,
        args.panel_type,
        args.output,
        args.threads,
        args.verbose_logging,
        args.compress,
        args.discard_snp,
    )
