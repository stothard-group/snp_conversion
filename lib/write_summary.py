#!/usr/bin/env python3
import pandas as pd
import lib.file_parsing as fp
import argparse
import sys
from collections import OrderedDict


# This module creates a summary file for each input file:
# File
# Sample	Total   input SNPs(%)	SNPs with data(%)	Same allele(%)  Different alleles(%)	Indels(%)


def indel_calc(row):
    indel_pc = round(row["indel"] / row["snp_gen"] * 100, 2)
    return indel_pc


def consisten_gen(row):
    consistent_gen = row["snp_gen"] - row["incon"]
    return consistent_gen


def get_col_widths(sample):
    length_dict = OrderedDict()
    sample_len = max(map(len, sample)) + 2
    # figure out if we have long or short
    if (
        len(sample[next(iter(sample))]) == 9
    ):  # finds the length of the value (list) for the first entry in the dict
        is_long = True
    else:
        is_long = False
    total_cols = [
        "tot_snp",
        "snp_gen",
        "no_gen",
        "hom",
        "het",
        "indel",
        "incon",
        "equiv",
        "inequiv",
    ]
    if is_long:
        columns = total_cols
    else:
        columns = total_cols.remove("equiv")
        columns = total_cols.remove("inequiv")
    sample_df = pd.DataFrame(sample, index=columns).transpose()
    sample_df["pc_indels"] = sample_df.apply(lambda row: indel_calc(row), axis=1)
    sample_df["con_gen"] = sample_df.apply(lambda row: indel_calc(row), axis=1)
    # get additional percentage values and remove columns we don't need, write to dict
    str_df = sample_df.applymap(str)
    max_tot_snp = str_df.tot_snp.map(len).max() + 2
    max_snp_gen = str_df.tot_snp.map(len).max() + 7
    max_hom = str_df.hom.map(len).max() + 7
    max_het = str_df.het.map(len).max() + 7
    # special indel calculations if <10 or 0
    max_pc_indels = str_df.pc_indels.map(len).max()
    if max_pc_indels == 0:
        indel_val = 5
    elif 0 < max_pc_indels < 10:
        indel_val = 6
    else:
        indel_val = 7
    max_indels = str_df.indel.map(len).max() + indel_val
    # work out consistent genotypes
    max_consistent = str_df.con_gen.map(len).max() + 2
    max_inconsistent = str_df.incon.map(len).max() + 2
    # Add all the lengths to the length_dict, except for the ones specific to Long format
    length_dict.update(
        {
            "max_sample": sample_len,
            "max_tot_snp": max_tot_snp,
            "max_snp_gen": max_snp_gen,
            "max_hom": max_hom,
            "max_het": max_het,
            "max_indels": max_indels,
            "max_consistent": max_consistent,
            "max_inconsistent": max_inconsistent,
        }
    )
    # get equivalencies if long
    if is_long:
        max_equiv = str_df.equiv.map(len).max() + 2
        max_inequiv = str_df.inequiv.map(len).max() + 2
        length_dict.update({"max_equiv": max_equiv, "max_inequiv": max_inequiv})
    # Figure out which are larger, header or values
    comparison_dict = {"sample_lengths": length_dict}
    head_dict = OrderedDict()
    head_dict.update(
        {
            "max_sample": 6,
            "max_tot_snp": 13,
            "max_snp_gen": 24,
            "max_hom": 13,
            "max_het": 15,
            "max_indels": 9,
            "max_consistent": 23,
            "max_inconsistent": 25,
        }
    )
    if is_long:
        head_dict.update({"max_equiv": 22, "max_inequiv": 24})
    comparison_dict.update({"header_lengths": head_dict})
    key_list = list(length_dict.keys())
    comparison_df = pd.DataFrame.from_dict(comparison_dict)
    sorted_comp_df = comparison_df.reindex(key_list, axis="index")
    sorted_comp_df["max_values"] = sorted_comp_df[
        ["sample_lengths", "header_lengths"]
    ].max(axis=1)
    column_widths = sorted_comp_df["max_values"].tolist()
    return column_widths


########################################################################################################################


def write_summary(
    input_df, file_type, f_name, inconsistency_val, equiv_cols, inequiv_cols, tab_format
):
    sample_dict = {}
    if file_type is None:
        exit("File type must be known to write a summary file")
    else:
        if file_type != "LONG":
            samples = list(input_df)
            samples.remove("Name")
            total_snps = len(input_df.index)
            for member in samples:
                none_list = ["--", "---"]
                number = input_df[~input_df[member].isin(none_list)].shape[
                    0
                ]  # make sure this works with affy data
                sub_df = input_df[["Name", member]].copy()
                nn_list = ["AA", "TT", "CC", "GG"]
                indel_list = ["II", "DD"]
                homozygous_snps = sub_df[sub_df[member].isin(nn_list)].copy()
                indels = sub_df[sub_df[member].isin(indel_list)].copy()
                num_dashes = total_snps - number
                num_homozygous = len(homozygous_snps.index)
                num_indels = len(indels.index)
                num_heterozygous = number - num_homozygous - num_indels
                # Total_snps, snps_with_data, snps_without_data, homozygous, heterozygous, indels
                sample_dict.update(
                    {
                        member: [
                            total_snps,
                            number,
                            num_dashes,
                            num_homozygous,
                            num_heterozygous,
                            num_indels,
                            inconsistency_val,
                            equiv_cols,
                            inequiv_cols,
                        ]
                    }
                )

        else:  # If the file is long
            samples_list = list(set(input_df["Sample ID"]))
            samples_index_dict = {}
            samples = []
            # reset index in input df
            input_df.reset_index(inplace=True)
            for n in samples_list:
                index = input_df.loc[input_df["Sample ID"] == n].index[0]
                samples_index_dict.update({index: n})
            for i in sorted(samples_index_dict):
                samples.append(samples_index_dict[i])
            for each in samples:
                sub_df = input_df[input_df["Sample ID"] == each]
                total_snps = len(sub_df.index)
                sub_df_no_gc = sub_df[
                    [
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
                ].copy()
                sub_df_no_gc.drop(["Sample ID"], axis=1, inplace=True)
                top_df = sub_df_no_gc[["Allele1 - Top", "Allele2 - Top"]].copy()
                top_df["Combine"] = top_df[["Allele1 - Top", "Allele2 - Top"]].apply(
                    lambda x: "".join(x), axis=1
                )
                nn_list = ["AA", "TT", "GG", "CC"]
                homozygous_snps = top_df[top_df["Combine"].isin(nn_list)].copy()
                indel_list = top_df[top_df["Combine"].isin(["II", "DD"])].copy()
                num_homozygous = len(homozygous_snps.index)
                num_indels = len(indel_list.index)
                nan_df = sub_df_no_gc[sub_df_no_gc == "-"].copy()
                nan_df.dropna(how="any", inplace=True)
                not_present = len(nan_df.index)  # aka num_dashes
                number = total_snps - not_present
                num_heterozygous = number - num_homozygous - num_indels
                # Total_snps, snps_with_data, snps_without_data, homozygous, heterozygous, indels
                sample_dict.update(
                    {
                        each: [
                            total_snps,
                            number,
                            not_present,
                            num_homozygous,
                            num_heterozygous,
                            num_indels,
                            inconsistency_val,
                            equiv_cols,
                            inequiv_cols,
                        ]
                    }
                )
        col_widths = get_col_widths(sample_dict)
        # sample, total_snp, snp_genotypes, hom, het, indels, consistent, inconsistent, [equiv, inequiv]

        # Write to file
        output_basename = f_name.replace(".txt", "")
        outfile_name = output_basename + "_summary.txt"
        with open(outfile_name, "w") as summary:
            summary.write(f_name + " Summary\n")
            if tab_format is True:
                summary.write(
                    "Sample\tTotal SNPs\tSNPs with genotypes\tHomozygous\tHeterozygous\tIndels\tConsistent genotypes\tInconsistent genotypes"
                )
                if file_type != "LONG":
                    summary.write("\n")
                else:
                    summary.write("\tEquivalent genotypes\tInequivalent genotypes\n")
            elif tab_format is False:
                # get max lengths of lines using col_widths dict

                summary_line = "{sample:{w1}}{tot_snp:{w2}}{snp_gen:{w3}}{hom:{w4}}{het:{w5}}{ind:{w6}}{con:{w7}}{incon:{w8}}".format(
                    sample="Sample",
                    tot_snp="Total SNPs",
                    snp_gen="SNPs with genotypes",
                    hom="Homozygous",
                    het="Heterozygous",
                    ind="Indels",
                    con="Consistent genotypes",
                    incon="Inconsistent genotypes",
                    w1=col_widths[0],
                    w2=col_widths[1],
                    w3=col_widths[2],
                    w4=col_widths[3],
                    w5=col_widths[4],
                    w6=col_widths[5],
                    w7=col_widths[6],
                    w8=col_widths[7],
                )
                summary.write(summary_line)
                if file_type != "LONG":
                    summary.write("\n")
                else:
                    summary_line = "{equiv:{w9}}{inequiv:{w10}}".format(
                        equiv="Equivalent genotypes",
                        inequiv="Inequivalent genotypes",
                        w9=col_widths[8],
                        w10=col_widths[9],
                    )
                    summary.write(summary_line + "\n")

            else:
                pass
            for items in sample_dict:
                percent_snps_with_data = (
                    str(sample_dict[items][1])
                    + " ("
                    + str(
                        round((sample_dict[items][1] / sample_dict[items][0]) * 100, 2)
                    )
                    + ")"
                )
                percent_same_allele = (
                    str(sample_dict[items][3])
                    + " ("
                    + str(
                        round((sample_dict[items][3] / sample_dict[items][1]) * 100, 2)
                    )
                    + ")"
                )
                percent_different_alelle = (
                    str(sample_dict[items][4])
                    + " ("
                    + str(
                        round((sample_dict[items][4] / sample_dict[items][1]) * 100, 2)
                    )
                    + ")"
                )
                percent_indels = (
                    str(sample_dict[items][5])
                    + " ("
                    + str(
                        round((sample_dict[items][5] / sample_dict[items][1]) * 100, 2)
                    )
                    + ")"
                )
                consistent_genotypes = sample_dict[items][1] - sample_dict[items][6]
                inconsistent_genotypes = sample_dict[items][6]
                if file_type != "LONG":
                    if tab_format:
                        output_line = "{sample}\t{tot_snp}\t{snp_gen}\t{hom}\t{het}\t{ind}\t{con}\t{incon}\n".format(
                            sample=items,
                            tot_snp=sample_dict[items][0],
                            snp_gen=percent_snps_with_data,
                            hom=percent_same_allele,
                            het=percent_different_alelle,
                            ind=percent_indels,
                            con=consistent_genotypes,
                            incon=inconsistent_genotypes,
                        )
                    else:
                        output_line = (
                            "{sample:<{w1}}{tot_snp:<{w2}}{snp_gen:<{w3}}{hom:<{w4}}"
                            "{het:<{w5}}{ind:<{w6}}{con:<{w7}}{incon:<{w8}}\n".format(
                                sample=items,
                                tot_snp=sample_dict[items][0],
                                snp_gen=percent_snps_with_data,
                                hom=percent_same_allele,
                                het=percent_different_alelle,
                                ind=percent_indels,
                                con=consistent_genotypes,
                                incon=inconsistent_genotypes,
                                w1=col_widths[0],
                                w2=col_widths[1],
                                w3=col_widths[2],
                                w4=col_widths[3],
                                w5=col_widths[4],
                                w6=col_widths[5],
                                w7=col_widths[6],
                                w8=col_widths[7],
                            )
                        )

                else:  # if file type IS Long
                    if tab_format:
                        output_line = (
                            "{sample}\t{tot_snp}\t{snp_gen}\t{hom}\t{het}\t{ind}\t{con}\t{incon}\t{equiv}\t"
                            "{inequiv}\n".format(
                                sample=items,
                                tot_snp=sample_dict[items][0],
                                snp_gen=percent_snps_with_data,
                                hom=percent_same_allele,
                                het=percent_different_alelle,
                                ind=percent_indels,
                                con=consistent_genotypes,
                                incon=inconsistent_genotypes,
                                equiv=equiv_cols,
                                inequiv=inequiv_cols,
                            )
                        )
                    else:
                        output_line = (
                            "{sample:<{w1}}{tot_snp:<{w2}}{snp_gen:<{w3}}{hom:<{w4}}"
                            "{het:<{w5}}{ind:<{w6}}{con:<{w7}}{incon:<{w8}}{equiv:<{w9}}{inequiv:<{w10}}\n".format(
                                sample=items,
                                tot_snp=sample_dict[items][0],
                                snp_gen=percent_snps_with_data,
                                hom=percent_same_allele,
                                het=percent_different_alelle,
                                ind=percent_indels,
                                con=consistent_genotypes,
                                incon=inconsistent_genotypes,
                                equiv=equiv_cols,
                                inequiv=inequiv_cols,
                                w1=col_widths[0],
                                w2=col_widths[1],
                                w3=col_widths[2],
                                w4=col_widths[3],
                                w5=col_widths[4],
                                w6=col_widths[5],
                                w7=col_widths[6],
                                w8=col_widths[7],
                                w9=col_widths[8],
                                w10=col_widths[9],
                            )
                        )
                summary.write(output_line)
            summary.close()


#  This is just for testing - do not run write_summary as a stand-alone program!
#  Use check_format with the --summarize option instead
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Writes summary of matrix and long format files - "
        "DO not run as a standalone program; run check_format with the "
        "--summary option instead. Standalone functionality is for testing "
        "purposes only and may give incorrect results."
    )
    parser.add_argument("--input", type=str, help="Input file to be summarized")
    parser.add_argument(
        "--input-format",
        type=str,
        choices=["TOP", "FWD", "DESIGN", "AB", "PLUS", "LONG"],
        help="Type of file(s) expected: 'TOP', 'FWD', 'DESIGN', 'AB', 'PLUS' or 'LONG'. \n"
        "Note that AB format does not distinguish between indels and variants, and will not report indels "
        "separately",
    )
    parser.add_argument(
        "-t",
        action="store_true",
        default=False,
        required=False,
        help="Output summary in tabular format. Default: False",
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    header_row, header_dict = fp.parse_header(args.input)
    # Make new dataframes
    with open(args.input, "r") as input_file:
        df = pd.read_csv(input_file, skiprows=header_row, sep="\t")
        df.rename(columns={"Unnamed: 0": "Name"}, inplace=True)
        inconsistency = 0
        inequivalent = 0
        equivalent = 0
        write_summary(
            df,
            args.input_format,
            args.input,
            inconsistency,
            equivalent,
            inequivalent,
            args.t,
        )
