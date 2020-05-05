#!/usr/bin/env python3


# This module creates PED and MAP files for use in PLINK. For more information on PED and MAP files, see
# http://zzz.bwh.harvard.edu/plink/data.shtml

# PED file specification - Individual ID and genotype only, everything is SPACE-delimited (not tab)
# snp_panel_genotypes: y-axis = marker, x-axis = individual


def create_ped_file(output_file_basename, snp_panel_genotypes, logfile):
    """
    Creates a basic PED file from a SNP panel - only contains Individual ID and genotypes, all whitespace delimited
    :param output_file_basename: basename of the output file for writing
    :param snp_panel_genotypes: snp panel data
    :param logfile: logfile
    :return: PED file written
    """
    ped_output_name = output_file_basename + ".ped"
    transposed_file = snp_panel_genotypes.set_index('Name').T
    genotype_columns = list(transposed_file.columns)
    for column in genotype_columns:
        a_column = column + "_A"
        b_column = column + "_B"
        transposed_file['{}'.format(a_column)], transposed_file['{}'.format(b_column)] = \
            zip(*transposed_file[column].apply(lambda x: list(x)))
    transposed_file.drop(columns=genotype_columns, inplace=True)
    transposed_file.replace({"-": "0"}, inplace=True)
    transposed_file.to_csv(ped_output_name, sep=" ", header=False)
    return transposed_file, logfile


def chromosome_replace(species_name, var_dataframe):
    """
    Generates a dict for replacement values for chromosomes X, Y, PAR, and MT according to PLINK specifications
    :param species_name: species name, either bos_taurus or sus_scrofa
    :param var_dataframe: variant dataframe
    :return: dict containing {old: new} values
    """
    chrom_list = var_dataframe["BLAST_chromosome"].tolist()
    if species_name == "bos_taurus":
        # Use non-standard chromosome IDs: use flag --cow or --chr-set 29 no-xy
        chr_dict = {}
        if "X" in chrom_list:
            chr_dict.update({"X": "30"})
        if "Y" in chrom_list:
            chr_dict.update({"Y": "31"})
        if "MT" in chrom_list:
            chr_dict.update({"MT": "32"})
        if "" in chrom_list:
            chr_dict.update({"": "0"})
    elif species_name == "sus_scrofa":
        # Use non-standard chromosome IDs: use flag --chr-set 18 no-xy
        chr_dict = {}
        if "X" in chrom_list:
            chr_dict.update({"X": "19"})
        if "Y" in chrom_list:
            chr_dict.update({"Y": "20"})
        if "MT" in chrom_list:
            chr_dict.update({"MT": "21"})
        if "" in chrom_list:
            chr_dict.update({"": "0"})
    else:
        chr_dict = {}
    return chr_dict


def create_map_file(output_file_basename, snp_panel_df, var_df, species, logfile):
    """
    Creates MAP file from a SNP panel file and var df - contains 4 columns - chromosome, SNP identifier (marker),
    genetic distance, and base pair position. Genetic distance is always 0.
    :param output_file_basename: basename of input file for writing output
    :param snp_panel_df: dataframe containing the SNP marker panel
    :param var_df: variant dataframe for getting chromosome position information
    :param species: species, required for setting chromosome number properly
    :param logfile: logfile
    :return: MAP file written
    """
    map_output_name = output_file_basename + ".map"
    # Get var df containing only panel data
    snp_panel_name = snp_panel_df["Name"].tolist()
    var_subset = var_df[var_df["Name"].isin(snp_panel_name)].copy()
    var_subset.drop(columns=["SNP", "BLAST_strand", "Reference_allele_forward_strand", "DESIGN_A", "DESIGN_B",
                             "FORWARD_A", "FORWARD_B", "PLUS_A", "PLUS_B", "TOP_A", "TOP_B"], inplace=True)
    # get new chromosome information
    replace_vals = chromosome_replace(species, var_subset)

    if bool(replace_vals):
        var_chr_replaced = var_subset[var_subset.BLAST_chromosome.replace(replace_vals)]
    else:
        var_chr_replaced = var_subset.copy()
    output_df = var_chr_replaced[["BLAST_chromosome", "Name", "BLAST_position"]].copy()
    output_df["Genetic_distance"] = "0"
    output_df = output_df[["BLAST_chromosome", "Name", "Genetic_distance", "BLAST_position"]]
    convert_dict = {"BLAST_position": "int64"}
    output_df2 = output_df.astype(convert_dict)
    output_df2.to_csv(map_output_name, sep=" ", header=False, index=False)
    return map_output_name, logfile
