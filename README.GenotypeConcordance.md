# Genotype Concordance

- FILE: genotype_concordance.py
- AUTH: Emily Herman (eherman@ualberta.ca)
- DATE: APR 7, 2020
- VERS: 2.0

This program determines the concordance of variants between VCF genotype files 
and SNP panel files (both Illumina and Affymetrix), using previously computed 
variant position files. 

The Genotype Concordance program is written in Python 3. See the INSTALL file 
for installation instructions.

README sections

[Quick Start](#quick-start)

[Usage](#usage)

[Program Options](#program-options)

[Input Files](#input-files)

[Variant Files](#genotype-conversion-key-files-variant-files)

[Filtering Options](#filtering-options)

[Program Output](#program-output)

* [Metastatistics](#metastatistics)
* [Statistics](#statistics)

[Extracting Discordant Positions](#extracting-discordant-positions)

[Special Options](#special-options)

[Example Commands](#example-commands)


## Quick Start

Make sure that you have installed the required python libraries listed in the 
INSTALL file.

1. Download the program using git clone: 
`git clone https://github.com/stothard-group/snp_conversion_2.git`
2. Create a directory to store the conversion files. By default, the program 
looks for these in the directory `variant_position_files`, however this can be 
changed with the ``--key-dir`` option.
3. Add the conversion files to this directory.
4. Change directories to the `concordance_sample_files` directory and execute 
the `sample_concordance.sh` script:
`./sample_concordance.sh`

The program is executed using several small datasets to test assess concordance 
of an Illumina (LONG format) and Affymetrix dataset with a sample VCF file. 
This script will move all output files into the directory 
`concordance_sample_files/test_output/`, and compare these files with 
the expected files in `concordance_sample_files/sample_output/`. If the files 
are the same, the output `No differences found` is printed to the screen for 
each comparison. A list of sample commands found in the script is given below in 
[Example Commands](#example-commands), with a short explanation of each.

## Usage

```
usage: ./genotype_concordance  [-h] --snp-panel SNP_PANEL --panel-type
                               {TOP,FWD,AB,PLUS,DESIGN,LONG,affymetrix}
                               --vcf-file VCF_FILE [--key-dir KEY_DIR]
                               --assembly ASSEMBLY
                               --species {bos_taurus, sus_scrofa}
                               [--filter-vcf]
                               [--qual QUAL] [--filter FILTER]
                               [--output OUTPUT]
                               [--output-type {basic,tabular,pretty}] [-v]
                               [--extract-discordant] [-t THREADS]
                              
The following arguments are required: --snp-panel, --panel-type, --vcf-file,
--assembly, --species
```

### Program Options

```
  --snp-panel SNP_PANEL
                        Name or path to the SNP panel file. Can be Illumina or
                        Affymetrix format
  --panel-type {TOP,FWD,AB,PLUS,DESIGN,LONG,affymetrix}
                        Type of panel file: 'TOP', 'FWD', 'AB', 'PLUS',
                        'DESIGN', 'LONG', 'affymetrix'
  --vcf-file VCF_FILE   Name or path to the VCF file containing genotype
                        information
  --key-dir KEY_DIR     Directory containing genotype conversion key files
                        (default directory: variant_position_files)
  --assembly ASSEMBLY   Assembly name - see README for full list of choices
  --species {bos_taurus,sus_scrofa}
                        Organism name
  --filter-vcf          [Optional] Use this flag to specify that the VCF file
                        should be filtered on QUAL and/or FILTER values. Must
                        be used in conjunction with --qual and/or --filter
                        parameters
  --qual QUAL           [Optional] Only perform concordance analysis using
                        variants with QUAL scores higher than this value. Must
                        be used in conjunction with the --filter-vcf flag
  --filter FILTER       [Optional] Only perform concordance analysis using
                        variants with these value(s) in the FILTER field. Can
                        be a single value or a comma-separated list. Must be
                        used in conjunction with --filter-vcf flag
  --output OUTPUT       [Optional] Output prefix to append to
                        '_metastatistics.txt' and '_statistics.txt' files
                        (default: 'concordance')
  --output-type {basic,tabular,pretty}
                        [Optional] Type of output for statistics file: basic
                        (tsv, stats only), tabular (tsv with extra info),
                        pretty (nice formatting) (default: tabular)
  -v, --verbose-logging
                        [Optional] Write program steps and errors to a log
                        file
  --extract-discordant  [Optional] Write discordant positions to a file called
                        [snp-panel]_discordant.txt
  -t THREADS, --threads THREADS
                        [optional] Number of threads to use during conversion
                        (default = 2)
```

## Input Files

#### SNP panel

The SNP panel file can be in either Illumina or Affymetrix format. For Illumina,
 the file can be in 'Long' format or any of the matrix formats. The program then
 converts the data to PLUS format (either Illumina or Affymetrix) prior to the 
 concordance analysis. While the program will do this for you with any type of 
 allowed input, **it is strongly recommended that you perform this 
 conversion separately prior to running genotype_concordance.py, by using the 
 convert_file utility within the snp_conversion program.** This is because 
 conversion can be a memory- and time-intensive process, and if the 
 concordance analysis is stopped for any reason, conversion will need to be 
 redone. 
 
 The SNP panel file may contain positions and samples not found in the VCF file.
 
#### VCF file

The VCF file can be compressed with zip, gzip, or bzip2, or can be in a
native, uncompressed format. The path to the VCF file must be specified with 
`--vcf-file`. The VCF file may contain positions and samples not found in the 
SNP panel file. 

## Genotype conversion key files (Variant Files)

Variant files contain strand format information for the markers in each panel, 
which is specific to the genomic assembly. These files are required for format 
checking and conversion prior to the concordance analysis. The assembly name 
and species information must therefore be specified when running `check_format` 
and `convert_file` utilities using the `--assembly` and `--species` options, 
respectively. SNP conversion utilities will detect the matching or best-matching 
variant file. The filename for all variant files will have the structure 
`[panel name].[assembly].[conversion|position].csv[.gz]`. Ideally, the user 
input file will contain all the markers present in a single variant file, and 
this file will be used in format checking. 

If there is not an exact match between the user input markers and those in the 
variant file, Genotype Concordance modules will select the file that contains 
all matching markers, with the fewest additional ones. If there are markers in 
the user input file that are not found in the best-matching variant file, SNP 
conversion will exit and print the orphan markers to the file 
`[snp panel basename]_problem_variants.txt`.

Variant files should be placed in the directory 
`variant_position_files/[species]/[assembly]/`. A different variant file 
directory can be specified, but it must preserve the 
`[dir]/[species]/[assembly]/` structure. Please store the variant files, and 
only the variant files, in their own directory.

The variant file headers (indicated by lines beginning with '#') contain more 
information about variant file structure. 

## Filtering Options

Positions in the VCF file can be filtered by values in the FILTER and QUAL 
fields. To filter positions, use the `--filter-vcf` flag, in conjunction with 
`--qual` and/or `--filter` options. By filtering on the QUAL value, the program 
will include only variants that have quality scores greater than this value. By 
 filtering on the FILTER value, the program will include only variants that have
 this term in the concordance analysis. The filter value may be a single value 
 or a list of values separated by commas (but no spaces).
 
 
## Program Output

The Genome Concordance program generates two main output files: a metastatistics 
report containing the overall concordance statistics of the two files, and a 
statistics report containing per-sample concordance statistics. The `--output` 
option can be used to specify the prefix for "_metastatistics.txt" and 
"_statistics.txt" output files (the default prefix is "concordance").

##### Metastatistics
The following meta-statistics are reported:

* Total samples in the SNP panel file
* Number and percentage of SNP panel samples in the VCF file
* A list of samples in the SNP panel file and VCF file
* A list of samples in the SNP panel file but *not* the VCF file
* Number of probes in the SNP panel file and positions in the VCF file
* Number of positions shared between the SNP panel file and VCF file, based on 
genomic location
* The percentage of SNP panel positions and VCF positions that the previous 
value represents

##### Statistics

There are three output formats for the per-sample concordance statistics: basic,
 tabular (default), and pretty, which can be specified with `--output-type`. 
Basic and tabular output formats are similar (both are tab-delimited), but 
tabular output has more detailed explanation of the output statistics. The basic
 format is provided for users who wish to extract statistics values easily. The 
pretty format is not tab-delimited; rather, it creates a nicely formatted output
 table.

For each sample, the following statistics are reported in this order with 
percentages following the raw values in brackets:

* Number of markers shared between SNP panel and VCF file

Of the shared markers:

* Number of informational markers (percentage of shared markers)
* Number of non-informational markers (percentage of shared markers)

Of the markers with genotype information (informational):

* Number of homozygous reference concordant (percentage of informational)
* Number of homozygous alternate concordant (percentage of informational)
* Number of heterozygous concordant (percentage of informational)

**Concordancy analysis**

Total concordant and discordant markers:

* Number of total concordant markers (percentage of informational) ***
* Number of true concordant markers (percentage of informational) ***
* Number of discordant markers (percentage of informational)

Concordance breakdown:

* Number of homozygous reference concordant markers (percentage of concordant)
* Number of homozygous alternate concordant markers (percentage of concordant)
* Number of discordant markers where the SNP panel and VCF file are both 
homozygous: one is the reference allele, one is the alternate allele (percentage
 of discordant)
* Number of discordant markers where the SNP panel is homozygous and the VCF 
file is heterozygous: panel has the ref or alt allele, and the VCF file is 
ref/alt (percentage of discordant)
* Number of discordant makers where the SNP panel is homozygous and the VCF file
 contains alleles that are neither ref nor alt (percentage of discordant)

* Number of heterozygous concordant markers: the heterozygous portion of 
**total** concordant markers (percentage of concordant) ***
* Number of true heterozygous concordant markers: the heterozygous portion of 
**true** concordant markers (percentage of concordant) ***
* Number of heterozygous concordant markers where the genotype is reversed (i.e.
 alt/ref) (percentage of concordant)
* Number of discordant markers where the SNP panel is heterozygous and the VCF 
file is homozygous (either ref or alt) (percentage of discordant)
* Number of discordant markers where the SNP panel is heterozygous and the VCF 
file contains alleles that are neither ref nor alt (percentage of discordant)

---
\***
A note on heterozygous concordance calculations

Sometimes heterozygous alleles appear "backwards" in the SNP panel compared to 
the VCF file, where the alternate allele is listed first. This program considers
 both ref/alt and alt/ref heterozygous alleles to be concordant, but the 
statistics distinguish between these cases. Heterozygous alleles where the SNP 
panel and VCF file alleles are both ref/alt are counted as True concordant, 
whereas Total concordant includes True concordant heterozygous alleles, all 
homozygous alleles, and alt/ref heterozygous alleles.

---

## Extracting discordant positions

Discordant positions can be extracted using the option `--extract-discordant`, 
which will write to a file prefixed with the SNP panel name, followed by 
"_discordant.txt". This file has the following columns:

| CHROM	| POS |	Name | Sample | VCF genotype | Panel genotype |
|-------|-----|------|--------|--------------|----------------|


1. Chromosome (CHROM) from VCF file
2. Position (POS) from VCF file
3. Marker name from SNP panel file
4. Name of the sample in which the position is discordant
5. Genotype in the VCF file
6. Genotype in the SNP panel

## Special Options

### Multi-threading

Some portions of the program can make use of multiple cores on Unix systems 
(Linux, Macintosh, and WSL - Windows Subsystem for Linux). This behaviour can be
controlled using the `--threads` option. On these systems, the default is 2 
threads, but more threads can be specified. Unfortunately, threading is not 
supported on Windows machines. 

### Verbose logging

The `--verbose-logging` flag will write progress and error messages to a log 
file, which takes the name of the SNP panel file as a prefix, followed by a 
date and timestamp, followed by ".log".

## Example commands

Below are some sample commands for genotype_concordance tasks. Sample data are 
provided to run each of these commands within the 
`concordance_sample_files/input_files/` directory. In 
`concordance_sample_files`, running the bash script `./sample_concordance.sh` 
will execute the commands below. The bash 
script will then move all output files to the directory 
`concordance_sample_files/test_output/` and compare these files to files 
previously generated in `concordance_sample_files/sample_output/`. Using these 
commands, messaging that is normally printed to the screen is piped to output 
files ending with "*_output.txt".
****
Note that you may have to make the script executable using 

```
chmod u+x sample_concordance.sh
``` 

Concordance analysis between an Illumina LONG file and a VCF file, with tabular 
output
```
python ../genotype_concordance \
--snp-panel input_files/G_CCGP_long_sample_input.txt --panel-type LONG \
--vcf-file input_files/SNPs_reduced_anon.vcf.gz --species bos_taurus \
--assembly ARS-UCD1_2_Btau5_0_1Y --output-type tabular \
--output concordance_test1 > test_output/long_vs_vcf_tab_concordance.txt
```

Concordance analysis between an Illumina LONG file and a VCF file, filtering on 
quality values, with pretty output
```
python ../genotype_concordance \
--snp-panel input_files/G_CCGP_long_sample_input.txt --panel-type LONG \
--vcf-file input_files/SNPs_reduced_anon.vcf.gz --species bos_taurus \
--assembly ARS-UCD1_2_Btau5_0_1Y --filter-vcf --qual 100 --output-type pretty \
--output concordance_q100_test2 \
> test_output/long_vs_vcf_q100_pretty_concordance.txt
```

Concordance analysis between an Affymetrix file and a VCF file, outputting a 
list of discordant positions
```
python ../genotype_concordance \
--snp-panel input_files/G_CCGP_affy_short_input.txt --panel-type affymetrix \
--vcf-file input_files/SNPs_reduced_anon.vcf.gz --species bos_taurus \
--assembly ARS-UCD1_2_Btau5_0_1Y --extract-discordant \
--output concordance_affy > test_output/affy_vs_vcf_concordance.txt
```