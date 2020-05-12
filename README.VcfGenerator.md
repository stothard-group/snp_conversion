# VCF Generator

- FILE: vcf_generator
- AUTH: Emily Herman (eherman@ualberta.ca)
- DATE: MAY 11, 2020
- VERS: 1.0

This program accepts Illumina genotype files in matrix format (AB, Forward, Top,
 Design, Plus) Long format, and Affymetrix files, and using previously computed 
 variant position files for SNP datasets, converts the data to VCF format.
 
vcf_generator is written in Python 3. Please read the INSTALL file prior to 
use.

README Sections

[Program Options](#program-options)

[Input Files](#input-files)

[Variant Files](#genotype-conversion-key-files-variant-files)

[Program Output](#program-output)

[Example Commands](#example-commands)

## Program Options
```
usage: vcf_generator [-h] --snp-panel SNP_PANEL
                     [--panel-type {TOP,FWD,AB,DESIGN,LONG,PLUS,affymetrix}]
                     [-v] [--key-dir KEY_DIR] --assembly ASSEMBLY --species
                     SPECIES [--output OUTPUT] [-t THREADS]
                     [-c] [--discard-snp]

Converts SNP panel file to VCF format file

optional arguments:
  -h, --help            show this help message and exit
  --snp-panel SNP_PANEL
                        Name or path to the SNP panel file. Can be Illumina or
                        Affymetrix format
  --panel-type {TOP,FWD,AB,DESIGN,LONG,PLUS,affymetrix}
                        Type of file(s) expected: 'TOP', 'FWD', 'AB',
                        'DESIGN', 'LONG', 'PLUS', or 'affymetrix' (Affymetrix)
  -v, --verbose-logging
                        [optional] Write output to both STDOUT and log file
  --key-dir KEY_DIR     Directory containing genotype conversion key files
                        (default = variant_position_files)
  --assembly ASSEMBLY   Assembly name (use conversion_list utility to see all 
  						available choices)
  --species 			Species name (use conversion_list utility to see all 
  						available choices)
  --output OUTPUT       Full output file name to store VCF results (default:
                        [input basename].vcf
  -t THREADS, --threads THREADS
                        [optional] Number of threads to use if conversion to
                        PLUS format is required (default = 2)
  -c, --compress        Compress VCF output file with gzip (not bgzip)
  --discard-snp         Set ID field values to '.' (default: ID field contains
                        SNP marker names)
```


## Input Files

### Acceptable file types as inputs

In general, this program accepts files of Illumina matrix, Illumina Long, and 
Affymetrix formats. See files in `sample_files/input_files` for examples of 
each format. The program then converts the data to PLUS format (either Illumina 
or Affymetrix) prior to writing the VCF file. While the program will do this 
for you with any type of allowed input, **it is strongly recommended that you 
perform this conversion separately prior to running vcf_generator, by using the 
 convert_file utility within the snp_conversion program.** This is because 
 conversion can be a memory- and time-intensive process, and if the 
 vcf_generator program is stopped for any reason, conversion will need to be 
 redone. 
 
The SNP panel file may contain indels, and may contain positions for which 
there is no associated information in the genotype conversion key file. In 
both cases, the genotype in the VCF file will be listed as `./.`. 


## Genotype conversion key files (Variant Files)

Variant files contain strand format information for the markers in each panel, 
which is specific to the genomic assembly. These files are required for format 
checking and conversion prior to vcf generation. The assembly name and species 
information must therefore be specified using the `--assembly` and `--species` 
options, respectively. The program will detect the matching or best-matching 
variant file. The filename for all variant files will have the structure 
`[panel name].[assembly].[conversion|position].csv[.gz]`. Ideally, the user 
input file will contain all the markers present in a single variant file, and 
this file will be used in format checking. These files can be generated 
using the Nextflow workflow found in the Genotype Conversion File Builder 
repository at 
https://github.com/stothard-group/genotype_conversion_file_builder

If there is not an exact match between the user input markers and those in the 
variant file, the program will select the file that contains 
all matching markers, with the fewest additional ones. If there are markers in 
the user input file that are not found in the best-matching variant file, the 
program will exit and print the orphan markers to the file 
`[snp panel basename]_problem_variants.txt`. 

Variant files should be placed in the directory 
`variant_position_files/[species]/[assembly]/`. A different variant file 
directory can be specified, but it must preserve the 
`[dir]/[species]/[assembly]/` structure. Please store the variant files, and 
only the variant files, in their own directory.

The variant file headers (indicated by lines beginning with '#') contain more 
information about variant file structure. 

## Program Output

The vcf_generator will create a VCF format file with the following metadata:
```
##fileformat=VCFv4.2
##filedate={1}
##source=vcf_generator.py
##reference={2}
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##vcf_generator_commandline.vcf_generator=<{3}>
```
1. Date (format: YearMonthDay)
2. Reference assembly name specified by `--assembly`
3. Parameters suppled to the vcf_generator by the user

### Main Data Fields

All missing data is indicated by `.`

CHROM and POS fields contain the chromosome and nucleotide position of the SNP, 
respectively. 

By default, the ID field contains the SNP panel marker names. To remove these 
values and replace them with a `.`, use the flag `--discard-snp`. 

REF and ALT fields contain the reference and alternative allele values defined 
in the genotype conversion key file. 

QUAL, FILTER, and INFO values are `.`

The FORMAT line is simply GT (genotype), and the genotype for each individual 
is listed in the subsequent columns. 

For more information on the VCF file format, see the 
[Samtools VCF specification document](https://github.com/samtools/hts-specs/blob/master/VCFv4.2.pdf).

## Example Commands

The following command uses a sample input file from `sample_files/input files/` 
and creates the file `50kv3_mFWD_14June2019_converted.vcf`. 

```
./vcf_generator --snp-panel sample_files/input_files/50kv3_mFWD_14June2019.txt \
--panel-type FWD --key-dir variant_position_files \
--assembly ARS-UCD1_2_Btau5_0_1Y --species bos_taurus \
--output 50kv3_mFWD_14June2019_converted.vcf 
```

