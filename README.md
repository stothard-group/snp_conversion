# SNP Conversion

- FILE: SNP_conversion.py
- AUTH: Emily Herman (eherman@ualberta.ca)
- DATE: APR 7, 2020
- VERS: 2.0

This program accepts Illumina genotype files in matrix format (AB, Forward, Top, Design, Plus) Long format, and Affymetrix
 files, and using previously computed variant position files for SNP datasets, checks the input file for formatting and 
 SNP accuracy. It can also convert between formats and merge files of the same format. 
 
SNP_conversion.py is written in Python 3. Please read the INSTALL file prior to use.

README Sections

[Program Utilities](#program-utilities)

[Program Options](#program-options)

[Input Files](#input-files)

[Variant Files](#variant-files)

[Notes on File Conversion](#notes-on-file-conversion)

[Notes on Options](#notes-on-options)

[Example Commands](#example-commands)


## Program Utilities

### Usage to check the format of input file(s) only:

```
usage: python SNP_conversion.py check_format [-h] [--input-dir INPUT_DIR]
                                      [--file-list FILE_LIST]
                                      [--input-format {TOP,FWD,AB,LONG,DESIGN,PLUS,mixed,affymetrix}]
                                      [--get-snp-panel] [--key-dir KEY_DIR]
                                      --species {bos_taurus, sus_scrofa}
                                      [-v] --assembly ASSEMBLY [-s]
                                      [--tabular] [--plink]


the following arguments are required: --assembly, --species
```
### Usage to convert input file(s) to new format (and check file formats):

```
usage: python SNP_conversion.py convert_file [-h] [--input-dir INPUT_DIR]
                                      [--file-list FILE_LIST]
                                      [--input-format {TOP,FWD,AB,LONG,DESIGN,PLUS,mixed,affymetrix}]
                                      --output-format
                                      {TOP,FWD,AB,PLUS,DESIGN,LONG,AFFY-PLUS}
                                      [--output-name OUTPUT_NAME] [-t THREADS]
                                      [-s] --assembly ASSEMBLY
                                      --species {bos_taurus, sus_scrofa}
                                      [--key-dir KEY_DIR] [--tabular] [-v] [--plink]


the following arguments are required: --output-format, --assembly, --species
```

### Usage to merge multiple files of the same format:
```
usage: python SNP_conversion.py merge_files [-h] --input-dir INPUT_DIR
                                     [--file-list FILE_LIST]
                                     [--input-format {TOP,FWD,AB,LONG,DESIGN,PLUS,affymetrix}]
                                     --output-name OUTPUT_NAME

the following arguments are required: --output-name, --input-format
```

## Program options

### General options:
```
  --input-dir INPUT_DIR
                        directory containing input file(s) (default: current directory)
  --file-list FILE_LIST
                        [optional] comma-separated list of files in the input
                        directory to be converted (no whitespace)
  --input-format {TOP,FWD,AB,LONG,DESIGN,PLUS,mixed,affymetrix}
                        Type of file(s) expected: 'TOP', 'FWD', 'AB', 'LONG',
                        'DESIGN', 'PLUS', 'mixed', or 'affymetrix'. 'mixed'
                        may not be used when merging files.
```
#### Check format options:
```
  --get-snp-panel       [optional] Will determine which genotype conversion
                        key files contain all SNPs in the input
  --key-dir KEY_DIR     Directory containing genotype conversion key files
                        (default directory: variant_position_files)
  -v, --verbose-logging
                        [optional] Write output to both STDOUT and log file
  --assembly ASSEMBLY   Assembly name - see README for full list of choices
  --species {bos_taurus,sus_scrofa}
                        Organism name
  -s, --summary         Summarize converted SNP file in *_summary.txt file
  --tabular TABULAR		Output summary file in tabular format (default: True)
  --plink				Creates PLINK flat files (PED and MAP) (default: False)
```

#### Convert File options:

```
  --output-format {TOP,FWD,AB,PLUS,DESIGN,LONG}
                        Type of file(s) to be created: 'TOP', 'FWD', 'AB',
                        'PLUS', 'DESIGN', 'LONG'. Only one type of output can
                        be specified at a time. PLUS converts the data to the
                        forward strand of the reference genome. LONG refers to
                        the long-format Illumina file.
  --output-name OUTPUT_NAME
                        Suffix of output file, without file extension. File
                        will be named [input_file].[output_name].txt (default
                     	: output)
  -t THREADS, --threads THREADS
                        [optional] Number of threads to use during conversion
                        (default = 2). This option will be ignored on Windows 
                        systems.
  -s, --summary         Summarize converted SNP file in *_summary.txt file
  --assembly ASSEMBLY   Assembly name - see README for full list of choices
  --species {bos_taurus,sus_scrofa}
                        Organism name
  --key-dir KEY_DIR     Directory containing genotype conversion key files
                        (default directory: variant_position_files)
  --tabular TABULAR		Output summary file in tabular format (default: True)
  --plink				Creates PLINK flat files (PED and MAP) (default: False)


```
#### Merge file options:
```
  --output-name OUTPUT_NAME
                        Name of merged file

```
## Note for Windows Users

This program can be used on Windows systems and on systems with Windows Subsystem for Linux (WSL). **If you have WSL 
installed**, it is strongly recommended that SNP Conversion is run in that partition (and the user input files are also 
in that partition). This will allow for parallel processing of input data and therefore faster program runtime.  **If 
WSL is not installed**, you will not be able to run the `sample_conversion.sh` Bash script within the `sample_files` 
directory. However, you will be able to run the sample commands within that file (also found at the end of this README). 
You can compare the output files from those commands with the contents of the corresponding files in 
`sample_files/sample_output/`, as these should be identical, with the exception of timestamps.

## Input Files

### Acceptable file types as inputs

In general, this program accepts files of Illumina matrix, Illumina Long, and Affymetrix formats. See files in 
`sample_files/input_files` for examples of each format.

For file conversion, the file type combinations are restricted to the following:

|  Input File | Output File  |
|---|---|
| Illumina Top, Forward, Design, Plus, Long | Illumina Top, Forward, Design, Plus, Long, **AB**  |
| Affymetrix native*, Affymetrix Plus  | Affymetrix Plus, Affymetrix native* |


*Equivalent to Illumina TOP format

**Note that AB format is not a valid input format for file conversion, as there is not enough information in the AB file 
to support conversion**

## Merge Files Inputs

**All input for the merge_files utility must be of the same file type. The merge_files utility does not check that 
this is the case, so the onus is on the user to use check_format, and ensure that only files of the same format are 
merged.**

### Internal file formatting

#### Illumina matrix and long formats

SNP conversion modules expect matrix and long files to have a header block and a data block:

[Header]

... header information ...

[Data]

... data ...

The following input format structures are accepted (with warnings):

- Files lacking the [Header] block (beginning with [Data])
- Files with an empty [Header] block ([Header]\n[Data])
- Files beginning with sample names in the case of matrix files, or Sample ID-Allele1-Allele2 etc. in the case of long
format files 

#### Affymetrix format

SNP conversion modules expect Affymetrix files to begin with the data, where the first line should contain sample names 
under the column name 'probeset_id', followed by two identical column names for each sample. If an Affymetrix file 
contains a header, remove the header information before submitting as an input file.


#### Notes on predicting the format of input files

If a file is submitted with the option `--input-format mixed` (or no input-format option), the script will try to 
determine the most likely input format. 

If the input file contains SNP markers determined to be potentially incorrect based on the predicted file
format, the script will NOT output the incorrect markers to a log file. To get a list of these markers, run check_format again 
using the predicted file format as the input format. 

#### Input file types

See [DNA strand designations](docs/dna_strand_designations.pdf) for an in-depth explanation of the Illumina Forward, Top, 
AB, and Design formats as they relate to the genome.

Briefly:
 - Design strand: the strand used by Illumina to design probes
 - Forward strand: Forward strand as designated by dbSNP (depends on NCBI genome build)
 - Plus strand: the 5' end of the plus strand is at the tip of the short (p) arm of the chromosome. This designation 
 also relies on genome build and is used by HapMap and 1000 Genomes Project
 - Top: nomenclature developed by Illumina using sequence-based context to assign strand designations that does not 
 change regardless of database or genome assembly used

The native Affymetrix format is the same as the Illumina Top format. SNP Conversion allows it to be converted to 
PLUS format and vice versa (but not to any other format).

**Important notes**:
 - Plus and Forward format are not necessarily the same
 - The AB format does not discriminate between SNPs and indels 


## Variant files

Variant files contain conversion information for the markers in each SNP panel, which is specific to the genomic 
assembly. Variant files can be found in the directory `variant_position_files/[species]/[assembly]/`. A different 
variant file directory can be specified, but it must preserve the `[dir]/[species]/[assembly]/` structure. Also, please 
store the variant files, and only the variant files, in their own directory.

Current acceptable assemblies listed by species are:
- bos_taurus
	- UMD3_1_chromosomes
	- ARS-UCD1_2_Btau5_0_1Y
- sus_scrofa
	- sscrofa11_1
	- sscrofa10_2
	- USMARCv1_0

 
 The assembly name and species must therefore be specified when running check_format and convert_file modules using 
the `--assembly` and `--species` options, respectively. The filename for all variant files will have the structure 
`[panel name].[assembly].[conversion|position].csv[.gz]`. SNP conversion modules will detect the matching or 
best-matching variant file. Ideally, the user input file will contain all the markers present in a single variant file, 
and this file will be selected for conversion. If there is not an exact match between the user input markers and those 
in the variant file, SNP conversion modules will select the file that contains all matching markers, with the fewest 
additional ones. If there are markers in the user input file that are not found in the best-matching variant file, SNP 
conversion will exit and print the orphan markers to the file `[snp panel basename]_problem_variants.txt`.

The option `--get-snp-panel` will return the name of the panel selected by the program.

Variant files must be in csv format, however they may be compressed using zip, gzip, or bzip2. Note that 
*tar.gz and *tar files are not permitted.  

The variant file headers (indicated by lines beginning with '#') contain more information about variant 
file structure. 

## Notes on file conversion

The convert_file module of SNP conversion runs check_format on all input files, and will not convert files with 
incorrect SNPs.

AB format matrix files cannot be converted to other file types.
 

## Notes on options

### Incorrect SNP threshold

By default, 95% of markers in these input files must be correct for the script to make a prediction. However, this 
value can be changed by editing the variable `minimum_correct_snp_fraction` at the top of `lib/check_format`.

### Creating a summary

The `--summary` flag can be used to create a summary of any Illumina or Affymetrix format file when running 
check_format or convert_file modules. For each sample, the summary file reports
 - Total number of variants
 - Number of variants with genotypes (i.e. not '--', '---', or 'NoCall')
 - Number of homozygous SNPs
 - Number of heterozygous SNPs
 - Number of indels
 - Number of consistent genotypes
 - Number of inconsistent genotypes
 
For Long format files, the summary file also reports:
 - Number of equivalent genotypes
 - Number of inequivalent genotypes
 
In Long format, genotype equivalency is a measure of whether all of the alleles for a marker are equivalent across 
 formats. For example, a set of alleles may be correctly formatted as "Top" and therefore considered consistent, but 
 are not equivalent with the alleles reported as "Forward" and "Design". 
 
The default summary output is in a "pretty" format. To output a tab-formatted summary file, use the `--tabular` flag in 
conjunction with `--summary`. 

A summary file cannot be generated when using the merge_files function. Instead, run check_format with the `--summary` 
flag. 

### Creating PLINK MAP and PED flat files

The option `--plink` can be used with check_format or convert_file to generate PED and MAP files for the input panel 
or converted panel, respectively. 

The PED output file is very simple; it contains only the Individual ID and the 
genotypes (in A/C/G/T format), separated by whitespace. To use this file with PLINK, specify the missing fields with 
`--no-fid --no-parents --no-sex --no-pheno`, or input this missing information. 

The MAP file contains the chromosome code, SNP identifier (the marker name), the genetic distance (set 
to 0),  and the base pair position. Chromosome codes are organism specific; to use this file in PLINK, specify `--cow` 
or `--chr-set 29 no-xy` for bos_taurus, and `--chr-set 18 no-xy` for sus_scrofa. All values are separated by whitespace.

For more information on these files, see the [PLINK data format page](http://zzz.bwh.harvard.edu/plink/data.shtml)

## Example commands

Below are some sample commands for common SNP_conversion tasks. Sample data are provided to run each of these commands 
within the `sanple_files` directory. In `sample_files`, running the script `./sample_conversions.sh` will run the 
commands below with the sample data and compare the outputs to files previously generated in `sample_output`. 
Note that you may have to make the script executable using 

```
chmod u+x sample_conversions.sh
``` 

### File format checking

Determine the format of all files in a directory, and return the variant file name (SNP panel) that most likely 
corresponds to each file

```
python SNP_conversion.py check_format --input-dir [input_dir] --assembly [assembly] --species [species] --get-snp-panel
```

Check whether a list of files in Illumina Forward format are correctly formatted and specify the directory
containing variant files

```
python SNP_conversion.py check_format --input-dir [input_dir] --file-list [file1,file2..fileN] --input-format FWD
 --assembly [assembly] --species [species] --key-dir [variant_file_dir]
```

Get a list of inconsistent markers in a file suspected to be in Top format, and write all output to a log file 

```
python SNP_conversion.py check_format --input-dir [input_dir] --file-list [file1] --input-format TOP
 --assembly [assembly] --species [species] --verbose
```

Output a tab-formatted summary file after checking the format of a Long-format file
 
 ```
python SNP_conversion.py check_format --input-dir [input_dir] --file-list [file1] --input-format Long
 --assembly [assembly] --species [species] --summary --tabular
```
The summary option can also be used with the convert_file module.

### File conversion


Convert an Illumina matrix file in Top format to a Long format file without specifying an output suffix
 
 ```
python SNP_conversion.py convert_file --input-dir [input_dir] --file-list [file1] --input-format TOP 
--output-format LONG --assembly [assembly] --species [species]
```

Convert a list of files of unknown or mixed formats to Forward format, specifying the output suffix 'FORWARD'

 ```
python SNP_conversion.py convert_file --input-dir [input_dir] --file-list [file1] --output-format FWD 
--output-name FORWARD --assembly [assembly] --species [species]
```

Convert a file from Affymetrix (native) to Affymetrix Plus format, specifying the number of threads and output 
suffix 'affy_plus'

```
python SNP_conversion.py convert_file --input-dir [input_dir] --file-list [file1] --input-format affymetrix 
--output-format AFFY-PLUS --output-name affy_plus --assembly [assembly] --species [species] --threads 4
```

### Merging files

Merge a list of files in Design format and output the file 'merged_design_files.txt'

```
python SNP_conversion.py merge_files --input-dir [input_dir] --file-list [file1,file2,file3] 
--input-format DESIGN --output-name merged_design_files.txt
```


