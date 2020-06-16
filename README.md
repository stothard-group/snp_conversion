# snp_conversion

- FILE: SNP_conversion.py
- AUTH: Emily Herman (eherman@ualberta.ca)
- DATE: APR 7, 2020
- VERS: 2.0

**Please see the [snp_conversion wiki](https://github.com/stothard-group/snp_conversion/wiki) for more information about the program.**

## General Description

snp_conversion is a collection of tools for working with genotype data from Illumina and Affymetrix SNP panels. The program can:
* Check that all SNPs in the panel are consistently formatted
* Determine the format of a panel
* Convert a consistently formatted panel to a different format
* Merge panels of the same format
* Create PLINK flat files from SNP panel genotype data
* Assess concordance of genotype data in the panel file with that of a user-specified VCF file
* Generate a VCF file from a SNP panel file

#### Usage information
```
usage: snp_conversion [-h] Tools ...

A selection of tools for working with genotype data from Illumina and
Affymetrix SNP panels To see the full usage for each Utility: ./snp_convserion
[Utility] --help

positional arguments:
  Tools
    check_format        Checks the format of the input file(s)
    convert_file        Converts input file(s) to user-specified format
    merge_files         Merges files of the same format
    conversion_list     Lists all species and assembly names available
    genotype_concordance
                        Determines concordance between SNP panel genotypes and
                        VCF file genotypes
    vcf_generator       Converts SNP panel file to VCF format file
```

There are five main tools within snp_conversion:
1. [check_format](https://github.com/stothard-group/snp_conversion/wiki/Format-Checking)
2. [convert_file](https://github.com/stothard-group/snp_conversion/wiki/File-Conversion)
3. [merge_files](https://github.com/stothard-group/snp_conversion/wiki/Merge-Files)
4. [genotype_concordance](https://github.com/stothard-group/snp_conversion/wiki/Genotype-Concordance)
5. [vcf_generator](https://github.com/stothard-group/snp_conversion/wiki/VCF-Generator)

See [Install page](https://github.com/stothard-group/snp_conversion/wiki/Installation) and [Quick Start Guide](https://github.com/stothard-group/snp_conversion/wiki/Quick-Start) to get started using the snp_conversion tools with a test dataset.

Most snp_conversion tools rely on genotype conversion files to check input data formatting. See the [Genotype Conversion Files](https://github.com/stothard-group/snp_conversion/wiki/Genotype-Conversion-Files) section for information about generating these files.
