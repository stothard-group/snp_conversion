#!/usr/bin/env python3
import argparse
import sys
import os
import file_format_checker as ffc
import file_conversion as convert
import merge_snp_files as merge

# Get cwd for --input-dir
cwd = os.getcwd()
# This script can call the file format checker and file conversion modules

parent_parser = argparse.ArgumentParser(
    add_help=False, description="Converts between Illumina matrix and long formats"
)
parser = argparse.ArgumentParser(
    add_help=True,
    description="Converts between Illumina matrix and long formats. To see "
    "the full usage for each function:  "
    "python SNP_conversion.py [Function] --help",
)


# Make general options
general_group = parent_parser.add_argument_group("General options")
general_group.add_argument(
    "--input-dir",
    type=str,
    required=False,
    default=cwd,
    help="directory containing input file(s) (default directory: ./)",
)
general_group.add_argument(
    "--file-list",
    type=str,
    help="[optional] comma-separated list of files in the input directory to be converted (no whitespace)",
)
general_group.add_argument(
    "--input-format",
    type=str,
    choices=["TOP", "FWD", "AB", "LONG", "DESIGN", "PLUS", "mixed", "affymetrix"],
    default="mixed",
    help="Type of file(s) expected: 'TOP', 'FWD', 'AB', 'LONG', 'DESIGN', 'PLUS', 'mixed', or 'affymetrix'. "
    "'mixed' may not be used when merging files.",
)

## Parsers for the 3 programs
subs = parser.add_subparsers(dest="command")
subs.metavar = "Functions"
sub_a = subs.add_parser(
    "check_format",
    parents=[parent_parser],
    help="Checks the format of the input file(s)",
)
sub_b = subs.add_parser(
    "convert_file",
    parents=[parent_parser],
    help="Converts input file(s) to user-specified format",
)
sub_c = subs.add_parser(
    "merge_files", parents=[parent_parser], help="Merges files of the same format"
)

## Function-specific options
# File conversion and file format checking

format_group = sub_a.add_argument_group("File Format Check Options")
format_group.add_argument(
    "--get-snp-panel",
    action="store_true",
    default=False,
    required=False,
    help="[optional] Will determine which genotype conversion key files contain all SNPs in the input",
)
format_group.add_argument(
    "--key-dir",
    type=str,
    default="variant_position_files",
    help="Directory containing genotype conversion key files (default directory: variant_position_files)",
)
format_group.add_argument(
    "--species",
    required=True,
    type=str,
    choices=["bos_taurus", "sus_scrofa"],
    help="Organism name",
)
format_group.add_argument(
    "-v",
    "--verbose-logging",
    action="store_true",
    default=False,
    required=False,
    help="[optional] Write output to both STDOUT and log file",
)
format_group.add_argument(
    "--assembly",
    type=str,
    required=True,
    help="Assembly name - see README for full list of choices",
)

format_group.add_argument(
    "-s",
    "--summary",
    action="store_true",
    default=False,
    required=False,
    help="Summarize converted SNP file in *_summary.txt file",
)
format_group.add_argument(
    "--tabular",
    action="store_true",
    default=False,
    required=False,
    help="Output summary file in tabular format (default: False)",
)
format_group.add_argument(
    "--plink",
    action="store_true",
    default=False,
    required=False,
    help="Creates PLINK flat files (PED and MAP) (default: False)",
)

# File conversion
conversion_group = sub_b.add_argument_group("File Conversion Options")
conversion_group.add_argument(
    "--output-format",
    type=str,
    required=True,
    choices=["TOP", "FWD", "AB", "PLUS", "DESIGN", "LONG", "AFFY-PLUS"],
    help="Type of file(s) to be created: 'TOP', 'FWD', 'AB', 'PLUS', 'DESIGN', 'LONG', 'AFFY-PLUS'. "
    "Only one type of output can be specified at a time. "
    "PLUS converts the data to the forward strand of the reference genome. "
    "LONG refers to the long-format Illumina file."
    "AFFY-TOP refers to an Affymetrix-format file with PLUS alleles instead of Affy's native Forward alleles",
)
conversion_group.add_argument(
    "--output-name",
    type=str,
    default="output",
    help="Suffix of output file, without file extension. File will be named [input_file].[output_name].txt "
    "(default = output)",
)
conversion_group.add_argument(
    "-t",
    "--threads",
    type=int,
    default=2,
    required=False,
    help="[optional] Number of threads to use during conversion (default = 2). This option will be ignored on Windows "
    "systems",
)
conversion_group.add_argument(
    "-s",
    "--summary",
    action="store_true",
    default=False,
    required=False,
    help="Summarize converted SNP file in *_summary.txt file",
)
conversion_group.add_argument(
    "--assembly",
    type=str,
    required=True,
    help="Assembly name - see README for full list of choices",
)
conversion_group.add_argument(
    "--key-dir",
    type=str,
    default="variant_position_files",
    help="Directory containing genotype conversion key files (default directory: variant_position_files)",
)
conversion_group.add_argument(
    "--species",
    required=True,
    type=str,
    choices=["bos_taurus", "sus_scrofa"],
    help="Organism name",
)
conversion_group.add_argument(
    "--tabular",
    action="store_true",
    default=False,
    required=False,
    help="Output summary file in tabular format (default: False)",
)
conversion_group.add_argument(
    "-v",
    "--verbose-logging",
    action="store_true",
    default=False,
    required=False,
    help="[optional] Write output to both STDOUT and log file",
)
conversion_group.add_argument(
    "--plink",
    action="store_true",
    default=False,
    required=False,
    help="Creates PLINK flat files (PED and MAP) (default: False)",
)


# Merge files
merge_group = sub_c.add_argument_group("Merge files Options")
merge_group.add_argument(
    "--output-name", type=str, required=True, help="Name of merged file"
)


if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()
input_dir = args.input_dir

if args.file_list is None:
    file_list = os.listdir(input_dir)
else:
    file_list = args.file_list

if args.command == "convert_file":
    if args.input_format == "AB":
        exit("Cannot convert AB format to another format")
    output = convert.convert_file(
        input_dir,
        args.file_list,
        args.input_format,
        args.output_format,
        args.output_name,
        args.threads,
        args.verbose_logging,
        args.key_dir,
        args.summary,
        args.assembly,
        args.tabular,
        args.species,
        args.plink,
    )
elif args.command == "check_format":
    # Send the appropriate commands to the format module
    log_file = None
    file_type, correct_format, return_log = ffc.file_format_check(
        input_dir,
        args.file_list,
        args.input_format,
        args.get_snp_panel,
        args.verbose_logging,
        args.key_dir,
        log_file,
        args.assembly,
        args.summary,
        args.tabular,
        args.species,
        args.plink,
    )
elif args.command == "merge_files":
    if not args.input_format:
        kill_message = "The following argument is required: --input-format"
        sys.exit(kill_message)
    else:
        if args.input_format == "mixed":
            kill_message = "'mixed' is not a valid option when merging files"
            sys.exit(kill_message)
        else:
            pass
        if args.input_format == "affymetrix":
            output = merge.merge_affy_files(
                input_dir, args.file_list, args.input_format, args.output_name
            )
        else:
            output = merge.merge_matrix_files(
                input_dir, args.file_list, args.input_format, args.output_name
            )

else:
    pass
