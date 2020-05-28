#!/bin/bash
set -e
if [[ ! -d test_output ]]; then
    mkdir test_output
fi

echo "Running test commands - please ignore system UserWarnings"


########################################################################################################################
# The commands in this block can be run individually within the test directory, or from within this script

## File format checking
# Determine the format of all files in the sample directory and retrieve SNP panel
python ../snp_conversion check_format --input-dir input_files/conversion --assembly UMD3_1_chromosomes --get-snp-panel --species bos_taurus > test_output/check_input_files_output.txt
# Check whether Illumina Forward format file is correctly formatted, specifying the key file directory. Output corresponding PED and MAP files
python ../snp_conversion check_format --input-dir input_files/conversion --file-list 50kv3_mFWD_14June2019.txt --plink --input-format FWD --assembly UMD3_1_chromosomes --conversion variant_position_files --species bos_taurus > test_output/check_forward_files_output.txt
# Get a list of inconsistent markers in a file suspected to be in Top format, and write all output to a log file
python ../snp_conversion check_format --input-dir input_files/conversion --file-list 50kv3_mTOP_inconsistent.txt --input-format TOP --assembly UMD3_1_chromosomes --verbose --species bos_taurus > test_output/check_inconsistent_files_output.txt

# Output a tab-formatted summary file after checking the format of a Long-format file
python ../snp_conversion check_format --input-dir input_files/conversion --file-list 50kv3_Long_14June2019.txt --input-format LONG --assembly UMD3_1_chromosomes --summary --tabular --species bos_taurus > test_output/check_long_format_file_output.txt

## File Conversion
# Convert an Illumina matrix file in Top format to a Long format file without specifying an output suffix
python ../snp_conversion convert_file --input-dir input_files/conversion --file-list 50kv3_mTOP_14June2019.txt --input-format TOP  --output-format LONG --assembly UMD3_1_chromosomes --species bos_taurus > test_output/convert_top_to_long_output.txt
# Convert a list of files of unknown or mixed formats to Forward format, specifying the output suffix 'FORWARD'
python ../snp_conversion convert_file --input-dir input_files/conversion --file-list 50kv3_Long_14June2019.txt,50kv3_mTOP_14June2019.txt --output-format FWD --output-name FORWARD --assembly UMD3_1_chromosomes --species bos_taurus > test_output/convert_mixed_to_forward_output.txt
# Convert a file from Affymetrix (native) to Affymetrix Plus format, specifying the number of threads and output suffix 'affy_plus'
python ../snp_conversion convert_file --input-dir input_files/conversion --file-list G_CCGP_Axiom_sample.txt --input-format AFFY --output-format AFFY-PLUS --output-name affy_plus --assembly UMD3_1_chromosomes --species bos_taurus --threads 2 > test_output/convert_affy_native_to_plus_output.txt

## Merging files
# Merge a list of files in Forward format and output the file 'merged_forward_files.txt'
python ../snp_conversion merge_files --input-dir input_files/conversion --file-list 50kv3_mFWD_part1.txt,50kv3_mFWD_part2.txt --input-format FWD --output merged_forward_files.txt > test_output/merge_forward_files_output.txt


########################################################################################################################
# The commands in this block can be run individually within the test directory, or from within this script

# Concordance analysis between an Illumina LONG file and a VCF file, with tabular output
python ../snp_conversion genotype_concordance  --snp-panel input_files/concordance/G_CCGP_long_sample_input.txt --panel-type LONG --vcf input_files/concordance/SNPs_reduced_anon.vcf.gz --species bos_taurus --assembly ARS-UCD1_2_Btau5_0_1Y --output-type tabular --output concordance_test1 > test_output/long_vs_vcf_tab_concordance.txt

# Concordance analysis between an Illumina LONG file and a VCF file, filtering on quality values, with pretty output
python ../snp_conversion genotype_concordance  --snp-panel input_files/concordance/G_CCGP_long_sample_input.txt --panel-type LONG --vcf input_files/concordance/SNPs_reduced_anon.vcf.gz --species bos_taurus --assembly ARS-UCD1_2_Btau5_0_1Y --filter-vcf --qual 100 --output-type pretty --output concordance_q100_test2 > test_output/long_vs_vcf_q100_pretty_concordance.txt

# Concordance analysis between an Affymetrix file and a VCF file, outputting a list of discordant positions
python ../snp_conversion genotype_concordance  --snp-panel input_files/concordance/G_CCGP_affy_short_input.txt --panel-type AFFY --vcf input_files/concordance/SNPs_reduced_anon.vcf.gz --species bos_taurus --assembly ARS-UCD1_2_Btau5_0_1Y --extract-discordant --output concordance_affy > test_output/affy_vs_vcf_concordance.txt

########################################################################################################################
# The commands in this block can be run individually within the test directory, or from within this script

# Generate a VCF file from a SNP panel file
python ../snp_conversion vcf_generator --snp-panel input_files/concordance/G_CCGP_affy_short_input.txt --species bos_taurus --assembly ARS-UCD1_2_Btau5_0_1Y --panel-type AFFY > test_output/snp_to_vcf.txt

########################################################################################################################

# move any stray output files
mv 50kv3_mTOP_inconsistent*.log test_output
mv 50kv3_Long_14June2019_summary.txt test_output
mv 50kv3_mTOP_14June2019.output.txt test_output
mv 50kv3_*.FORWARD.txt test_output
mv G_CCGP_Axiom_sample.affy_plus.txt test_output
mv merged_forward_files.txt test_output
mv 50kv3_mFWD_14June2019.map test_output
mv 50kv3_mFWD_14June2019.ped test_output
mv concordance_test1_metastatistics.txt test_output/concordance_test1_metastatistics.txt
mv concordance_test1_statistics.txt test_output/concordance_test1_statistics.txt
mv concordance_q100_test2_metastatistics.txt test_output/concordance_q100_test2_metastatistics.txt
mv concordance_q100_test2_statistics.txt test_output/concordance_q100_test2_statistics.txt
mv concordance_affy_metastatistics.txt test_output/concordance_affy_metastatistics.txt
mv concordance_affy_statistics.txt test_output/concordance_affy_statistics.txt
mv G_CCGP_affy_short_input_discordant.txt test_output/G_CCGP_affy_short_input_discordant.txt
mv G_CCGP_affy_short_input.vcf test_output/G_CCGP_affy_short_input.vcf



#compare new output to sample output
new_output=test_output
old_output=sample_output
new_files=($( find ${new_output} -type f -print0 | perl -ne 'my @files = split(/\0/, $_); foreach(@files) { if (!($_ =~ m/\.svn/)) {print "$_\n";}}'))
#make sure there is a line for diff to ignore
ignore_line='Conversion time for file.+'
ignore_line2='.+Checking the format of the SNP panel'
ignore_line3='.+Finding the matching variant file'

for (( i=0; i<${#new_files[@]}; i++ ));
do
    old_file=${old_output}`echo "${new_files[$i]}" | perl -nl -e 's/^[^\/]+//;' -e 'print $_'`
    # make sure this can handle the log files
    if [[ ${new_files[$i]: -4} == ".log" ]];
    then
    	new_file_fix=`echo "${new_files[$i]}" | perl -nl -e  's/-[0-9]{14}.log//;' -e 'print$_'`
    	new_file_basename=`echo "${new_file_fix}" | perl -nl -e 's/.+[\/]//;' -e 'print $_'`
		old_file_fix=($( find ${old_output} -type f -print0 | perl -ne 'my @files = split(/\0/, $_); foreach(@files) { if ($_ =~ m/\.log/)  {print "$_\n";}}'))
		for (( m=0; m<${#old_file_fix[@]}; m++ ));
		do
			if [[ ${old_file_fix[$m]} =~ $new_file_basename ]];
			then
				old_file_match=${old_file_fix[$m]}
			fi
		done
		echo "Comparing ${old_file_match} to ${new_files[$i]}"
    	set +e
    	diff -u ${old_file_match} ${new_files[$i]}
    	if [[ $? -eq 0 ]]; then
			echo "No differences found"
    	fi
    	set -e
    else
		# compare non-log files
		echo "Comparing ${old_file} to ${new_files[$i]}"
		set +e
		diff -u -I='[Conversion time for file *]' ${old_file} ${new_files[$i]}
		if [[ $? -eq 0 ]]; then
			echo "No differences found"
		fi
    fi
    set -e
done