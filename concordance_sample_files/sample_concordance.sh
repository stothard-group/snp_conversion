#!/bin/bash

set -e
if [[ ! -d test_output ]]; then
    mkdir test_output
fi

echo "Running test commands - please ignore system UserWarnings"

########################################################################################################################
# The commands in this block can be run individually within the sample_concordance directory, or from within this script

# Concordance analysis between an Illumina LONG file and a VCF file, with tabular output
python ../genotype_concordance --snp-panel input_files/G_CCGP_long_sample_input.txt --panel-type LONG --vcf-file input_files/SNPs_reduced_anon.vcf.gz --species bos_taurus --assembly ARS-UCD1_2_Btau5_0_1Y --output-type tabular --output concordance_test1 > test_output/long_vs_vcf_tab_concordance.txt

# Concordance analysis between an Illumina LONG file and a VCF file, filtering on quality values, with pretty output
python ../genotype_concordance --snp-panel input_files/G_CCGP_long_sample_input.txt --panel-type LONG --vcf-file input_files/SNPs_reduced_anon.vcf.gz --species bos_taurus --assembly ARS-UCD1_2_Btau5_0_1Y --filter-vcf --qual 100 --output-type pretty --output concordance_q100_test2 > test_output/long_vs_vcf_q100_pretty_concordance.txt

# Concordance analysis between an Affymetrix file and a VCF file, outputting a list of discordant positions
python ../genotype_concordance --snp-panel input_files/G_CCGP_affy_short_input.txt --panel-type affymetrix --vcf-file input_files/SNPs_reduced_anon.vcf.gz --species bos_taurus --assembly ARS-UCD1_2_Btau5_0_1Y --extract-discordant --output concordance_affy > test_output/affy_vs_vcf_concordance.txt

########################################################################################################################


# move any stray output files
mv concordance_test1_metastatistics.txt test_output/concordance_test1_metastatistics.txt
mv concordance_test1_statistics.txt test_output/concordance_test1_statistics.txt
mv concordance_q100_test2_metastatistics.txt test_output/concordance_q100_test2_metastatistics.txt
mv concordance_q100_test2_statistics.txt test_output/concordance_q100_test2_statistics.txt
mv concordance_affy_metastatistics.txt test_output/concordance_affy_metastatistics.txt
mv concordance_affy_statistics.txt test_output/concordance_affy_statistics.txt
mv G_CCGP_affy_short_input_discordant.txt test_output/G_CCGP_affy_short_input_discordant.txt


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
    	diff -u -I='.+Checking the format of the SNP panel' -I='.+Finding the matching variant file' ${old_file_match} ${new_files[$i]}
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