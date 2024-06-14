#! /bin/bash

output_file="combined_file.txt"

# Create an empty file or clear the existing one

> "$output_file"

cut -f1 counts/SRR25880304/featurecounts_results.txt | grep -v ^# | paste > $output_file


# Iterations through each directory
for file in counts/*
do

doc=${file}/featurecounts_results.txt

# Using paste to combine the files side by side
grep -v "^#" $doc | cut -f7 > file1.txt
paste "$output_file" file1.txt > temp_file.txt

# Replace the original combined file with the temporary file
mv temp_file.txt "$output_file"


done

# Clean up the column headers
sed s/_sorted.bam//g combined_file.txt > combined_final.txt
sed -i s+/gscratch/cnyam/soybean/output/bam/++g combined_final.txt
sed -i s/.Wm82.a2.v1//g combined_final.txt
