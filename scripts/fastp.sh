#!/bin/bash

#fastp.sh

#SBATCH --job-name=fastp.sh
#SBATCH --account=dh_lab_shrec

#SBATCH -o /gscratch/cnyam/soybean/data/cleaned_data/fastp_output.out
#SBATCH -e /gscratch/cnyam/soybean/data/cleaned_data/fastp_error.out

#SBATCH --mail-user=cnyam@uwyo.edu
#SBATCH --mail-type=ALL

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12

#SBATCH --time=24:00:00


module load miniconda3/23.11.0

conda activate biodata


for file in /gscratch/cnyam/soybean/data/*1.fastq.gz

do

filename=$(basename -- $file) 
read=${filename%_1.fastq.gz} 

read1=${read}_1.fastq.gz
read2=${read}_2.fastq.gz
input=/gscratch/cnyam/soybean/data/
output=/gscratch/cnyam/soybean/data/cleaned_data/

echo $read1
echo $read2 
fastp \
  --in1 ${input}${read1} \
  --out1 ${output}${read1} \
  --in2 ${input}${read2} \
  --out2 ${output}${read2} \
  --json ${output}${read}.json \
  --html ${output}${read}.html

done

echo "job complete"
