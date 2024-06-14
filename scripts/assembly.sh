#!/bin/bash

#assembly.sh

#SBATCH --job-name=assembly.sh
#SBATCH --account=dh_lab_shrec

#SBATCH -o /gscratch/cnyam/soybean/docs/assembly_output.out
#SBATCH -e /gscratch/cnyam/soybean/data/assembly_error.out

#SBATCH --mail-user=cnyam@uwyo.edu
#SBATCH --mail-type=ALL

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12

#SBATCH --time=24:00:00


module load miniconda3/23.11.0

conda activate biodata


for file in /gscratch/cnyam/soybean/data/cleaned_data/*1.fastq.gz

do

# Gets only the file name without all the extensions
filename=$(basename -- $file)
read=${filename%_1.fastq.gz}

# Assigning the variables needed to perform the alignment \
# read1, read2, the reference, and output directory
read1=${read}_1.fastq.gz
read2=${read}_2.fastq.gz
input=/gscratch/cnyam/soybean/data/cleaned_data/
reference=/gscratch/cnyam/soybean/data/reference/assembly/Gmax_275_v2.0_hisat2.idx
output_dir=/gscratch/cnyam/soybean/output/bam/

echo "Working on ${read}"

# Aligning using hisat2, converting from sam to bam, \
# sorting and indexing using samtools
hisat2 -x ${reference} \
-1 ${input}${read1} \
-2 ${input}${read2} | \
samtools view -S -b - | \
samtools sort -o ${output_dir}${read}_sorted.bam - &&
samtools index ${output_dir}${read}_sorted.bam

done

echo "Job Completed"

