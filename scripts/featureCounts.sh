#!/bin/bash

#featureCounts.sh

#SBATCH --job-name=featureCounts.sh
#SBATCH --account=dh_lab_shrec

#SBATCH -o /gscratch/cnyam/soybean/docs/featureCounts_output.out
#SBATCH -e /gscratch/cnyam/soybean/docs/featureCounts_error.out

#SBATCH --mail-user=cnyam@uwyo.edu
#SBATCH --mail-type=ALL

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12

#SBATCH --time=24:00:00


module load miniconda3/23.11.0

conda activate biodata



# Loop to iterate through the fastq files
for file in /gscratch/cnyam/soybean/output/bam/*.bam
do

# Gets only the file name without all the extensions
filename=$(basename -- $file)
read=${filename%_sorted.bam}

# Assigning the variables needed to perform the alignment \
# read1, read2, the reference, and output directory
bam=${read}_sorted.bam
annotation=/gscratch/cnyam/soybean/data/reference/annotation/Gmax_275_Wm82.a2.v1.gene.gtf
mkdir -p /gscratch/cnyam/soybean/output/counts/${read}
output_dir=/gscratch/cnyam/soybean/output/counts/${read}
input=/gscratch/cnyam/soybean/output/bam/


echo "Working on ${read}"

featureCounts -p -T 10 --countReadPairs -C -B \
-a $annotation \
-o ${output_dir}/featurecounts_results.txt \
${input}${bam}

done

echo "***COMPLETED***"

