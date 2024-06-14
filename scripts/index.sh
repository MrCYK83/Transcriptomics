#!/bin/bash

#index.sh

#SBATCH --job-name=index.sh
#SBATCH --account=dh_lab_shrec

#SBATCH -o /gscratch/cnyam/soybean/docs/assembly_index_output.out
#SBATCH -e /gscratch/cnyam/soybean/docs/assemble_index_error.out

#SBATCH --mail-user=cnyam@uwyo.edu
#SBATCH --mail-type=ALL

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12

#SBATCH --time=24:00:00


module load miniconda3/23.11.0

conda activate biodata

hisat2-build \
/gscratch/cnyam/soybean/data/reference/assembly/Gmax_275_v2.0.fa \
/gscratch/cnyam/soybean/data/reference/assembly/Gmax_275_v2.0_hisat2.idx

echo "Job complete"
