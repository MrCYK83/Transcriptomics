#!/bin/bash

#test_fastp.sh

#SBATCH --job-name=test_fastp.sh
#SBATCH --account=dh_lab_shrec

#SBATCH -o /gscratch/cnyam/soybean/test_output.out
#SBATCH -e /gscratch/cnyam/soybean/test_error.out

#SBATCH --mail-user=cnyam@uwyo.edu
#SBATCH --mail-type=ALL

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12

#SBATCH --time=24:00:00


module load miniconda3/23.11.0

conda activate biodata

fastp --help
