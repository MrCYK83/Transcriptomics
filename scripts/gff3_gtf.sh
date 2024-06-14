#!/bin/bash

#gff3_gtf.sh

#SBATCH --job-name=gff3_gtf.sh
#SBATCH --account=dh_lab_shrec

#SBATCH -o /gscratch/cnyam/soybean/docs/gff3_gtf_output.out
#SBATCH -e /gscratch/cnyam/soybean/docs/gff3_gtf_error.out

#SBATCH --mail-user=cnyam@uwyo.edu
#SBATCH --mail-type=ALL

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12

#SBATCH --time=24:00:00


module load miniconda3/23.11.0

conda activate biodata

gffread \
/gscratch/cnyam/soybean/data/reference/annotation/Gmax_275_Wm82.a2.v1.gene.gff3 \
-T -o \
/gscratch/cnyam/soybean/data/reference/annotation/Gmax_275_Wm82.a2.v1.gene.gtf


echo "Job Complete"
