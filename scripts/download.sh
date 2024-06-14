#!/bin/bash

#download.sh PRJNA1011827

#SBATCH --job-name=download.sh
#SBATCH --account=dh_lab_shrec

#SBATCH -o /gscratch/cnyam/Transcriptomics/fastp_output.out
#SBATCH -e /gscratch/cnyam/Transcriptomics/fastp_error.out

#SBATCH --mail-user=cnyam@uwyo.edu
#SBATCH --mail-type=ALL

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12

#SBATCH --time=24:00:00



curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/021/SRR25880321/SRR25880321_1.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880321_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/021/SRR25880321/SRR25880321_2.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880321_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/017/SRR25880317/SRR25880317_1.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880317_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/017/SRR25880317/SRR25880317_2.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880317_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/018/SRR25880318/SRR25880318_1.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880318_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/018/SRR25880318/SRR25880318_2.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880318_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR25880304/SRR25880304_1.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880304_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR25880304/SRR25880304_2.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880304_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/010/SRR25880310/SRR25880310_1.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880310_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/010/SRR25880310/SRR25880310_2.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880310_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/008/SRR25880308/SRR25880308_1.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880308_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/008/SRR25880308/SRR25880308_2.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880308_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/011/SRR25880311/SRR25880311_1.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880311_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/011/SRR25880311/SRR25880311_2.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880311_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/016/SRR25880316/SRR25880316_1.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880316_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/016/SRR25880316/SRR25880316_2.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880316_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/007/SRR25880307/SRR25880307_1.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880307_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/007/SRR25880307/SRR25880307_2.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880307_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/005/SRR25880305/SRR25880305_1.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880305_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/005/SRR25880305/SRR25880305_2.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880305_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/009/SRR25880309/SRR25880309_1.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880309_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/009/SRR25880309/SRR25880309_2.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880309_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR25880306/SRR25880306_1.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880306_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR25880306/SRR25880306_2.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880306_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/019/SRR25880319/SRR25880319_1.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880319_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/019/SRR25880319/SRR25880319_2.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880319_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/014/SRR25880314/SRR25880314_1.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880314_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/014/SRR25880314/SRR25880314_2.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880314_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/020/SRR25880320/SRR25880320_1.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880320_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/020/SRR25880320/SRR25880320_2.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880320_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/013/SRR25880313/SRR25880313_1.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880313_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/013/SRR25880313/SRR25880313_2.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880313_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/015/SRR25880315/SRR25880315_1.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880315_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/015/SRR25880315/SRR25880315_2.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880315_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/012/SRR25880312/SRR25880312_1.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880312_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/012/SRR25880312/SRR25880312_2.fastq.gz -o /gscratch/cnyam/soybean/data/SRR25880312_2.fastq.gz
