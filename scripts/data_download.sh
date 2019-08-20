#!/bin/bash 
#SBATCH --job-name=data_download
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load sratoolkit/2.8.2
cd ../raw_data

fastq-dump --split-files SRR1517848
