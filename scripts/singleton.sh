#!/bin/bash 
#SBATCH --job-name=rm_singleton
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=40G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date


cd ../align

module load samtools

samtools view -@ 8 -F 0x04 -b SRR1517848.bam > SRR1517848_filtered.bam

