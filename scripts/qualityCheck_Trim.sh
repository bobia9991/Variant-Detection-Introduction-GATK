#!/bin/bash 
#SBATCH --job-name=quality_control
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=10G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

mkdir ../{fastqc,trimmed}

module load sickle
module load fastqc


#Quality check of Reads

fastqc -t 4 -o ../fastqc ../raw_data/SRR1517848_1.fastq ../raw_data/SRR1517848_2.fastq


#Trimming of reads

sickle pe \
-t sanger \
-f ../raw_data/SRR1517848_1.fastq \
-r ../raw_data/SRR1517848_2.fastq \
-o ../trimmed/trimmed_SRR1517848_1.fastq \
-p ../trimmed/trimmed_SRR1517848_2.fastq \
-l 45 -q 25 \
-s ../trimmed/singles_SRR1517848.fastq

#Quality check of Reads post trimming
fastqc -t 4 -o ../fastqc ../trimmed/trimmed_SRR1517848_1.fastq ../trimmed/trimmed_SRR1517848_2.fastq




