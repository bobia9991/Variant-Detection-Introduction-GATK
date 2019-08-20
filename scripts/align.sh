#!/bin/bash 
#SBATCH --job-name=align
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=40G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

mkdir -p ../align

cd ../align

module load bwa/0.7.17
bwa mem -t 4 ../resources/chr20.fasta ../trimmed/trimmed_SRR1517848_1.fastq ../trimmed/trimmed_SRR1517848_2.fastq  -o SRR1517848.sam

#SAM to BAM CONVERSION

module load samtools/1.7
samtools view -@ 4 -bS SRR1517848.sam >> SRR1517848.bam


