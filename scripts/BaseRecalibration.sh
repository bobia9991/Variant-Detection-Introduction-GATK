#!/bin/bash 
#SBATCH --job-name=BaseCalbiration
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

mkdir ../baserecalbn
cd ../baserecalbn

module load GATK/4.0
export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratch

# Baserecalibration takes place in two step

# STEP1
gatk BaseRecalibrator \
        -I ../reorder/SRR1517848_karyotype.bam \
        -R ../resources/chr20.fasta \
        --known-sites ../resources/chr20_dbSNPref.vcf \
        -O recal_data.table

# STEP2
gatk ApplyBQSR \
        -R ../resources/chr20.fasta \
        -I ../reorder/SRR1517848_karyotype.bam \
        --bqsr-recal-file recal_data.table \
        -O SRR1517848_recalb.bam



