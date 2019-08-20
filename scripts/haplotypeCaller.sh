#!/bin/bash 
#SBATCH --job-name=hapCaller
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

mkdir ../haplotypes
cd ../haplotypes

module load GATK/4.0
export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratch

gatk HaplotypeCaller \
        --reference ../resources/Homo_sapiens_assembly38.fasta \
        --input ../baserecalbn/SRR1517848_recalb.bam \
        --output SRR1517848_haplotype.vcf

