#!/bin/bash 
#SBATCH --job-name=readgroup
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

mkdir ../readgroup
cd ../readgroup


module load picard/2.9.2
export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratch

java -jar $PICARD AddOrReplaceReadGroups \
        INPUT=../nonduplicates/SRR1517848_nodup.bam \
        OUTPUT=SRR1517848_rg.bam \
        RGID=group1 \
        RGSM=SRR1517848 \
        RGPL=illumina \
        RGLB=1 \
        RGPU=barcode \
        CREATE_INDEX=True


