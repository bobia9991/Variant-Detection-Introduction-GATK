#!/bin/bash 
#SBATCH --job-name=deduplication
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

mkdir ../nonduplicates

cd ../nonduplicates

module load picard/2.9.2
export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratch

java -jar $PICARD MarkDuplicates \
        INPUT=../align/SRR1517848_filtered_sort.bam \
        OUTPUT=SRR1517848_nodup.bam \
        REMOVE_DUPLICATES=Ture \
        METRICS_FILE=SRR1517848_metrics.txt \
        CREATE_INDEX=True

