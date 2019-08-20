#!/bin/bash 
#SBATCH --job-name=bam_sort
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

module load picard/2.9.2
export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratch

java -jar $PICARD SortSam \
        INPUT=SRR1517848_filtered.bam \
        OUTPUT=SRR1517848_filtered_sort.bam \
        SORT_ORDER=coordinate \
        CREATE_INDEX=True


