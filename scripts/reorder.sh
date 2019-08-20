#!/bin/bash 
#SBATCH --job-name=reorder
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

mkdir ../reorder
cd ../reorder


module load picard/2.9.2
export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratch

java -jar $PICARD ReorderSam \
        INPUT=../readgroup/SRR1517848_rg.bam \
        OUTPUT=SRR1517848_karyotype.bam \
        REFERENCE=../resources/Homo_sapiens_assembly38.fasta \
        CREATE_INDEX=True

