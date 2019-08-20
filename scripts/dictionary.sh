#!/bin/bash 
#SBATCH --job-name=fasta_index
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

cd ../resources

module load picard/2.9.2
export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratch

java -jar $PICARD CreateSequenceDictionary \
        REFERENCE=Homo_sapiens_assembly38.fasta \
        OUTPUT=Homo_sapiens_assembly38.dict \
        CREATE_INDEX=True
