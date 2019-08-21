# Variant Detection Introduction using GATK 

This repository is a usable, publicly available tutorial for analyzing differential expression data and creating topological gene networks. All steps have been provided for the UConn CBC Xanadu cluster here with appropriate headers for the Slurm scheduler that can be modified simply to run.  Commands should never be executed on the submit nodes of any HPC machine.  If working on the Xanadu cluster, you should use sbatch scriptname after modifying the script for each stage.  Basic editing of all scripts can be performed on the server with tools such as nano, vim, or emacs.  If you are new to Linux, please use <a href="https://bioinformatics.uconn.edu/unix-basics">this</a> handy guide for the operating system commands.  In this guide, you will be working with common bio Informatic file formats, such as <a href="https://en.wikipedia.org/wiki/FASTA_format">FASTA</a>, <a href="https://en.wikipedia.org/wiki/FASTQ_format">FASTQ</a>, <a href="https://en.wikipedia.org/wiki/SAM_(file_format)">SAM/BAM</a>, and <a href="https://en.wikipedia.org/wiki/General_feature_format">GFF3/GTF</a>. You can learn even more about each file format <a href="https://bioinformatics.uconn.edu/resources-and-events/tutorials/file-formats-tutorial/">here</a>. If you do not have a Xanadu account and are an affiliate of UConn/UCHC, please apply for one <a href="https://bioinformatics.uconn.edu/contact-us/">here</a>.

<div id="toc_container">
<p class="toc_title">Contents</p>
<ul class="toc_list">
    <li><a href="#Header_1"> 1. Introduction to Variant Detection using Whole Exome Sequencing</>
    <li><a href="#Header_2"> 2. Sample Data Download </>
    <li><a href="#Header_3"> 3. Quality Check and filtering of the reads</>
    <li><a href="#Header_4"> 4. Preparing the Reference Sequence</>
    <li><a href="#Header_5"> 5. Aligning of Reads</>
    <li><a href="#Header_6"> 6. SAM to BAM Conversion and Remove Singletons</>
    <li><a href="#Header_7"> 7. Sort BAM files using PICARD</>
    <li><a href="#Header_8"> 8. Remove PCR Duplicates using PICARD</>
    <li><a href="#Header_9"> 9. Add Read Group Information</>
    <li><a href="#Header_10"> 10. Reorder BAM file</>
    <li><a href="#Header_11"> 11. Base Recalibration</>
    <li><a href="#Header_12"> 12. Variant Calling</>
</ul>
</div> 



### Introduction to Variant Detection using Whole Exome Sequencing 
  


Genome-wide sequencing methods have provided a method of deep understanding in the  sequence variations between genotype and phenotype. The [1000 genome project](http://www.1000genomes.org)  have indeed made inroads into the study of population genetics, including the investigation of causal variants of genes for various human syndromes. Next Generation Sequencing (NGS) technology is evolving at a rapid rate and new sequencing platforms are been developed. Whole Genome Sequencing (WGS) have been used comprehensively to detecting genomic variations such as single nucleotide variants (SNVs), Copy number variants (CNVs), insertions and deletions (InDels) and chromosomal rearrangements.  However at low cost, whole exome sequencing (WES) platforms, have performed well in high sequencing coverage and readily interpreting protein coding exons associated compared to WGS platforms. This technique has also led to more samples to be analyzed and can generate targeted DNA sequences and identify substantially more genetic variations.   


The exome is the part of the genome composed of exons, the sequences which, when transcribed, remain within the mature RNA after introns are removed by RNA splicing and contribute to the final protein product encoded by that gene. It consists of all DNA that is transcribed into mature RNA in cells of any type, as distinct from the transcriptome, which is the RNA that has been transcribed only in a specific cell population. The exome of the human genome consists of roughly 180,000 exons constituting about 1% of the total genome.


### Sample Data Download 
For the GATK pipe line we will be using SRR1517848 dataset belongs to PRJNA59853 project.  The project was designed by Broad institute as a part of generating haplotype map for human genome.  SRR1517848 dataset is random exon sequencing of genomic DNA (paired-end library) on Illumina HiSeq 2000 platform with read length of 75bps.  The first step will be to download the data from SRA database using `fastq-dump` application of `sratoolkit`.While downloading we will be splitting the SRA format file into two fastq files; forward and reverse (R1 and R2)strand by using `--split-files` option in the command

<pre>
#!/bin/bash
#SBATCH --job-name=data_download
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

mkdir -p ../raw_data

module load sratoolkit/2.8.2
cd ../raw_data

fastq-dump --split-files SRR1517848</pre>
The full script for slurm shedular can be found in the scripts folder by the name <a href="/scripts/data_download.sh">data_download.sh</a>.

The succesful execution of the script will output 2 fastq files in `raw_data` folder.

<pre>
raw_data/
├── SRR1517848_1.fastq
├── SRR1517848_2.fastq
</pre>  


### Quality Check and filtering of the reads

In this step we will check the quality of the reads using `fastqc` and trim reads to remove low quality bases with help of `sickle`.  A detailed account of these tools is explained in the <a href="https://github.com/CBC-UCONN/Variant-Calling-using-freebayes-and-Annotation">Variant-Calling-using-freebayes-and-Annotation</a> tutorial.
<pre>
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

mkdir -p ../{fastqc,trimmed}

module load sickle
module load fastqc

</pre>   
 



#### Quality check of Reads

```bash
fastqc -t 4 -o ../fastqc ../raw_data/SRR1517848_1.fastq ../raw_data/SRR1517848_2.fastq
```   


#### Trimming of reads

```bash
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
```  

The full slurm script is called [qualityCheck_Trim.sh](/scripts/qualityCheck_Trim.sh).  

In this script `qualityCheck_Trim.sh` the dataset is passed through 3 different applications.  First we use`fastqc` to check the quality of the raw reads, second reads are trimmed of low quality bases using `sickle` and lastly check the qualirty of reads post trimming using `fastqc`. Fastqc produces html reports that can be downloaded to your local machine to visualise the results. 

Now these reads are ready for downstream processing. The next step in the workflow is to map reads to human genome using `bwa`.  The human genome is available as fasta file but it requires some pre-processing before we can use it to map reads. The next few steps are aimed at preparing the reference sequence (human genome fasta file).

### Pre-processing steps of Reference Sequence
* Generate BWA index
* Generate Fasta File Index 
* Create a Sequence Dictionary

In here we are preparing these files upfront so the GATK will be able to use the human genome FASTA file as a reference.

### Generate the BWA index  

First run the following bwa command to create the index, given that you have reference chr20.fasta file already downloaded in to the folder `resources` . 

<pre style="color: silver; background: black;">
#!/bin/bash
#SBATCH --job-name=bwa_index
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

module load bwa

cd ../resources

bwa index -p chr20 chr20.fasta

</pre>
The full slurm script for creating the index can be found at scripts folder by the name, <a href="/scripts/bwa_index.sh">bwa_index.sh</a> .

This will create the files listed below. These files will be used by `BWA` while mapping reads to reference. 
<pre>
resources/
├── chr20.fa
├── chr20.amb
├── chr20.ann
├── chr20.bwt
├── chr20.pac
└── chr20.sa
</pre>


### Generate Fasta File Index  
Using `samtools` we will create a index of the reference fasta file.  
<pre>
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

module load samtools

cd ../resources

samtools faidx chr20.fasta</pre>

The full slurm script for creating the index can be found at scripts folder by the name, <a href="/scripts/fasta_index.sh">fasta_index.sh</a>.
This will create:
<pre>
resources/
└──  chr20.fasta.fai
</pre>

It will consist of one record per line for each of the contigs in the fasta file. Where each record is composed of 
* contig name
* size
* location
* bases per lane
* bytes per lane


### Create Sequence Dictionary
Use the picard tools to create the dictionary by the following command:
<pre>
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
        REFERENCE=chr20.fasta \
        OUTPUT=chr20.dict \
        CREATE_INDEX=True
</pre>
The full slurm script for creating the dictionary can be found at scripts folder by the name, <a href="/scripts/dictionary.sh">dictionary.sh</a>.

This will create:
<pre>
resources/
└── chr20.dict
</pre>

This is formated like a SAM file header and when running GATK it automatically looks for these files.

### Aligning of Reads

Using BWA aligner we are going to align the reads to the reference fasta file. As explained in freebayes tutorial, BWA is generally slower than Bowtie2 with similar sensitivity and both tools can perform gapped alignment for the identification of indels and can effectively map paired-end reads. However, BWA is a bit more accurate and provides information on which alignments are trustworthy. Small numbers of bad alignments can result in many false variant calls, so accuracy is paramount, and is the basis for choosing BWA.

Since we have paired-end reads command will look like:
<pre>
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

</pre>
The full scrip is called <a href="/scripts/align.sh">align.sh</a> and can be found in the scripts folder.

<pre>
Usage: bwa mem [options] reference.fasta read1.fa read2.fa -o outputname.sam
mem    The BWA-MEM algorithm performs local alignment  
-t     Number of threads  
</pre>

Alignment will create SAM files, once the alignment is done for all the samples we will end up with:
<pre>
<strong>align</strong>/
├── SRR1517848.sam
</pre>

### SAM to BAM Conversion
The BWA aligner will create a `.sam` file (Sequence Alignment/Map format) with mapping information of the reads. SAM files are human readable and can be large files. For the easy processing through the programs we will convert these files to binary format which is the BAM format using SAMtools.  

So the code will look like:
<pre>module load samtools/1.7
samtools view -@ 8 -bS INPUT.sam > OUTPUT.bam </pre> 

<pre>Useage: samtools view [options] 
-@    number of treads
-b    output in BAM format
-S    input format auto detected
</pre>

This will create BAM format files:
<pre>
<strong>align/</strong>
├── SRR1517848.sam
├── SRR1517848.bam
</pre>


### Remove Singletons
The unmapped reads are called singletons. The samtools flags can be used to remove these tagged reads.  
Using the following command we will remove the singletons:
<pre>
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
</pre>

This will produce BAM files:
<pre>
<strong>align/</strong>
├── SRR1517848.sam
├── SRR1517848.bam
├── SRR1517848_filtered.bam

</pre>

The full slurm script for SAM to BAM conversion and singletons removal can be found at scripts folder by the name <a href="/scripts/singletons.sh">singletons.sh</a>


### Sort BAM files using PICARD

Once the singletons removed the BAM files are sorted using PICARD tools. For GATK analysis the BAM files need to be correctly formatted as well. The correct formatting includes:
* It must be aligned
* It must be sorted in coordinate order
* It must list the read groups with sample names in the header
* Every read must belong to a read group.
* The BAM file must pass Picard ValidateSamFile validation


In order to pass the files into Picard for processing in this step we will sort the aligned reads in the coordinate order.

The command is as follows:
<pre>
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
        CREATE_INDEX=True</pre>

The full script is called <a href="/scripts/picard_sort.sh">picard_sort.sh</a> and can be found at scripts folder.

This will create sorted BAM format files:
<pre>
<strong>align/</strong>
├── SRR1517848.sam
├── SRR1517848.bam
├── SRR1517848_filtered.bam
├── SRR1517848_filtered_sort.bam
</pre>


#### Tip
After sorting your reads you can check whether the reads are sorted accordingly and does it have the `SO: coordinate` flag to satisfy the GATK requirements by using `samtools` command to check the Header:  
`samtools view -H SRR1517848_filtered_sort.bam`   

which will produce: 
<pre>
@HD	VN:1.5	SO:coordinate
@SQ	SN:chr1	LN:249250621
@SQ	SN:chr2	LN:243199373
@SQ	SN:chr3	LN:198022430
@SQ	SN:chr4	LN:191154276
@SQ	SN:chr5	LN:180915260
@SQ	SN:chr6	LN:171115067
@SQ	SN:chr7	LN:159138663
@SQ	SN:chrX	LN:155270560
@SQ	SN:chr8	LN:146364022
@SQ	SN:chr9	LN:141213431
@SQ	SN:chr10	LN:135534747
@SQ	SN:chr11	LN:135006516
@SQ	SN:chr12	LN:133851895
@SQ	SN:chr13	LN:115169878
@SQ	SN:chr14	LN:107349540
@SQ	SN:chr15	LN:102531392
@SQ	SN:chr16	LN:90354753
@SQ	SN:chr17	LN:81195210
@SQ	SN:chr18	LN:78077248
@SQ	SN:chr20	LN:63025520
@SQ	SN:chrY	LN:59373566
@SQ	SN:chr19	LN:59128983
@SQ	SN:chr22	LN:51304566
@SQ	SN:chr21	LN:48129895
@SQ	SN:chr6_ssto_hap7	LN:4928567
@SQ	SN:chr6_mcf_hap5	LN:4833398
. 
.
</pre>


### Remove PCR Duplicates using PICARD

During the sequencing the same DNA molecules can be sequenced multiple times resulting in duplicates. These reads should not be counted as information in variant detection. In this step we will mark the duplicate reads and will remove them.  

Following command will remove the duplicate reads from each sample file.
<pre>
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
</pre>
The full script is called <a href="/scripts/markduplicates.sh">markduplicates.sh</a> and can be found at scripts folder.

This will result in duplicates removed BAM files which will be:
<pre>
<strong>noduplicates/</strong>
├── SRR1517848_metrics.txt
├── SRR1517848_nodup.bai
├── SRR1517848_nodup.bam
</pre>


### Add Read Group Information

In this section we will be adding meta data about the sample. Adding meta data is very important is downstream analysis of your data, and these information is visible to GATK analysis tools. In here we use the minimal read group information for the samples and some are important tags.

In the SAM/BAM file the read group information is indicated in @RG tag which signify the "read group".

* ID : globally unique string which identify this run. Usually this linked to the lane where the data was run.

* SM : associated name in the DNA sample. This will be the sample identifier and it is the most important tag. In GATK all the analysis is done by sample, and this will selects which sample group it will belong to.

* PL : platform used. eg: "Illumina", "Pacbio", "iontorrent"

* LB : an identifier of the library from this DNA was sequenced. This field is important for future reference and quality control. In the case of errors associated with DNA preparation, this will link the data to the laboratory preparation step.

* PU : platform unit identifier for the run. The generic identifier will allow to go back to the machine, time and where it was run. Usually this is a flowcell-barcode-lane unique identifier.

To learn more about SAM tools tags please refer the [SAM tools format](http://samtools.github.io/hts-specs/SAMv1.pdf).

The following Picard tools command will add the read group information to each sample.
<pre>
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
</pre>

The full script is called <a href="/scripts/add_readgroup.sh">add_readgroup.sh</a> and can be found in scripts folder.
The above command will add reads groups to each sample and will created BAM files:
<pre>
readgroup/
├── SRR1517848_rg.bai
├── SRR1517848_rg.bam
</pre>

#### How to check the reads have read group information ?
You can do this by quick samtools and unix commands using:  
`samtools view -H SRR1517848_rg.bam | grep '^@RG'`  
which will give you:
<pre>@RG	ID:group1	LB:1	PL:illumina	SM:SRR1517848	PU:barcode</pre>

The presence of the `@RG` tags indicate the presence of read groups. Each read group has a `SM` tag, indicating the sample from which the reads belonging to that read group originate.
q
In addition to the presence of a read group in the header, each read must belong to one and only one read group. Given the following example reads.


### Reorder BAM file

Samtools `sort` command sorts the chromosomes/contigs in lexicographical order, however GATK downstream analysis requires the chromosmes/contigs to be sorted based on their order of occurance in fasta file.  This steps acieves that.  Reads which are mapped to contigs which are absent in the new reference file are rejected or dropped. This step can run faster if we provide the indexed BAM file.

The following command will reorder the BAM file using PICARD tools:
<pre>
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
        REFERENCE=../resources/chr20.fasta \
        CREATE_INDEX=True
</pre>

The full script is called <a href="/scripts/reorder.sh">reorder.sh</a> and can be found in scripts folder.
This will create karyotype BAM files:
<pre>
<strong>reorder/</strong>
├── SRR1517848_karyotype.bai
├── SRR1517848_karyotype.bam
</pre>

### Base Recalibration
Variant calling algorithms rely heavily on the quality scores assigned to the individual base calls in each sequence read. These scores are per-base estimates of error emitted by the sequencing machines. Unfortunately the scores produced by the machines are subject to various sources of systematic technical error, leading to over- or under-estimated base quality scores in the data. Base quality score recalibration (BQSR) is a process in which we apply machine learning to model these errors empirically and adjust the quality scores accordingly. This allows us to get more accurate base qualities, which in turn improves the accuracy of our variant calls. The base recalibration process involves two key steps: first the program builds a model of covariation based on the data and a set of known variants (which you can bootstrap if there is none available for your organism), then it adjusts the base quality scores in the data based on the model. There is an optional but highly recommended step that involves building a second model and generating before/after plots to visualize the effects of the recalibration process. This is useful for quality control purposes. This tool performs the first step described above: it builds the model of covariation and produces the recalibration table. It operates only at sites that are not in dbSNP; we assume that all reference mismatches we see are therefore errors and indicative of poor base quality. This tool generates tables based on various user-specified covariates (such as read group, reported quality score, cycle, and context). Assuming we are working with a large amount of data, we can then calculate an empirical probability of error given the particular covariates seen at this site, where p(error) = num mismatches / num observations. The output file is a table (of the several covariate values, number of observations, number of mismatches, empirical quality score).
The process take place in two steps In first step the recalibration statistics are calculated and then applied to bam file in second step.

For this step we will require a vcf file of known varaints that are know about the human geneome `chr20.dbsnp138.vcf`.

<pre>

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
export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratchi

# Baserecalibration takes place in two step

# STEP1

gatk BaseRecalibrator \
   -I ../reorder/SRR1517848_karyotype.bam \
   -R ../resources/chr20.fasta \
   --known-sites /UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Detection-Introduction-GATK/resources/Homo_sapiens_assembly38.dbsnp138.vcf \
   -O recal_data.table

# STEP2

 gatk ApplyBQSR \
   -R ../resources/chr20.fasta \
   -I ../reorder/SRR1517848_karyotype.bam \
   --bqsr-recal-file recal_data.table \
   -O SRR1517848_recalb.bam

</pre>

The full script is called <a href="/scripts/reorder.sh">BaseRecalibration.sh</a> and can be found in scripts folder.


### Variant Calling 

In this step we will call the variants using HaplotypeCaller in GATK software. 

<pre>
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
        --reference ../resources/chr20.fasta \
        --input ../baserecalbn/SRR1517848_recalb.bam \
        --output SRR1517848_haplotype.vcf
</pre>

The full script is called <a href="/scripts/haplotypeCaller.sh">haplotypeCaller.sh</a> and can be found in scripts folder.
This creates a VCF file called fileName_haplotype.vcf, containing all the variant sites the HaplotypeCaller evaluate including both SNPs and Indels.
<pre>
<strong>haplotypes/</strong>
├── SRR1517848_haplotype.vcf
├── SRR1517848_haplotype.vcf.idx
</pre>

A preview of the vcf file is shown here
<pre>

##fileformat=VCFv4.2
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##GATKCommandLine=<ID=Hapl
..
..
..
##contig=<ID=HLA-DRB1*16:02:01,length=11005,assembly=38>
##source=HaplotypeCaller
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SRR1517848
chr20   158721  .       G       A       17.81   .       AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=8.91;SOR=0.693  GT:AD:DP:GQ:PL  1/1:0,2:2:6:45,6,0
chr20   174042  .       C       T       15.65   .       AC=2;AF=1.00;AN=2;DP=1;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=58.00;QD=15.65;SOR=1.609 GT:AD:DP:GQ:PL  1/1:0,1:1:3:42,3,0
chr20   174047  .       G       C       15.65   .       AC=2;AF=1.00;AN=2;DP=1;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=58.00;QD=15.65;SOR=1.609 GT:AD:DP:GQ:PL  1/1:0,1:1:3:42,3,0
chr20   174049  .       C       G       15.65   .       AC=2;AF=1.00;AN=2;DP=1;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=58.00;QD=15.65;SOR=1.609 GT:AD:DP:GQ:PL  1/1:0,1:1:3:42,3,0
chr20   174065  .       T       A       18.59   .       AC=2;AF=1.00;AN=2;DP=1;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=58.00;QD=18.59;SOR=1.609 GT:AD:DP:GQ:PL  1/1:0,1:1:3:45,3,0
chr20   174068  .       C       T       18.59   .       AC=2;AF=1.00;AN=2;DP=1;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=58.00;QD=18.59;SOR=1.609 GT:AD:DP:GQ:PL  1/1:0,1:1:3:45,3,0
chr20   174075  .       C       A       18.59   .       AC=2;AF=1.00;AN=2;DP=1;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=58.00;QD=18.59;SOR=1.609 GT:AD:DP:GQ:PL  1/1:0,1:1:3:45,3,0
</pre>


The variants identified can be furter filtered or sorted on different parameters that suits the downstream analysis.  <a href="https://software.broadinstitute.org/gatk/">GATK</a> has a range of tools available to achieve that.  A more detailed account of various tools can be found in the "Best Practices" section on the <a href="https://software.broadinstitute.org/gatk/">GATK</a> website.

