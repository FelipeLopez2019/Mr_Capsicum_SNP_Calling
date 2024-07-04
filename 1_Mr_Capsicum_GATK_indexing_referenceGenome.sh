#!/bin/bash
#$ -N script_Allium
#$ -V
#$ -l h_vmem=20G
#$ -l h_data=20G
#$ -S /bin/bash
#$ -pe PE 1
#$ -cwd


# cd /home/llopez/Allium/GATK/reference
# Genome indexing for BWA
bwa index /home/frodaf/vg-v1.46.0/input/Capsicum/reads/Trimming_Felipe/Genomes/GCA_002878395.3_UCD10Xv1.1_genomic.fasta
# Creation of 2 new index files for reference (in formats used by Samtools and GATK)
samtools faidx /home/frodaf/vg-v1.46.0/input/Capsicum/reads/Trimming_Felipe/Genomes/GCA_002878395.3_UCD10Xv1.1_genomic.fasta
java -jar /home/frodaf/vg-v1.46.0/input/Capsicum/reads/Trimming_Felipe/GATK/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar CreateSequenceDictionary -R /home/frodaf/vg-v1.46.0/input/Capsicum/reads/Trimming_Felipe/Genomes/GCA_002878395.3_UCD10Xv1.1_genomic.fasta -O /home/frodaf/vg-v1.46.0/input/Capsicum/reads/Trimming_Felipe/Genomes/GCA_002878395.3_UCD10Xv1.1_genomic.dict
