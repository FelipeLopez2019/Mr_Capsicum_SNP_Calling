#!/bin/sh
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=2-00:00:00
#SBATCH --job-name=Allium01
#SBATCH -o %x_%j.out      # File to which STDOUT will be written
#SBATCH -e %x_%j.err      # File to which STDERR will be written

module load python/3.9_miniconda-4.10.3
source activate gatk4

cd /home/serverxxxx-lflopezh/genome/Espeletia
# Genome indexing for BWA
bwa index /home/serverxxxx-lflopezh/genome/GCA_002878395.3_UCD10Xv1.1_genomic.fasta
# Creation of 2 new index files for reference (in formats used by Samtools and GATK)
samtools faidx /home/serverxxxx-lflopezh/genome/GCA_002878395.3_UCD10Xv1.1_genomic.fasta
gatk CreateSequenceDictionary -R /home/serverxxxx-lflopezh/genome/GCA_002878395.3_UCD10Xv1.1_genomic.fasta -O /home/serverxxxx-lflopezh/genome/GCA_002878395.3_UCD10Xv1.1_genomic.dict