#!/bin/sh
#SBATCH --partition=longjobs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=13-00:00:00
#SBATCH --job-name=Allium01
#SBATCH -o %x_%j.out      # File to which STDOUT will be written
#SBATCH -e %x_%j.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL


module load python/3.9_miniconda-4.10.3
source activate gatk4

# Define el archivo de referencia
reference="/home/serverxxxx-lflopezh/genome/GCA_002878395.3_UCD10Xv1.1_genomic.fasta"
gatk  GenotypeGVCFs \
   -R $reference \
   -V combined.g.vcf.gz \
   -L CHR7 \
   -O VCF_final_CHR7.vcf