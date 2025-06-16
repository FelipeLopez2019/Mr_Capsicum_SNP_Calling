#!/bin/sh
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=13-00:00:00
#SBATCH --job-name=Allium01
#SBATCH -o %x_%j.out      # File to which STDOUT will be written
#SBATCH -e %x_%j.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL


module load python/3.9_miniconda-4.10.3
source activate gatk4

gatk GatherVcfs \
-I /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR1.depured.vcf.recode.vcf \
-I /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR2.depured.vcf.recode.vcf \
-I /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR3.depured.vcf.recode.vcf \
-I /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR4.depured.vcf.recode.vcf \
-I /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR5.depured.vcf.recode.vcf \
-I /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR6.depured.vcf.recode.vcf \
-I /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR7.depured.vcf.recode.vcf \
-I /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR8.depured.vcf.recode.vcf \
-I /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR9.depured.vcf.recode.vcf \
-I /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR10.depured.vcf.recode.vcf \
-I /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR11.depured.vcf.recode.vcf \
-I /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR12.depured.vcf.recode.vcf \
-O /home/serverxxxx-lflopezh/Capsicum_GATK4/Mr_Capsicum_final.vcf.gz
