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


# Define los parámetros de filtrado
missing_individual_threshold=0.3  # Porcentaje máximo permitido de datos perdidos por individuo
missing_loci_threshold=0.1      # Porcentaje máximo permitido de datos perdidos por loci
maf_threshold=0.05               # Frecuencia alélica mínima
min_depth=10                     # Profundidad mínima de cobertura
max_depth=100                    # Profundidad máxima de cobertura


cftools --gzvcf /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR1.vcf --max-missing 0.7 --minDP 3 --recode --out /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR1.depured.vcf
vcftools --gzvcf /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR2.vcf --max-missing 0.7 --minDP 3 --recode --out /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR2.depured.vcf
vcftools --gzvcf /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR3.vcf --max-missing 0.7 --minDP 3 --recode --out /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR3.depured.vcf
vcftools --gzvcf /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR4.vcf --max-missing 0.7 --minDP 3 --recode --out /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR4.depured.vcf
vcftools --gzvcf /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR5.vcf --max-missing 0.7 --minDP 3 --recode --out /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR5.depured.vcf
vcftools --gzvcf /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR6.vcf --max-missing 0.7 --minDP 3 --recode --out /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR6.depured.vcf
vcftools --gzvcf /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR7.vcf --max-missing 0.7 --minDP 3 --recode --out /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR7.depured.vcf
vcftools --gzvcf /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR8.vcf --max-missing 0.7 --minDP 3 --recode --out /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR8.depured.vcf
vcftools --gzvcf /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR9.vcf --max-missing 0.7 --minDP 3 --recode --out /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR9.depured.vcf
vcftools --gzvcf /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR10.vcf --max-missing 0.7 --minDP 3 --recode --out /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR10.depured.vcf
vcftools --gzvcf /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR11.vcf --max-missing 0.7 --minDP 3 --recode --out /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR11.depured.vcf
vcftools --gzvcf /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR12.vcf --max-missing 0.7 --minDP 3 --recode --out /home/serverxxxx-lflopezh/Capsicum_GATK4/VCF_final_CHR12.depured.vcf