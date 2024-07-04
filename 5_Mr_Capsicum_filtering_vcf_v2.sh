#!/bin/bash
#$ -N script_Allium
#$ --variant
#$ -l h_vmem=20G
#$ -l h_data=20G
#$ -S /bin/bash
#$ -pe PE 1
#$ -cwd



# Define los par�metros de filtrado
missing_individual_threshold=0.3  # Porcentaje m�ximo permitido de datos perdidos por individuo
missing_loci_threshold=0.1      # Porcentaje m�ximo permitido de datos perdidos por loci
maf_threshold=0.05               # Frecuencia al�lica m�nima
min_depth=10                     # Profundidad m�nima de cobertura
max_depth=100                    # Profundidad m�xima de cobertura


# Define el archivo de entrada y salida
input_vcf="/home/frodaf/vg-v1.46.0/input/Capsicum/reads/Trimming_Felipe/GATK/finalisimo3.vcf"
output_vcf="/home/frodaf/vg-v1.46.0/input/Capsicum/reads/Trimming_Felipe/GATK/finalisimo3_clean.vcf"

# Paso 1: Filtrar por porcentaje de datos faltantes
vcftools --vcf $input_vcf --max-missing missing_loci --recode --out filtered_na

# Paso 2: Filtrar por profundidad de cobertura
vcftools --vcf filtered_na.recode.vcf --minDP min_depth --recode --out filtered_na_depth

# Filtra por frecuencia al�lica m�nima
vcftools --vcf filtered_na_depth.recode.vcf --maf $maf_threshold --recode --out $output_vcf

# Limpieza de archivos temporales
rm filtered_na.recode.vcf filtered_na_depth.recode.vcf