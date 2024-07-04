#!/bin/bash
#$ -N script_Allium
#$ --variant
#$ -l h_vmem=20G
#$ -l h_data=20G
#$ -S /bin/bash
#$ -pe PE 1
#$ -cwd

# Define los archivos de entrada .g.vcf
gvcfs=("100.g.vcf"
"101.g.vcf"
"102.g.vcf"
"103.g.vcf"
"104.g.vcf"
"105.g.vcf"
"108.g.vcf"
"109.g.vcf"
"10.g.vcf"
"110.g.vcf"
"111.g.vcf"
"112.g.vcf"
"113.g.vcf"
"114.g.vcf"
"115.g.vcf"
"116.g.vcf"
"117.g.vcf"
"118.g.vcf"
"119.g.vcf"
"11.g.vcf"
"120.g.vcf"
"121.g.vcf"
"122A.g.vcf"
"122B.g.vcf"
"123.g.vcf"
"124.g.vcf")

# Define el archivo de referencia
reference="/home/frodaf/vg-v1.46.0/input/Capsicum/reads/Trimming_Felipe/Genomes/GCA_002878395.3_UCD10Xv1.1_genomic.fasta"


# Define el archivo de salida combinado
combined_gvcf="/home/frodaf/vg-v1.46.0/input/Capsicum/reads/Trimming_Felipe/GATK/combined.g.vcf"

# Construye el comando con todos los archivos de entrada
input_args=""
for gvcf in "${gvcfs[@]}"; do
    input_args+=" --variant $gvcf"
done

# Ejecuta CombineGVCFs para combinar los archivos g.vcf
java -jar /home/frodaf/vg-v1.46.0/input/Capsicum/reads/Trimming_Felipe/GATK/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar CombineGVCFs \
    -R $reference \
    $input_args \
    -O $combined_gvcf