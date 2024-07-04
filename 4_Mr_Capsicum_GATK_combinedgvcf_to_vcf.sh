#!/bin/bash
#$ -N script_Allium
#$ --variant
#$ -l h_vmem=20G
#$ -l h_data=20G
#$ -S /bin/bash
#$ -pe PE 1
#$ -cwd

java -jar /home/frodaf/vg-v1.46.0/input/Capsicum/reads/Trimming_Felipe/GATK/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar GenotypeGVCFs \
   -R /home/frodaf/vg-v1.46.0/input/Capsicum/reads/Trimming_Felipe/Genomes/GCA_002878395.3_UCD10Xv1.1_genomic.fasta \
   -V combined.g.vcf \
   -O finalisimo3.vcf