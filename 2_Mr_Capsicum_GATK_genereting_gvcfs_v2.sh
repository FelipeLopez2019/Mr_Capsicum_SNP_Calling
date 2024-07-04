#!/bin/bash
#$ -N script_Allium
#$ -V
#$ -l h_vmem=20G
#$ -l h_data=20G
#$ -S /bin/bash
#$ -pe PE 1
#$ -cwd

## Building g.VCF files of each sample per GATK4.0 protocol

for R1_file in /home/frodaf/vg-v1.46.0/input/Capsicum/reads/Trimming_Felipe/Trimmomatic/*_R1_pareado.fastq.gz
do
    # Locating the destination of the process
    R2_file="${R1_file/_R1/_R2}"
    filename=$(basename "$R1_file")
    # Extraer el nombre de la muestra
    SUBSTRING=$(echo "$filename" | cut -d'_' -f 1)
    echo `date`
    echo "using $SUBSTRING"

    # Alignment by bwa
    bwa mem -t 4 -R "@RG\tID:${SUBSTRING}\tSM:${SUBSTRING}\tPL:Illumina" /home/frodaf/vg-v1.46.0/input/Capsicum/reads/Trimming_Felipe/Genomes/GCA_002878395.3_UCD10Xv1.1_genomic.fasta $R1_file $R2_file > ${SUBSTRING}.sam

    # Clean SAM files
    java -jar /home/frodaf/vg-v1.46.0/input/Capsicum/reads/Trimming_Felipe/GATK/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar CleanSam -INPUT ${SUBSTRING}.sam -VALIDATION_STRINGENCY SILENT -OUTPUT ${SUBSTRING}.clean.sam
    rm ${SUBSTRING}.sam

    # Convert from SAM to BAM
    samtools view -S ${SUBSTRING}.clean.sam -b -o ${SUBSTRING}.clean.bam

    # Ensure all mate information is in sync between each read and its mate
    java -Xms2g -Xmx15g -jar /home/frodaf/vg-v1.46.0/input/Capsicum/reads/Trimming_Felipe/GATK/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar FixMateInformation -INPUT ${SUBSTRING}.clean.bam -VALIDATION_STRINGENCY SILENT -OUTPUT ${SUBSTRING}.clean.fix.bam
    rm ${SUBSTRING}.clean.sam ${SUBSTRING}.clean.bam

    # Sort BAM file
    java -jar /home/frodaf/vg-v1.46.0/input/Capsicum/reads/Trimming_Felipe/GATK/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar SortSam -INPUT ${SUBSTRING}.clean.fix.bam -OUTPUT ${SUBSTRING}.sort.bam -SORT_ORDER coordinate -VALIDATION_STRINGENCY SILENT
    rm ${SUBSTRING}.clean.fix.bam

    # Remove duplicates
    java -jar /home/frodaf/vg-v1.46.0/input/Capsicum/reads/Trimming_Felipe/GATK/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar MarkDuplicates -I ${SUBSTRING}.sort.bam -O ${SUBSTRING}.dup.bam -REMOVE_DUPLICATES TRUE -METRICS_FILE ${SUBSTRING}.dup.metrics.txt
    rm ${SUBSTRING}.sort.bam

    # Index BAM file
    samtools index ${SUBSTRING}.dup.bam

    # Get variants (gVCF) for a sample
    java -jar /home/frodaf/vg-v1.46.0/input/Capsicum/reads/Trimming_Felipe/GATK/gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar HaplotypeCaller -I ${SUBSTRING}.dup.bam -R /home/frodaf/vg-v1.46.0/input/Capsicum/reads/Trimming_Felipe/Genomes/GCA_002878395.3_UCD10Xv1.1_genomic.fasta -O ${SUBSTRING}.g.vcf -ERC GVCF -ploidy 2
    rm ${SUBSTRING}.dup.bam

    echo `date`
    echo "finished obtaining g.VCF of each sample"
done