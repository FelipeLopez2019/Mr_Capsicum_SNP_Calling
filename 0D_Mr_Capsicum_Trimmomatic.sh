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

adaptadores="/home/apolo-lflopezh/Espel_GATK4/adapters"

OUTPUT_DIR="/home/apolo-lflopezh/Reads/Trimmomatic"
# Especifica la ruta completa al directorio de entrada de bbduk_trimming
INPUT_DIR="/home/apolo-lflopezh/Reads"
# Cambia el numero de hilos (-t) segun sea necesario
for R1_file in "$INPUT_DIR"/*_1.fastq.gz
do
	R2_file="${R1_file/_1/_2}"
	name1=$(basename "$R1_file" .fastq.gz)
	name2=$(basename "$R2_file" .fastq.gz)
	name=$(basename "$R1_file" _1.fastq.gz)
	echo $name1
	echo $name2
	trimmomatic PE -phred33 $R1_file $R2_file \
  	"$OUTPUT_DIR"/"$name"_R1_pareado.fastq.gz "$OUTPUT_DIR"/"$name"_R1_no_pareado.fastq.gz \
  	"$OUTPUT_DIR"/"$name"_R2_pareado.fastq.gz "$OUTPUT_DIR"/"$name"_R2_no_pareado.fastq.gz \
  	ILLUMINACLIP:$adaptadores/TruSeq3-PE-2.fa:2:30:10:8:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50 HEADCROP:12 
done
echo $(date)
echo "finalizado"