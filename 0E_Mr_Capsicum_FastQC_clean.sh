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

fastqc -o /home/serverxxxx-lflopezh/FastQC --nogroup -t 2 /home/serverxxxx-lflopezh/Reads/Trimmomatic/*_pareado.fastq.gz

echo `date`
echo "finalizado"
