#!/bin/sh
#SBATCH --partition=longjobs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=5-00:00:00
#SBATCH --job-name=Allium01
#SBATCH -o %x_%j.out      # File to which STDOUT will be written
#SBATCH -e %x_%j.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL


module load python/3.9_miniconda-4.10.3
source activate gatk4

## Building g.VCF files of each sample per GATK4.5.0.0 protocol
# alignment by bwa

OUTPUT_DIR="/home/serverxxxx-lflopezh/Capsicum_GATK4"
# Especifica la ruta completa al directorio de entrada de bbduk_trimming
INPUT_DIR="/home/serverxxxx-lflopezh/Reads/Trimmomatic"
# Cambia el numero de hilos (-t) segun sea necesario
for R1_file in "$INPUT_DIR"/*_R1_pareado.fastq.gz
do
    R2_file="${R1_file/_R1/_R2}"
    name1=$(basename "$R1_file" _pareado.fastq.gz)
    name2=$(basename "$R2_file" _pareado.fastq.gz)
    SUBSTRING=$(basename "$R1_file" _R1_pareado.fastq.gz)
    # alignment by bwa
    ## example shows: ${SUBSTRING}.fastq.gz
    bwa mem -R '@RG\tID:${SUBSTRING}\tSM:${SUBSTRING}\tPL:Illumina' /home/serverxxxx-lflopezh/genome/GCA_002878395.3_UCD10Xv1.1_genomic.fasta \
    $R1_file $R2_file > ${SUBSTRING}.sam
    #clean sam files
    gatk CleanSam -INPUT ${SUBSTRING}.sam -VALIDATION_STRINGENCY SILENT -OUTPUT ${SUBSTRING}.clean.sam
    rm ${SUBSTRING}.sam
    #This tool ensures that all mate information is in sync between each read and its mate.
    gatk FixMateInformation -INPUT ${SUBSTRING}.clean.sam -VALIDATION_STRINGENCY SILENT -OUTPUT ${SUBSTRING}.clean.fix.sam
    rm ${SUBSTRING}.clean.sam
    #Convert from sam to bam
    samtools view -S ${SUBSTRING}.clean.fix.sam -b -o ${SUBSTRING}.clean.fix.bam
    rm ${SUBSTRING}.clean.fix.sam
    #order alignment
    #via samtools
    #samtools sort -o G1_1.sort.bam -O BAM  G1_1.bam
    #via gatk
    gatk SortSam -INPUT ${SUBSTRING}.clean.fix.bam -OUTPUT ${SUBSTRING}.sort.bam -SORT_ORDER coordinate -VALIDATION_STRINGENCY SILENT
    rm ${SUBSTRING}.clean.fix.bam
    #statistics
    #chmod 0777 Bamtools/1.0
    #bamtools stats -in ${SUBSTRING}.sort.bam > ${SUBSTRING}.sort.bam.STATS_CON_DUP
    #Remove duplicates
    gatk MarkDuplicates -I ${SUBSTRING}.sort.bam -O ${SUBSTRING}.dup.bam -REMOVE_DUPLICATES TRUE  -METRICS_FILE ${SUBSTRING}.dup.metrics.txt
    rm ${SUBSTRING}.sort.bam
    #Statistics after cleaning
    bamtools stats -in ${SUBSTRING}.dup.bam > ${SUBSTRING}.dup.bam.STATS_SIN_DUP
    # index for gatk
    samtools index ${SUBSTRING}.dup.bam
    # second file is dictionary reference
    #via samtools
    #rm GCF_000499845.1_PhaVulg1_0_genomic.dict
    #samtools dict -o /home/serverxxxx-lflopezh/genome/AlliumCepa_v1.2._onlychr.dict /home/serverxxxx-lflopezh/genome/GCA_002878395.3_UCD10Xv1.1_genomic.fasta
    #via gatk
    #rm GCF_000499845.1_PhaVulg1_0_genomic.dict
    # Get variants (gvcf) for a sample
    gatk HaplotypeCaller -I ${SUBSTRING}.dup.bam -R /home/serverxxxx-lflopezh/genome/GCA_002878395.3_UCD10Xv1.1_genomic.fasta -O ${SUBSTRING}.g.vcf -ERC GVCF -ploidy 2
    rm ${SUBSTRING}.dup.bam
    echo `date`
    echo "finished obtaining g.VCF of each sample"
    #rm ${SUBSTRING}.g.vcf.gz
done

# After the loop in each sample continue to collapse the gvcf to vcf
# gatk GenotypeGVCFs -R /home/llopez/Allium/GATK/reference/GenomeNPU2023.fasta -V Sample1.g.vcf -V Sample2.g.vcf -V Sample86.g.vcf -O All_samples.vcf
# echo `date`
# echo "completed obtaining molecular variants"
sed -i 's/$SUBSTRING/91/g'  91.g.vcf
sed -i 's/$SUBSTRING/263/g' 263.g.vcf
sed -i 's/$SUBSTRING/163/g' 163.g.vcf
sed -i 's/$SUBSTRING/119/g' 119.g.vcf
sed -i 's/$SUBSTRING/21/g'  21.g.vcf
sed -i 's/$SUBSTRING/148/g' 148.g.vcf
sed -i 's/$SUBSTRING/244/g' 244.g.vcf
sed -i 's/$SUBSTRING/153/g' 153.g.vcf
sed -i 's/$SUBSTRING/126/g' 126.g.vcf
sed -i 's/$SUBSTRING/125/g' 125.g.vcf
sed -i 's/$SUBSTRING/131/g' 131.g.vcf
sed -i 's/$SUBSTRING/185/g' 185.g.vcf
sed -i 's/$SUBSTRING/17/g'  17.g.vcf
sed -i 's/$SUBSTRING/254/g' 254.g.vcf
sed -i 's/$SUBSTRING/212/g' 212.g.vcf
sed -i 's/$SUBSTRING/99/g'  99.g.vcf
sed -i 's/$SUBSTRING/50/g'  50.g.vcf
sed -i 's/$SUBSTRING/277/g' 277.g.vcf
sed -i 's/$SUBSTRING/232/g' 232.g.vcf
sed -i 's/$SUBSTRING/105/g' 105.g.vcf
sed -i 's/$SUBSTRING/71/g'  71.g.vcf
sed -i 's/$SUBSTRING/251/g' 251.g.vcf
sed -i 's/$SUBSTRING/228/g' 228.g.vcf
sed -i 's/$SUBSTRING/143/g' 143.g.vcf
sed -i 's/$SUBSTRING/112/g' 112.g.vcf
sed -i 's/$SUBSTRING/24/g'  24.g.vcf
sed -i 's/$SUBSTRING/250/g' 250.g.vcf
sed -i 's/$SUBSTRING/230/g' 230.g.vcf
sed -i 's/$SUBSTRING/281/g' 281.g.vcf
sed -i 's/$SUBSTRING/179/g' 179.g.vcf
sed -i 's/$SUBSTRING/178/g' 178.g.vcf
sed -i 's/$SUBSTRING/102/g' 102.g.vcf
sed -i 's/$SUBSTRING/223/g' 223.g.vcf
sed -i 's/$SUBSTRING/76/g'  76.g.vcf
sed -i 's/$SUBSTRING/133/g' 133.g.vcf
sed -i 's/$SUBSTRING/264/g' 264.g.vcf
sed -i 's/$SUBSTRING/70/g'  70.g.vcf
sed -i 's/$SUBSTRING/150/g' 150.g.vcf
sed -i 's/$SUBSTRING/142/g' 142.g.vcf
sed -i 's/$SUBSTRING/144/g' 144.g.vcf
sed -i 's/$SUBSTRING/86/g'  86.g.vcf
sed -i 's/$SUBSTRING/61/g'  61.g.vcf
sed -i 's/$SUBSTRING/106/g' 106.g.vcf
sed -i 's/$SUBSTRING/180/g' 180.g.vcf
sed -i 's/$SUBSTRING/245/g' 245.g.vcf
sed -i 's/$SUBSTRING/98/g'  98.g.vcf
sed -i 's/$SUBSTRING/224/g' 224.g.vcf
sed -i 's/$SUBSTRING/27/g'  27.g.vcf
sed -i 's/$SUBSTRING/26/g'  26.g.vcf
sed -i 's/$SUBSTRING/124/g' 124.g.vcf
sed -i 's/$SUBSTRING/140/g' 140.g.vcf
sed -i 's/$SUBSTRING/201/g' 201.g.vcf
sed -i 's/$SUBSTRING/166/g' 166.g.vcf
sed -i 's/$SUBSTRING/139/g' 139.g.vcf
sed -i 's/$SUBSTRING/65/g'  65.g.vcf
sed -i 's/$SUBSTRING/229/g' 229.g.vcf
sed -i 's/$SUBSTRING/33/g'  33.g.vcf
sed -i 's/$SUBSTRING/97/g'  97.g.vcf
sed -i 's/$SUBSTRING/38/g'  38.g.vcf
sed -i 's/$SUBSTRING/118/g' 118.g.vcf
sed -i 's/$SUBSTRING/5/g'   5.g.vcf
sed -i 's/$SUBSTRING/243/g' 243.g.vcf
sed -i 's/$SUBSTRING/259/g' 259.g.vcf
sed -i 's/$SUBSTRING/3/g'   3.g.vcf
sed -i 's/$SUBSTRING/141/g' 141.g.vcf
sed -i 's/$SUBSTRING/108/g' 108.g.vcf
sed -i 's/$SUBSTRING/89/g'  89.g.vcf
sed -i 's/$SUBSTRING/165A/g' 165A.g.vcf
sed -i 's/$SUBSTRING/248/g' 248.g.vcf
sed -i 's/$SUBSTRING/35/g'  35.g.vcf
sed -i 's/$SUBSTRING/15/g'  15.g.vcf
sed -i 's/$SUBSTRING/100/g' 100.g.vcf
sed -i 's/$SUBSTRING/195/g' 195.g.vcf
sed -i 's/$SUBSTRING/270/g' 270.g.vcf
sed -i 's/$SUBSTRING/176/g' 176.g.vcf
sed -i 's/$SUBSTRING/117/g' 117.g.vcf
sed -i 's/$SUBSTRING/48/g'  48.g.vcf
sed -i 's/$SUBSTRING/4/g'   4.g.vcf
sed -i 's/$SUBSTRING/113/g' 113.g.vcf
sed -i 's/$SUBSTRING/165B/g' 165B.g.vcf
sed -i 's/$SUBSTRING/20/g'  20.g.vcf
sed -i 's/$SUBSTRING/147/g' 147.g.vcf
sed -i 's/$SUBSTRING/289/g' 289.g.vcf
sed -i 's/$SUBSTRING/181/g' 181.g.vcf
sed -i 's/$SUBSTRING/156/g' 156.g.vcf
sed -i 's/$SUBSTRING/279/g' 279.g.vcf
sed -i 's/$SUBSTRING/116/g' 116.g.vcf
sed -i 's/$SUBSTRING/32/g'  32.g.vcf
sed -i 's/$SUBSTRING/187/g' 187.g.vcf
sed -i 's/$SUBSTRING/67/g'  67.g.vcf
sed -i 's/$SUBSTRING/184/g' 184.g.vcf
sed -i 's/$SUBSTRING/41/g'  41.g.vcf
sed -i 's/$SUBSTRING/128/g' 128.g.vcf
sed -i 's/$SUBSTRING/273/g' 273.g.vcf
sed -i 's/$SUBSTRING/193/g' 193.g.vcf
sed -i 's/$SUBSTRING/18/g'  18.g.vcf
sed -i 's/$SUBSTRING/57/g'  57.g.vcf
sed -i 's/$SUBSTRING/92/g'  92.g.vcf
sed -i 's/$SUBSTRING/28/g'  28.g.vcf
sed -i 's/$SUBSTRING/171/g' 171.g.vcf
sed -i 's/$SUBSTRING/84/g'  84.g.vcf
sed -i 's/$SUBSTRING/60/g'  60.g.vcf
sed -i 's/$SUBSTRING/217/g' 217.g.vcf
sed -i 's/$SUBSTRING/69/g'  69.g.vcf
sed -i 's/$SUBSTRING/227/g' 227.g.vcf
sed -i 's/$SUBSTRING/111/g' 111.g.vcf
sed -i 's/$SUBSTRING/278/g' 278.g.vcf
sed -i 's/$SUBSTRING/154/g' 154.g.vcf
sed -i 's/$SUBSTRING/137/g' 137.g.vcf
sed -i 's/$SUBSTRING/51/g'  51.g.vcf
sed -i 's/$SUBSTRING/25/g'  25.g.vcf
sed -i 's/$SUBSTRING/132/g' 132.g.vcf
sed -i 's/$SUBSTRING/10/g'  10.g.vcf
sed -i 's/$SUBSTRING/87/g'  87.g.vcf
sed -i 's/$SUBSTRING/19/g'  19.g.vcf
sed -i 's/$SUBSTRING/276/g' 276.g.vcf
sed -i 's/$SUBSTRING/88/g'  88.g.vcf
sed -i 's/$SUBSTRING/240/g' 240.g.vcf
sed -i 's/$SUBSTRING/136/g' 136.g.vcf
sed -i 's/$SUBSTRING/206/g' 206.g.vcf
sed -i 's/$SUBSTRING/266/g' 266.g.vcf
sed -i 's/$SUBSTRING/31/g'  31.g.vcf
sed -i 's/$SUBSTRING/267/g' 267.g.vcf
sed -i 's/$SUBSTRING/107/g' 107.g.vcf
sed -i 's/$SUBSTRING/127/g' 127.g.vcf
sed -i 's/$SUBSTRING/194/g' 194.g.vcf
sed -i 's/$SUBSTRING/85/g'  85.g.vcf
sed -i 's/$SUBSTRING/16/g'  16.g.vcf
sed -i 's/$SUBSTRING/78/g'  78.g.vcf
sed -i 's/$SUBSTRING/253/g' 253.g.vcf
sed -i 's/$SUBSTRING/157/g' 157.g.vcf
sed -i 's/$SUBSTRING/93/g'  93.g.vcf
sed -i 's/$SUBSTRING/149/g' 149.g.vcf
sed -i 's/$SUBSTRING/151/g' 151.g.vcf
sed -i 's/$SUBSTRING/231/g' 231.g.vcf
sed -i 's/$SUBSTRING/160/g' 160.g.vcf
sed -i 's/$SUBSTRING/34/g'  34.g.vcf
sed -i 's/$SUBSTRING/129/g' 129.g.vcf
sed -i 's/$SUBSTRING/115/g' 115.g.vcf
sed -i 's/$SUBSTRING/205/g' 205.g.vcf
sed -i 's/$SUBSTRING/222/g' 222.g.vcf
sed -i 's/$SUBSTRING/260/g' 260.g.vcf
sed -i 's/$SUBSTRING/55/g'  55.g.vcf
sed -i 's/$SUBSTRING/202/g' 202.g.vcf
sed -i 's/$SUBSTRING/161/g' 161.g.vcf
sed -i 's/$SUBSTRING/274/g' 274.g.vcf
sed -i 's/$SUBSTRING/255/g' 255.g.vcf
sed -i 's/$SUBSTRING/72/g'  72.g.vcf
sed -i 's/$SUBSTRING/122A/g' 122A.g.vcf
sed -i 's/$SUBSTRING/164/g' 164.g.vcf
sed -i 's/$SUBSTRING/268/g' 268.g.vcf
sed -i 's/$SUBSTRING/12/g'  12.g.vcf
sed -i 's/$SUBSTRING/130/g' 130.g.vcf
sed -i 's/$SUBSTRING/158/g' 158.g.vcf
sed -i 's/$SUBSTRING/200/g' 200.g.vcf
sed -i 's/$SUBSTRING/265/g' 265.g.vcf
sed -i 's/$SUBSTRING/210/g' 210.g.vcf
sed -i 's/$SUBSTRING/1/g'   1.g.vcf
sed -i 's/$SUBSTRING/155/g' 155.g.vcf
sed -i 's/$SUBSTRING/174/g' 174.g.vcf
sed -i 's/$SUBSTRING/103/g' 103.g.vcf
sed -i 's/$SUBSTRING/196/g' 196.g.vcf
sed -i 's/$SUBSTRING/114/g' 114.g.vcf
sed -i 's/$SUBSTRING/177/g' 177.g.vcf
sed -i 's/$SUBSTRING/123/g' 123.g.vcf
sed -i 's/$SUBSTRING/233/g' 233.g.vcf
sed -i 's/$SUBSTRING/192/g' 192.g.vcf
sed -i 's/$SUBSTRING/110/g' 110.g.vcf
sed -i 's/$SUBSTRING/246/g' 246.g.vcf
sed -i 's/$SUBSTRING/29/g'  29.g.vcf
sed -i 's/$SUBSTRING/215/g' 215.g.vcf
sed -i 's/$SUBSTRING/234/g' 234.g.vcf
sed -i 's/$SUBSTRING/242/g' 242.g.vcf
sed -i 's/$SUBSTRING/46/g'  46.g.vcf
sed -i 's/$SUBSTRING/168/g' 168.g.vcf
sed -i 's/$SUBSTRING/182/g' 182.g.vcf
sed -i 's/$SUBSTRING/44/g'  44.g.vcf
sed -i 's/$SUBSTRING/96/g'  96.g.vcf
sed -i 's/$SUBSTRING/241/g' 241.g.vcf
sed -i 's/$SUBSTRING/252/g' 252.g.vcf
sed -i 's/$SUBSTRING/53/g'  53.g.vcf
sed -i 's/$SUBSTRING/220/g' 220.g.vcf
sed -i 's/$SUBSTRING/95/g'  95.g.vcf
sed -i 's/$SUBSTRING/288/g' 288.g.vcf
sed -i 's/$SUBSTRING/45/g'  45.g.vcf
sed -i 's/$SUBSTRING/214/g' 214.g.vcf
sed -i 's/$SUBSTRING/198/g' 198.g.vcf
sed -i 's/$SUBSTRING/101/g' 101.g.vcf
sed -i 's/$SUBSTRING/37/g'  37.g.vcf
sed -i 's/$SUBSTRING/282/g' 282.g.vcf
sed -i 's/$SUBSTRING/262/g' 262.g.vcf
sed -i 's/$SUBSTRING/63/g'  63.g.vcf
sed -i 's/$SUBSTRING/146/g' 146.g.vcf
sed -i 's/$SUBSTRING/209/g' 209.g.vcf
sed -i 's/$SUBSTRING/191/g' 191.g.vcf
sed -i 's/$SUBSTRING/261/g' 261.g.vcf
sed -i 's/$SUBSTRING/297/g' 297.g.vcf
sed -i 's/$SUBSTRING/236/g' 236.g.vcf
sed -i 's/$SUBSTRING/52/g'  52.g.vcf
sed -i 's/$SUBSTRING/145/g' 145.g.vcf
sed -i 's/$SUBSTRING/249/g' 249.g.vcf
sed -i 's/$SUBSTRING/291/g' 291.g.vcf
sed -i 's/$SUBSTRING/293/g' 293.g.vcf
sed -i 's/$SUBSTRING/294/g' 294.g.vcf
sed -i 's/$SUBSTRING/54/g'  54.g.vcf
sed -i 's/$SUBSTRING/219/g' 219.g.vcf
sed -i 's/$SUBSTRING/30/g'  30.g.vcf
sed -i 's/$SUBSTRING/109/g' 109.g.vcf
sed -i 's/$SUBSTRING/73/g'  73.g.vcf
sed -i 's/$SUBSTRING/290/g' 290.g.vcf
sed -i 's/$SUBSTRING/216/g' 216.g.vcf
sed -i 's/$SUBSTRING/221/g' 221.g.vcf
sed -i 's/$SUBSTRING/122B/g' 122B.g.vcf
sed -i 's/$SUBSTRING/64/g'  64.g.vcf
sed -i 's/$SUBSTRING/22/g'  22.g.vcf
sed -i 's/$SUBSTRING/272/g' 272.g.vcf
sed -i 's/$SUBSTRING/58/g'  58.g.vcf
sed -i 's/$SUBSTRING/190/g' 190.g.vcf
sed -i 's/$SUBSTRING/134/g' 134.g.vcf
sed -i 's/$SUBSTRING/189/g' 189.g.vcf
sed -i 's/$SUBSTRING/300/g' 300.g.vcf
sed -i 's/$SUBSTRING/104/g' 104.g.vcf
sed -i 's/$SUBSTRING/188/g' 188.g.vcf
sed -i 's/$SUBSTRING/283/g' 283.g.vcf
sed -i 's/$SUBSTRING/23/g'  23.g.vcf
sed -i 's/$SUBSTRING/295/g' 295.g.vcf
sed -i 's/$SUBSTRING/138/g' 138.g.vcf
sed -i 's/$SUBSTRING/208/g' 208.g.vcf
sed -i 's/$SUBSTRING/66/g'  66.g.vcf
sed -i 's/$SUBSTRING/284/g' 284.g.vcf
sed -i 's/$SUBSTRING/49/g'  49.g.vcf
sed -i 's/$SUBSTRING/257/g' 257.g.vcf
sed -i 's/$SUBSTRING/6/g'   6.g.vcf
sed -i 's/$SUBSTRING/258/g' 258.g.vcf
sed -i 's/$SUBSTRING/287/g' 287.g.vcf
sed -i 's/$SUBSTRING/36/g'  36.g.vcf
sed -i 's/$SUBSTRING/169/g' 169.g.vcf
sed -i 's/$SUBSTRING/121/g' 121.g.vcf
sed -i 's/$SUBSTRING/203/g' 203.g.vcf
sed -i 's/$SUBSTRING/256/g' 256.g.vcf
sed -i 's/$SUBSTRING/81/g'  81.g.vcf
sed -i 's/$SUBSTRING/275/g' 275.g.vcf
sed -i 's/$SUBSTRING/59/g'  59.g.vcf
sed -i 's/$SUBSTRING/286/g' 286.g.vcf
sed -i 's/$SUBSTRING/11/g'  11.g.vcf
sed -i 's/$SUBSTRING/199/g' 199.g.vcf
sed -i 's/$SUBSTRING/40/g'  40.g.vcf
sed -i 's/$SUBSTRING/152/g' 152.g.vcf
sed -i 's/$SUBSTRING/285/g' 285.g.vcf
sed -i 's/$SUBSTRING/292/g' 292.g.vcf
sed -i 's/$SUBSTRING/298/g' 298.g.vcf
sed -i 's/$SUBSTRING/296/g' 296.g.vcf
sed -i 's/$SUBSTRING/204B2/g' 204B2.g.vcf
sed -i 's/$SUBSTRING/120/g' 120.g.vcf
sed -i 's/$SUBSTRING/39/g'  39.g.vcf
sed -i 's/$SUBSTRING/237/g' 237.g.vcf
sed -i 's/$SUBSTRING/170/g' 170.g.vcf
sed -i 's/$SUBSTRING/68/g'  68.g.vcf
sed -i 's/$SUBSTRING/239/g' 239.g.vcf
sed -i 's/$SUBSTRING/74/g'  74.g.vcf
sed -i 's/$SUBSTRING/94/g'  94.g.vcf
sed -i 's/$SUBSTRING/47/g'  47.g.vcf
sed -i 's/$SUBSTRING/218A/g' 218A.g.vcf
sed -i 's/$SUBSTRING/280/g' 280.g.vcf
sed -i 's/$SUBSTRING/235/g' 235.g.vcf
sed -i 's/$SUBSTRING/269/g' 269.g.vcf
sed -i 's/$SUBSTRING/175/g' 175.g.vcf
sed -i 's/$SUBSTRING/77/g'  77.g.vcf
sed -i 's/$SUBSTRING/219B/g' 219B.g.vcf
sed -i 's/$SUBSTRING/211/g' 211.g.vcf
sed -i 's/$SUBSTRING/204B1/g' 204B1.g.vcf
sed -i 's/$SUBSTRING/226/g' 226.g.vcf
sed -i 's/$SUBSTRING/299/g' 299.g.vcf
sed -i 's/$SUBSTRING/90/g'  90.g.vcf
sed -i 's/$SUBSTRING/213/g' 213.g.vcf
sed -i 's/$SUBSTRING/82/g'  82.g.vcf
sed -i 's/$SUBSTRING/238/g' 238.g.vcf
sed -i 's/$SUBSTRING/204A/g' 204A.g.vcf
sed -i 's/$SUBSTRING/62/g'  62.g.vcf
sed -i 's/$SUBSTRING/167/g' 167.g.vcf
sed -i 's/$SUBSTRING/197/g' 197.g.vcf
sed -i 's/$SUBSTRING/225/g' 225.g.vcf
sed -i 's/$SUBSTRING/75/g'  75.g.vcf