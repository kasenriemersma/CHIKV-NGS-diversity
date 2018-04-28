#!/bin/sh

echo "What is the path to the directory for the output files? e.g. /Path/To/Output"

read outputdir

echo "What CHIKV variant is this? WT, HIFI, DM HIFI, LOFI, or LOFI_MUT?"

read variant

echo "What is the name of the raw read 1 file? (fdsaf_r1)"

read read1

echo "What is the name of the raw read 2 file? (fdsaf_r1)"

read read2

echo "What is the name of the sample?"

read sample

cd $outputdir

# QA raw reads

fastqc $read1.fastq

fastqc $read2.fastq

# Quality trim and adapter trim raw reads

java -jar ~/Bioinformatics/Trimmomatic0_36/trimmomatic-0.36.jar PE $read1.fastq $read2.fastq paired_$read1.fastq unpaired_$read1.fastq paired_$read2.fastq unpaired_$read2.fastq ILLUMINACLIP:/Users/kasenriemersma/Bioinformatics/Trimmomatic0_36/adapters/TruSeq3-PE.fa:2:30:10 HEADCROP:15 LEADING:35 SLIDINGWINDOW:4:35 MINLEN:50

# Primer trim reads

cutadapt -g ^CACACGTAGCCTACCAGTTTC -g ^GAGGACTAGAATCAAATGG -g ^GTGCGGCTTCTTCAATATG -g ^ACTGAATGCAGCCTTCGTAG -g ^AGCGTTGATGGCGAGATAC \
-g ^CACTACAGGAAGTACCAATGG -g ^GCAGAAGCCGACAGCAAGTA -g ^CAACGAGCCGTATAAGTATTGG -g ^TTCAGCCTGGACACCTTTC -g ^CCCGCTCTGTCCTCAAGCTG \
-g ^GATGCGGCTGCTGTCATGAC -g ^CACTGTTCTTAAAGGACTC -g ^GAACGTCTCGCTTGATGCT -g ^CCACATAGGTATGCTGTCGCC -g ^GTTTGGGTTGGGATGAACT \
-g ^GCCATACCCACCATCGACAG -g ^CGCTCTTACCGGGTTTGTTG -g ^GGAAGAGTTCGGTATGCTATG -a GAAACTGGTAGGCTACGTGTG$ -a CCATTTGATTCTAGTCCTC$ \
-a CATATTGAAGAAGCCGCAC$ -a CTACGAAGGCTGCATTCAGT$ -a GTATCTCGCCATCAACGCT$ -a CCATTGGTACTTCCTGTAGTG$ -a TACTTGCTGTCGGCTTCTGC$ \
-a CCAATACTTATACGGCTCGTTG$ -a GAAAGGTGTCCAGGCTGAA$ -a CAGCTTGAGGACAGAGCGGG$ -a GTCATGACAGCAGCCGCATC$ -a GAGTCCTTTAAGAACAGTG$ \
-a AGCATCAAGCGAGACGTTC$ -a GGCGACAGCATACCTATGTGG$ -a AGTTCATCCCAACCCAAAC$ -a CTGTCGATGGTGGGTATGGC$ -a CAACAAACCCGGTAAGAGCG$ \
-a CATAGCATACCGAACTCTTCC$ -G ^CACACGTAGCCTACCAGTTTC -G ^GAGGACTAGAATCAAATGG -G ^GTGCGGCTTCTTCAATATG -G ^ACTGAATGCAGCCTTCGTAG \
-G ^AGCGTTGATGGCGAGATAC -G ^CACTACAGGAAGTACCAATGG -G ^GCAGAAGCCGACAGCAAGTA -G ^CAACGAGCCGTATAAGTATTGG -G ^TTCAGCCTGGACACCTTTC \
-G ^CCCGCTCTGTCCTCAAGCTG -G ^GATGCGGCTGCTGTCATGAC -G ^CACTGTTCTTAAAGGACTC -G ^GAACGTCTCGCTTGATGCT -G ^CCACATAGGTATGCTGTCGCC \
-G ^GTTTGGGTTGGGATGAACT -G ^GCCATACCCACCATCGACAG -G ^CGCTCTTACCGGGTTTGTTG -G ^GGAAGAGTTCGGTATGCTATG -A GAAACTGGTAGGCTACGTGTG$ \
-A CCATTTGATTCTAGTCCTC$ -A CATATTGAAGAAGCCGCAC$ -A CTACGAAGGCTGCATTCAGT$ -A GTATCTCGCCATCAACGCT$ -A CCATTGGTACTTCCTGTAGTG$ \
-A TACTTGCTGTCGGCTTCTGC$ -A CCAATACTTATACGGCTCGTTG$ -A GAAAGGTGTCCAGGCTGAA$ -A CAGCTTGAGGACAGAGCGGG$ -A GTCATGACAGCAGCCGCATC$ \
-A GAGTCCTTTAAGAACAGTG$ -A AGCATCAAGCGAGACGTTC$ -A GGCGACAGCATACCTATGTGG$ -A AGTTCATCCCAACCCAAAC$ -A CTGTCGATGGTGGGTATGGC$ \
-A CAACAAACCCGGTAAGAGCG$ -A CATAGCATACCGAACTCTTCC$ -e 0.06 -O 6 -o qtatptlf_$read1.fastq -p qtatptlf_$read2.fastq ./paired_$read1.fastq ./paired_$read2.fastq

# QA quality-trimmed/adapter-trimmed/primer-trimmed reads

fastqc qtatptlf_$read1.fastq

fastqc qtatptlf_$read2.fastq

if [[ "$variant" = "WT" ]]; then
	reference="CHIKVFLIC06049_WT"
elif [[ "$variant" = "HIFI" ]]; then
		reference="CHIKVFLIC06049_Y"
	elif [[ "$variant" = "DM HIFI" ]]; then
			reference="CHIKVFLIC06049_DY"
		elif [[ "$variant" = "LOFI" ]]; then
				reference="CHIKVFLIC06049_G"
fi

# Merge paired end reads with Flash
flash -m 100 -M 137 -x 0.019 -o flashmerged_$sample $outputdir/qtatptlf_$read1.fastq $outputdir/qtatptlf_$read2.fastq

fastqc ./flashmerged_$sample.extendedFrags.fastq

# Crop head and tail of merged reads
java -jar ~/Bioinformatics/Trimmomatic0_36/trimmomatic-0.36.jar SE ./flashmerged_$sample.extendedFrags.fastq flash_crop_$sample.fastq CROP:150 HEADCROP:25

# Normalize depth with BBnorm
~/Bioinformatics/bbmap/bbnorm.sh in=$outputdir/flash_crop_$sample.fastq out=flash_crop_norm_$sample.fastq target=2000

mkdir "QC"

mv ./*fastqc* ./QC/

# Align merged reads with BWA MEM
bwa index /Users/kasenriemersma/Bioinformatics/Reference/$reference.fasta

bwa mem -t 4 /Users/kasenriemersma/Bioinformatics/Reference/$reference.fasta $outputdir/flash_crop_norm_$sample.fastq > flash_crop_norm_$sample.bam

# Sort sam by coordinate
java -jar ~/Bioinformatics/picard_2.18.1/picard.jar SortSam I=$outputdir/flash_crop_norm_$sample.bam O=sorted_flash_crop_norm_$sample.bam SORT_ORDER=coordinate

# Use LoFreq* to identify true variants
lofreq call -f /Volumes/Coffey_Lab/Kasen/Reference_CHIKV_Sequences/$reference.fasta -o flash_crop_norm_variants_$sample.vcf ./sorted_flash_crop_norm_$sample.bam

# Use SNPdat to characterize variants
perl /Users/kasenriemersma/snpdat/SNPdat_v1.0.5.pl -i ./flash_crop_norm_variants_$sample.vcf -g /Volumes/Coffey_Lab/Kasen/GTF_files/$reference.txt -f /Volumes/Coffey_Lab/Kasen/Reference_CHIKV_Sequences/$reference.fasta -o flash_crop_norm_snpdat_$sample.txt

# Generate nucleotide counts per position
samtools index -b ./sorted_flash_crop_norm_$sample.bam

pysamstats -f /Users/kasenriemersma/Bioinformatics/Reference/$reference.fasta -t variation -D 1000000 --format csv --output $outputdir/ntcounts_$sample.csv ./sorted_flash_crop_norm_$sample.bam

# Calculate diversity metrics with R
Rscript /Users/kasenriemersma/Bioinformatics/DiversityMetrics_v1_1.R $outputdir $sample $variant

exit 1