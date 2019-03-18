#!/bin/bash

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

# Merge paired end reads with BBtools
~/Bioinformatics/bbmap/bbmerge.sh in1=$outputdir/paired_$read1.fastq in2=$outputdir/paired_$read2.fastq out=bbmerge_$sample.fastq qtrim=t trimq=35 minq=35 tno=t minoverlap=50 pfilter=1

fastqc ./bbmerge_$sample.fastq

#Length filter and read quality filter
~/Bioinformatics/bbmap/bbduk.sh in=$outputdir/bbmerge_$sample.fastq out=bbmerge_lfrq_$sample.fastq maxlength=150 mbq=40

# Normalize depth with BBnorm
~/Bioinformatics/bbmap/bbnorm.sh in=$outputdir/bbmerge_lfrq_$sample.fastq out=bbmerge_norm_$sample.fastq target=2000

mkdir "QC"

mv ./*fastqc* ./QC/

# Align merged reads with BWA MEM
bwa index /Users/kasenriemersma/Bioinformatics/Reference/$reference.fasta

bwa mem -t 4 -B 15 -O 135 /Users/kasenriemersma/Bioinformatics/Reference/$reference.fasta $outputdir/bbmerge_norm_$sample.fastq > bbmerge_norm_$sample.bam

# Sort sam by coordinate
java -jar ~/Bioinformatics/picard_2.18.1/picard.jar SortSam I=$outputdir/bbmerge_norm_$sample.bam O=sorted_bbmerge_norm_$sample.bam SORT_ORDER=coordinate

# Generate nucleotide counts per position
samtools index -b ./sorted_bbmerge_norm_$sample.bam

pysamstats -f /Users/kasenriemersma/Bioinformatics/Reference/$reference.fasta -t variation -D 1000000 --format csv --output $outputdir/ntcounts_$sample.csv ./sorted_bbmerge_norm_$sample.bam

# Calculate diversity metrics with R
Rscript /Users/kasenriemersma/Bioinformatics/DiversityMetrics_v1_1.R $outputdir $sample $variant

# Use LoFreq* to identify true variants
lofreq call -f /Users/kasenriemersma/Bioinformatics/Reference/$reference.fasta -o bbmerge_norm_variants_$sample.vcf ./sorted_bbmerge_norm_$sample.bam

# Use SNPdat to characterize variants
perl /Users/kasenriemersma/Bioinformatics/snpdat/SNPdat_v1.0.5.pl -i ./bbmerge_norm_variants_$sample.vcf -g /Users/kasenriemersma/Bioinformatics/GTF_files/$reference.txt -f /Users/kasenriemersma/Bioinformatics/Reference/$reference.fasta -o bbmerge_norm_snpdat_$sample.txt

# Use SNPgenie to characterize selection
perl /Users/kasenriemersma/Bioinformatics/SNPGenie-master/snpgenie.pl --vcfformat=2 --snpreport=bbmerge_norm_variants_$sample.vcf --fastafile=/Users/kasenriemersma/Bioinformatics/Reference/$reference.fasta --gtffile=/Users/kasenriemersma/Bioinformatics/GTF_files/$reference.txt

exit 1
