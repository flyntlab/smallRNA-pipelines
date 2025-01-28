#!/bin/bash

#SBATCH --partition=node
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=1000:00:00
#SBATCH --job-name=all_Andre
#SBATCH --output=myoutput-%j.txt

cd "$SLURM_SUBMIT_DIR"

module load bowtie
module load bedtools
module load samtools

#^^^^^^^^^^^^^^^^^^^^^^^^
#Optional Content for SLURM Users

#Concatenate all raw reads into one fastq file
cat ../*.fq > all.raw.fastq

#Build reference index for target organism's genome and record chromosome info
bowtie-build ../*genome.fasta genoBwt1
samtools faidx ../*genome.fasta
awk '{print $1"\t"$2}' < ../*geno.fasta.fai > chrSizes.txt

#Separate all RNAs by size
#Typically piRNAs (24-33nt) are largest and mi/siRNAs are smaller (18-24nt)

awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; if (length(seq) >= 15 && length(seq) <= 32) {print header, seq}}' \
< all.raw.fasq > all.fastq

#rm all.raw.fastq

awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; if (length(seq) >= 21 && length(seq) <= 24) {print header, seq}}' \
< all.fastq > small.fastq &

awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; if (length(seq) >= 21 && length(seq) <= 21) {print header, seq}}' \
< all.fastq > small21.fastq &

awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; if (length(seq) >= 26 && length(seq) <= 30) {print header, seq}}' \
< all.fastq > large.fastq &

#Map all fastq files to bam
bowtie -q -p 20 -a -m50 --no-unal -x genoBwt1 all.fastq -S |
samtools view -@ 20 -q 10 -b |
samtools sort -@ 20 -m 6G > all.bam
samtools index all.bam

bowtie -q -p 20 -a -m50 --no-unal -x genoBwt1 small.fastq -S |
samtools view -@ 20 -q 10 -b |
samtools sort -@ 20 -m 6G > small.bam
samtools index small.bam

bowtie -q -p 20 -a -m50 --no-unal -x genoBwt1 small21.fastq -S |
samtools view -@ 20 -q 10 -b |
samtools sort -@ 20 -m 6G > small21.bam
samtools index small21.bam

bowtie -q -p 20 -a -m50 --no-unal -x genoBwt1 large.fastq -S |
samtools view -@ 20 -q 10 -b |
samtools sort -@ 20 -m 6G > large.bam
samtools index large.bam

python overlapping_reads.py --input small.bam \
--minquery 19 --maxquery 19 --mintarget 19 --maxtarget 19 --overlap 17 --output siRNA2-19.fasta
python overlapping_reads.py --input small.bam \
--minquery 20 --maxquery 20 --mintarget 20 --maxtarget 20 --overlap 18 --output siRNA20.fasta
python overlapping_reads.py --input small.bam \
--minquery 22 --maxquery 22 --mintarget 22 --maxtarget 22 --overlap 20 --output siRNA22.fasta
python overlapping_reads.py --input small.bam \
--minquery 21 --maxquery 21 --mintarget 21 --maxtarget 21 --overlap 19 --output siRNA21.fasta
python overlapping_reads.py --input small.bam \
--minquery 23 --maxquery 23 --mintarget 23 --maxtarget 23 --overlap 21 --output siRNA23.fasta

cat siRNA2* > dicer.fa && rm siRNA2*

bowtie -p 20 -f -a -m100 --no-unal -x genoBwt1 dicer.fa -S |
samtools view -@ 20 -q 10 -b |
samtools sort -@ 20 -m 6G > dicer.bam
samtools index dicer.bam

python overlapping_reads.py --input all.bam \
--minquery 25 --maxquery 30 --mintarget 25 --maxtarget 30 --overlap 10 --output pingpong.fa

bowtie -p 20 -f -a -m100 --no-unal -x genoBwt1 pingpong.fa -S |
samtools view -@ 20 -q 10 -b |
samtools sort -@ 20 -m 6G > pingpong.bam
samtools index pingpong.bam

bedtools genomecov -bg -ibam all.bam  | awk '$4 > 10' > HE.tmp.bedgraph
bedtools merge -d 500 -i HE.tmp.bedgraph > HE.tmp.merge.bed

awk '{n=$2; x=$3; print $1"\t"$2"\t"$3"\t"x-n}' < HE.tmp.merge.bed |
awk '$4 > 40' > HE.all.bed

bedtools multicov -bams small.bam all.bam -bed HE.all.bed |
awk '{n=$5; x=$6; print $1"\t"$2"\t"$3"\t"n/x}' |
awk '$4 > 0.5' > HE.small.bed

bedtools multicov -bams large.bam all.bam -bed HE.all.bed |
awk '{n=$5; x=$6; print $1"\t"$2"\t"$3"\t"n/x}' |
awk '$4 > 0.5' > HE.large.bed

bedtools bamtobed -i all.bam | awk '{print $3-$2 }' | sort | uniq -c > all_read_counts_Whole.txt

bedtools intersect -u -a HE.all.bed -b pingpong.bam > HE.pingpong.bed

bedtools intersect -u -a HE.all.bed -b dicer.bam > HE.dicer.bed
