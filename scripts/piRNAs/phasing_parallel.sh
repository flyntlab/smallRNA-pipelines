#!/bin/bash

#SBATCH --partition=node
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=1000:00:00
#SBATCH --job-name=cram_speed
#SBATCH --output=myoutput-%j.txt

cd "$SLURM_SUBMIT_DIR"

# samtools index ../all.bam

module load gcc
module load R
#module load samtools
#module load bedtools



chmod +x phasingPlotsParallel.R
#split -a 15 -l 7 HE.pingpong.bed

#split -a 10 -l 1 HE.all.bed


#split -a 10 -l 50 HE.all.bed

split -a 10 -l 1 HE.large.bed

#When splitting these files , the second integer dictates how many lines are pasted into each file
#The trade off is higher variance for better processing speed.





# Ensure the parallel tool is available
module load parallel

process_file() {
    local file=$1
  #Make sure to use genome.fasta.fai file

#If using CRAM
samtools view -b -t GWHACBE00000000.genome.fasta -ML "$file" all.cram|
   # samtools view -b -ML "$file" |
bedtools bamtobed |
    awk 'BEGIN{OFS="\t"} { for (i=1;i<=1;++i) { if ($6=="-") { print $1,$2-30,$2+30+1,$4,$5,$6 } else { print $1,$3-30-1,$3+30,$4,$5,$6 }}}' |
    bedtools getfasta -fi GWHACBE00000000.genome.fasta -bed stdin -fo stdout -s -name -tab |
    python piPipes_nuc_percentage.py 20 | cat -n > "phasing_${file}.txt"

    ./phasingPlotsParallel.R "phasing_${file}.txt" "phasZscore_${file}.txt"

    awk '{print $5}' < "phasZscore_${file}.txt" > "${file}.TZ.txt"
    awk '{print $2}' < "phasZscore_${file}.txt" > "${file}.AZ.txt"
}

export -f process_file

# Run the process_file function in parallel
parallel process_file ::: x*

paste -s *TZ.txt > allTphase_parallel.txt
paste -s *AZ.txt > allAphase_parallel.txt

# Optionally remove the split files if no longer needed
# rm x*
