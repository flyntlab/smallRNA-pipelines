#!/bin/bash

#SBATCH --partition=node
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=1000:00:00
#SBATCH --job-name=cram_conv
#SBATCH --output=myoutput-%j.txt

cd "$SLURM_SUBMIT_DIR"


#samtools view -C -h -T GWHACBE00000000.genome.fasta -o all.cram all.bam
samtools index all.cram
