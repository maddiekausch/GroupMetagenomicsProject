#!/bin/bash
#SBATCH --mail-type=END,FAIL --mail-user=mrk143@georgetown.edu
#SBATCH --job-name="project_trim"
#SBATCH --output="%x.o%j"

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=03:00:00
#SBATCH --mem=10G

# Load Trimmomatic module
shopt -s expand_aliases
module load trimmomatic

# Define paths and variables
R1=/home/mrk143/fastq_raw/SRR6996009/SRR6996009.sra_1.fastq.gz
R2=/home/mrk143/fastq_raw/SRR6996009/SRR6996009.sra_2.fastq.gz
OUTDIR=/home/mrk143/project_reads_trimmed
mkdir -p $OUTDIR

trimmomatic PE -threads $SLURM_CPUS_PER_TASK \
  $R1 $R2 \
  $OUTDIR/R1_paired.fq.gz  $OUTDIR/R1_unpaired.fq.gz \
  $OUTDIR/R2_paired.fq.gz  $OUTDIR/R2_unpaired.fq.gz \
  ILLUMINACLIP:/home/mrk143/TruSeq3-PE.fa:2:30:10 \
  SLIDINGWINDOW:4:20 MINLEN:50

