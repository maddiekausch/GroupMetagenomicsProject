
## Master Notes for project
3/12/2026
GOAL: Bog frozen rep A (SAMN08784152), quality check and trim reads
# Download fastq files...
# Go to home directory
cd ~

# Load anaconda
module load anaconda3

# Create an environment with SRA tools (only need to do this once)
conda create -n sra_env -c bioconda sra-tools
# Type "yes" when prompted

# Activate the environment
conda activate sra_env

# Make a folder for  raw data and go into it
mkdir fastq_raw
cd fastq_raw

# Download sample
prefetch SAMN08784152

# Go into the downloaded folder
cd SRR6996009

# Convert to fastq format, then compress
fasterq-dump *.sra
gzip *.fastq

# ORGANIZE FILES
# Go back to your main project folder
cd ~

# Make a folder for cleaned reads (for later)
mkdir cleaned_reads

# Upload raw data to the class bucket
gsutil cp raw/*.fastq.gz gs://gu-biology-dept-class/group3/raw/

# RUNFASTQ ON RAW FILES
# Enter interactive mode on a compute node
srun --pty bash

# Load FastQC
module load fastqc

# Make an output folder
mkdir -p fastqc_out

# Run FastQC on your files
fastqc -o fastqc_out raw/*.fastq.gz

# Check outputs appeared
ls fastqc_out

# RUN TRIMMOMATIC
#!/bin/bash
#SBATCH --mail-type=END,FAIL --mail-user=mrk143@georgetown.edu # send this to the email when the slurm is complete 
#SBATCH --job-name= "project_trim"
#SBATCH --output="%x.o%j"

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=03:00:00
#SBATCH --mem=10G

#Load Trimmomatic module
shopt -s expand_aliases
module load trimmomatic

#Define paths and variables
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

  
# RUN FASTQC ON TRIMMED FILES -> Team 2
#mkdir -p fastqc_cleaned
#fastqc -o fastqc_cleaned project_reads_trimmed/R*_paired.fq.gz
```
Compare these reports to the raw ones — quality should be better.

# INSTALL MEGAHIT
#!/bin/bash
#SBATCH --job-name=megahit_SRR37587558   	# how job appears in the queue
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8                 
#SBATCH --mem=32G                         
#SBATCH --time=03:00:00                   
#SBATCH --output=/home/hbw18/bioinfoproject/logs/megahit_test1.%j.out      
#SBATCH --error=/home/hbw18/bioinfoproject/logs/megahit_test1.%j.err       

#note %j = job ID

# ==== Load mamba/conda module (students: no need to change) ====
module load mamba
source $(mamba info --base)/etc/profile.d/conda.sh

# Activate the environment where you had MEGAHIT installed
conda activate megahit-env

# ==== Set paths and filenames (students: edit this block!) ====

# Directory where the cleaned reads live
READDIR=/home/hbw18/trimmed/trim1 #can input where the cleaned reads are 

# Input read files (paired-end)
READ1=${READDIR}/SRR37587558_R1_paired.fastq.gz
READ2=${READDIR}/SRR37587558_R2_paired.fastq.gz

# Output directory (give it a name, it will be created by MEGAHIT)
OUTDIR=/home/hbw18/bioinfoproject/megahit/hbw18_megahit_out

# ==== Run MEGAHIT ====

megahit \
  -1 ${READ1} \
  -2 ${READ2} \
  -t ${SLURM_CPUS_PER_TASK} \
  -o ${OUTDIR}

echo "Done. Contigs should be in: ${OUTDIR}/final.contigs.fa"

nano megahit.slurm

#COUNTING THE SEQUENCES OF THE FASTQC -Team 2
#sbatch megahit1.sbatch
#ls # this was the output so the slurm was succesful: hbw18_megahit_out  megahit1.sbatch
#cd hbw18_megahit_out/
#ls
#grep -c "^>" final.contigs.fa #this is the number of sequences in the file 


#CHECKING THE QUALITY
#mamba activate megahit-env
#mamba install -c bioconda seqkit
#install:\seqkit stats -a final.contigs.fa


