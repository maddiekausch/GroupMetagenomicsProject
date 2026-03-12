
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
trimmomatic PE \
  raw/sample_R1.fastq.gz raw/sample_R2.fastq.gz \
  cleaned_reads/sample_R1_paired.fastq.gz cleaned_reads/sample_R1_unpaired.fastq.gz \
  cleaned_reads/sample_R2_paired.fastq.gz cleaned_reads/sample_R2_unpaired.fastq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
  SLIDINGWINDOW:4:20 \
  MINLEN:50
  
# RUN FASTQC ON TRIMMED FILES
mkdir -p fastqc_trimmed_out
fastqc -o fastqc_trimmed_out cleaned_reads/*_paired.fastq.gz
```
Compare these reports to the raw ones — quality should be better.

