
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

 # 3/17/26 GROUP 1
# GOAL: Bog frozen rep A (SAMN08784152), assemble contigs from trimmed reads
# INSTALL MEGAHIT
module load mamba/
Create environment with megahit 
mamba create -y -n megahit-env -c conda-forge -c bioconda megahit

# WRITE SLURM SCRIPT FOR MEGAHIT
slurm script:
#!/bin/bash
#SBATCH --job-name=megahit_sample_name   	# how job appears in the queue
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8                 
#SBATCH --mem=32G                         
#SBATCH --time=03:00:00                   
#SBATCH --output=/home/sb1949/bioinform/classproject/logs/megahit_test1.%j.out      
#SBATCH --error=/home/sb1949/bioinform/classproject/logs/megahit_test1.%j.err       

#note %j = job ID

#==== Load mamba/conda module ====
module load mamba
source $(mamba info --base)/etc/profile.d/conda.sh

#Activate the environment where you had MEGAHIT installed
conda activate megahit-env

#==== Set paths and filenames ====

#Directory where the cleaned reads live
READDIR=/home/sb1949/bioinform/classproject/trimmedreads

#Input read files (paired-end)
READ1=${READDIR}/mrk143_group_project_files_project_reads_trimmed_R1_paired.fq.gz
READ2=${READDIR}/mrk143_group_project_files_project_reads_trimmed_R2_paired.fq.gz

#Output directory (give it a name, it will be created by MEGAHIT)
OUTDIR=/home/sb1949/bioinform/classproject/megahit/mrk143_megahit_out

#==== Run MEGAHIT ====

megahit \
  -1 ${READ1} \
  -2 ${READ2} \
  -t ${SLURM_CPUS_PER_TASK} \
  -o ${OUTDIR}

echo "Done. Contigs should be in: ${OUTDIR}/final.contigs.fa"


# CHECK OUTPUT
Navigate to output folder
cd /home/mrk143/project/megahit/bog_frozen_A_megahit_out
Should see: final.contigs.fa, log, and intermediate files
ls
Count how many contigs assembled (lines starting with >)
grep -c ">" final.contigs.fa
Peek at first few contigs
head -20 final.contigs.fa

# INSTALL SEQKIT (run on login node, only need to do once)
module load mamba/
mamba activate megahit-env
mamba install -c bioconda seqkit

# RUN SEQKIT STATS ON ASSEMBLY
seqkit stats -a final.contigs.fa

#  UPLOAD ASSEMBLY TO CLASS BUCKET
gsutil cp megahit/bog_frozen_A_megahit_out/final.contigs.fa gs://gu-biology-dept-class/group3/megahit/


# RUN FASTQC ON TRIMMED FILES -> Team 2
#mkdir -p fastqc_cleaned
#fastqc -o fastqc_cleaned project_reads_trimmed/R*_paired.fq.gz

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

#==== Load mamba/conda module (students: no need to change) ====
module load mamba
source $(mamba info --base)/etc/profile.d/conda.sh

#Activate the environment where you had MEGAHIT installed
conda activate megahit-env

#==== Set paths and filenames (students: edit this block!) ====

#Directory where the cleaned reads live
READDIR=/home/hbw18/trimmed/trim1 #can input where the cleaned reads are 

#Input read files (paired-end)
READ1=${READDIR}/SRR37587558_R1_paired.fastq.gz
READ2=${READDIR}/SRR37587558_R2_paired.fastq.gz

#Output directory (give it a name, it will be created by MEGAHIT)
OUTDIR=/home/hbw18/bioinfoproject/megahit/hbw18_megahit_out

#==== Run MEGAHIT ====

megahit \
  -1 ${READ1} \
  -2 ${READ2} \
  -t ${SLURM_CPUS_PER_TASK} \
  -o ${OUTDIR}

echo "Done. Contigs should be in: ${OUTDIR}/final.contigs.fa"

nano megahit.slurm

# COUNTING THE SEQUENCES OF THE FASTQC -Team 2
#sbatch megahit1.sbatch
#ls # this was the output so the slurm was succesful: hbw18_megahit_out  megahit1.sbatch
#cd hbw18_megahit_out/
#ls
#grep -c "^>" final.contigs.fa #this is the number of sequences in the file 


# CHECKING THE QUALITY
#mamba activate megahit-env
#mamba install -c bioconda seqkit
#install:\seqkit stats -a final.contigs.fa

# Results from Seqkit

file              format  type  num_seqs     sum_len  min_len  avg_len  max_len   Q1   Q2   Q3  sum_gap  N50  N50_num  
final.contigs.fa  FASTA   DNA     73,883  41,566,611      200    562.6   97,639  335  387  488        0  497    2,449       

Q20(%)  Q30(%)  AvgQual  GC(%)  sum_n
0       0        0   47.5      0


# SLURM SCRIPT FOR VIRSORTER2

#!/bin/bash
#SBATCH --job-name=virsorter_sample3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=03:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mrk143@georgetown.edu
#SBATCH --output=/home/mrk143/group_project_files/logs/virsorter.%j.out
#SBATCH --error=/home/mrk143/group_project_files/logs/virsorter.%j.err

#==== Load mamba (students: no need to change) ====
module load mamba
source $(mamba info --base)/etc/profile.d/conda.sh

#Activate the environment where you had VirSorter2 installed
mamba activate vs2-env

#==== Set paths and filenames ====
#set up directories
INDIR=/home/mrk143/group_project_files/megahit/mrk143_megahit_out        #directory where input will come from
OUTROOT=/home/mrk143/group_project_files/virsorter       #directory output will go

SAMPLE_ID=sample3                                #just the basic sample name
INPUT="${INDIR}/final.contigs.fa"                #contig file name/location
OUTDIR="${OUTROOT}/vs2-${SAMPLE_ID}"
mkdir -p "${OUTDIR}"

#==== Run virsorter2 with >5kb cutoff and DNA virus categories first
echo "Running VirSorter2 on ${INPUT}"
virsorter run \
  -w "${OUTDIR}" \
  -i "${INPUT}" \
  --keep-original-seq \
  --include-groups dsDNAphage,NCLDV,ssDNA \
  --min-length 5000

echo "Done."

# Find your results.
(vs2-env) [hbw18@m12-controller ~]$ cd bioinfoproject
(vs2-env) [hbw18@m12-controller bioinfoproject]$ cd virosorter
(vs2-env) [hbw18@m12-controller virosorter]$ cd vs2-
(vs2-env) [hbw18@m12-controller vs2-]$ ls
config.yaml               final-viral-combined.fa  iter-0
final-viral-boundary.tsv  final-viral-score.tsv
(vs2-env) [hbw18@m12-controller vs2-]$ grep -c "^>" final-viral-combined.fa
41
(vs2-env) [hbw18@m12-controller vs2-]$ module load mamba/
(vs2-env) [hbw18@m12-controller vs2-]$ mamba activate megahit-env
(megahit-env) [hbw18@m12-controller vs2-]$ seqkit seq -m 5000 final-viral-combined.fa | grep -c “>”
(megahit-env) [hbw18@m12-controller vs2-]$ seqkit seq -m 5000 final-viral-combined.fa > final-viral-combined_min5kb.fa
[WARN] you may switch on flag -g/--remove-gaps to remove spaces
(megahit-env) [hbw18@m12-controller vs2-]$ grep -c "^>" final-viral-combined_min5kb.fa
41

# Install vclust
(megahit-env) [hbw18@m12-controller vs2-]$ module load mamba
(megahit-env) [hbw18@m12-controller vs2-]$ mamba create -n votu-env -c bioconda -c conda-forge vclust
(megahit-env) [hbw18@m12-controller vs2-]$ mamba activate votu-env

# Prefilter similar genome sequence pairs before conducting pairwise alignments.
(votu-env) [hbw18@m12-controller vs2-]$ vclust prefilter -i final-viral-combined_min5kb.fa -o fltr.txt

# Align similar genome sequence pairs and calculate pairwise ANI measures.
(votu-env) [hbw18@m12-controller vs2-]$ vclust align -i final-viral-combined_min5kb.fa -o ani.tsv --filter fltr.txt

# Cluster genome sequences based on given ANI measure and minimum threshold (these files were generated in the previous steps)
(votu-env) [hbw18@m12-controller vs2-]$ vclust cluster -i ani.tsv -o clusters.tsv --ids ani.ids.tsv --metric ani --ani 0.95 --out-repr

# Make a list of the vOTU headers
(votu-env) [hbw18@m12-controller vs2-]$ awk '{print $2}' clusters.tsv | sort -u > votu_seeds.txt

# Put these vOTU “seed” sequences into a new file and deactivate mamba.
(votu-env) [hbw18@m12-controller vs2-]$ mamba deactivate
(megahit-env) [hbw18@m12-controller vs2-]$ seqkit grep -f votu_seeds.txt final-viral-combined_min5kb.fa > votus_final.fna

# Check Numbers
(megahit-env) [hbw18@m12-controller vs2-]$ wc -l votu_seeds.txt
42 votu_seeds.txt
(megahit-env) [hbw18@m12-controller vs2-]$ grep -c ">" votus_final.fna
41

# Edit the text to take out header (which was messing up the count)
(megahit-env) [hbw18@m12-controller vs2-]$ nano votu_seeds.txt
(megahit-env) [hbw18@m12-controller vs2-]$ wc -l votu_seeds.txt
41 votu_seeds.txt

# 3/24/26
# VirSorter2
bash#!/bin/bash
#SBATCH --job-name=virsorter_sample3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=03:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mrk143@georgetown.edu
#SBATCH --output=/home/mrk143/group_project_files/logs/virsorter.%j.out
#SBATCH --error=/home/mrk143/group_project_files/logs/virsorter.%j.err

module load mamba
source $(mamba info --base)/etc/profile.d/conda.sh
mamba activate vs2-env

INDIR=/home/mrk143/group_project_files/megahit/mrk143_megahit_out
OUTROOT=/home/mrk143/group_project_files/virsorter
SAMPLE_ID=sample3
INPUT="${INDIR}/final.contigs.fa"
OUTDIR="${OUTROOT}/vs2-${SAMPLE_ID}"
mkdir -p "${OUTDIR}"

echo "Running VirSorter2 on ${INPUT}"
virsorter run \
  -w "${OUTDIR}" \
  -i "${INPUT}" \
  --keep-original-seq \
  --include-groups dsDNAphage,NCLDV,ssDNA \
  --min-length 5000

echo "Done."
VirSorter2 Output
config.yaml               final-viral-combined.fa
final-viral-boundary.tsv  final-viral-score.tsv
iter-0/
bash# Count viral contigs
grep -c "^>" final-viral-combined.fa
# 41

# Filter for >5kb and check count
seqkit seq -m 5000 final-viral-combined.fa > final-viral-combined_min5kb.fa
grep -c "^>" final-viral-combined_min5kb.fa
# 41

# Step 7: vclust — Cluster vOTUs
bash# Install vclust (only once)
module load mamba
mamba create -n votu-env -c bioconda -c conda-forge vclust
mamba activate votu-env
bash# 1. Prefilter similar genome sequence pairs before pairwise alignments
vclust prefilter -i final-viral-combined_min5kb.fa -o fltr.txt

# 2. Align and calculate pairwise ANI
vclust align -i final-viral-combined_min5kb.fa -o ani.tsv --filter fltr.txt

# 3. Cluster at 95% ANI
vclust cluster -i ani.tsv -o clusters.tsv --ids ani.ids.tsv \
  --metric ani --ani 0.95 --out-repr

# 4. Make a list of vOTU representative sequence headers
awk '{print $2}' clusters.tsv | sort -u > votu_seeds.txt

# 5. Extract representative sequences
mamba deactivate
mamba activate megahit-env
seqkit grep -f votu_seeds.txt final-viral-combined_min5kb.fa > votus_final.fna
Results:
bashwc -l votu_seeds.txt

# 42 votu_seeds.txt  ← had a stray header line

nano votu_seeds.txt  # removed header

wc -l votu_seeds.txt
# 41 votu_seeds.txt

grep -c ">" votus_final.fna
# 41

# Step 8: CheckV — Quality Metrics
Purpose: Evaluate completeness, contamination, and trim host regions from proviruses.
Input: vOTUs.fna (after VirSorter2 + vclust)
What CheckV does:
Identifies and trims host regions (ORF annotation, gene/GC content)
Estimates completeness vs. a database of complete viral genomes
Identifies closed genomes (terminal repeats)
Assigns quality tiers: Complete / High / Medium / Low / Undetermined

Output columns: completeness (%), quality tier, estimated genome length, host contamination estimate, cleaned sequences.
Why trim host segments for AMG studies? Host-derived genes are annotated as metabolic genes more often than viral ones. Even minor contamination inflates AMG detection. Only regions confidently identified as viral should be used downstream, and candidate AMGs must be flanked by viral-like genes and not at scaffold edges.

# Setup
bashmkdir checkv && cd checkv
module load checkv
checkv download_database ./

# Slurm Script
bash#!/bin/bash
#SBATCH --job-name=checkv
#SBATCH --output=/home/NETID/bioinform/classproject/logs/checkv-%j.out
#SBATCH --error=/home/NETID/bioinform/classproject/logs/checkv-%j.err
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=NETID@georgetown.edu

module load checkv

CHECKVDB="/home/NETID/iceland/virusid/checkv/checkv-db-v1.5"
SAMPLE_ID="vOTUs"
INPUT="/home/NETID/bioinform/classproject/votus/votus_10kb_6samples.fna"
OUTDIR="/home/NETID/bioinform/classproject/checkv/${SAMPLE_ID}"
mkdir -p "${OUTDIR}"

echo "Running CheckV on ${INPUT}"
checkv end_to_end "${INPUT}" "${OUTDIR}" -d "${CHECKVDB}" -t ${SLURM_CPUS_PER_TASK}
echo "Done."
Find Results
bash# Results are in:
quality_summary_votus.tsv
Note how many contigs are Complete / High / Medium / Low quality.
Get Pooled Class vOTUs
bashgcloud storage cp gs://gu-biology-dept-class/ClassProject/votus_10kb_6samples.fna [your_destination]

Why use pooled class vOTUs? Pooling all samples before mapping ensures reads from every sample are aligned to the same reference, making cross-sample abundance comparisons valid. Rare viruses present in only one sample can still be detected across others.


# Step 9: Bowtie2 — Read Alignment
Purpose: Map trimmed reads back to the vOTU reference to measure abundance (coverage).
Setup
bashmkdir bowtie2 && cd bowtie2

# Copy vOTU reference file into this directory
cp [path/to/votus_10kb_6samples.fna] .

module load bowtie2

# Build the index
bowtie2-build votus_10kb_6samples.fna votu_index
Slurm Script
bash#!/bin/bash
#SBATCH --job-name=bowtie2_vOTUs
#SBATCH --output=/home/NETID/bioinform/classproject/logs/bowtie-%j.out
#SBATCH --error=/home/NETID/bioinform/classproject/logs/bowtie-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=NETID@georgetown.edu
#SBATCH --time=8:00:00
#SBATCH --mem=16G

# --------- SET UP ----------
SAMPLE="sampleN_netid"   # EDIT: your sample number + netid
INDEX="/home/NETID/bioinform/classproject/bowtie2/index/votu_index"
OUTPUTDIR="/home/NETID/bioinform/classproject/bowtie2/${SAMPLE}"
R1="/home/NETID/.../trimmedreads/NETID_cleaned_reads_R1_paired.fq.gz"   # EDIT
R2="/home/NETID/.../trimmedreads/NETID_cleaned_reads_R2_paired.fq.gz"   # EDIT

# --------- LOAD MODULES ----------
module purge
module load bowtie2/2.5.4

# --------- RUN ----------
mkdir -p "${OUTPUTDIR}"
cd "${OUTPUTDIR}"
mkdir -p logs

echo "Running bowtie2 on sample ${SAMPLE}"
bowtie2 -p 8 -x "${INDEX}" -1 "${R1}" -2 "${R2}" \
  | samtools view -bS - > "${SAMPLE}.bam"
echo "Finished alignment"

# --------- SORT & INDEX ----------
echo "Sorting..."
samtools sort "${SAMPLE}.bam" > "${SAMPLE}_sorted.bam"
echo "Indexing..."
samtools index "${SAMPLE}_sorted.bam"
echo "Finished ${SAMPLE}"
Output files:
FileDescriptionSAMPLE.bamCompressed binary alignmentSAMPLE_sorted.bamSorted BAM (required for indexing)SAMPLE_sorted.bam.baiBAM index (enables fast lookup)

# Upload to Class Bucket
gcloud storage cp SAMPLE_sorted.bam gs://gu-biology-dept-class/ClassProject/bam
gcloud storage cp SAMPLE_sorted.bam.bai gs://gu-biology-dept-class/ClassProject/bam

