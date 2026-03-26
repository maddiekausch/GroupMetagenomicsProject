
## Master Notes for project
3/12/2026

GOAL: We have the Bog frozen rep A (SAMN08784152), and we want to quality check and trim the reads. The product will be used for future analysis. 
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

We successfully trimmed the data, and we can now use the trimmed data for assembly. 

# 3/17/26 
# GOAL: Bog frozen rep A (SAMN08784152), assemble contigs from trimmed reads
# INSTALL MEGAHIT
module load mamba/
Create environment with megahit 
mamba create -y -n megahit-env -c conda-forge -c bioconda megahit

# WRITE SLURM SCRIPT FOR MEGAHIT
slurm script:
#!/bin/bash
#SBATCH --job-name= megahit_SRR37587558	# how job appears in the queue
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8                 
#SBATCH --mem=32G                         
#SBATCH --time=03:00:00                   
#SBATCH --output=/home/mrk143/bioinfoproject/logs/megahit_test1.%j.out         
#SBATCH --error=/home/mrk143/bioinfoproject/logs/megahit_test1.%j.err       
     

#note %j = job ID

#==== Load mamba/conda module ====
module load mamba
source $(mamba info --base)/etc/profile.d/conda.sh

#Activate the environment where you had MEGAHIT installed
conda activate megahit-env

#==== Set paths and filenames ====

#Directory where the cleaned reads live
READDIR=/home/mrk143/trimmed/trim1

#Input read files (paired-end)
READ1=${READDIR}/mrk143_group_project_files_project_reads_trimmed_R1_paired.fq.gz
READ2=${READDIR}/mrk143_group_project_files_project_reads_trimmed_R2_paired.fq.gz

#Output directory (give it a name, it will be created by MEGAHIT)
OUTDIR=/home/mrk143/bioinfoproject/megahit/mrk143_megahit_out

# Run MegaHit -> assembles contigs from overlap patterns of the short reads 

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

# Results from the SeqKit -> we found a lot of contigs from this sample, the low N50 shows that we had shorter contigs
Number of Contigs: 73,883
Total Length: Around 41.5 MB
Average Length: 562 BP 
N50: 497 

#  UPLOAD ASSEMBLY TO CLASS BUCKET
gsutil cp megahit/bog_frozen_A_megahit_out/final.contigs.fa gs://gu-biology-dept-class/group3/megahit/

# 3/19/26 
# GOAL: Identify viral sequences and group them into different viral populations, want to run Virosorter 

# Load mamba 
module load mamba
source $(mamba info --base)/etc/profile.d/conda.sh

# Activate the environment where you had VirSorter2 installed
mamba activate vs2-env
rm -rf db
visorter setup -d db -j 4

# Create Slurm Script for running virsorter 
#!/bin/bash
#SBATCH --job-name=virsorter_mrk143
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=03:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mrk143@georgetown.edu
#SBATCH --output=/home/mrk143/group_project_files/logs/virsorter.%j.out
#SBATCH --error=/home/mrk143/group_project_files/logs/virsorter.%j.err

# ==== Load mamba ====
module load mamba
source $(mamba info --base)/etc/profile.d/conda.sh

# Activate the environment where VirSorter2 was installed
mamba activate vs2-env

# ==== Set paths and filenames ====
# directory where MEGAHIT output is
INDIR=/home/mrk143/group_project_files/megahit/mrk143_megahit_out

# directory where VirSorter2 output will go
OUTROOT=/home/mrk143/group_project_files/virsorter
mkdir -p "${OUTROOT}"

# basic sample name
SAMPLE_ID=sample3_mrk143

# contig file from MEGAHIT
INPUT="${INDIR}/final.contigs.fa"

# output folder
OUTDIR="${OUTROOT}/vs2-${SAMPLE_ID}"
mkdir -p "${OUTDIR}"

# ==== Run VirSorter2 ====
echo "Running VirSorter2 on ${INPUT}"

virsorter run \
  -w "${OUTDIR}" \
  -i "${INPUT}" \
  --keep-original-seq \
  --include-groups dsDNAphage,NCLDV,ssDNA \
  --min-length 5000 \
  -j ${SLURM_CPUS_PER_TASK}

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

# Checkv.
- completeness, contamination, and provirus trimming
- Viral contigs can be partial, chimeric, or include host DNA.  CheckV evaluates the quality of each viral contig and removes host segments from proviruses.

Input: viral contigs (single‑contig genomes or fragments), after VirSorter2 and clustering → your vOTUs

Identify and trim host regions (annotates ORFs and uses gene content, GC content, and other features to detect host‑like segments)
Estimate completeness (It compares each contig to a large database of complete viral genomes using amino‑acid identity and alignment coverage - (e.g., 40% complete vs 98% complete).
Identify closed genomes (It looks for terminal repeats)
Assign quality tiers: it classifies each sequence into quality categories

# Set up checkv

You’ll need your file with the vOTUs. (note its full path)
Set up a directory called “checkv”, move into that directory to download the database in the next step.

module load checkv						#its available as a module on the HPC
checkv download_database ./				#make sure you’re in your checkv folder!

# Slurm script for checkv

#!/bin/bash
#SBATCH --job-name=checkv
#SBATCH --output=/home/mrk143/group_project_files/logs/checkv-%j.out
#SBATCH --error=/home/mrk143/group_project_files/logs/checkv-%j.err
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mrk143@georgetown.edu


#==== Load checkv program module (students: no need to change) ====

module load checkv


#==== Set variables, paths, and filenames (students: edit this block!) ====

CHECKVDB="/home/mrk143/group_project_files/checkv/checkv-db-v1.5"

SAMPLE_ID="vOTUs"
INPUT="/home/mrk143/group_project_files/votus/votus_final.fna"
OUTDIR="/home/mrk143/group_project_files/checkv/${SAMPLE_ID}"

mkdir -p "${OUTDIR}"


#==== run checkv (students: no need to change. The second line is the command) ====
echo "Running CheckV on ${INPUT}"
checkv end_to_end "${INPUT}" "${OUTDIR}" -d "${CHECKVDB}" -t ${SLURM_CPUS_PER_TASK}
echo "Done."


# Find your results.
After checkv completes, you’ll see a few files. The information you want to look at is in: 
quality_summary_votus.tsv.  

Have the note-taker check it out and mention how many “complete” “high quality” or “low quality” you have.
32 low quality reads
8 non determined reads
1 complete read

Next - grab the pooled vOTUs from the buckett:  

$ gcloud storage cp gs://gu-biology-dept-class/ClassProject/votus_10kb_6samples.fna [location]
These consist of the class’s virsorter contigs, filtered for >10kb, and clustered.

#  Running bowtie2 + samtools
Purpose: Map trimmed reads back to the vOTU reference to measure abundance (coverage).

Create a bowtie2 directory for yourself.

Copy the votu file (votus_10kb_6samples.fna) into this directory. 

Find your “trimmed reads” and note the full path names for R1 and R2.

# Build your index with bowtie2 

Load bowtie2. 

$ srun --pty bash
$ module load bowtie2
$ bowtie2-build votus_10kb_6samples.fna votu_index
$ exit

# slurm script for bowtie2
#!/bin/bash
#SBATCH --job name=bowtie2_vOTUs                 
#SBATCH --output=/home/mrk143/group_project_files/logs/bowtie-%j.out
#SBATCH --error=/home/mrk143/group_project_files/logs/bowtie-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mrk143@georgetown.edu
#SBATCH --time=8:00:00                                
#SBATCH --mem=16G                        

#---------SET UP----------
SAMPLE=”sample3_mrk143”  	#whatever your sample number is!
INDEX="home/mrk143/group_project_files/bowtie2/votu_index"
OUTPUTDIR="home/mrk143/group_project_files/bowtie2/${SAMPLE}"

#--------- LOAD MODULES ----------
module purge
module load bowtie2/2.5.4

#--------- RUN BOWTIE2 AND PIPE TO SAMTOOLS ----------

#First make output and log directories; move into OUTPUTDIR
mkdir -p "${OUTPUTDIR}"
cd "${OUTPUTDIR}"
mkdir -p logs

echo "Running bowtie2 on sample ${SAMPLE}"

bowtie2 -p 8 -x "${INDEX}" -1 "/home/mrk143/group_project_files/project_reads_trimmed/R1_paired_clean.fq.gz" -2 "/home/mrk143/group_project_files/project_reads_trimmed/R2_paired_clean.fq.gz" \
| samtools view -bS - > "${SAMPLE}.bam"

echo "Finished running bowtie2 and performing compression"

#---------sort and index files
echo "Sorting"
samtools sort "${SAMPLE}.bam" > "${SAMPLE}_sorted.bam"

echo "Indexing"
samtools index "${SAMPLE}_sorted.bam"

echo "Finished ${SAMPLE}"




# Important: Upload all your bowtie output files to bucket!!! 

Put into class project directory and make sure they’re labeled (e.g., “sample2_mrk_sorted.bam”)!!!!

$ gcloud storage cp [file] gs://gu-biology-dept-class/ClassProject/bam
   

