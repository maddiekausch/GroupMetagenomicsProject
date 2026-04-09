
## Master Notes for project
3/12/2026

GOAL: We have the Bog frozen rep A (SAMN08784152), and we want to quality check and trim the reads. The product will be used for future analysis. 

This section covers downloading raw sequencing reads from NCBI's SRA (Sequence Read Archive), converting them to FASTQ format, running quality control with FastQC, and trimming low-quality bases with Trimmomatic. These are standard first steps in any metagenomics or genomics workflow. Raw reads from the sequencer often contain adapter sequences, low-quality tails, and other artifacts that need to be removed before assembly. 

# Download fastq files...
We use NCBI's SRA Toolkit (prefetch and fasterq-dump) to download and convert reads. prefetch downloads the compresssed .sra file, and fasterq-dump converts it to the paired-end FASTQ files we need. We gzip compress them afterward to save disk space. 

```bash
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

```

# RUNFASTQ ON RAW FILES
FastQC generates an HTML quality report for each FASTQ file. Key things we're checking: per-base quality scores (Phred > or equal to 20 across most of the read), adapter content, and sequence length distribution. These metrics tell us what Trimmmomatic settings to use. 

```bash

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

```

# RUN TRIMMOMATIC
Trimmomatic trims adapter sequences and low-quality bases from paired-end reads. Key parameters we used:
- PE: paired end mode to proccess R1 and R2 together so pairs stay in sync
- Illuminaclip: removes TruSeq3 PE adapter sequences (2 allowed mismatches, palindrome clip threshold 30, simple clip threshold 10)
- sliding window: 4:20: slides a 4-base window across the read and clips when average quality drops below Phred 20
- Minlen:50: discards reads shorter than 50 bp after trimming. Too short reads don't assemble well
Output is four files: paired R1, unpaired R1, paired R2, unpaired R2. We'll only use the paired files for assembly.

```bash
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
```

We successfully trimmed the data, and we can now use the trimmed data for assembly. 

# 3/17/26 
# GOAL: Bog frozen rep A (SAMN08784152), assemble contigs from trimmed reads
Assembly takes the trimmed short reads and reconstructs longer genomic sequences (contigs) by finding overlapping k-mer patterns. We're using MEGAHIT which is well suited for metagenomes because it uses a succinct de Bruijn graph and is memory efficient. This is important when sample contains DNA from many different organisms at varying abundances. 

```bash 
# INSTALL MEGAHIT
module load mamba/
# Create environment with megahit 
mamba create -y -n megahit-env -c conda-forge -c bioconda megahit
```
# WRITE SLURM SCRIPT FOR MEGAHIT
We request 8 CPUs and 32 GB RAM. MEGAHIT can be memory intensive for large metagenomes, but 32G should be sufficient for this sample size. The job runs for up to 3 hours, which is conservative for this data volume. 

slurm script:
```bash
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
```

# CHECK OUTPUT
After assembly, we check that the output file exists and get basic stats. grep -c ">" counts FASTA headers = number of contigs. SeqKit gives us more detailed stats including N50, which is a standard assembly quality metric: the length at which 50% of all assembled bases are contained in contigs of that size or loneger. Higher is better. 

```bash
#Navigate to output folder
cd /home/mrk143/project/megahit/bog_frozen_A_megahit_out
#Should see: final.contigs.fa, log, and intermediate files
ls
#Count how many contigs assembled (lines starting with >)
grep -c ">" final.contigs.fa
#Peek at first few contigs
head -20 final.contigs.fa

# INSTALL SEQKIT (run on login node, only need to do once)
module load mamba/
mamba activate megahit-env
mamba install -c bioconda seqkit

# RUN SEQKIT STATS ON ASSEMBLY
# SeqKit is a fast toolkit for FASTA/FASTQ manipulation. the -a flag gives all statistics including N50, min/max length, and GC content. 

seqkit stats -a final.contigs.fa
```
# Results from the SeqKit -> we found a lot of contigs from this sample, the low N50 shows that we had shorter contigs
Number of Contigs: 73,883
Total Length: Around 41.5 MB
Average Length: 562 BP 
N50: 497

The low N50 (497bp) tells us most contigs are quite short. This is expected for environmental metagenomes, where reads from many organisms at varying abundances are assembled together. Short, fragmented contigs are common when coverage is uneven for the Viral identification step, we will filter to contigs greater than or equal to 5,000 bp since VirSorter2 requires longer sequences to reliably identify viral signals. 

```bash

#  UPLOAD ASSEMBLY TO CLASS BUCKET
gsutil cp megahit/bog_frozen_A_megahit_out/final.contigs.fa gs://gu-biology-dept-class/group3/megahit/
```
# 3/19/26 
# GOAL: Identify viral sequences and group them into different viral populations, want to run Virosorter 
VirSorter2 is a machine-learning-based tool that classifies contigs as viral or non-viral. It uses hallmark viral genes, sequence features, and gene content to score each contig. We then cluster the viral sequences into vOTUs (viral Operational Taxonomic Units) using vclust. We use 95% ANI (Average Nucleotide Identity) as the species-level cutoff for viruses. 
```bash
# Load mamba 
module load mamba
source $(mamba info --base)/etc/profile.d/conda.sh

# Activate the environment where you had VirSorter2 installed
mamba activate vs2-env
rm -rf db
visorter setup -d db -j 4
```
# Create Slurm Script for running virsorter 
Key parameters: --include-groups dsDNAphage, NCLDV, ssDNA are for targeting double-stranded DNA phages, nucleocytoplasmic large DNA viruses, and single-stranded DNA viruses, which are the most relevant groups for this bog metagenome. --min-length 5000 filters out short contigs that don't have enough gene content to reliably classify. --keep-original-seq preserves full contig sequences rather than trimming predicted host regions -- we keep this on because CheckV will do host trimming later. 
```bash
nano virosorter.batch

# slurm script

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
```
```bash
# Find your results.
(vs2-env) [mrk143@m12-controller ~]$ cd group_project_files
(vs2-env) [mrk143@m12-controller group_project_files]$ cd virsorter
(vs2-env) [mrk143@m12-controller virsorter]$ cd vs2-sample3_mrk143
(vs2-env) [mrk143@m12-controller vs2-sample3_mrk143]$ ls
config.yaml               final-viral-combined.fa  iter-0
final-viral-boundary.tsv  final-viral-score.tsv
(vs2-env) [mrk143@m12-controller vs2-sample3_mrk143]$ grep -c "^>" final-viral-combined.fa
41
(vs2-env) [mrk143@m12-controller vs2-sample3_mrk143]$ module load mamba/
(vs2-env) [mrk143@m12-controller vs2-sample3_mrk143]$ mamba activate megahit-env
(megahit-env) [mrk143@m12-controller vs2-sample3_mrk143]$ seqkit seq -m 5000 final-viral-combined.fa | grep -c “>”
(megahit-env) [mrk143@m12-controller vs2-sample3_mrk143]$ seqkit seq -m 5000 final-viral-combined.fa > final-viral-combined_min5kb.fa
[WARN] you may switch on flag -g/--remove-gaps to remove spaces
(megahit-env) [mrk143@m12-controller vs2-sample3_mrk143]$ grep -c "^>" final-viral-combined_min5kb.fa
41
```
# Cluster into vOTUs with vclust
vclust clusters viral genomes by Average Nucleotide Identity (ANI). We use 95% ANI as the species level cutoff (this is standard for viral taxonomy). Workflow--> prefilter (fast k-mer screening to find candidate pairs) --> align (compute pairwise ANI) --> cluster. Output seeds are the representative sequences for each cluster. i.e., one sequence per vOTU. 

```bash
# Install vclust
(megahit-env) [mrk143@m12-controller vs2-sample3_mrk143]$ module load mamba
(megahit-env) [mrk143@m12-controller vs2-sample3_mrk143]$ mamba create -n votu-env -c bioconda -c conda-forge vclust
(megahit-env) [mrk143@m12-controller vs2-sample3_mrk143]$ mamba activate votu-env

# Prefilter similar genome sequence pairs before conducting pairwise alignments.
(votu-env) [mrk143@m12-controller vs2-sample3_mrk143]$ vclust prefilter -i final-viral-combined_min5kb.fa -o fltr.txt

# Align similar genome sequence pairs and calculate pairwise ANI measures.
(votu-env) [mrk143@m12-controller vs2-sample3_mrk143]$ vclust align -i final-viral-combined_min5kb.fa -o ani.tsv --filter fltr.txt

# Cluster genome sequences based on given ANI measure and minimum threshold (these files were generated in the previous steps)
(votu-env) [mrk143@m12-controller vs2-sample3_mrk143]$ vclust cluster -i ani.tsv -o clusters.tsv --ids ani.ids.tsv --metric ani --ani 0.95 --out-repr

# Make a list of the vOTU headers
(votu-env) [mrk143@m12-controller vs2-sample3_mrk143]$ awk '{print $2}' clusters.tsv | sort -u > votu_seeds.txt

# Put these vOTU “seed” sequences into a new file and deactivate mamba.
(votu-env) [mrk143@m12-controller vs2-sample3_mrk143]$ mamba deactivate
(megahit-env) [mrk143@m12-controller vs2-sample3_mrk143]$ seqkit grep -f votu_seeds.txt final-viral-combined_min5kb.fa > votus_final.fna

# Check Numbers
(megahit-env) [mrk143@m12-controller vs2-sample3_mrk143]$ wc -l votu_seeds.txt
42 votu_seeds.txt
(megahit-env) [mrk143@m12-controller vs2-sample3_mrk143]$ grep -c ">" votus_final.fna
41

# Edit the text to take out header (which was messing up the count)
(megahit-env) [mrk143@m12-controller vs2-sample3_mrk143]$ nano votu_seeds.txt
(megahit-env) [mrk143@m12-controller vs2-sample3_mrk143]$ wc -l votu_seeds.txt
41 votu_seeds.txt
```
Note-- wc -l votu_seeds.txt initially returned 42 because the file had a header line. We edited the file with nano to remove the header. 

# 3/24/26
# Goal: Evaluate genome quality and measure abundance 
Viral contigs in environmental metagenomes are often incomplete, chimeric, or contain host DNA (proviruses or phages integrated into bacterial chromosomes). CheckV addressed this by estimating completeness by comparing each contig to a database of complete viral genomes using amino-acid identity and alignment coverage. It detected and trimmed host regions from proviruses (using gene content, GC content, and other signals). It also identified closed/complete genomes by looking for terminal repeats and assigned quality tiers: complete, high quality (over 90% complete), medium quality (50-90%), low quality (less than 50%) or undetermined.
Input: our vOTU representative sequences from the clustering step. 

# Checkv
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
```bash
module load checkv						#its available as a module on the HPC
checkv download_database ./				#make sure you’re in your checkv folder!
```
# Slurm script for checkv
```bash
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

# Load checkv program module 

module load checkv

# Set variables, paths, and filenames 

CHECKVDB="/home/mrk143/group_project_files/checkv/checkv-db-v1.5"

SAMPLE_ID="vOTUs"
INPUT="/home/mrk143/group_project_files/votus/votus_final.fna"
OUTDIR="/home/mrk143/group_project_files/checkv/${SAMPLE_ID}"

mkdir -p "${OUTDIR}"

#  run checkv
echo "Running CheckV on ${INPUT}"
checkv end_to_end "${INPUT}" "${OUTDIR}" -d "${CHECKVDB}" -t ${SLURM_CPUS_PER_TASK}
echo "Done."

```

# Our results.
32 low quality reads
8 non determined reads
1 complete read

Having only 1 complete genome out of 41 vOTUs is expected because the environmental samples are sequenced at mixed, uneven coverage and most viral genomes are only partially recovered. Low-quality reads can still be used for taxonomy and abundance analysis. The 1 complete genome is particularly valuable. 

# Grab votus and put them in class bucket
```bash
$ gcloud storage cp gs://gu-biology-dept-class/ClassProject/votus_10kb_6samples.fna [location]
```
These consist of the class’s virsorter contigs, filtered for >10kb, and clustered.

The class combined vOTUs from all 6 samples and filtered for greater than 10 kb. Now map our trimmed reads back to this combined reference to measure the abundance of each vOTU in our specific sample. Tells us which viruses are actually present and at what relative abundance. Bowtie2 does the read alignment, samtools converts, sorts, and indexes the output for downstream analysis. 

#  Running bowtie2 + samtools
Goal: Map trimmed reads back to the vOTU reference to measure abundance and coverage.

# Build your index with bowtie2 
Before mapping reads, Bowtie2 requires a pre-built index of the reference sequences. Allows fast lookups instead of scanning every sequence. Only needs to be done once per reference file. 

Load bowtie2. 
```bash
$ srun --pty bash
$ module load bowtie2
$ bowtie2-build votus_10kb_6samples.fna votu_index
$ exit
```
# slurm script for bowtie2
We pipe Bowtie2 output directly into samtools to convert from SAM to BAM (compressed binary) format. This avoids writing a large intermediate SAM file to disk. We then sort and index the BAM file, which is required for any downstream coverage or abundance calculations. 
```bash
#!/bin/bash
#SBATCH --job name=bowtie2_vOTUs                 
#SBATCH --output=/home/mrk143/group_project_files/logs/bowtie-%j.out
#SBATCH --error=/home/mrk143/group_project_files/logs/bowtie-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mrk143@georgetown.edu
#SBATCH --time=8:00:00                                
#SBATCH --mem=16G                        

# Set Up
SAMPLE=”sample3_mrk143”  	#whatever your sample number is!
INDEX="home/mrk143/group_project_files/bowtie2/votu_index"
OUTPUTDIR="home/mrk143/group_project_files/bowtie2/${SAMPLE}"

# Load Modules
module purge
module load bowtie2/2.5.4

# Run BOWTIE2 and Pipe to Samtools

# First make output and log directories; move into OUTPUTDIR
mkdir -p "${OUTPUTDIR}"
cd "${OUTPUTDIR}"
mkdir -p logs

echo "Running bowtie2 on sample ${SAMPLE}"

bowtie2 -p 8 -x "${INDEX}" -1 "/home/mrk143/group_project_files/project_reads_trimmed/R1_paired_clean.fq.gz" -2 "/home/mrk143/group_project_files/project_reads_trimmed/R2_paired_clean.fq.gz" \
| samtools view -bS - > "${SAMPLE}.bam"

echo "Finished running bowtie2 and performing compression"

# Sort and index files

echo "Sorting"
samtools sort "${SAMPLE}.bam" > "${SAMPLE}_sorted.bam"

echo "Indexing"
samtools index "${SAMPLE}_sorted.bam"

echo "Finished ${SAMPLE}"
```
# Upload all  bowtie output files to bucket
```bash
$ gcloud storage cp home/mrk143/group_project_files/bowtie2/sample3_mrk143/sample3_mrk143.bam gs://gu-biology-dept-class/ClassProject/bam
```
# 4/7/26
Goal: Create ecological data visualizations to better understand the diversity of our samples. These graphs were created using R and R studio. 

# Bar graphs for Richness and Shannon Diversity
The votu richness graph helps us compare how many unique vOTUs present are in each sample. The shannon diversity measures both richness + evenness in each sample, which shows how evenly taxa are distributed within each sample.

```bash
# Install packages. These are required R packages for diversity calculations and plotting.
install.packages(c("vegan", "pheatmap", "RColorBrewer"))

# Load libraries in order that you can use the packages you just downloaded for data visualization. 
library(vegan)
library(pheatmap)
library(RColorBrewer)

# Read input file. It also import the vOTU abundance table. check.names preserves the original sample names. Important make sure you have set your working directory to where the csv is located. You can change your working directory by going into the session tab on R and then clicking the folder you want as your directory. 
abund_raw <- read.csv(
  "votus_12367_coverm_TPM_3votusremoved.csv",
  header = TRUE,
  check.names = FALSE)

# Format data into a usable form. This sets the vOTU IDs as row names and remove the non-numeric column.
rownames(abund_raw) <- abund_raw$Contig
abund <- abund_raw[, -1]

# Sanity check. This confirms the dimensions and correct labeling of rows (vOTUs) and columns (samples).
dim(abund)
head(rownames(abund))
colnames(abund)

# Transpose matrix. For the graphs we need samples as rows and taxa as columns.
abund_t <- t(abund)

# Remove empty samples. Exclude samples with zero total abundance to avoid errors.
abund_t <- abund_t[rowSums(abund_t) > 0, , drop = FALSE]


# Calculating Richness. Counts how many taxa are present per sample. 
richness <- specnumber(abund_t)

# Calculating Shannon Diversity, which is richness and evenness. 
shannon <- diversity(abund_t, index = "shannon")

# Combine results into one table.
alpha_div <- data.frame(
  Sample   = rownames(abund_t),
  Richness = richness,
  Shannon  = shannon)
alpha_div

# Plot richness per sample
barplot(
  alpha_div$Richness,
  names.arg = alpha_div$Sample,
  las = 2,
  ylab = "vOTU richness (count)",
  main = "Per-sample viral richness")

# Plot Shannon diversity per sample
barplot(
  alpha_div$Shannon,
  names.arg = alpha_div$Sample,
  las = 2,
  ylab = "Shannon diversity",
  main = "Per-sample viral diversity")

```
<img width="1410" height="1312" alt="Richness Bargraph" src="https://github.com/user-attachments/assets/b41e3062-2dde-4fc9-bf9e-85dcfd981341" />
<img width="1410" height="1312" alt="Shannon Graph" src="https://github.com/user-attachments/assets/1106b5d1-12bd-4b93-b5ec-b4fe4230f8e9" />

Some interpretation of the graphs: 
- Richness varied substantially across samples
- Large richness does not always mean high Shannon Diversity. This could be beacause despite high richness, the community can be dominated by a few taxa. - Samples with medium richness exhibited higher Shannon values (more even distribution)

# Figure results of heatmap
Heatmaps visualize which taxa are present and how abundant they are across samples. 
```bash

# Create variable for data, so you don't have to keep on writing the full excel file. 
filename <- "ClassProject_votus_12367_coverm_TPM.xlsx"  # Excel file to read

# Keep vOTUs with greater TPMs. This is so our data is not hinging on stray outliers and focusing on meaningful patterns.  
tpm_threshold <- 10

# Color gradient for heatmap                                    
heatmap_colors <- c("#440154", "#31688e", "#35b779", "#fde725") 

# Load required libraries
library(readxl)    # for reading Excel files
library(pheatmap)  # for making heatmaps

# Read the Excel file
cov <- read_xlsx(filename)

# Keep only the Contig column and columns with TPM values from the Excel File. 
tpm_cols <- grepl("TPM$", names(cov))  # find all column names ending in 'TPM'
cov_tpm <- cov[ , c("Contig", names(cov)[tpm_cols])]  # subset only relevant columns

# Remove S1k141_26921||full. Some vOTUs may be outliers and can mess up the heatmap
cov_tpm <- subset(cov_tpm, Contig != "S1k141_26921||full")

# Drop very low-abundance vOTUs. This cleans your data more. 
cov_tpm$max_tpm <- apply(cov_tpm[ , -1], 1, max, na.rm = TRUE) 
cov_tpm <- subset(cov_tpm, max_tpm > tpm_threshold)              
cov_tpm$max_tpm <- NULL  

# Stop if no vOTUs passed the filter. This prevents downstream errors and informs user that threshold may be too high
if (nrow(cov_tpm) == 0) {
  stop("No vOTUs passed the TPM threshold. Try lowering tpm_threshold.")}

# Prepare matrix for heatmap. Rows = vOTUs, Cols = samples.
mat <- as.matrix(cov_tpm[ , -1])
rownames(mat) <- cov_tpm$Contig

# Log-transform for color scaling. This makes differences easier to see.
mat_log <- log10(mat + 1)  # +1 to avoid log(0)

# Create the Heat Map. 
pheatmap(
  mat_log,
  cluster_rows = TRUE,  # cluster vOTUs
  cluster_cols = TRUE,  # cluster samples
  scale = "none",       # don't scale rows or columns
  color = colorRampPalette(heatmap_colors)(100),  # gradient colors
  fontsize_row = 8,
  fontsize_col = 10,
  main = "vOTU relative abundance (log10 TPM + 1)")
```
<img width="549" height="857" alt="heatmap bioinformatics (1)" src="https://github.com/user-attachments/assets/4dcd3379-c5e5-4d3c-b820-b991d5ecbd08" />

Some interpretation of the graphs:
- Some vOTUs are in high abundance in specific samples while being absent or rare in others, which indicates strong compositional differences between samples
-  Some samples share pretty similar viral communities, while others are more distinct. A possible explanation could be due to key dominant taxa affecting others prescence




