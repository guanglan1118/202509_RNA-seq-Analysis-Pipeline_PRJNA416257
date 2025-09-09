# RNA-seq Analysis Pipeline / 202509_PRJNA416257
## Datasets
- Dataset PRJNA416257
- Nat Commun.2019/PMID: 30655535
- ONECUT2 is a driver of neuroendocrine prostate cancer

## Folder layout
~~~
# bash
project_PRJNA416257/
├─ sra/              # sra
├─ sra_fastq/        # FASTQs
├─ ref/              # reference (STAR) STAR_index_gencodev44
├─ mapping/          # raw_counts.csv
├─ r/                # R scripts
└─ results/          # DE tables, plots, GSEA
~~~
## 0) Get the SRR runs & metadata (SRA → FASTQ)
<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106305>

<https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA416257&o=acc_s%3Aa>

On the PRJNA416257 page, click “Send to” → “Run Selector” → “Run Selector” → "Download Metadata/Accession List". 
Save as **metadata.csv**
Edit a minimal metadata.csv with columns:

![metadata](figures/metadata.png)

## 1) Download FASTQs
### 1.1) intsall fasterq-dump 

*fasterq-dump is a tool from SRA-tools (the NCBI Sequence Read Archive toolkit)*

*fasterq-dump does not support --gzip* 

~~~
# bash
# Create a dedicated environment
conda create -n sra -c bioconda -c conda-forge sra-tools pigz
conda activate sra

# Verify installation
which fasterq-dump
fasterq-dump --version  #fasterq-dump : 3.2.1
~~~

### 1.2) Downloads the sequencing run from NCBI SRA
### 1.2.1) Downloads single cases
*download .sra files* 
~~~
bash
pwd # project_PRJNA416257/raw/
prefetch SRR7179504
~~~
This will produce files like:
- raw/SRR7179504.sra

*fastq-dump to extract .sra files into a FASTQ file*
~~~
bash
pwd # project_PRJNA416257/raw_fastq
mkdir -p raw_fastq
# Extract into raw_fastq/
fastq-dump --split-files -O raw_fastq raw/SRR7179504.sra
~~~
This will produce files like:
- raw_fastq/SRR7179504_1.fastq
- which means it might single-end sequencing

check single or paired end sequence
~~~
# bash
# check if SRA file is complete
vdb-validate raw/SRR7179504.sra

# check basic info of SRA file 
vdb-dump --info raw/SRR7179504.sra
~~~

~~~
# check single or paired end sequence
vdb-dump -f tab -C READ_LEN,READ_TYPE -R 1-3 raw/SRR7179504.sra
~~~
This will produce files like:
- 76, 0   SRA_READ_TYPE_BIOLOGICAL, SRA_READ_TYPE_TECHNICAL
- 76, 0   SRA_READ_TYPE_BIOLOGICAL, SRA_READ_TYPE_TECHNICAL
- 76, 0   SRA_READ_TYPE_BIOLOGICAL, SRA_READ_TYPE_TECHNICAL

### 1.2.2) Downloads batch cases
~~~
# bash
# (1) download .sra files into ./sra directory
mkdir -p sra

for sra_id in SRR7179504 SRR7179505 SRR7179506 SRR7179507 \
              SRR7179508 SRR7179509 SRR7179510 SRR7179511 \
              SRR7179520 SRR7179521 SRR7179522 SRR7179523 \
              SRR7179524 SRR7179525 SRR7179526 SRR7179527 \
              SRR7179536 SRR7179537 SRR7179540 SRR7179541
do
    echo "Currently downloading: $sra_id"
    prefetch $sra_id --output-directory ./sra
done
~~~

~~~
# bash
# (2) convert to fastq.gz
mkdir -p fastq
for sra_id in SRR7179504 SRR7179505 SRR7179506 SRR7179507 \
              SRR7179508 SRR7179509 SRR7179510 SRR7179511 \
              SRR7179520 SRR7179521 SRR7179522 SRR7179523 \
              SRR7179524 SRR7179525 SRR7179526 SRR7179527 \
              SRR7179536 SRR7179537 SRR7179540 SRR7179541
do
    echo "Generating fastq for: $sra_id"
    fastq-dump --outdir fastq --gzip --skip-technical \
               --readids --read-filter pass --dumpbase --split-3 --clip \
                ./$sra_id/$sra_id.sra
done       
~~~
This will produce files like:
pwd #sra/fastq
- SRR7179504_pass.fastq.gz  SRR7179508_pass.fastq.gz  SRR7179520_pass.fastq.gz  SRR7179524_pass.fastq.gz  SRR7179536_pass.fastq.gz
- SRR7179505_pass.fastq.gz  SRR7179509_pass.fastq.gz  SRR7179521_pass.fastq.gz  SRR7179525_pass.fastq.gz  SRR7179537_pass.fastq.gz
- SRR7179506_pass.fastq.gz  SRR7179510_pass.fastq.gz  SRR7179522_pass.fastq.gz  SRR7179526_pass.fastq.gz  SRR7179540_pass.fastq.gz
- SRR7179507_pass.fastq.gz  SRR7179511_pass.fastq.gz  SRR7179523_pass.fastq.gz  SRR7179527_pass.fastq.gz  SRR7179541_pass.fastq.gz

### 1.2.3) Concatenating FASTQ files
~~~
pwd #sra/fastq
cat SRR7179504_pass.fastq.gz SRR7179505_pass.fastq.gz SRR7179506_pass.fastq.gz SRR7179507_pass.fastq.gz  > LNCAP_Normoxia_S1.fastq.gz
cat SRR7179508_pass.fastq.gz SRR7179509_pass.fastq.gz SRR7179510_pass.fastq.gz SRR7179511_pass.fastq.gz  > LNCAP_Normoxia_S2.fastq.gz
cat SRR7179520_pass.fastq.gz SRR7179521_pass.fastq.gz SRR7179522_pass.fastq.gz SRR7179523_pass.fastq.gz  > LNCAP_Hypoxia_S1.fastq.gz
cat SRR7179524_pass.fastq.gz SRR7179525_pass.fastq.gz SRR7179526_pass.fastq.gz SRR7179527_pass.fastq.gz  > LNCAP_Hypoxia_S2.fastq.gz
~~~
This will produce files like:
~~~
LNCAP_Hypoxia_S1.fastq.gz   SRR7179504_pass.fastq.gz  SRR7179508_pass.fastq.gz  SRR7179520_pass.fastq.gz  SRR7179524_pass.fastq.gz  SRR7179536_pass.fastq.gz
LNCAP_Hypoxia_S2.fastq.gz   SRR7179505_pass.fastq.gz  SRR7179509_pass.fastq.gz  SRR7179521_pass.fastq.gz  SRR7179525_pass.fastq.gz  SRR7179537_pass.fastq.gz
LNCAP_Normoxia_S1.fastq.gz  SRR7179506_pass.fastq.gz  SRR7179510_pass.fastq.gz  SRR7179522_pass.fastq.gz  SRR7179526_pass.fastq.gz  SRR7179540_pass.fastq.gz
LNCAP_Normoxia_S2.fastq.gz  SRR7179507_pass.fastq.gz  SRR7179511_pass.fastq.gz  SRR7179523_pass.fastq.gz  SRR7179527_pass.fastq.gz  SRR7179541_pass.fastq.gz
~~~
~~~
mv SRR7179536_pass.fastq.gz PC3_Normoxia_S1.fastq.gz
mv SRR7179537_pass.fastq.gz PC3_Normoxia_S2.fastq.gz
mv SRR7179540_pass.fastq.gz PC3_Hypoxia_S1.fastq.gz
mv SRR7179541_pass.fastq.gz PC3_Hypoxia_S2.fastq.gz
~~~

This will produce files like:
~~~
LNCAP_Hypoxia_S1.fastq.gz   PC3_Hypoxia_S1.fastq.gz   SRR7179504_pass.fastq.gz  SRR7179508_pass.fastq.gz  SRR7179520_pass.fastq.gz  SRR7179524_pass.fastq.gz
LNCAP_Hypoxia_S2.fastq.gz   PC3_Hypoxia_S2.fastq.gz   SRR7179505_pass.fastq.gz  SRR7179509_pass.fastq.gz  SRR7179521_pass.fastq.gz  SRR7179525_pass.fastq.gz
LNCAP_Normoxia_S1.fastq.gz  PC3_Normoxia_S1.fastq.gz  SRR7179506_pass.fastq.gz  SRR7179510_pass.fastq.gz  SRR7179522_pass.fastq.gz  SRR7179526_pass.fastq.gz
LNCAP_Normoxia_S2.fastq.gz  PC3_Normoxia_S2.fastq.gz  SRR7179507_pass.fastq.gz  SRR7179511_pass.fastq.gz  SRR7179523_pass.fastq.gz  SRR7179527_pass.fastq.gz
~~~
We won’t need the individual SRA runs anymore, so we can remove them all using the command rm SRR*, which removes all the files in the folder that begin with “SRR”. 
~~~
rm SRR*
~~~
~~~
LNCAP_Hypoxia_S1.fastq.gz  LNCAP_Normoxia_S1.fastq.gz  PC3_Hypoxia_S1.fastq.gz  PC3_Normoxia_S1.fastq.gz
LNCAP_Hypoxia_S2.fastq.gz  LNCAP_Normoxia_S2.fastq.gz  PC3_Hypoxia_S2.fastq.gz  PC3_Normoxia_S2.fastq.gz
~~~

## 2) Mapping reads using STAR
~~~
zcat LNCAP_Hypoxia_S1.fastq.gz | head -4

@SRR7179520.1.1 1 length=76
GTGAANATAGGCCTTAGAGCACTTGANGTGNTAGNGCANGTNGNNCCGGAACGNNNNNNNNAGGTNGNNNGNGTTG
+SRR7179520.1.1 1 length=76
AAAAA#EEEEEEEEEEEEEEEEEEEE#EEE#EEE#EEE#EE#E##EEEEEEEE########EEEE#E###E#EAEA
~~~

**Building the genome index**
~~~
mkdir ref/
wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz
gunzip *.gz
~~~

**Building STAR index**






**Mapping reads**

mapping.sh 
~~~
#!/usr/bin/env bash
#BSUB -J star_map_gencodev44
#BSUB -o logs/star_map_gencodev44.%J.out
#BSUB -e logs/star_map_gencodev44.%J.err
#BSUB -n 16
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=4000]"     # ~4 GB/core -> ~64 GB total
#BSUB -M 64000                  # hard memory limit (MB)
#BSUB -W 24:00                  # walltime hh:mm
# (no -q line; uses default queue)

set -euo pipefail

# =========================
# Project-specific settings
# =========================
PROJECT="/research/groups/yanggrp/home/glin/work_2025/Sep/project_PRJNA416257"
FASTQ_DIR="${PROJECT}/sra/fastq"   # your FASTQs are here (single-end)
GENOME_DIR="${PROJECT}/ref/STAR_index_gencodev44"
GTF="/research/groups/yanggrp/home/glin/work_2025/Sep/project_PRJNA1014743/ref/gencode.v44.annotation.gtf"
OUTROOT="${PROJECT}/mapping"       # all results go here
THREADS=${LSB_DJOB_NUMPROC:-16}    # LSF core count (default 16)

# Optional: tidy small intermediates after each sample
CLEAN_INTERMEDIATES=true

# Optional: activate your environment
# source ~/.bashrc
# conda activate sra
# or: module load STAR

# =========================
# Pre-flight checks
# =========================
mkdir -p "${OUTROOT}" logs

if ! command -v STAR >/dev/null 2>&1; then
  echo "ERROR: STAR not found in PATH." >&2
  exit 1
fi

if [[ ! -s "${GENOME_DIR}/Genome" ]]; then
  echo "ERROR: STAR index appears incomplete at ${GENOME_DIR} (missing Genome file)." >&2
  exit 1
fi

if [[ ! -s "${GTF}" ]]; then
  echo "ERROR: GTF not found at ${GTF}" >&2
  exit 1
fi

shopt -s nullglob
FASTQS=("${FASTQ_DIR}"/*.fastq.gz)
if (( ${#FASTQS[@]} == 0 )); then
  echo "No FASTQ files found in ${FASTQ_DIR}" >&2
  exit 1
fi

echo "Found ${#FASTQS[@]} FASTQ files in ${FASTQ_DIR}"
echo "Output directory: ${OUTROOT}"
echo "Using index: ${GENOME_DIR}"
echo "Using GTF: ${GTF}"
echo "Threads: ${THREADS}"

# =========================
# Alignment loop (single-end)
# =========================
for fq in "${FASTQS[@]}"; do
  sample=$(basename "${fq}" .fastq.gz)
  outdir="${OUTROOT}/${sample}"
  mkdir -p "${outdir}"

  echo "[$(date)] Mapping ${sample}"

  STAR \
    --runThreadN "${THREADS}" \
    --genomeDir "${GENOME_DIR}" \
    --readFilesIn "${fq}" \
    --readFilesCommand zcat \
    --sjdbGTFfile "${GTF}" \
    --twopassMode Basic \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "${outdir}/${sample}_" \
    --quantMode GeneCounts

  # Index BAM if samtools is available
  if command -v samtools >/dev/null 2>&1; then
    samtools index -@ "${THREADS}" "${outdir}/${sample}_Aligned.sortedByCoord.out.bam"
  fi

  # Optional cleanup: keep BAM + BAI + counts + final log
  if [[ "${CLEAN_INTERMEDIATES}" == "true" ]]; then
    rm -f "${outdir}/${sample}_Aligned.out.bam" 2>/dev/null || true
    rm -f "${outdir}/${sample}_SJ.out.tab" 2>/dev/null || true
    rm -f "${outdir}/${sample}_Log.out" 2>/dev/null || true
    rm -f "${outdir}/${sample}_Log.progress.out" 2>/dev/null || true
    rm -f "${outdir}/${sample}_Chimeric.out.junction" 2>/dev/null || true
  fi
done

# =========================
# Summarize mapping metrics
# =========================
summary="${OUTROOT}/mapping_summary.tsv"
echo -e "sample\tUniquely_mapped%\tReads_in_genes_Unstranded\tReads_in_genes_FirstStrand\tReads_in_genes_SecondStrand" > "${summary}"

for log in "${OUTROOT}"/*/*_Log.final.out; do
  s=$(basename "${log}" _Log.final.out)
  # Parse "Uniquely mapped reads %" from STAR Log.final.out (robust to spacing)
  uniq_pct=$(awk -F '|' '/Uniquely mapped reads %/ {gsub(/%/,"",$2); gsub(/^[ \t]+|[ \t]+$/,"",$2); print $2}' "${log}")
  rpg_dir="${OUTROOT}/${s%_*}"
  rpg="${rpg_dir}/${s}_ReadsPerGene.out.tab"
  if [[ -f "${rpg}" ]]; then
    # Sum counts across genes (skip header rows N_*)
    unstr=$(awk 'BEGIN{sum=0} !/^N_/ {sum+=$2} END{print sum}' "${rpg}")
    firsts=$(awk 'BEGIN{sum=0} !/^N_/ {sum+=$3} END{print sum}' "${rpg}")
    seconds=$(awk 'BEGIN{sum=0} !/^N_/ {sum+=$4} END{print sum}' "${rpg}")
  else
    unstr=NA; firsts=NA; seconds=NA
  fi
  echo -e "${s}\t${uniq_pct}\t${unstr}\t${firsts}\t${seconds}" >> "${summary}"
done

echo "[$(date)] Done."
echo "Results per sample are in: ${OUTROOT}/<sample>/"
echo "Project summary written to: ${summary}"
~~~

~~~
mkdir -p logs mapping
bsub < mapping.sh
~~~

This will produce files like:
mapping/mapping_summary.tsv
...

**select the read counts by sample**
mapping/make_counts.py
~~~
#bash 
python3 make_counts.py
~~~
This will produce files like:
  Counts: raw_counts.csv
  QC:     qc.csv

**select the read counts by sample**clear 

















