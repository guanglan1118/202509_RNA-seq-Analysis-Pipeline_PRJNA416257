# RNA-seq Analysis Pipeline / 202509_PRJNA416257
## Datasets
- Dataset PRJNA416257
- Nat Commun.2019/PMID: 30655535
- ONECUT2 is a driver of neuroendocrine prostate cancer

## Folder layout
~~~
# bash
project_PRJNA416257/
├─ raw/              # sra
├─ raw_fastq/        # FASTQs
├─ ref/              # reference (Salmon index, GTF/FA)
├─ qc/               # FastQC & MultiQC
├─ quant/            # Salmon outputs per sample
├─ meta/             # metadata.csv, tx2gene.csv
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
- raw/SRR26030905/SRR26030905.sra
- raw/SRR26030906/SRR26030906.sra
- raw/SRR26030907/SRR26030907.sra
- raw/SRR26030908/SRR26030908.sra
- raw/SRR26030909/SRR26030909.sra
- raw/SRR26030910/SRR26030910.sra

### 1.3) Converts .sra archive files into plain FASTQ files
*convert single file* 
~~~
#bash
mkdir -p raw_fastq
fasterq-dump raw/SRR26030905/SRR26030905.sra -e 8 -p --split-files -t tmp -O raw_fastq/
~~~

*convert batch file* 
~~~
#bash
pwd #project_PRJNA1014743
for sra in raw/*/*.sra; do
  SRR=$(basename "$sra" .sra)
  echo "==> Converting $SRR"
  fasterq-dump "$sra" -e 8 -p --split-files -t tmp -O raw_fastq/
done
~~~
This will produce files like:
- raw_fastq/SRR26030905_1.fastq
- raw_fastq/SRR26030905_2.fastq
- raw_fastq/SRR26030906_1.fastq
- raw_fastq/SRR26030906_2.fastq
- raw_fastq/SRR26030907_1.fastq
- raw_fastq/SRR26030907_2.fastq
- raw_fastq/SRR26030908_1.fastq
- raw_fastq/SRR26030908_2.fastq
- raw_fastq/SRR26030909_1.fastq
- raw_fastq/SRR26030909_2.fastq
- raw_fastq/SRR26030910_1.fastq
- raw_fastq/SRR26030910_2.fastq

## 2) Raw data QC
### 2.1) intsall fastqc
~~~
# bash
conda install -c bioconda fastqc
fasterq-dump --version  #FastQC v0.12.1

conda install -c bioconda multiqc
multiqc --version  #multiqc, version 1.30     
~~~

### 2.2) run QC
~~~
# bash
mkdir -p qc/fastqc qc/multiqc
# Run FastQC (6 threads)
fastqc -t 6 -o qc/fastqc raw_fastq/*.fastq

# Summarize reports
multiqc -o qc/multiqc qc/fastqc
~~~
This will produce files like:
- qc/fastqc/SRR26030905_1.fastqc.html
- qc/fastqc/SRR26030905_1.fastqc.zip
- qc/fastqc/SRR26030905_2.fastqc.html
- qc/fastqc/SRR26030905_2.fastqc.zip
...

- qc/multiqc/multiqc_report.html
1. ~25M reads per sample
2. ~50% GC content (normal for human)
3. Low adapter contamination
4. <1% overrepresented sequences
5. Quality scores are consistently high


## 3) Quantification 
### 3.1) Build a decoy-aware index
~~~
# bash
# ref/
mkdir -p ref

# Download (use HTTPS instead of FTP to avoid firewall issues)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz

gunzip *.gz
~~~

This will produce files like:
- gencode.v44.annotation.gtf  
- gencode.v44.transcripts.fa  
- GRCh38.primary_assembly.genome.fa

~~~
# bash
# Create decoys list (chromosome headers from the genome FASTA)
grep "^>" GRCh38.primary_assembly.genome.fa | cut -d " " -f1 | sed 's/>//g' > decoys.txt
~~~
This will produce files like:
- decoys.txt 
 
~~~
# bash
# Make gentrome (transcripts + genome)
cat gencode.v44.transcripts.fa GRCh38.primary_assembly.genome.fa > gencode.v44.gentrome.fa
~~~
This will produce files like:
- gencode.v44.gentrome.fa 

~~~
conda activate sra
conda install -c bioconda salmon=1.10.3
salmon --version  #version : 1.10.3
~~~

~~~
# Build Salmon index (decoy-aware)
salmon index \
  -t gencode.v44.gentrome.fa \
  -d decoys.txt \
  -i salmon_gencode_v44_decoy \
  --gencode \
  -p 8
~~~

This will produce files like:
[info] Building perfect hash
[info] Index built successfully


~~~
quant/
 ├── CTRL1/
 ├── CTRL2/
 ├── CTRL3/
 ├── TRT1/
 ├── TRT2/
 └── TRT3/
~~~
 
### 3.2) Quantify with recommended flags
Bias correction and selective alignment generally improve estimates.
~~~
# adjust sample names to match your files (yours looked like SRR* earlier)
mkdir -p ../quant
for S in CTRL1 CTRL2 CTRL3 TRT1 TRT2 TRT3
do
  salmon quant \
    -i refs/salmon_gencode_v44_decoy \
    -l A \
    -1 raw_fastq/${S}_R1.fastq \
    -2 raw_fastq/${S}_R2.fastq \
    --validateMappings \
    --gcBias \
    --seqBias \
    --numBootstraps 100 \
    -p 8 \
    -o quant/$S
done
~~~

### 3) Prepare tx2gene for gene-level summarization
You’ll need this for tximport → DESeq2/edgeR.
~~~
# From the GENCODE GTF
awk '$3=="transcript" { 
  tx=""; gene=""; 
  for(i=9;i<=NF;i++){
    if($i~/^transcript_id/) {tx=$(i+1); gsub(/"|;/, "", tx)}
    if($i~/^gene_id/)       {gene=$(i+1); gsub(/"|;/, "", gene)}
  }
  if(tx!="" && gene!="") print tx"\t"gene
}' gencode.v44.annotation.gtf > tx2gene_gencode_v44.tsv
~~~




## 4) Import counts into R

~~~
library(tximport)
library(DESeq2)
library(readr)

# Metadata
coldata <- read.csv("metadata.csv", row.names=1)

# Files
files <- file.path("quant", rownames(coldata), "quant.sf")
names(files) <- rownames(coldata)

# Transcript → gene mapping
tx2gene <- read.csv("tx2gene_gencode_v44.csv")  # transcript_id,gene_id

# Import
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

dds <- DESeqDataSetFromTximport(txi, colData=coldata, design=~ batch + condition)

# Prefilter
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep,]
~~~

## 5) QC & Visualization
~~~
vsd <- vst(dds)

plotPCA(vsd, intgroup=c("condition","batch"))

library(pheatmap)
dists <- dist(t(assay(vsd)))
pheatmap(as.matrix(dists))
~~~

## 6) Differential Expression

~~~
dds <- DESeq(dds)

res <- results(dds, contrast=c("condition","treated","control"))
res <- lfcShrink(dds, contrast=c("condition","treated","control"), type="apeglm")

summary(res)

res_sig <- subset(as.data.frame(res), padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(res_sig, "DE_genes.csv", row.names=TRUE)
~~~


## 7) Visualization (Volcano & Heatmap)
~~~
library(ggplot2)

res_df <- as.data.frame(res)
ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(alpha=0.5) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", color="red") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="blue")

topgenes <- head(order(res$padj), 30)
pheatmap(assay(vsd)[topgenes,], cluster_rows=TRUE, cluster_cols=TRUE,
         annotation_col=coldata)

~~~


## 8) Pathway Enrichment (GSEA)
~~~
library(fgsea)
library(msigdbr)

msig <- msigdbr(species="Homo sapiens", category="H") |>
        split(~gene_symbol)

ranks <- res_df$log2FoldChange
names(ranks) <- rownames(res_df)

fg <- fgsea(msig, stats=ranks, minSize=15, maxSize=500, nperm=10000)
write.csv(fg[order(fg$padj),], "GSEA_results.csv")
~~~

## 9) Deliverables

At the end you would have:

- QC reports (FastQC, MultiQC, PCA, heatmaps).

- Counts matrix and DESeq2 objects (dds.rds, vsd.rds).

- DE gene table (DE_genes.csv).

- Figures (Volcano, Heatmap).

- Pathway results (GSEA_results.csv).




















