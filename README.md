This project is related with RNA-seq analysis for Lucia Clemens-Daxinger lab.
It includes:
- Differentian gene expressiona analysis
- Differential repetitive classes analysis
- Differential exon expression analysis
- QC
- Contamination checks

Requeriments
-
- Linux OS/ MacOS X
- Bioconda (https://bioconda.github.io/)
- Snakemake (http://snakemake.readthedocs.io/en/latest/getting_started/installation.html)



# Installation
-
Go to the directory where you want to make the analysis, eg. /exports/humgen/USERDIR (replace USERDIR by your directory in Shark)
```
git clone https://git.lumc.nl/dsanleongranado/RNA-seq-snakemake.git
```
This will create a directory named RNA-seq-snakemake and inside there are the scripts and files.

An alternative way to do it is download the project (download button on https://git.lumc.nl/dsanleongranado/RNA-seq-snakemake web page) and unzip on the location you will execute the pipeline
Preparation

# Run

## Modify samples file 
Go inside the directory

```
cd RNA-seq-snakemake
```
Create a file called samples.tsv, it should be a text file with tab separated columns. The first row should include the column names.

If it is to analyze a GEO data:
* Download the Table of the SRA project. It is a text file with information about the samples.
* Rename the column where the run id appears by **SampleID**
* Rename the column where there is the condition information by **condition**. 
* If there are samples you don't want, just delete the rows related to those samples
* Rename the text file to samples.tsv and copy it inside the RNA-seq-snakemake

If the analysis is over fastq files, create text file (samples.tsv) with the following columns:
- SampleID: sample name
- Condition
- R1: fastq file path for forward reads
- R2: fastq file path for reverse reads


## Execution
-
If the cluster you are executing the pipeline is DMRAA compatible:
```
sh snakemakeLauncher.sh
```
# Software included
This is the software included in this pipeline
* [MultiQC](https://multiqc.info/)
* [SRA tools](https://github.com/ncbi/sra-tools)
* [centrifuge](https://ccb.jhu.edu/software/centrifuge/index.shtml)
* [Krona](https://github.com/marbl/Krona/wiki)
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
* [STAR aligner](https://github.com/alexdobin/STAR)
* [Samtools](http://www.htslib.org/)
* [Picard tools](https://broadinstitute.github.io/picard/)
* [rnaseqc](https://github.com/broadinstitute/rnaseqc)
* [HTseq-count](https://htseq.readthedocs.io/en/release_0.11.1/count.html)
* [Bedtools](https://bedtools.readthedocs.io/en/latest/)
* [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
* [DEXSEQ](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html)
* [Glimma](https://bioconductor.org/packages/release/bioc/html/Glimma.html)
* [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)


## Modifying experiment_params.yaml
The main cofiguration of the pipeline is inside the file experiment_params.yalm

This file contains the variables related to the experiment:
- GENOME_ZIP_FASTA_URL: genome sequence path
- GENOME_ZIP_GFF_UR: genome annotation path
- GENOME_ZIP_GTF_URL: genome annotation in GTF format
- GENOME_REFFLAT: genome annotation in refflat format
- GENOME_CHROM_SIZE: chromosomes sizes
- REPCLASSES_ZIP_bed_URL: bed file with repetitive elements positions
- REPCLASSES_ZIP_TABLE_URL: table with the repetitve element information. It contains 3 columns: repetitive name - repetitive class - repetitive family
- centrifuge: contains the index id and the path to the centrifuge index
- samples: if you have your samples file in other place or other name, here you can change it
- reference: the path where you want store the index of the genome or an existing genome
- contrast: contrast to check (NOTE: for now only accept two conditions)
- sraDatabase: Path to the GEO data when it is downloaded
- adapter5 and adapter3: adapters of the sequences if they have. By default are automatically detected but if it is not possible, the values of those variables will be selected.
- stranded: it indicates if the RNA-seq is strand specific or not. It can be "no" "RF" or "FR".
- preprocess: This flag is to indicate if the filtering, quality controls etc is previously done. and the pipeline will start with the alignment. If the value is "no" the pipeline will start from the begining.
- multiplemapping: it activates the multiple mapping for the differential expression.
- deseq2: it activates the gene differential expression analysis
- DEXEQ: It activates the exon differential expression analysis
- repetitiveAnalysis: It activates the repetitive classs differential expression analysis
- QC: It activates the quality controls.

By default, the genome is mm10 in mouse and the annotation is Gencode

## Modifying cluster_config.yml

Some times it is necessary to add more memory for a specific step (e.q. the number of samples is very big), and in this file it is possible to allow more memory.
To do that, look for the name of the rule that needs memory and change the value.

plotKrona:
  mem: 10G
this step can use 10Gb of memory and we need 30 Gb to run our pipeline because we have a lot of samples.

plotKrona:
  mem: 30G

## Output

The results will be inside the folder analysis. Inside this folder there are several folders with the outputs:

```
+-- experiment_params.yml - Configuration file with variables like, genome fasta file, annotation file, steps to run, (see section Modifying experiment_params.yaml)
+-- samples.tsv - Text file with the samples information, sampleID, path to the RAW data, samples metadata.
+-- annotation - 
|   +-- DEXSEQ.gff - Annotation created by the pipeline to run the differential exon analysis
+-- counts
|   +-- all.tsv - Table with read counts per gene for every sample. The rows are the genes and the columns are the samples.
|   +-- sample1.countsPerGeneperExon.csv - File with read counts per exon and gene. There are one file per sample.
|   +-- ...
|   +-- repetitive_elements.tsv - Table with read counts per repetitive element in the genome. The rows are the repetitive elements and the colums are the samples. The order of the samples is the same than the order in samples.tsv
+-- dataDir - Folder with a copy of original data, in the case that the original samples are divided in multiple files. This folder will contain the concatenated files.
|   +-- sample1_R1.fastq.gz - forward reads . There are 1 file per sample
|   +-- sample1_R2.fastq.gz - reverse reads. There are 1 file per sample
|   +-- ...
+-- filteredDataDir - Folder with the reads without the adaptors
|   +-- sample1_R1_val_1.fq.gz
|   +-- sample1_R2_val_1.fq.gz
|   +-- ...
+-- finalDataDir - Folder with the reads without the adaptors and low quality reads
|   +-- sample1_R1.fastq.gz
|   +-- sample1_R2.fastq.gz
|   +-- ...
+-- mappedDataDir - Folder with the alignments
|   +-- sample1
|       +-- Sample1.dedup.bam - Deduplicated alignments
|       +-- Sample1.dedup.bai - Index of the deduplicated alignment
|       +-- Sample1S214Aligned.sortedByCoord.out.bam - Alignment (it includes replicates)
|       +-- Sample1S214Aligned.sortedByCoord.out.bam.bai - Index of the alignment
|       +-- Sample1_STARgenome - Temporal data created by STAR . See [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) 
|       +-- Sample1_STARpass1 - Temporal data created by STAR
|       +-- Sample1_STARtmp - Temporal data created by STAR
|       +-- Sample1Log.final.out - Temporal data created by STAR
|       +-- Sample1Log.out - Temporal data created by STAR
|   +-- ... There are one folder per sample
+-- log - Folder with software logs
|   +-- centrifuge - Logs of centrifuge software
|   +-- mapping - Logs of STAR software
|   +-- trimGalorePE - Logs of trimGalore software
+-- QC - Folder with quality control information
|   +-- finalDataDir - Folder with FastQC output of final reads
|   +-- dataDir - Folder with FastQC output of original reads
|   +-- filterDataDir - Folder with FastQC output of filtered reads
|   +-- picard - Folder with picard output of final reads
|   +-- rnaseqc - Folder with rnaseqc output of final reads
|   +-- samtools - Folder with samtools output of final reads
+-- results
|   +-- centrifuge - Folder with centrifuge output. It contains the files with the number of reads per bacteria/virus in our sample.
|   +-- DEXSEQ - Folder with the table with differential exon expression analysis. It contains a pdf with some figures like, PCA plot, clustering plot and MA plot. It contains a webpage to navegate in the results.
|   +-- diffexp - Folder with the table with differential gene expression  analysis. It contains a pdf with some figures like, PCA plot, clustering plot and MA plot. It contains a webpage to navegate in the results.
|   +-- krona - Folder with a webpage with the information of centrifuge folder. It is the final contamination report.
|   +-- report - Folder with the final QC report.
+-- toUCSC : folder with the bigwig files. They can be used to see the data in UCSC browser.
+-- toUCSC_norm : folder with the bigwig files. This files are normalized by RPKM. They can be used to see the data in UCSC browser.

```
-
## Final tables
Inside the results, if the pipeline finished properly it si possible to find excel files with the results.
The files contains gene annotations like position in the genome, description, ENSEMBL ID, gene name, gene type.
They contains some important information like the logFC (log2 foldchange), logCPM ( log2 normalized averaged expression), p-value, adj-pvalue (FDR, BH adjustment), the counts per gene and sample. 

There is one table inside results/diffexp with the differential gene expression.
There is one table inside results/DEXSEQ with the differential exon expression.
There is one table inside results/repClasses with the differential gene expression.

## Final report

The final report can be found in analysis/results/report/multiqc_report.html, this report it was generated with [MultiQC](https://multiqc.info/). In the MultiQC manual it is possible to find more details.

There are multiple sections that they can be accessed in the left pannel.

### Sections

1.  General Statistics : Table with importan information like:
    
  *  % rRNA
  *  % mRNA
  *  % Aligned reads
  *  averaged insert size
  *  % of duplicates. Value from alignment
  *  % of aligned reads versus the total size of the sample
  *  Number of aligned reads (in millions)
  *  % of trimmed reads. The trimmed reads are the reads that they contains the adaptor or low quality at the end of the sequence
  *  % of duplicates . Value from the fastq file. This value don't take in account the other pair, so, it is less accurated than the % of duplicates from the alignment.
  *  % GC percentaje of Gs and Cs
  *  Number of sequences. In milllions.

The table will contain several rows per sample because it has information from the alignment and from the individual files from the sample (R1 reads and R2 reads).
There are extra rows with the sufix _val_1 and _val_2 because this is the information from the filtered reads. And there is per sample an extra row with the sufix STARpass1 because the pipeline align 2 times the reads to the genome to find new splicing events. And this line contains the first alignment.

2.  HTSeq Count:

HTSeq Count is part of the HTSeq Python package - it takes a file with aligned sequencing reads, plus a list of genomic features and counts how many reads map to each feature.
   
3. Picard

  *  Alignment Summary: Information about hoy many reads are alignmed. Numbers obtained by Picard tools. Plase note that Picard's read counts are divided by two for paired-end data.
  *  Base Distribution: Plot shows the distribution of bases by cycle.
  *  Insert size: Plot shows the number of reads at a given insert size. Reads with different orientations are summed.
  *  Mark Duplicates: Plots shows the amount of duplicates in the alignment
  *  RnaSeqMetrics Assignment : Number of bases in primary alignments that align to regions in the reference genome. They are divided by Coding regions, UTR, Intronic, Intergenic and not aligned.

4. STAR: Statistics generated by the aligner.
  * Aligment Scores. Information about the amount of reads uniquely mapped, mappled to multiple loci or discarded.
  * Gene Counts. Amout of reads mapping versus genes.
5. Cutadapt: This plot shows the number of reads with certain lengths of adapter trimmed. Obs/Exp shows the raw counts divided by the number expected due to sequencing errors. A defined peak may be related to adapter length. See the cutadapt documentation for more information on how these numbers are generated.
6. FastQC
  *  Sequence Counts: Sequence counts for each sample. Duplicate read counts are an estimate only.
  *  Sequence Quality Histograms: The mean quality value across each base position in the read.
  *   Per Sequence Quality Scores: The number of reads with average quality scores. Shows if a subset of reads has poor quality.
  *  Per Base Sequence Content: The proportion of each base position for which each of the four normal DNA bases has been called. NOTE: In Illumina always the first 9 nucleotides has a bias.. it is not a problem
  *  Per Sequence GC Content: The average GC content of reads. Normal random library typically have a roughly normal distribution of GC content.
  *  Per Base N Content: The percentage of base calls at each position for which an N was called.
  *   Sequence Length Distribution: The distribution of fragment sizes (read lengths) found.
  *   Sequence Duplication Levels: The relative level of duplication found for every sequence.
  *  Overrepresented sequences: The total amount of overrepresented sequences found in each library.
  *   Adapter Content: The cumulative percentage count of the proportion of your library which has seen each of the adapter sequences at each position.

## Developing version

-
```
git clone https://git.lumc.nl/dsanleongranado/RNA-seq-snakemake.git --branch DEXSEQ --single-branch DEXSEQ_snakemake
```
## M & M
Quality assessment of the raw sequencing reads was done using FastQC  v0.11.2. Adapters were removed by TrimGalore v0.4.5 using default parameters for paired-end Illumina reads, after which, quality filtering was performed by the same software. Reads smaller than 20bps and those with an error rate (TrimGalore option "-e") higher than 0.1 were discarded, after which a final quality assessment of the filtered reads was done with FastQC to identify possible biases left after filtering. The remaining reads were mapped to  the mouse reference genome (build mm10) using the STAR aligner v2.5.1 using default parameters with the following exceptions: “–outputMultimapperOrder random” and  “–twopassMode basic”. The quantification was done by HTSeq-count v0.91 using the GENCODE MV23 annotation with the option “–stranded no”. The statistical analysis was done using DESeq2 v1.2.0 (R package). The final list of differential expressed genes contains genes for which the adjusted p-value (Benjamini-Hochberg correction) < 0.05.
