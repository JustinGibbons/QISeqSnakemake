# QISeqSnakemake

QISeqSnakemake is an implementation of the [QISeq workflow](https://github.com/ThomasDOtto/QISeq) in snakemake. The original publication by [Bronner et al 2016](https://pubmed.ncbi.nlm.nih.gov/27197223/) was published in Genome Research.

QISeq is a forward genetics approach that allows the identification of genes involved in a phenotype without the need to isolate clones. The QISeq methodology is a modification of the [TraDIS method](https://pubmed.ncbi.nlm.nih.gov/26794317/) to work with _Plasmodium_ parasites that uses a _piggyBac_ transposase that randomly inserts into TTAA sites and a modified Illumina sequencing protocol [(Bronner et al 2016)](https://pubmed.ncbi.nlm.nih.gov/27197223/). 

[QISeq has been used to identify](https://www.ncbi.nlm.nih.gov/sites/myncbi/justin.gibbons.1/collections/66719328/public/) the essential genomes of 2 _Plasmodium_ species as well as identifying genes important for drug response, fever response, oxidative stress, infection of sickle-trait cells and gametocyte development. 

### What's new?
  1. _Snakemake_ makes the workflow more user friendly and modifiable.
  2. _Snakemake_ allows the workflow to efficiently scale to larger sample sizes and sequencing depths                  through better management of job submissions and handeling of temporary files.

## Installation

QISeqSnakemake can be installed and run on Linux and Mac operating systems.

### Installing dependencies
A conda environment will be created to install all of the dependencies.

```Insert code for creating conda environment```

### Snakefiles and Slurm submission

Two snakefiles are provided in this GitHub repository: snakefile_forward and snakefile_reverse. The one you use depends on the library preparation method used to create the QISeq Illumina library. The snakefile_forward file is used to run jobs when the transposon sequence is in the R1 read and snakefile_reverse is used when the transpson sequence is in the R2 read as illustrated in Figure 1:
![figure1](https://github.com/user-attachments/assets/d7ffb9d1-d400-429b-b2b8-b43cb8b4a7d5)

A script named run_qiseq_pipeline.sh is also provided to assist with submitting the snakemake workflows to the Slurm workload manager. For submitting the workflows using other workload managers please see the [snakemake documentation](https://snakemake.readthedocs.io/en/v7.7.0/executing/cluster.html).

## Workflow overview

There are 4 general stages to QISeq data processing and analysis:
  1. Preprocessing
      - Transposon sequence trimming.
      - Indexing of reference genome.
  2. Sequence alignment and alignment processing
     - Align reads to genome using Bowtie.
     - Sort alignments, add read groups and convert to bam format.
     - Mark duplicates to allow counting unique insertions
  3. Insertion identification and quantification
      - Count number of unique reads at each site.
      - Count total number of reads mapped to each of the identified sites.
      - Combine insertion site counts from each sample into a single file. 
  4. QC, Annotation and analysis
     - Identifying insertion distances from TTAA sites and select a maximum distance therhold.
     - Choose read count cutoff.
     - Keep mutants identified in multiple replicates.
     - Annotate the quality controled insertion site to the nearest gene.
     - Analyze the resulting annotated counts table.
     

The QISeqSnakemake workflow performs the first 3 stages. 
![figure2](https://github.com/user-attachments/assets/d0e2de6c-ac64-4a8e-bf9d-b7288f3eba51)

## Inputs and Outputs

### Input files
- Sample Paired-end fastq files
- Reference genome as a fasta file

### Input parameters
- **Offset**: Used to calculate the insertion site position. Varies by sequencing platform:
    - NovaSeq: -2
    - HiSeq: 4
    - Miseq: 5
- **insertSize**: The insert size selected for during library prepartion. Used by bowtie2 to search for valid alignments.
- **Min_Unique**: Minimum number of unique reads required to identify an insertion site.
- **Merge_distance**: Merge insertion calls between samples if they are $\leq$ Merge_distance bp from each other.

## Output files
- Insertion site counts
- Sample bam files
- Samtools stats reports







