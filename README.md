# snakeMake

## Overview

This repository contains a collection of Snakemake workflows and supporting scripts for genomics and transcriptomics analysis, designed for use on SGE or SLURM clusters (e.g., VSC). The pipelines cover a range of tasks including genome assembly, annotation, RNA-seq analysis, chimeric RNA detection, and more.

## Features

- **Multiple Genome Alignment**: Automated workflows for aligning multiple genomes and extracting conserved elements.
- **De Novo Genome Assembly**: Pipelines for assembly using Newbler and MIRA.
- **Genome Annotation**: Automated annotation using evidence-based and ab initio methods (in progress).
- **Chimeric RNA Extraction**: Identify chimeric RNAs from STAR-mapped reads.
- **Functional Annotation**: Annotate predicted proteins with InterProScan, MEROPS, CAZy, TCDB, and more.
- **Differential Expression**: Identify differentially expressed genes from RNA-seq count data.
- **lncRNA Discovery**: Pipeline for long non-coding RNA identification and quantification.

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- Python 3.x
- R (with edgeR, limma, optparse, etc.)
- Perl (for some utility scripts)
- Cluster environment (SGE or SLURM recommended)
- Bioinformatics tools: STAR, Trinity, Cufflinks, PASA, Augustus, SNAP, InterProScan, BLAST, bedtools, featureCounts, RepeatMasker, etc.

See each workflow for specific dependencies.

## Usage

Each workflow is defined in a `.Snakemake` file. To run a workflow:

```bash
snakemake -d $(pwd) -s $(pwd)/<workflow>.Snakemake --stats snakemake.stats -j <jobs> --cluster 'qsub <cluster_params>'
```

Replace `<workflow>` with the desired pipeline (e.g., `genomeAnnotation`, `differentialExpression`, `chimericRNA`, etc.).

## Workflow Descriptions

Below is a summary of each main workflow, with typical inputs and outputs. See the top of each `.Snakemake` file for more details and required file formats.

### genomeAnnotation.Snakemake
- **Purpose:** Automated genome annotation using ab initio and evidence-based methods (e.g., PASA, Trinity, Cufflinks, Augustus, SNAP, Scipio, CEGMA).
- **Input:** Genome FASTA, RNA-seq reads, protein databases.
- **Output:** Annotated GFF/GTF files, intermediate files for each tool.

### differentialExpression.Snakemake
- **Purpose:** Differential gene expression analysis from RNA-seq count data.
- **Input:** FeatureCounts or similar count matrix, sample metadata.
- **Output:** Lists of differentially expressed genes, summary statistics.

### functionalAnnotation.Snakemake
- **Purpose:** Functional annotation of predicted proteins using InterProScan, MEROPS, CAZy, TCDB, and BLAST.
- **Input:** Protein FASTA files, annotation databases.
- **Output:** Annotated TSV/merged files, summary tables.

### chimericRNA.Snakemake
- **Purpose:** Detect and extract chimeric RNAs from STAR-mapped RNA-seq data.
- **Input:** STAR chimeric junction files, genome FASTA.
- **Output:** Chimeric RNA lists, GFF/TSV files.

### lncRNA.Snakemake
- **Purpose:** Identify and quantify long non-coding RNAs (lncRNAs) from transcriptome assemblies.
- **Input:** Genome annotation, transcript assemblies, protein databases, mapped reads.
- **Output:** lncRNA GFF files, count tables, differential expression results.

### multipleGenomeAlignments.Snakemake
- **Purpose:** Perform multiple genome alignments and extract conserved elements.
- **Input:** Genome FASTA files for all species.
- **Output:** Alignment files (MAF, GFF), summary statistics.

### Other scripts
- Supporting scripts in Perl and R are used for data processing, statistics, and format conversion. See comments in each script for usage.

### Example Commands

#### Genome Annotation
```bash
snakemake -d $(pwd) -s $(pwd)/genomeAnnotation.Snakemake --stats snakemake.stats -j 100 --cluster 'qsub {params.cluster}'
```

#### Differential Expression
```bash
snakemake -d $(pwd) -s $(pwd)/differentialExpression.Snakemake --stats snakemake.stats -j 100 --cluster 'qsub {params.cluster}'
```

#### Functional Annotation
```bash
snakemake -d $(pwd) -s $(pwd)/functionalAnnotation.Snakemake --stats snakemake.stats -j 100 --cluster 'qsub {params.cluster}'
```

#### Chimeric RNA Extraction
```bash
snakemake -d $(pwd) -s $(pwd)/chimericRNA.Snakemake --stats snakemake.stats -j 100 --cluster 'qsub {params.cluster}'
```

#### lncRNA Discovery
```bash
snakemake -d $(pwd) -s $(pwd)/lncRNA.Snakemake --stats snakemake.stats -j 100 --cluster 'qsub {params.cluster}'
```

#### Multiple Genome Alignment
```bash
snakemake -d $(pwd) -s $(pwd)/multipleGenomeAlignments.Snakemake --stats snakemake.stats -j 100 --cluster 'qsub {params.cluster}'
```

## Directory Structure

- `*.Snakemake` — Main Snakemake workflow files for each analysis type
- `*.pl`, `*.R` — Supporting scripts (Perl, R)
- `InterProscanDictionary.dat` — Example annotation dictionary
- `README.md` — This file

## Notes

- Some workflows require editing paths and environment variables to match your cluster setup.
- See comments in each workflow for details on required input files and expected outputs.
- For questions or contributions, please open an issue or pull request.





