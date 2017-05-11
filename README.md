# seqcheck
Simple mini-pipeline to sanity check new sequencing data

This [Nextflow](https://www.nextflow.io/) pipeline processes a set of FASTQ files (completely independently &mdash; i.e., treating individual read files from a paired end experiment as separate samples) using the following tools:

+ [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - performs basic QC on raws reads in the FASTQ file
+ [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) - aligns FASTQ to a reference genome
+ [Salmon](https://salmon.readthedocs.io/en/latest/) - maps FASTQ to a reference transcriptome

The pipeline was designed with the intention of using a simple series of *fast* tools to obtain basic QC info about raw and aligned/mapped reads for each file.

**Note:** no adaptor or quality trimming of reads is performed in the current version of the pipeline, which could affect alignment/mapping rates.

Outputs from all tools for all samples are combined into a single summary HTML report using [MultiQC](http://multiqc.info/).

## Usage

For a folder of FASTQ files on Synapse:
```
./nextflow run jaeddy/seqcheck --syndir syn8262420 --build GRCh38
```

&nbsp;

or, if FASTQ files are already stored locally:
```
./nextflow run jaeddy/seqcheck --indir <path-to-fastq-folder> --build GRCh38
```

## Resources

The pipeline currently behaves best when run on a single machine (either locally or on Amazon EC2). Further testing/development is needed to enable distributed processing across a cluster of compute notes (including with Nextflow's autoscaling capabilities).

## Outputs

Some assumptions for interpreting overall results in MultiQC:

+ Samples from whole genome or exome sequencing should generally have higher rates of alignment/mapping (**% Aligned**) to the genome (with HISAT2) than to the transcriptome (with Salmon) &mdash; this difference should be more pronounced for whole genome FASTQs.
+ Samples from RNA sequencing should have roughly equal alignment/mapping rates to the genome and transcriptome.

Other QC metrics (duplication rates, quality score distributions, GC content, etc.) can be interpreted based on the [documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for FastQC.
