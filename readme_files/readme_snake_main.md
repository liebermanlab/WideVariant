# Snakemake pipeline


## Overview

This is a tutorial for the Snakemake pipeline portion of Lieberman Lab's SNV-calling pipeline. The core function of the Snakemake pipeline is to process raw sequencing reads (remove adapters and trim reads based on quality), align reads to a reference genome, identify candidate SNV positions, and save data arrays with stats on the alignment of each sample at each candidate SNV position. These data arrays are used downstream in the [local analysis](readme_local_main) step. 


### Alternative modes of the Snakemake pipeline

This tutorial focuses on using the Snakemake pipeline to align reads to a reference genome and identify candidate SNVs. However, it should be noted that the Snakemake pipeline can run in alternative modes, including the "assembly" mode which uses the processed sequencing reads from each sample to create a sample-specific genome assembly and the "bracken" mode which uses processed sequencing reads to estimate taxonomic abundance. For more information, see [How to run the Snakemake pipeline](readme_snake_run.md).


## Philosphy and organization of WideVariant

The philosophy and organization of the pipeline, as well as detailed descriptions of the data structures used and steps executed in the Snakemake pipeline can be found [here](readme_snake_organization.md). This organization is fairly different than standard different bioinformatic workflows because it designed to retain raw data, and to be easy to rerun on different subsets of the data. 


## How to run the Snakemake pipeline

Here are instructions for [how to run](readme_snake_run.md) the Snakemake pipeline.

The inputs to the Snakemake pipeline include: raw sequencing data files (forward and reverse fastq files for each bacterial isolate) and a reference genome (fasta file).


## Detailed descriptions of the Snakemake

Detailed descriptions of the steps executed in the Snakemake pipeline can be found [here](readme_snake_rules.md). Many steps use well-established bioinformatics tools.


## Helpful hints for using the command line

Users who are not already familiar with the linux terminal may find [this brief overview](readme_snake_basics.md) helpful.


## Wishlist

Here is a [wishlist](readme_snake_wishlist.md) of upgrades to the current version of the Snakemake pipeline.


## Table of contents

[Main Lieberman Lab pipeline README](../README.md)
* [Snakemake pipeline](readme_snake_main.md)
	* [Philosophy and technical details about the snakemake pipeline](readme_snake_rules.md)
	* [How to run the snakemake pipeline](readme_snake_run.md)
	* [Wishlist for snakemake pipeline upgrades](readme_snake_wishlist.md)
	* [Helpful hints for using the command line](readme_snake_basics.md)
* [Local analysis](readme_local_main.md)
	* [How to run the local analysis script](readme_local_run.md)
	* [Wishlist for local analysis upgrades](readme_local_wishlist.md)
	* [Python best practices](readme_local_best.md)
