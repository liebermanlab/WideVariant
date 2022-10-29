# Lieberman Lab pipeline tutorial

This is a tutorial for running the Lieberman Lab's core SNV-calling pipeline.

## Overview

The purpose of this pipeline is to detect single nucleotide differences between closely related bacterial isolates. 

* Inputs: 
	* short-read sequencing data of closely related bacterial isolates
	* an annotated reference genome
* Outputs: 
	* table of high-quality SNVs that differentiate isolates from each other
	* parsimony tree that visualizes how the isolates are related to each other 

The pipeline is split into two main components, as described below.


### 1. Snakemake pipeline

The first portion of the pipeline aligns raw sequencing data from bacterial isolates to a reference genome and identifies candidate SNV positions. This step is implemented in a workflow management system called [Snakemake](http://snakemake.readthedocs.io) and is executed on a [SLURM cluster](https://slurm.schedmd.com/documentation.html). More information is available [here](readme_files/readme_snake_main.md).


## 2. Local python analysis

The second portion of the pipeline filters candidate SNVs based on data arrays generated in the first portion and generates a high-quality SNV table and a parsimony tree. This step is implemented with a custom python script. More information can be found [here](readme_files/readme_local_main.md).



## Example use cases

Previous iterations of this pipeline have been used to study:
* [_C. acnes_ biogeography in the human skin microbiome](https://www.sciencedirect.com/science/article/pii/S1931312821005783)
* [Adaptive evolution of _S. aureus_ on patients with atopic dermatitis](https://www.biorxiv.org/content/10.1101/2021.03.24.436824v3)
* [Adaptive evolution of _B. fragilis_ on healthy people](https://www.sciencedirect.com/science/article/pii/S1931312819301593)



## Table of contents

[Main Lieberman Lab pipeline README](README.md)
* [Snakemake pipeline](readme_files/readme_snake_main.md)
	* [How to run the snakemake pipeline](readme_files/readme_snake_run.md)
	* [Technical details about the snakemake pipeline](readme_files/readme_snake_rules.md)
	* [Wishlist for snakemake pipeline upgrades](readme_files/readme_snake_wishlist.md)
	* [Helpful hints for using the command line](readme_files/readme_snake_basics.md)
* [Local analysis](readme_files/readme_local_main.md)
	* [How to run the local analysis script](readme_files/readme_local_run.md)
	* [Wishlist for local analysis upgrades](readme_files/readme_local_wishlist.md)
	* [Python best practices](readme_files/readme_local_best.md)
