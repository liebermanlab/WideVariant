# Wishlist for snakemake pipeline upgrades


## General

* ...


## Mapping step

* Fix outstanding issues with de-duplicating reads
* ...


## Case step

* Harmonize data objects coming out of GUS to match the input data objects for the local python scripts
	* Switch around axes of data arrays such that the 0th axis represents samples and the 1st axis represents position on the genome
	* Make sure all arrays are numpy arrays
* Compute and save a simplified coverage matrix that includes median coverage per sample over each contig
	* Code for computing this already exists in the local analysis script
	* Inputs to `build_candidate_mutation_table.py` will need to be udpated to include the reference genome (so that it is possible to know where the contig boundaries are)
* Make it possible to align a group of samples to multiple reference genomes
* ...


## Assembly step

* Start using bakta instead of prokka for annotations
* Add ability to make lineage co-assemblies with a specified number of reads per sample
* ...


## Bracken step

* ...


## Table of contents

[Main Lieberman Lab pipeline README](../README.md)
* [Snakemake pipeline](readme_snake_main.md)
	* [How to run the snakemake pipeline](readme_snake_run.md)
	* [Technical details about the snakemake pipeline](readme_snake_rules.md)
	* [Wishlist for snakemake pipeline upgrades](readme_snake_wishlist.md)
	* [Helpful hints for using the command line](readme_snake_basics.md)
* [Local analysis](readme_local_main.md)
	* [How to run the local analysis script](readme_local_run.md)
	* [Wishlist for local analysis upgrades](readme_local_wishlist.md)
	* [Python best practices](readme_local_best.md)
