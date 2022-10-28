# How to run the local analysis script


## Outline

This page describes the steps you should take to run the local analysis, which filters the candidate SNVs and generates a parsimony tree. 

The local analysis script uses custom python classes and functions to facilitate your analysis of candidate SNVs. For example, there are classes for a candidate mutation table, basecalls across SNV positions, and a reference genome. All of these classes have attributes and methods, many of which are already used in the analysis script. If you aren't familiar with object-oriented programming, it might be helpful to read a little about what classes, methods, and attributes are. There are also functions that operate on these classes. And, of course, if you want to do something beyond the basics covered in this script, you can write your own new functions too!


## 1. Generate necessary input files

The code needed to run the local analysis can be found in the directory `local_analysis/`:
* `local_analysis/new_snv_script.py` - main script you will run interactively to analyze your data
* `local_analysis/modules/snv_module_recoded.py` - defines the classes and functions used in `new_snv_script.py`

You will also need to supply information on your reference genome:
* `/path/to/reference/genome/genome.fasta` - a FASTA file of the reference genome to which you aligned the reads from your samples; this file must be named `genome.fasta`
* `/path/to/reference/genome/annotations.gff` - a GFF file with annotations of your reference genome, which must be in the same directory as your FASTA file

Finally, you will need to have already generated a candidate mutation table and coverage matrix using the snakemake pipeline:
* `/path/to/data/candidate_mutation_table.pickle.gz` - candidate SNV positions along with statistics about the read alignments from each sample at these positions
* `/path/to/data/cov_raw_sparsecsr_mat.npz` - raw coverage matrix containing information on the number of reads from each sample aligned to each position on the genome


## 2. Set up your python environment and IDE

### Creating a conda environment

If you have not done so already, you should [install conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) on your local machine. Conda is a package management system and environment management system. This helps ensure that the python packages (and their dependencies) that you need to run the local analysis scripts are installed correctly and that the lab's code is reproducible.

Once conda is installed, you can create a python coding environment: ```conda env create -f environment.yaml```.

Then you can activate your new environment (called `pythoncoding`, as defined in the yaml file you used to create the environment): ```conda activate pythoncoding```.

### Launching Spyder

You can then launch [Spyder](https://www.spyder-ide.org) (a python IDE) from the terminal: ```spyder new_snv_script.py```. Once Spyder launches, check to make sure you in the right working directory (the working directory is displayed at top right of the window). You should be in the same directory as your main python script. 

You can run the script cell-by-cell by pressing ctrl+enter. 


## 3. Set paths to your reference genome and your data

At the top of the script, you can set the paths to your reference genome and data by changing the variables called `dir_ref_genome_data` and `dir_cmt_data` respectively. 


## 4. Set SNV filters

When you get to the filtering section, you will need to adjust the filters to be appropriate for your dataset. SNV filtering is an iterative process, and you should experiment with looser/stricter filters and with different kinds of filters. 

### Types of filters

This script already includes a few different kinds of filters that target:
* low quality samples;
* low quality basecalls;
* low quality positions on the genome (across samples);
* recombination regions.

### Evaluating filter thresholds

Histograms are a great way to check if your filter thresholds are justified. Make lots of them!

### Evaluating SNV quality

The function ```plot_interactive_scatter_barplots``` makes an interactive plot, where each dot represents a SNV that passed your filters. Clicking one of the dots will generate a bar chart, where you'll be able to see all basecalls on forward and reverse reads across all of your samples.

The purpose of the interactive plots is to 

Signatures of high-quality SNVs include:
* high read coverage over the SNV position
* consistent basecalls in forward and reverse reads
* consitent number of reads aligned in the forward and reverse directions

Signs of questionable SNVs include:
* impure alleles (alignment issues or contamination issues)
* differences in the basecalls in forward vs reverse reads (alignment issues)
* differences in the number of forward vs reverse reads (alignment issues)
* lots of SNVs close together on the genome (likely recombination)
* and more!


## 5. Make a tree

The last section of the script prepares data for a parsimony tree. Please note that you will need to have the files ```dnapars.app``` and ```dnapars``` in your working directory ([dnapars](https://evolution.genetics.washington.edu/phylip/doc/dnapars.html) is a parsimony program).

The output is a Newick-formatted tree (```.tree``` file)

### Rooting your tree

To see your tree, you can open the file in a program like [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) or [ITOL](https://itol.embl.de).

If you have outgroup samples, then you can root the tree on those. If you do not have outgroup samples, then your tree will be unrooted.


## Table of contents

[Main Lieberman Lab pipeline README](../README.md)
* [Snakemake pipeline](readme_snake_main.md)
	* [How to run the snakemake pipeline](readme_snake_run.md)
	* [Technical information on the snakemake pipeline](readme_snake_rules.md)
	* [Helpful hints for using the command line](readme_snake_basics.md)
	* [Wishlist for snakemake pipeline upgrades](readme_snake_wishlist.md)
* [Local analysis](readme_local_main.md)
	* [How to run the local analysis script](readme_local_run.md)
	* [Wishlist for local analysis upgrades](readme_local_wishlist.md)