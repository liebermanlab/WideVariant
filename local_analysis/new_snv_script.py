# -*- coding: utf-8 -*-
"""

SUMMARY:

    This script demonstrates how to go from a candidate mutation table and a 
    reference genome to a high-quality SNV table and a parsimony tree.

    This script assumes that your input is whole genome sequencing data from
    isolates that has been processed through the Lieberman Lab's standard 
    Snakemake pipeline.
    
    This script is intended to be run interactively rather than as a single 
    function to enable dynamic adjustment of filters to meet the particulars of 
    your sequencing runs, including sample depth and amount of 
    cross-contamination. In our experience, universal filters do not exist, so 
    manual inspection of filtering is essential to SNV calling.


WARNINGS!

    1. The SNV filters illustrated in this script may not be appropriate for 
    your dataset. It is critical that you use the built-in tools for manual 
    inspection of SNV quality. 
    
    2. This script uses the new version of the Lieberman Lab's SNV module, 
    which as of October 2022 is NOT compatible with old versions. If you add 
    functionality from previous versions, you must ensure that data structure
    and indexing are updated appropriately. See documentation for more 
    information on the standards implemented here.


VERSION HISTORY:

    YEAR.MONTH; Name: Add new revisions here!
    
    2022.10; Arolyn: Major overhaul of main script and python module. Not 
    compatible with previous versions. Introduced classes/methods for candidate 
    mutation table data, basecalls, and reference genomes. Implemented 
    consistency in which python packages are used (for example, all numerical
    arrays are now numpy arrays and heterogeneous arrays are now pandas 
    dataframes) and in indexing of genomes and nucleotides. Added many 
    functions for data visualization and quality control. 
    
    2022.04; Tami, Delphine, Laura: Lieberman Lab Hackathon
    
    Additional notes: This module is based on a previous MATLAB version.


@author: Lieberman Lab at MIT. Authors include: Tami Lieberman, Idan Yelin, 
Felix Key, Arolyn Conwill, A. Delphine Tripp, Evan Qu, Laura Markey

"""


#%%#####################
## SET UP ENVIRONMENT ##
########################

# Import python packages
import sys
import os
import numpy as np
import pandas as pd
from scipy import stats
import datetime, time

# Import Lieberman Lab SNV-calling python package
dir_py_scripts = "/Users/arolyn/Dropbox (MIT)/Postdoc/Research/CLABSI/analysis_v2/modules" 
sys.path.insert(0, dir_py_scripts)
import snv_module_recoded as snv # SNV calling module

# Get timestamp
ts = time.time() 
timestamp = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S') # used in output files



#%%#############################
## DATA IMPORT AND PROCESSING ##
################################


#%% Dataset info

dataset_name = 'Escherichia_coli_P-11'
data_file_cmt = 'data/group_1_candidate_mutation_table.pickle.gz'
data_file_cov = 'data/group_1_coverage_matrix_raw.pickle.gz'
dir_ref_genome = 'ref_genomes/Escherichia_coli_P-11-001'
ref_genome_name = 'Escherichia_coli_P-11-001'
samples_to_exclude = [] # option to exclude specific samples manually

# Make subdirectory for this dataset
dir_output = 'output'
os.system( "mkdir -p " + dir_output );



#%% Generate candidate mutation table object

# Import candidate mutation table data generated in Snakemake

# Use this version for updated candidate mutation table matrices
# [quals,p,counts,in_outgroup,sampleNames,indel_counter] = \
#     snv.read_candidate_mutation_table_pickle_gzip(cmtFile) 


# Use this version for old candidate mutation table matrices
[quals,p,counts,in_outgroup,sample_names,indel_counter] = \
    snv.read_old_candidate_mutation_table_pickle_gzip( data_file_cmt ) 


# Create instance of candidate mutation table class
my_cmt = snv.cmt_data_object( sample_names, 
                             in_outgroup, 
                             p, 
                             counts, 
                             quals, 
                             indel_counter 
                             )


#%% Import reference genome information

# Create instance of reference genome class
my_rg = snv.reference_genome_object( dir_ref_genome )
my_rg_annot = my_rg.annotations
my_rg_annot_0 = my_rg_annot[0]



#%% Process raw coverage matrix

# Create instance of simple coverage class (just a summary of coverage matrix data to save memory)
# my_cov = snv.cov_data_object_simple( snv.read_cov_mat_gzip( covFile ), 
#                                     my_cmt.sample_names, 
#                                     my_rg.genome_length, 
#                                     my_rg.contig_starts, 
#                                     my_rg.contig_names 
#                                     )

# Create instance of full coverage matrix class (includes full coverage matrix as attribute)
my_cov = snv.cov_data_object( snv.read_cov_mat_gzip( data_file_cov ), \
                             my_cmt.sample_names, \
                             my_rg.genome_length, \
                             my_rg.contig_starts, \
                             my_rg.contig_names \
                             )



#%% Exclude any samples listed above as bad

samples_to_exclude_bool = np.array( [x in samples_to_exclude for x in my_cmt.sample_names] )

my_cmt.filter_samples( ~samples_to_exclude_bool )
my_cov.filter_samples( ~samples_to_exclude_bool )



#%%###################
## FILTER BASECALLS ##
######################


#%% FILTER BASECALLS

# Create instance of basecalls class for initial calls
my_calls = snv.calls_object( my_cmt )


#%% Filter parameters

# Remove samples that are not high quality
filter_parameter_sample_across_sites = {
                                        'min_average_coverage_to_include_sample': 10, # remove samples that have low coverage
                                        'max_frac_Ns_to_include_sample': 0.2 # remove samples that have too many undefined base (ie. N).
                                        }

# Remove sites within samples that are not high quality
filter_parameter_site_per_sample = {
                                    'min_major_nt_freq_for_call' : 0.85, # on individual samples, a calls' major allele must have at least this freq
                                    'min_cov_per_strand_for_call' : 5,  # on individual samples, calls must have at least this many reads on the fwd/rev strands individually
                                    'min_qual_for_call' : 30, # on individual samples, calls must have this minimum quality score
                                    'max_frac_reads_supporting_indel' : 0.33 # on individual samples, no more than this fraction of reads can support an indel at any given position
                                    }

# Remove sites across samples that are not high quality
filter_parameter_site_across_samples = {
                                        'max_fraction_ambigious_samples' : 0.25, # across samples per position, the fraction of samples that can have undefined bases
                                        'min_median_coverage_position' : 5, # across samples per position, the median coverage
                                        'max_mean_copynum' : 4, # mean copy number at a positions across all samples
                                        'max_max_copynum' : 7 # max maximum copynumber that a site can have across all samples
                                        }


#%% Filter samples based on coverage

# Identify samples with low coverage and make a histogram
[ low_cov_samples, goodsamples_coverage ] = snv.filter_samples_by_coverage( 
    my_cov.get_median_cov_of_chromosome(), 
    filter_parameter_sample_across_sites['min_average_coverage_to_include_sample'], 
    my_cov.sample_names, 
    True, 
    dir_output 
    )

# Filter candidate mutation table, coverage, and calls objects
my_cmt.filter_samples( goodsamples_coverage )
my_cov.filter_samples( goodsamples_coverage )
my_calls.filter_samples( goodsamples_coverage )


#%% Filter calls per position per sample

# Filter based on coverage

my_calls.filter_calls_by_element( 
    my_cmt.fwd_cov < filter_parameter_site_per_sample['min_cov_per_strand_for_call'] 
    ) # forward strand coverage too low

my_calls.filter_calls_by_element( 
    my_cmt.rev_cov < filter_parameter_site_per_sample['min_cov_per_strand_for_call'] 
    ) # reverse strand coverage too low

# Filter based on quality

my_calls.filter_calls_by_element( 
    my_cmt.quals < filter_parameter_site_per_sample['min_qual_for_call'] 
    ) # quality too low

# Filter based on major allele frequency

my_calls.filter_calls_by_element( 
    my_cmt.major_nt_freq < filter_parameter_site_per_sample['min_major_nt_freq_for_call'] 
    ) # major allele frequency too low


#%% Filter positions with indels

with np.errstate(divide='ignore',invalid='ignore'):
    # compute number of reads supporting an indel
    frac_reads_supporting_indel = np.sum(my_cmt.indel_stats,axis=2)/my_cmt.coverage # sum reads supporting insertion plus reads supporting deletion
    frac_reads_supporting_indel[ ~np.isfinite(frac_reads_supporting_indel) ] = 0
    # note: this fraction can be above zero beacuse the number of reads supporting an indel includes a +/-3 bp window around a given position on the genome
my_calls.filter_calls_by_element( 
    frac_reads_supporting_indel > filter_parameter_site_per_sample['max_frac_reads_supporting_indel'] 
    ) # too many reads supporting indels


#%% Filter positions that look iffy across samples

my_calls.filter_calls_by_position( 
    my_calls.get_frac_Ns_by_position() > filter_parameter_site_across_samples['max_fraction_ambigious_samples'] 
    ) # too many samples with ambiuguous calls at this position
snv.filter_histogram( 
    my_calls.get_frac_Ns_by_position(), 
    filter_parameter_site_across_samples['max_fraction_ambigious_samples'], 
    'Fraction Ns by position'
    )

my_calls.filter_calls_by_position( 
    np.median( my_cmt.coverage, axis=0 ) < filter_parameter_site_across_samples['min_median_coverage_position'] 
    ) # insufficient median coverage across samples at this position
# # Optional code to make a histogram
# snv.filter_histogram( 
#     np.median( my_cmt.coverage, axis=0 ), 
#     filter_parameter_site_across_samples['min_median_coverage_position'], 
#     'Median coverage by position'
#     )

copy_number_per_sample_per_pos = my_cmt.coverage / np.expand_dims( my_cov.get_median_cov_of_chromosome(), 1) # compute copy number
copy_number_avg_per_pos = np.mean( copy_number_per_sample_per_pos, axis=0 ) # mean copy number at each position
my_calls.filter_calls_by_position( 
    copy_number_avg_per_pos > filter_parameter_site_across_samples['max_mean_copynum'] 
    ) # average copy number too high
copy_number_max_per_pos = np.max( copy_number_per_sample_per_pos, axis=0 ) # mean copy number at each position
my_calls.filter_calls_by_position( 
    copy_number_max_per_pos > filter_parameter_site_across_samples['max_max_copynum'] 
    ) # average copy number too high
# # Optional code to make a histogram
# snv.filter_histogram( 
#     copy_number_avg_per_pos, 
#     filter_parameter_site_across_samples['max_mean_copynum'] , 
#     'Average copy number by position'
#     )


#%% Filter samples that have too many ambiguous calls

# Identify samples with many ambiguous basecalls and make a histogram
pos_to_consider = my_calls.p[ np.any(  my_calls.calls, axis=0 ) ] # mask positions with no basecalls in any samples
[ samples_with_toomanyNs, goodsamples_nonambig ] = snv.filter_samples_by_ambiguous_basecalls( 
    my_calls.get_frac_Ns_by_sample( pos_to_consider ), 
    filter_parameter_sample_across_sites['max_frac_Ns_to_include_sample'], 
    my_calls.sample_names, 
    my_calls.in_outgroup, # does not filter outgroup samples!!!
    True, 
    dir_output 
    )

# Filter candidate mutation table, coverage, and calls objects
my_cmt.filter_samples( goodsamples_nonambig )
my_cov.filter_samples( goodsamples_nonambig )
my_calls.filter_samples( goodsamples_nonambig )



#%%##########################
## INFER ANCESTRAL ALLELES ##
#############################


#%% Part 1: Get ancestral nucleotide from outgroup

# Ancestral alleles should be inferred from an outgroup.

# Filtered calls for outgroup samples only
calls_outgroup = my_calls.get_calls_in_outgroup_only()
# Switch N's (0's) to NaNs
calls_outgroup_N_as_NaN = calls_outgroup.astype('float') # init ()
calls_outgroup_N_as_NaN[ calls_outgroup_N_as_NaN==0 ] = np.nan

# Infer ancestral allele as the most common allele among outgroup samples (could be N)
calls_ancestral = np.zeros( my_calls.num_pos, dtype='int') # init as N's
outgroup_pos_with_calls = np.any(calls_outgroup,axis=0) # positions where the outgroup has calls
calls_ancestral[outgroup_pos_with_calls] = stats.mode( calls_outgroup_N_as_NaN[:,outgroup_pos_with_calls], axis=0, nan_policy='omit' ).mode.squeeze()

# Report number of ancestral alleles inferred from outgroup
print('Number of candidate SNVs with outgroup alleles: ' + str(sum(outgroup_pos_with_calls)) + '.')
print('Number of candidate SNVs missing outgroup alleles: ' + str(sum(calls_ancestral==0)) + '.')


#%% Part 2: Fill in any missing data with nucleotide from reference

# WARNING! Rely on this method with caution especially if the reference genome 
# was derived from one of your ingroup samples.

# # Pull alleles from reference genome across p
# calls_reference = my_rg.get_ref_NTs_as_ints( my_rg.p2contigpos(p) )

# # Update ancestral alleles
# pos_to_update = ( calls_ancestral==0 )
# calls_ancestral[ pos_to_update ] = calls_reference[ pos_to_update ]



#%%#########################
## IDENTIFY SNV POSITIONS ##
############################


#%% Compute mutation quality

# Grab filtered calls from ingroup samples only
calls_ingroup = my_calls.get_calls_in_sample_subset( np.logical_not( my_calls.in_outgroup ) )
quals_ingroup = my_cmt.quals[ np.logical_not( my_calls.in_outgroup ),: ]
num_samples_ingroup = sum( np.logical_not( my_calls.in_outgroup ) )
# Note: Here we are only looking for SNV differences among ingroup samples. If
# you also want to find SNV differences between the ingroup and the outgroup
# samples (eg mutations that have fixed across the ingroup), then you need to
# use calls and quals matrices that include outgroup samples.

# Compute quality
[ mut_qual, mut_qual_samples ] = snv.compute_mutation_quality( calls_ingroup, quals_ingroup ) 
# note: returns NaN if there is only one type of non-N call


#%% Identify suspected recombination positions

# Filter params
filter_parameter_recombination = {
                                    'distance_for_nonsnp' : 1000, # region in bp on either side of goodpos that is considered for recombination
                                    'corr_threshold_recombination' : 0.75 # minimum threshold for correlation
                                    }

# Find SNVs that are are likely in recombinant regions
[ p_recombo, recombo_bool ] = snv.find_recombination_positions( \
    my_calls, my_cmt, calls_ancestral, mut_qual, my_rg, \
    filter_parameter_recombination['distance_for_nonsnp'], \
    filter_parameter_recombination['corr_threshold_recombination'], \
    True, dir_output \
    )
    
# Save positions with likely recombination
if len(p_recombo)>0:
    with open( dir_output + '/snvs_from_recombo.csv', 'w') as f:
        for p in p_recombo:
            f.write(str(p)+'\n')



#%% Determine which positions have high-quality SNVs

# Filters
filter_SNVs_not_N = ( calls_ingroup != snv.nts2ints('N') ) # mutations must have a basecall (not N)
filter_SNVs_not_ancestral_allele = ( calls_ingroup != np.tile( calls_ancestral, (num_samples_ingroup,1) ) ) # mutations must differ from the ancestral allele
filter_SNVs_quals_not_NaN = ( np.tile( mut_qual, (num_samples_ingroup,1) ) >= 1) # alleles must have strong support 
filter_SNVs_not_recombo = np.tile( np.logical_not(recombo_bool), (num_samples_ingroup,1) ) # mutations must not be due to suspected recombination

# Fixed mutations per sample per position
fixedmutation = \
    filter_SNVs_not_N \
    & filter_SNVs_not_ancestral_allele \
    & filter_SNVs_quals_not_NaN \
    & filter_SNVs_not_recombo # boolean

hasmutation = np.any( fixedmutation, axis=0) # boolean over positions (true if at lest one sample has a mutation at this position)

goodpos_bool = np.any( fixedmutation, axis=0 )
goodpos_idx = np.where( goodpos_bool )[0]
num_goodpos = len(goodpos_idx)

print('Num mutations: '+str(num_goodpos))



#%% Make and annotate a SNV table

# Prepare data for SNV table

my_cmt_goodpos = my_cmt.copy()
my_cmt_goodpos.filter_positions( goodpos_bool )
my_cmt_goodpos_ingroup = my_cmt_goodpos.copy()
my_cmt_goodpos_ingroup.filter_samples( np.logical_not( my_cmt_goodpos_ingroup.in_outgroup ) )

my_calls_goodpos = my_calls.copy()
my_calls_goodpos.filter_positions( goodpos_bool )
calls_goodpos = my_calls_goodpos.calls
calls_goodpos_ingroup = calls_goodpos[ np.logical_not( my_calls_goodpos.in_outgroup ),: ]

p_goodpos = my_calls_goodpos.p

calls_ancestral_goodpos = calls_ancestral[ goodpos_bool ]

# Generate SNV table

# Parameters
promotersize = 250; # how far upstream of the nearest gene to annotate something a promoter mutation (not used if no annotation)

# Make a table (pandas dataframe) of SNV positions and relevant annotations
mutations_annotated = snv.annotate_mutations( \
    my_rg, \
    p_goodpos, \
    np.tile( calls_ancestral[goodpos_idx], (my_cmt_goodpos_ingroup.num_samples,1) ), \
    calls_goodpos_ingroup, \
    my_cmt_goodpos_ingroup, \
    fixedmutation[:,goodpos_idx], \
    mut_qual[:,goodpos_bool].flatten(), \
    promotersize \
    ) 


#%% SNV quality control plots

# Note: These data visualizations are intended to help you evaluate if your SNV
# filters are appropriate. Do not proceed to the next step until you are 
# convinced that your filters are strict enough to filter low-quality SNVs, but
# not so strict that good SNVs are eliminated as well. 

# Clickable bar charts for each SNV position
snv.plot_interactive_scatter_barplots( \
    p_goodpos, \
    mut_qual[0,goodpos_idx], \
    'pos', \
    'qual', \
    my_cmt_goodpos.sample_names, \
    mutations_annotated, \
    my_cmt_goodpos.counts )

# Heatmaps of basecalls, coverage, and quality over SNV positions
if num_goodpos>0:
    snv.make_calls_qc_heatmaps( my_cmt_goodpos, my_calls_goodpos, True, dir_output )



#%%###########################
## PARSIMONY AND TREEMAKING ##
##############################

#%% Filter calls for tree

# Note: Here we are using looser filters than before

# Choose subset of samples or positions to use in the tree by idx
samplestoplot = np.arange(my_cmt_goodpos.num_samples) # default is to use all samples 
goodpos4tree = np.arange(num_goodpos) # default is to use all positions 

# Get calls for tree
my_calls_tree = snv.calls_object( my_cmt_goodpos ) # re-initialize calls

# Apply looser filters than befor (want as many alleles as possible)
filter_parameter_calls_for_tree = {
                                    'min_cov_for_call' : 5, # on individual samples, calls must have at least this many fwd+rev reads
                                    'min_qual_for_call' : 30, # on individual samples, calls must have this minimum quality score
                                    'min_major_nt_freq_for_call' : 0.8,  # on individual samples, a call's major allele must have at least this freq
                                    }

my_calls_tree.filter_calls_by_element( 
    my_cmt_goodpos.coverage < filter_parameter_calls_for_tree['min_cov_for_call'] 
    ) # forward strand coverage too low

my_calls_tree.filter_calls_by_element( 
    my_cmt_goodpos.quals < filter_parameter_calls_for_tree['min_qual_for_call'] 
    ) # quality too low

my_calls_tree.filter_calls_by_element( 
    my_cmt_goodpos.major_nt_freq < filter_parameter_calls_for_tree['min_major_nt_freq_for_call'] 
    ) # major allele frequency too low

# Make QC plots again using calls for the tree
if num_goodpos>0:
    snv.make_calls_qc_heatmaps( my_cmt_goodpos, my_calls_tree, False, dir_output ) 



#%% Make a tree

if num_goodpos>0:

    # Get nucleotides of basecalls (ints to NTs)
    calls_for_treei = my_calls_tree.calls[ np.ix_(samplestoplot,goodpos4tree) ] # numpy broadcasting of row_array and col_array requires np.ix_()
    calls_for_tree = snv.ints2nts(calls_for_treei) # NATCG translation
    
    # Sample names for tree
    treesampleNamesLong = my_cmt_goodpos.sample_names
    for i, samplename in enumerate(treesampleNamesLong):
        if not samplename[0].isalpha():
            treesampleNamesLong[i]='S'+treesampleNamesLong[i] # sample names are modified to make parsing easier downstream
    sampleNamesDnapars = [ "{:010d}".format(i) for i in range(my_cmt_goodpos.num_samples)]
    
    # Add inferred ancestor and reference 
    calls_ancestral_for_tree =  np.expand_dims( snv.ints2nts(calls_ancestral_goodpos), axis=0 )
    calls_reference_for_tree = np.expand_dims( my_rg.get_ref_NTs( my_calls_tree.p ), axis=0 )
    calls_for_tree_all = np.concatenate( (calls_ancestral_for_tree,calls_reference_for_tree,calls_for_tree) ,axis=0) # first column now outgroup_nts; outgroup_nts[:, None] to make ndims (2) same for both
    sampleNamesDnapars_all = np.append([ 'Sanc','Sref' ],sampleNamesDnapars)
    treesampleNamesLong_all = np.append([ 'inferred_ancestor','reference_genome' ],treesampleNamesLong) 
    
    # Build tree
    snv.generate_tree( \
        calls_for_tree_all.transpose(), \
        treesampleNamesLong_all, \
        sampleNamesDnapars_all, \
        ref_genome_name, \
        dir_output, \
        "snv_tree_"+ref_genome_name, \
        buildTree='PS' \
        )



#%%###################
## SAVE DATA TABLES ##
######################


#%% Write SNV table to a tsv file

# Note: important to use calls for tree (rather than calls for finding SNVs)

# Make table
if num_goodpos>0:
    output_tsv_filename = dir_output + '/' + 'snv_table_mutations_annotations.tsv'
    snv.write_mutation_table_as_tsv( \
        p_goodpos, \
        mut_qual[0,goodpos_idx], \
        my_cmt_goodpos.sample_names, \
        mutations_annotated, \
        calls_for_tree_all, \
        treesampleNamesLong, \
        output_tsv_filename \
        )

    
#%% Write table of tree distances to a csv file

if num_goodpos>0:
    
    # Compute number of SNVs separating each sample from the inferred ancestor
    fixedmutation_tree = ( calls_for_tree != calls_ancestral_for_tree ) & ( calls_for_tree != 'N' ) # boolean
    
    dist_to_anc_by_sample = np.sum( fixedmutation_tree, axis=1 )
    # Save to a file
    with open( dir_output + '/snv_table_tree_distances.csv', 'w') as f:
        f.write('sample_name,num_SNVs_to_ancestor\n')
        for i,name in enumerate(treesampleNamesLong):
            f.write(name + ',' + str(dist_to_anc_by_sample[i]) + '\n' )


#%% Write a table of summary stats

print( 'Number of samples in ingroup: ' + str(num_samples_ingroup) + '.') 
print( 'Number of good SNVs found: ' + str(num_goodpos) + '.') 
print('Number of good SNVs with outgroup alleles: ' + str(sum(calls_ancestral_goodpos>0)) + '.')
print('Number of good SNVs missing outgroup alleles: ' + str(sum(calls_ancestral_goodpos==0)) + '.')
if num_goodpos>0:
    dist_to_anc_by_sample_ingroup = dist_to_anc_by_sample[ np.logical_not(my_calls_tree.in_outgroup) ]
    dmrca_median = np.median( dist_to_anc_by_sample_ingroup )
    dmrca_min = np.min( dist_to_anc_by_sample_ingroup )
    dmrca_max = np.max( dist_to_anc_by_sample_ingroup )
else:
    dmrca_median = 0
    dmrca_min = 0
    dmrca_max = 0
print( 'Median dMRCA: ' + str(dmrca_median) + '.') 
print( 'Min dMRCA: ' + str(dmrca_min) + '.') 
print( 'Max dMRCA: ' + str(dmrca_max) + '.') 

with open( dir_output + '/snv_table_simple_stats.csv', 'w') as f:
    f.write('dataset,genome,num_samples_ingroup,num_snvs,num_snvs_with_outgroup_allele,dmrca_median,dmrca_min,dmrca_max,\n')
    f.write(dataset_name+','+ref_genome_name+','+str(num_samples_ingroup)+','+str(num_goodpos)+','+str(sum(calls_ancestral_goodpos>0))+','+str(dmrca_median)+','+str(dmrca_min)+','+str(dmrca_max)+',\n')



#%%##################
## CONTIG ANALYSIS ##
#####################

# This section examines contig coverage.


#%% Make a plot of presence/absence of each contig

my_cov.plot_heatmap_contig_copy_num( dir_output )

# Record data in files
snv.write_generic_csv( my_cov.median_coverage_by_contig, my_cov.contig_names, my_cov.sample_names, dir_output+'/contig_table_median_coverage.csv' )
snv.write_generic_csv( my_cov.copy_number_by_contig, my_cov.contig_names, my_cov.sample_names, dir_output+'/contig_table_copy_number.csv' )


#%% Make plots of coverage traces of each contig

if type(my_cov) == snv.cov_data_object:
    for contig_num in np.linspace(1,my_cov.num_contigs,my_cov.num_contigs).astype(int):
        print('Generating copy number traces for contig ' + str(contig_num) + '...')
        my_cov.make_coverage_trace(contig_num,100,dir_output);
elif type(my_cov) == snv.cov_data_object_simple:
    print('Coverage matrix object type is cov_data_obj_simple, not cov_data_obj. Raw coverage matrix not available for copy number traces.')

