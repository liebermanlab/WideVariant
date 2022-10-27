#!/bin/bash


### Launch snakemake to run jobs via SLURM


# ### Check for clean
#
# # https://snakemake.readthedocs.io/en/stable/project_info/faq.html#how-do-i-remove-all-files-created-by-snakemake-i-e-like-make-clean
# if [ "$1" == "clean" ]; then
#     echo 'rm $(snakemake --summary | tail -n+2 | cut -f1)'
#     snakemake --summary | tail -n+2 | cut -f1
#     rm -f $(snakemake --summary | tail -n+2 | cut -f1)
#     exit 0
# fi


### Information for the SLURM scheduler (parameters pulled from json)

SM_PARAMS="job-name ntasks partition exclude time mail-user mail-type error output"
SM_ARGS=" --parsable --cpus-per-task {cluster.cpus-per-task} --mem-per-cpu {cluster.mem-per-cpu-mb}"

for P in ${SM_PARAMS}; do SM_ARGS="$SM_ARGS --$P {cluster.$P}"; done
echo "SM_ARGS: ${SM_ARGS}"


### Make subdirectory for logs files
mkdir -p logs


### Run snakemake pipeline

# -p : print to standard output
# -j : total number of jobs executed in parallel on cluster
# --restart-times : how many times snakemake should re-batch a failed job
# -n : for a dryrun (will not actually batch jobs)
# --unlock : to unlock the directory

snakemake -p \
    $* \
    --latency-wait 10 \
    -j 1000 \
    --cluster-config $(dirname $0)/cluster.slurm.json \
    --cluster "sbatch $SM_ARGS" \
    --cluster-status /scratch/mit_lieberman/scripts/slurm_status.py \
    --reason \
    --use-conda \
    --rerun-incomplete \
    --restart-times 1\ 
    -s Snakefile.py \
    --keep-going \
    --conda-prefix /scratch/mit_lieberman/tools/conda_snakemake/ \
    #-n 
    # --unlock \
