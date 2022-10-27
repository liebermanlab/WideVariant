#!/bin/bash
#SBATCH --job-name gus.main
#SBATCH -p defq,sched_mem1TB,sched_mem1TB_centos7
#SBATCH -n 1
#SBATCH --time=1-00:00:00
#SBATCH --mem=8GB
#SBATCH -o mainout.txt
#SBATCH -e mainout.txt
#SBATCH --mail-user=YOUR_EMAIL_HERE
#SBATCH --mail-type=ALL
#SBATCH --exclude=node327

# Activate conda environment (may need to change name of env)
source activate snakemake

# Run snakemake pipeline
bash snakemakeslurm.sh

# Print "Done!!!" at end of main log file
echo Done!!!
