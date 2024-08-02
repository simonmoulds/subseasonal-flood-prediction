#!/usr/bin/env bash
#SBATCH --partition=General
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=smoulds@ed.ac.uk
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=8000
#SBATCH --job-name=seasonal-flood-prediction
#SBATCH --time=12:00:00

source /opt/conda/etc/profile.d/conda.sh
conda activate snakemake

# # PYTHONPATH sometimes causes issues
# export PYTHONPATH=

# # Make sure R scripts can access local libraries
# export R_LIBS_USER=$HOME/local/rlibs

# Increase the limit on number of files
ulimit -n 12288

# Unset R_LIBS env var
unset R_LIBS

# Unlock, just in case
snakemake \
    --quiet \
    --snakefile workflow/Snakefile \
    --cores 1 \
    --configfile config/config.yml \
    --use-conda \
    --conda-frontend conda \
    --unlock

# Now run the workflow
snakemake \
    --quiet \
    --snakefile workflow/Snakefile_envs \
    --cores 1 \
    --configfile config/config.yml \
    --use-conda \
    --conda-frontend conda \
    --rerun-incomplete

snakemake \
    --quiet \
    --snakefile workflow/Snakefile \
    --cores 24 \
    --configfile config/config.yml \
    --use-conda \
    --conda-frontend conda \
    --rerun-incomplete

