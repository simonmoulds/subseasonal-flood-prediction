#!/usr/bin/env bash
#SBATCH --partition=short
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=simon.moulds@ouce.ox.ac.uk
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem-per-cpu=16000
#SBATCH --job-name=seasonal-flood-prediction
#SBATCH --time=12:00:00

# cd $SCRATCH || exit 1

module load Anaconda3/2022.05
# module load Mamba/4.14.0-0
module load R/4.1.0-foss-2021a
module unload GDAL/3.3.0-foss-2021a
module unload GEOS/3.9.1-GCC-10.3.0
module unload PROJ/8.0.1-GCC-10.3.0

source activate $DATA/envs/snakemake

# PYTHONPATH sometimes causes issues
export PYTHONPATH=

# # Make sure R scripts can access local libraries
# export R_LIBS_USER=$HOME/local/rlibs

# Increase the limit on number of files
ulimit -n 12288

# Unset R_LIBS env var
unset R_LIBS

# Unlock, just in case
snakemake \
    --quiet \
    --snakefile workflow/Snakefile_preprocessing \
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
    --snakefile workflow/Snakefile_preprocessing \
    --cores 24 \
    --configfile config/config.yml \
    --use-conda \
    --conda-frontend conda \
    --rerun-incomplete

