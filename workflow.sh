#!/bin/bash -l

#SBATCH -A naiss2025-5-219
#SBATCH -p memory
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=880GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=krzysztof.jurdzins@scilifelab.se
#SBATCH -t 4-00:00:00
#SBATCH -J EnvPredict

# unzip data/seq_files.zip
# unzip data/env_files.zip

snakemake --rerun-incomplete --cores all --unlock
snakemake --dag --cores all | dot -Tpng > dag.png
snakemake --rulegraph --cores all | dot -Tpng > rulegraph.png
snakemake --rerun-incomplete --cores all --latency-wait 60 --verbose