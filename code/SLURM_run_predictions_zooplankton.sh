#!/bin/bash -l

#SBATCH -A naiss2024-5-183
#SBATCH -p node
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=krzysztof.jurdzins@scilifelab.se
#SBATCH -t 10-00:00:00
#SBATCH -C mem1TB
#SBATCH -J EnvPredict-zooplankton

conda activate EnvPredict
Rscript --vanilla RF_XGB_ab.R -b "../env_data/combined/zooplankton_filtered.tsv" -m "" -w "/crex/proj/snic2020-6-126/projects/plankton_monitoring/EnvPredict/code/"