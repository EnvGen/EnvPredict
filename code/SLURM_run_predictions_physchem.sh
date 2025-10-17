#!/bin/bash -l

#SBATCH -A naiss2024-5-183
#SBATCH -p core
#SBATCH -n 20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=krzysztof.jurdzins@scilifelab.se
#SBATCH -t 6-00:00:00
#SBATCH -J EnvPredict


conda activate EnvPredict
Rscript --vanilla RF_XGB_ab.R -w "/crex/proj/snic2020-6-126/projects/plankton_monitoring/EnvPredict/code/"
