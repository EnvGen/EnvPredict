#!/bin/bash -l

#SBATCH -A naiss2024-5-183
#SBATCH -p node
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=krzysztof.jurdzins@scilifelab.se
#SBATCH -t 4-00:00:00
#SBATCH -C mem1TB
#SBATCH -J EnvPredict-zooplankton


# module add R/4.1.1
# module add R_packages/4.1.1
conda activate EnvPredict
# Rscript --vanilla filter_data.R -w "/crex/proj/snic2020-6-126/projects/plankton_monitoring/EnvPredict/code/"
# Rscript --vanilla RF_XGB_ab.R -w "/crex/proj/snic2020-6-126/projects/plankton_monitoring/EnvPredict/code/"
Rscript --vanilla RF_XGB_ab.R -b "../env_data/combined/zooplankton_filtered.tsv" -w "/crex/proj/snic2020-6-126/projects/plankton_monitoring/EnvPredict/code/"
# Rscript --vanilla RF_XGB_ab.R -b "../env_data/combined/phytoplankton_filtered.tsv" -w "/crex/proj/snic2020-6-126/projects/plankton_monitoring/EnvPredict/code/"
# Rscript --vanilla RF_XGB_ab.R -b "../env_data/combined/picoplankton_filtered.tsv" -w "/crex/proj/snic2020-6-126/projects/plankton_monitoring/EnvPredict/code/"
# Rscript --vanilla RF_XGB_ab.R -b "../env_data/combined/bacterioplankton_filtered.tsv" -w "/crex/proj/snic2020-6-126/projects/plankton_monitoring/EnvPredict/code/"
