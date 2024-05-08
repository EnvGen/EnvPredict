#!/bin/bash -l
 
#SBATCH -A naiss2024-5-183
#SBATCH -p core
#SBATCH -n 20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=krzysztof.jurdzins@scilifelab.se
#SBATCH -t 01:00:00
#SBATCH -J creating_EnvPredict_env

conda env create -f environment.yml
conda activate EnvPredict
Rscript --vanilla Install_Rpackages.R 20
