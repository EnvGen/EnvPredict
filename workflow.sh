#!/bin/bash -l

unzip data/seq_files.zip
unzip data/env_files.zip

snakemake --rerun-incomplete --cores all --unlock
snakemake --dag --cores all | dot -Tpng > dag.png
snakemake --rulegraph --cores all | dot -Tpng > rulegraph.png
snakemake --forceall --cores all --latency-wait 60 --verbose
# snakemake --rerun-incomplete --cores all --latency-wait 60 --verbose

# sbatch --j=envpredict --time=24:00:00 -A naiss2025-5-219 -p shared --nodes=1 --ntasks-per-node=128 --mail-type=ALL --mail-user=krzysztof.jurdzins@scilifelab.se workflow.sh