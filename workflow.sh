#!/bin/bash -l

unzip data/seq_files.zip
unzip data/env_files.zip

snakemake --rerun-incomplete --cores all --unlock
snakemake --dag --cores all | dot -Tpng > dag.png
snakemake --rulegraph --cores all | dot -Tpng > rulegraph.png
snakemake --rerun-incomplete --cores all --latency-wait 60 --verbose