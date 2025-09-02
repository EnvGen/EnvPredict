#!/bin/bash


# This script uses the DM.py provided in DeepMicro: https://github.com/minoh0201/DeepMicro
# Create log directory
mkdir -p logs
echo "==== DeepMicro run started $(date) ====" > master_log_april_1.txt

# Function to run a model and log output
run_model() {
  echo "[$(date)] Starting: $*" >> master_log_april_1.txt

  # Generate log file name from method + file + dims
  model_tag=$(echo "$*" | sed -E 's/.*--(ae|vae|cae).*--dm ([^ ]*).*-cd ([^ ]*).*/\1_\3_dm\2/' | tr '/' '_' | tr ' ' '_')
  log_file="logs/${model_tag}.log"

  $* > "$log_file" 2>&1
  status=$?

  if [ $status -eq 0 ]; then
    echo "[$(date)] ✅ SUCCESS: $model_tag" >> master_log_april_1.txt
  else
    echo "[$(date)] ❌ FAILED: $model_tag (Exit code $status)" >> master_log_april_1.txt
  fi
}

# SAE, 16S
for dm in 32 64 128 256 512; do
  run_model python DM.py -r 1 --no_clf -cd filtered_norm_seqtab_16S_pivoted_march_2025.csv --ae -dm $dm --save_rep
done

# SAE, 18S
for dm in 32 64 128 256 512; do
  run_model python DM.py -r 1 --no_clf -cd filtered_norm_seqtab_18S_pivoted_march_2025.csv --ae -dm $dm --save_rep
done

# VAE, 16S
for dm in 32,4 32,8 32,16 64,4 64,8 64,16 128,4 128,8 128,16 256,4 256,8 256,16 512,4 512,8 512,16; do
  run_model python DM.py -r 1 --no_clf -cd filtered_norm_seqtab_16S_pivoted_march_2025.csv --vae -dm $dm --save_rep
done

# VAE, 18S
for dm in 32,4 32,8 32,16 64,4 64,8 64,16 128,4 128,8 128,16 256,4 256,8 256,16 512,4 512,8 512,16; do
  run_model python DM.py -r 1 --no_clf -cd filtered_norm_seqtab_18S_pivoted_march_2025.csv --vae -dm $dm --save_rep
done

# CAE, 16S
for dm in 64,32 32,16 16,8 8,4 4,2 64,32,16 32,16,8 16,8,4 8,4,2 4,2,1; do
  run_model python DM.py -r 1 --no_clf -cd filtered_norm_seqtab_16S_pivoted_march_2025.csv --cae -dm $dm --save_rep
done

# CAE, 18S
for dm in 64,32 32,16 16,8 8,4 4,2 64,32,16 32,16,8 16,8,4 8,4,2 4,2,1; do
  run_model python DM.py -r 1 --no_clf -cd filtered_norm_seqtab_18S_pivoted_march_2025.csv --cae -dm $dm --save_rep
done

# DAE, 16S
for dm in 128,64,32 256,128,64 512,256,128 1024,512,256 2048,1024,512 64,32 128,64 256,128 512,256 1024,512; do
  run_model python DM.py -r 1 --no_clf -cd filtered_norm_seqtab_16S_pivoted_march_2025.csv --ae -dm $dm --save_rep
done

# DAE, 18S
for dm in 128,64,32 256,128,64 512,256,128 1024,512,256 2048,1024,512 64,32 128,64 256,128 512,256 1024,512; do
  run_model python DM.py -r 1 --no_clf -cd filtered_norm_seqtab_18S_pivoted_march_2025.csv --ae -dm $dm --save_rep
done

echo "==== DeepMicro run completed $(date) ====" >> master_log_april_1.txt
