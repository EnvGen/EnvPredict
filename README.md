# EnvPredict

## Required Data Folders

This repository does not include the data folders:
1. [env_data](https://drive.google.com/drive/folders/1TcN7MbLRw7bE0u1Z7SBG0QSNSt92K_3-?usp=drive_link) containing all the data from SHARKweb (physico-chemical data, phytoplankton & zooplankton microscopy counts, bacterioplankton & picoplankton epifluorescence counts + biovolume + carbon sequestration), as well as metadata for sequencing samples with integrated physico-chemical parameters averaged over sampling depths (physical_chemical_processed_translation.tsv).
2. [seq_data](https://drive.google.com/drive/folders/1yWm7nvBzqGdLw2qnjFn0anjrHqwQmtKG?usp=drive_link) containing the count tables, asv sequences, and taxonomic annotation, both raw and after [barnapp](https://github.com/tseemann/barrnap)https://github.com/tseemann/barrnap filtartion.

## Scripts

The scripts can be found in the /code/ subfolder

### Modelling/predictions scripts

envpredict_core.R - creates and runs all the XGBoost and Random Forest models.\
tabpnf_predictions.py - runs the TabPNF predictions

### Evaluation scripts

envpredict_eval.R - compares predictions of physicochemical data based on metabarcoding and microscopy data, as well as predictions of phyto- and zoo-plankton based on different data and approaches (Fig. 3 & 4).\
interannual_comparison.Rmd - compares predictions of physicochemical parameters for the 2015-2017 dataset based on models trained on the 2019-2020 dataset (different dataset), 2015-2017 dataset (same dataset), or both datasets (FIg. 5, Suppplementary Fig. S1 & S2).\
interannual_comparison_linear_regression.Rmd - runs an analysis analogous to interannual_comparison.Rmd, but checks if the actual and predicted values correlate with each other, not if they are the same.

## Running the EQRS analysis with HEAT
1. Unzip HEAT.zip
2. Run HEAT.R
3. Run Plotting_heat.R
