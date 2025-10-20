# EnvPredict

## Required Data Folders

This repository includes zipped data folders:
1. [env_data](https://drive.google.com/drive/folders/1TcN7MbLRw7bE0u1Z7SBG0QSNSt92K_3-?usp=drive_link) containing all the data from SHARKweb (physico-chemical data, phytoplankton & zooplankton microscopy counts, bacterioplankton & picoplankton epifluorescence counts + biovolume + carbon sequestration), as well as metadata for sequencing samples with integrated physico-chemical parameters averaged over sampling depths (physical_chemical_processed_translation.tsv).
2. [seq_data](https://drive.google.com/drive/folders/1yWm7nvBzqGdLw2qnjFn0anjrHqwQmtKG?usp=drive_link) containing the count tables, asv sequences, and taxonomic annotation, both raw and after [barnapp](https://github.com/tseemann/barrnap)https://github.com/tseemann/barrnap filtartion.


## Scripts

The scripts can be found in the /code/ subfolder

### Modelling/predictions scripts

**envpredict_core.R** - creates and runs all the XGBoost and Random Forest models.\
**tabpnf_predictions.py** - runs the TabPNF predictions. \
Prep_data_for_deep_micro.py - prepares data to be used as input for autoencoders to obtain Deep Representations \
createDeepRepresentations.sh - creatses Deep Representations with different autencoders. \
Prep_data_for_classifiers_after_deep_micro.py - prerpares the Deep Representations for downstream analyses.


### Evaluation scripts

**Plot_Figure_2.Rmd** - compares different prediction algorithms and predicitions made on different types of data (16S vs 18S, different taxonomic levels)\
**envpredict_eval.R** - compares predictions of physicochemical data based on metabarcoding and microscopy data, as well as predictions of phyto- and zoo-plankton based on different data and approaches (Fig. 3 & 4).\
**interannual_comparison.Rmd** - compares predictions of physicochemical parameters for the 2015-2017 dataset based on models trained on the 2019-2020 dataset (different dataset), 2015-2017 dataset (same dataset), or both datasets (FIg. 5, Suppplementary Fig. S1 & S2).\
interannual_comparison_linear_regression.Rmd - runs an analysis analogous to interannual_comparison.Rmd, but checks if the actual and predicted values correlate with each other, not if they are the same.\

## Pipeline

A Snakefile is available in the main directory of the repository, which allows to rerun most of the analysis with snakemake. The scripts that are run by the pipeline are highlighted in the section above, and are visualized in the running order in dag.png and rulegraph.png. We decided not to include the analysis of Deep Representations in the pipeline, since they have been obtained using a GPU and a considerable amount of computation, yielding poor results. We have not included the Ecological Quality Ratios (EQRs) analysis in the pipeline, and required extra data preparation steps.

By adjusting the reading data section of envpredict_core.R, this script can be adjusted to other datasets. Consider skipping the TabPNF part (by commenting out the rule and its output from the Snakefile and the parts that depend on this analysis from Plot_Figure_2.Rmd), since it takes a considerable amount of time to run (~two days on one node, 256 cores, of Dardel, the KTH HPC cluster).

## Running the Ecological Quality Ratios (EQRs) analysis with HEAT

1. Unzip code/HEAT.zip
2. Run HEAT.R
3. Run Plotting_heat.R
