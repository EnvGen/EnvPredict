
rule all:
    input:
        "plots/Figure_2/FigS1.pdf",
        "plots/plots_microscopy/predicting_zooplankton.pdf",
        "plots/interannual_comparison_plots/phys_chem_interannual_comparison.pdf"

rule core:
    input:
        "seq_files/18S/filtered_seqtab_18S.tsv",
        "seq_files/18S/filtered_taxa_18S.tsv",
        "seq_files/18S/seqtab_18S.tsv",
        "seq_files/18S/taxa_18S.tsv",
        "seq_files/16S/filtered_seqtab_16S.tsv",
        "seq_files/16S/filtered_taxa_16S.tsv",
        "env_files/physical_chemical_processed_translation.tsv",
        "env_files/phytoplankton_processed.tsv",
        "env_files/zooplankton_processed.tsv"
    output:
        "output/DifferentTaxonomicLevels/norm_clade_counts_16S_2_RF10fold_Predictions.tsv",
        "output/zooplankton_predicted/norm_asv_counts_16S_zoo_plan_genus_RF10fold_Predictions.tsv",
        "output/predict_2015_2017/norm_asv_counts_16S_phys_chem_2015_2017_different_years_Predictions.tsv"
    shell:
        """
        cd code && Rscript --verbose envpredict_core.R
        """

rule TabPFN:
    input:
        "output/DifferentTaxonomicLevels/norm_clade_counts_16S_2_RF10fold_Predictions.tsv",
        "env_files/physical_chemical_processed_translation.tsv"
    output:
        "output/TabPFN/norm_clade_counts_16S_2_RF10fold_Predictions.tsv"
    shell:
        """
        cd code && python tabpnf_predictions.py
        """


rule Fig2:
    input:
        "output/DifferentTaxonomicLevels/norm_clade_counts_16S_2_RF10fold_Predictions.tsv",
        "output/TabPFN/norm_clade_counts_16S_2_RF10fold_Predictions.tsv"
    output:
        "plots/Figure_2/FigS1.pdf"
    shell:
        """
        cd code && Rscript --verbose -e "rmarkdown::render('Plot_Figure_2.Rmd')"
        """

rule Fig34:
    input:
        "output/zooplankton_predicted/norm_asv_counts_16S_zoo_plan_genus_RF10fold_Predictions.tsv"
    output:
        "plots/plots_microscopy/predicting_zooplankton.pdf"
    shell:
        """
        cd code && Rscript --verbose envpredict_eval.R
        """

rule Fig5:
    input:
        "output/predict_2015_2017/norm_asv_counts_16S_phys_chem_2015_2017_different_years_Predictions.tsv"
    output:
        "plots/interannual_comparison_plots/phys_chem_interannual_comparison.pdf"
    shell:
        """
        cd code && Rscript --verbose -e "rmarkdown::render('interannual_comparison.Rmd')"
        """
