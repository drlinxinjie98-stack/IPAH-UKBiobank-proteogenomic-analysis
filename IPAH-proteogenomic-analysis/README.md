# Proteogenomic profiling of idiopathic pulmonary arterial hypertension

This repository contains the R scripts used for the analyses reported in the study, organized in their intended execution order.

## Repository contents

| Script | Purpose |
|---|---|
| `01_data_cleaning.R` | UK Biobank phenotype, medication, procedure, mortality, physical-activity, imaging, and proteomic data preparation |
| `02_discovery_individual_data.R` | Cohort construction, matching, and individual-level discovery analyses |
| `03_mr_cross_epidemiology.R` | Integration of epidemiological associations with Mendelian randomization results |
| `04_cross_ancestry.R` | Exploratory cross-ancestry effect-direction concordance analysis |
| `05_sensitivity_analysis.R` | Sensitivity analyses |
| `06_clustering.R` | Protein-based clustering, survival comparisons, heatmaps, and cluster interaction analyses |
| `07_manhattan_plots.R` | Manhattan-style plots of protein association results |
| `08_enrichment.R` | GO, KEGG, and Reactome enrichment analyses |
| `09_predictive_models.R` | Predictive modeling, ROC comparisons, and SHAP-based interpretation |

## Data availability and governance

No individual-level participant data are included in this repository. UK Biobank data must be accessed by eligible researchers under an approved application and analyzed in accordance with the UK Biobank data-use agreement and applicable platform requirements. Other input data are subject to the access conditions of their respective sources.

The scripts use project-relative paths beginning with `data/` and write outputs beneath `data/` or `Results/`. Users must provide appropriately authorized input data in the expected directory structure. File names document the required inputs; the repository does not create or distribute restricted datasets.

## Software environment

The code was prepared under R 4.5.2. Versions of the directly used R packages are listed in `package_versions.tsv`. Package versions should be checked against the environment used for the final archived analysis.

## Reproducibility notes

- Run the scripts in numerical order when reproducing the complete workflow.
- Review all input and output paths before execution.
- Several analyses use parallel processing; adjust the number of cores to the available computing environment.
- Key statistical specifications, covariates, thresholds, matching settings, and modeling parameters are defined in the scripts and described in the manuscript Methods.
- No participant identifiers, derived individual-level data, cached workspaces, credentials, or access tokens should be committed to this repository.

## Code availability

The analysis code supporting the findings of this study is provided in this repository. Questions about the code may be directed to the corresponding author. Access to the underlying individual-level data is governed separately by the relevant data custodians and is not provided through this repository.

