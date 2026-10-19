# Spatial Proteomics and Protein Subcellular Localization Prediction

A reproducible R-based spatial proteomics workflow using the `hyperLOPIT2015` dataset to explore protein abundance profiles, visualize localization patterns, predict protein subcellular localization with TAGM-MAP, and compare predictions with stored dataset results.

## Project Overview

Spatial proteomics investigates the subcellular localization of proteins by measuring their abundance across cellular fractions. Proteins localized to the same organelle tend to exhibit similar abundance profiles.

This project uses the published `hyperLOPIT2015` mouse embryonic stem-cell dataset, provided by the Bioconductor package `pRolocdata`, as a training and reproducibility example.

The workflow is inspired by spatial proteomics approaches used in published research, including work on *Eimeria tenella*. The current analysis uses the hyperLOPIT2015 dataset rather than experimental *Eimeria* data.

## Objectives

- Inspect quantitative protein abundance profiles and localization marker annotations.
- Explore marker profiles and mean abundance profiles across subcellular compartments.
- Visualize protein distributions using Principal Component Analysis (PCA) and t-distributed Stochastic Neighbor Embedding (t-SNE).
- Train a TAGM-MAP model using proteins with known localization labels.
- Predict localization assignments for proteins initially labelled as unknown.
- Filter predictions using a localization probability threshold of 0.99.
- Compare newly generated predictions with the TAGM-MAP results stored in the dataset.

## Dataset

**Dataset:** `hyperLOPIT2015`  
**Source:** `pRolocdata` Bioconductor package  
**Organism:** Mouse (*Mus musculus*), embryonic stem cells  
**Publication PMID:** [26754106](https://pubmed.ncbi.nlm.nih.gov/26754106/)

The dataset contains 5,032 protein features and 20 quantitative samples across two replicates. Protein annotations include known subcellular localization markers and proteins labelled as unknown.

The quantitative data are provided by the package, so a separate download of the raw dataset is not required for this workflow.

## Analysis Pipeline

1. **Load data:** Import the dataset and inspect protein, sample, and annotation metadata.
2. **Explore marker profiles:** Compare abundance distributions for selected localization classes and calculate mean profiles.
3. **Dimensionality reduction:** Use PCA and t-SNE to explore similarities and patterns among protein profiles.
4. **Model training:** Fit a TAGM-MAP model using known localization marker classes.
5. **Localization prediction:** Estimate class assignments and associated probabilities for proteins.
6. **Confidence filtering:** Calculate localization probability as:

   `localization probability = MAP allocation probability × (1 − outlier probability)`

   Retain predictions with a localization probability of at least 0.99.
7. **Visualization:** Plot marker classes and final confident assignments.
8. **Validation:** Compare newly generated predictions with the TAGM-MAP results stored in the dataset.

## Results From the Existing Analysis

The current notebook reports the following results:

| Metric | Reported result |
|---|---:|
| Total proteins | 5,032 |
| Initially unknown proteins | 4,106 |
| Unknown proteins passing the 0.99 threshold | 1,983 |
| Percentage of unknown proteins passing the threshold | 48.3% |
| Agreement in localization assignments across all proteins | 99.44% |
| Agreement among initially unknown proteins | 99.32% |
| Pearson correlation between stored and newly calculated MAP probabilities | 0.9982 |
| Proteins with differing localization assignments | 28 |

These values describe the existing analysis run and should be confirmed by rerunning the notebook in the target R environment.

A prediction passing the probability threshold is a model-based result, not independent experimental confirmation of a protein's biological localization.

## Repository Structure

```text
spatial-proteomics-hyperLOPIT/
├── README.md
├── analysis/
│   ├── spatial_proteomics_hyperLOPIT2015.Rmd
│   └── spatial_proteomics_hyperLOPIT2015.html
└── results/
    ├── figures/
```

## Requirements

- R
- RStudio (recommended)
- Bioconductor
- `MSnbase`
- `pRoloc`
- `pRolocdata`
- `Rtsne`

## Installation

Install the required packages in R:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "MSnbase",
  "pRoloc",
  "pRolocdata"
))

install.packages("Rtsne")
```

## Run the Analysis

1. Clone or download this repository.
2. Open the `.Rmd` file in RStudio.
3. Install the required packages.
4. Knit the notebook to HTML, or run the code chunks sequentially.
5. Review the figures, prediction tables, and validation results.

The dataset is loaded directly from `pRolocdata` using `data(hyperLOPIT2015)`.

## Outputs

The intended outputs include:

- PCA and t-SNE visualizations.
- Marker abundance profile plots.
- A visualization of final localization assignments.
- A table of confident predicted protein localizations.
- A summary of predicted localization classes.
- A comparison of newly generated and stored TAGM-MAP results.
- An HTML report containing the analysis and interpretation.

Only outputs that have actually been generated and checked should be described as available files in the repository.

## Limitations

- Model-based predictions require independent evidence before being treated as confirmed biological localizations.
- Localization classes have different numbers of marker proteins, which can affect model performance.
- Results may vary with package versions, model fitting, and computational settings.

## References

- hyperLOPIT2015 dataset publication: [PubMed PMID 26754106](https://pubmed.ncbi.nlm.nih.gov/26754106/)
- Bioconductor package: [pRoloc](https://bioconductor.org/packages/pRoloc/)
- Bioconductor package: [pRolocdata](https://bioconductor.org/packages/pRolocdata/)

## Author

**Shrouk Hany**

This repository is intended for educational, exploratory, and reproducibility purposes.
