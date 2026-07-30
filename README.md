# Distribution patterns of polyploid plants

**Author:** Teresa Zeni  
**Last updated:** July 2026

## Overview

This repository contains the R scripts used for the analyses and figure generation associated with a manuscript submitted to *Nature Communications*.

## Code and citation

The R scripts associated with this study are archived on Zenodo: https://doi.org/10.5281/zenodo.21701400 

If you use these scripts, please cite the Zenodo record.

The repository contains the following R scripts:

### 1. Cost-Distance Calculation

**Script:** `00_cost_distance.R`

Calculates conductance rasters and cost distances for three categories of plants based on substrate preference.

**Inputs:**

-   Lithology data (not included in this repository; see below)
-   Conversion tables derived from Supplementary Table 1
-   Site point data

**Outputs:**

-   Conductance rasters
-   Processed point shapefiles
-   Distance matrices

------------------------------------------------------------------------

### 2. Statistical Analyses and Figures

**Script:** `01_figures.R`

Performs descriptive statistics and selected GLM/GLMM analyses and generates figures presented in the manuscript.

**Inputs:**

-   Supplementary Data 1
-   Optional spatial layers for map generation

**Outputs:**

-   Manuscript figures

------------------------------------------------------------------------

### 3. GLMM/GLM Analyses

**Script:** `02_models_loop.R`

Fits GLMMs and GLMs to model the probability of polyploidy as a function of predictor variables for the full dataset and for mixed-ploidy genera and species.

**Inputs:**

-   Supplementary Data 1

**Outputs:**

-   Model summaries for GLMM analyses across response × distance × model × random-effect structures
-   GLMs for mixed-ploidy species/aggregates
-   GLMs for mixed-ploidy genera
-   Figures showing model predictions
-   Figures showing model effects

------------------------------------------------------------------------

### 4. Phylogenetic Correction

**Script:** `03_phylosig_pglm.R`

Calculates phylogenetic signal of ploidy and fits phylogenetic generalized linear models (`phyloglm`).

**Inputs:**

-   Phylogenetic tree
-   Supplementary Data 1

**Outputs:**

-   Phylogenetic signal estimates (Pagel's lambda)
-   Phylogenetic GLM results

------------------------------------------------------------------------

## Required Files and External Resources

The following files are required to reproduce the analyses.

### Included in the repository

-   R scripts
-   Conversion tables derived from Supplementary Table 1

### Required but not included

-   **Supplementary Data 1** — required for `01_figures.R`, `02_models_loop.R`, and `03_phylosig_pglm.R`
-   **Lithology data:** `full_geo_lim_3035_V2.gpkg` — published in Donnini et al. (2019)
-   **Phylogenetic tree:** Source Data associated with Fig. 2 in Wootton et al. (2025), *Nature Plants*

Please refer to the original publications and associated data repositories for access to these external resources.

------------------------------------------------------------------------

## 1. System Requirements

### Operating System

-   Windows 11 (tested)

### R Version

-   R ≥ 4.2.0
-   Analyses conducted using R 4.4.1

### Required R Packages

The analyses require the following R packages:

-   `phytools` 2.4-4
-   `maps` 3.4.2.1
-   `phylolm` 2.6.5
-   `phyr` 1.1.3
-   `here` 1.0.1
-   `ape` 5.8-1
-   `broom.mixed` 0.2.9.6
-   `broom` 1.0.7
-   `MuMIn` 1.48.11
-   `performance` 0.15.2
-   `visreg` 2.7.0
-   `dotwhisker` 0.8.3
-   `sjlabelled` 1.2.0
-   `sjmisc` 2.8.10
-   `car` 3.1-3
-   `carData` 3.0-5
-   `DiagrammeR` 1.0.11
-   `ggeffects` 1.7.2
-   `sjPlot` 2.8.16
-   `piecewiseSEM` 2.3.0
-   `lme4` 1.1-35.5
-   `ggpubr` 0.6.0
-   `rstatix` 0.7.2
-   `janitor` 2.2.0
-   `patchwork` 1.3.2
-   `tidyterra` 0.6.1
-   `gdistance` 1.6.4
-   `Matrix` 1.7-0
-   `igraph` 2.1.1
-   `raster` 3.6-30
-   `sp` 2.1-4
-   `terra` 1.7-83
-   `sf` 1.0-19
-   `rio` 1.2.3
-   `lubridate` 1.9.3
-   `forcats` 1.0.0
-   `stringr` 1.5.1
-   `dplyr` 1.1.4
-   `purrr` 1.0.2
-   `readr` 2.1.5
-   `tidyr` 1.3.1
-   `tibble` 3.2.1
-   `ggplot2` 4.0.2
-   `tidyverse` 2.0.0

No non-standard hardware is required.

------------------------------------------------------------------------

## 2. Installation

1.  Install R from CRAN and, optionally, RStudio.
2.  Install the required R packages.
3.  Place the R scripts in the `scripts/` directory.
4.  Place the required input files in the appropriate directories:

``` text
data_in/
├── conversion tables derived from Supplementary Table 1

data/
└── pro/
    └── Supplementary_Data_1.xlsx

GIS_data/
└── required shapefiles and GeoPackages
```

------------------------------------------------------------------------

## 3. Reproducing the Analyses

### Cost-Distance Calculation

1.  Obtain `full_geo_lim_3035_V2.gpkg` from the original source (Donnini et al., 2019).
2.  Place the file in the appropriate `GIS_data/` directory.
3.  Run `00_cost_distance.R`.

**Outputs:**

-   Conductance rasters
-   Processed point shapefiles
-   Distance matrices for three distance versions and three plant categories

**Expected runtime:** approximately 3 hours.

------------------------------------------------------------------------

### Statistical Analyses and Figures

1.  Ensure that Supplementary Data 1 is available in the appropriate directory.
2.  Run `01_figures.R`.

**Outputs:**

-   PDF/PNG figures

**Expected runtime:** approximately 5 minutes.

------------------------------------------------------------------------

### GLMM/GLM Analyses

1.  Ensure that Supplementary Data 1 is available in the appropriate directory.
2.  Run `02_models_loop.R`.

**Outputs:**

-   Model summaries for GLMM analyses across response × distance × model × random-effect structures
-   GLMs for mixed-ploidy species
-   GLMs for mixed-ploidy genera
-   Figures showing model predictions
-   Figures showing model effects

**Expected runtime:** approximately 5 minutes.

------------------------------------------------------------------------

### Phylogenetic Correction

1.  Ensure that the required phylogenetic tree and Supplementary Data 1 are available.
2.  Run `03_phylosig_pglm.R`.

**Outputs:**

-   Pagel's lambda estimates
-   Phylogenetic GLM (`phyloglm`) results

**Expected runtime:** approximately 5 minutes.

------------------------------------------------------------------------

## 4. Analysis Options

The scripts can be adapted to explore alternative analytical settings.

Users may:

-   Adjust filtering options
-   Modify model variables and predictor selection
-   Subset analyses by light preference
-   Subset analyses by dispersal distance
-   Subset analyses by taxonomic group
-   Change distance and ploidy versions
-   Skip computationally intensive GLMM or pSEM analyses for faster execution

------------------------------------------------------------------------

## 5. Citation

If you use the R scripts in this repository, please cite the software repository and the associated publication when available.

A citation for the archived version of the software will be added after the repository is archived on Zenodo.

The associated manuscript is currently submitted to *Nature Communications*.

------------------------------------------------------------------------

## 6. License

The R scripts in this repository are licensed under the **GNU General Public License v3.0 (GPL-3.0)**.

See the `LICENSE` file for the full license text.

------------------------------------------------------------------------

## 7. Contact

For questions regarding the R scripts or associated datasets, please contact:

Teresa Zeni  
**Email:** [tzzteresa\@gmail.com](mailto:tzzteresa@gmail.com)\
**ORCID:** <https://orcid.org/0009-0001-8982-7973>
