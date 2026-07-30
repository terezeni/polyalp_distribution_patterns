# Distribution patterns of polyploid plants

Author: TZ  

Contact: tzzteresa@gmail.com  

Last Update: 07/2026


## Overview:

This repository contains the following R scripts:


1. **Cost-Distance Calculation**

- Script: 00_cost_distance.R

- Optional: results already included in Source Data 

- Calculates conductance rasters and cost distances for three categories of plants based on substrate preference.  

- Inputs: lithology shapefile (not provided), conversion tables (Supplementary Table 1), site points   

- Outputs: raster files, processed point shapefiles, distance matrices.



3. **Statistical Analyses and Plots**  

- Script: 01_figures.R  

- Performs descriptive statistics, some GLM/GLMM analyses, and produces manuscript figures  

- Inputs: Supplementary Data 1, optional layers for map  

- Outputs: figures



4. **Large GLMM/GLM Loops**  

- Script: 02_models_loop.R  

- Fits GLMMs and GLMs to model the probability of polyploidy as a function of predictors for whole dataset and for mixed-ploidy genera/species.

- Inputs: Supplementary Data 1

- Outputs: model summaries for all GLMM results across response × distance × model × random structure, GLMs for mixed-ploidy species / aggregates, and GLMs for mixed-ploidy genera; figures for model prefictions; figures for effects



5. **Phylogenetic correction**

- Script: 03_phylosig_pglm.R

- Calculate phylogenetic signal of ploidy and fit phyloglm

- Inputs: phyloalps tree, Supplementary Data 1



## Required files:

- R scripts

- Supplementary Data 1 (for 01_figures.R, 02_models_loop.R, 03_phylosig_pglm.R)

- Not provided for 00_cost_distance.R: full_geo_lim_3035_V2.gpkg (published in Donnini et al. 2019)

- Conversion tables for  00_cost_distance.R: from Supplementary Table 1
  
- Phylogeny (Source Data Fig. 2 in [Wootton et al., Nature Plants (2025)](https://www.nature.com/articles/s41477-025-02001-0) )



## 1. System requirements

Operating Systems: tested on Windows 11



R Version: R >= 4.2.0 (analyses conducted in R 4.4.1)



Required R packages:  
phytools_2.4-4
maps_3.4.2.1
phylolm_2.6.5
phyr_1.1.3
here_1.0.1
ape_5.8-1 
broom.mixed 0.2.9.6  
broom 1.0.7  
MuMIn 1.48.11  
performance 0.15.2  
visreg 2.7.0  
dotwhisker 0.8.3  
sjlabelled 1.2.0  
sjmisc 2.8.10  
car 3.1-3  
carData 3.0-5  
DiagrammeR 1.0.11  
ggeffects 1.7.2  
sjPlot 2.8.16  
piecewiseSEM 2.3.0  
lme4 1.1-35.5  
ggpubr 0.6.0  
rstatix 0.7.2  
janitor 2.2.0  
patchwork 1.3.2  
tidyterra 0.6.1  
gdistance 1.6.4  
Matrix 1.7-0  
igraph 2.1.1  
raster 3.6-30  
sp 2.1-4  
terra 1.7-83  
sf 1.0-19  
rio 1.2.3  
lubridate 1.9.3  
forcats 1.0.0  
stringr 1.5.1  
dplyr 1.1.4  
purrr 1.0.2  
readr 2.1.5  
tidyr 1.3.1  
tibble 3.2.1  
ggplot2 4.0.2  
tidyverse 2.0.0

Non-standard hardware: None required

## 2. Installation

1. Install R (CRAN), optionally RStudio, and required R packages  (typical installation times, here)

2. Place scripts in appropriate folder: scripts/ 

3. Place input files in appropriate folders:

- data_in/: conversion tables (Supplementary Table 1)

- data/pro/: Supplementary Data 1 (Supplementary_Data_1.xlsx)

- GIS_data/: shapefiles, GeoPackages



## 3. Demo

**Cost-Distance Calculation**

- Prepare full_geo_lim_3035_V2.gpkg (published in Donnini et al. 2019)

- Run script 00_cost_distance.R

- Output: conductance rasters, processed point shapefiles, distance matrices for 3 distance versions and 3 plant categories

- Expected run time: ca. 3 h



**Statistical Analyses and Plots**

- Run script 01_figures.R  

- Output: PDF/PNG figures

- Expected run time: 5 min



**Large GLMM/GLM Loops**

- Run script 02_models_loop.R  

- Output: model summaries for all GLMM results across response × distance × model × random structure, GLMs for mixed-ploidy species, and GLMs for mixed-ploidy genera; figures for model predictions; figures for effects

- Expected run time: 5 min



**Phylogenetic correction**

- Run script 03_phylosig_pglm.R

- Output: Pagel's Lambda and phylogenetic GLM summary

- Expected run time: 5 min





## 4. Instructions

1. Adjust filtering options, model variables, and predictor selections as needed.  

2. Run scripts in R or RStudio.  



Optional:

- Subset analyses by light preference, dispersal distance, or taxon

- Change distance and ploidy versions

- Skip GLMMs or pSEM for faster runs



## 5. Contact

For questions regarding scripts or datasets, contact the author TZ.  



