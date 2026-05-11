\# README for submitted MS



Author: TZ  

Contact: tzzteresa@gmail.com  

Last Update: 06/2026



\## Overview:

This repository contains the following R scripts:



1\. \*\*Cost-Distance Calculation\*\*

&#x20;  - Script: 00\_cost\_distance.R

&#x20;  - Optional: results already included in Source Data 1 

&#x20;  - Calculates conductance rasters and cost distances for three categories of plants based on substrate preference.  

&#x20;  - Inputs: lithology shapefile (not provided), conversion tables (Supplementary Table 1), site points   

&#x20;  - Outputs: raster files, processed point shapefiles, distance matrices.



3\. \*\*Statistical Analyses and Plots\*\*  

&#x20;  - Script: 01\_figures.R  

&#x20;  - Performs descriptive statistics, some GLM/GLMM analyses, and produces manuscript figures  

&#x20;  - Inputs: Source Data 1, optional layers for map  

&#x20;  - Outputs: figures



4\. \*\*Large GLMM/GLM Loops\*\*  

&#x20;  - Script: 02\_models\_loop.R  

&#x20;  - Fits GLMMs and GLMs to model the probability of polyploidy as a function of predictors for whole dataset and for mixed-ploidy genera/species.

&#x20;  - Inputs: Source Data 1

&#x20;  - Outputs: model summaries for all GLMM results across response × distance × model × random structure, GLMs for mixed-ploidy species / aggregates, and GLMs for mixed-ploidy genera; figures for model prefictions; figures for effects



5\. \*\*Phylogenetic correction\*\*

&#x20;  - Script: 03\_phylosig\_pglm.R

&#x20;  - Calculate phylogenetic signal of ploidy and fit phyloglm

&#x20;  - Inputs: phyloalps tree, Source Data 1



\## Required files:

\- R scripts

\- Source Data 1 (for 01\_figures.R, 02\_models\_loop.R, 03\_phylosig\_pglm.R)

\- Not provided for 00\_cost\_distance.R: full\_geo\_lim\_3035\_V2.gpkg (published in Donnini et al. 2019)

\- Conversion tables for  00\_cost\_distance.R: from Supplementary Table 1 



\## 1. System requirements

Operating Systems: tested on Windows 11



R Version: R >= 4.2.0 (analyses conducted in R 4.4.1)



Required R packages:  

phytools\_2.4-4

maps\_3.4.2.1

phylolm\_2.6.5

phyr\_1.1.3

here\_1.0.1

ape\_5.8-1 

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



\## 2. Installation

1\. Install R (CRAN), optionally RStudio, and required R packages  (typical installation times, here)

2\. Place scripts in appropriate folder: scripts/ 

3\. Place input files in appropriate folders:

&#x20;  - data\_in/: conversion tables (Supplementary Table 1)

&#x20;  - data/pro/: Source Data 1 (Source\_Data.xlsx)

&#x20;  - GIS\_data/: shapefiles, GeoPackages



\## 3. Demo

\*\*Cost-Distance Calculation\*\*

\- Prepare full\_geo\_lim\_3035\_V2.gpkg (published in Donnini et al. 2019)

\- Run script 00\_cost\_distance.R

\- Output: conductance rasters, processed point shapefiles, distance matrices for 3 distance versions and 3 plant categories

\- Expected run time: ca. 3 h



\*\*Statistical Analyses and Plots\*\*

\- Run script 01\_figures.R  

\- Output: PDF/PNG figures

\- Expected run time: 5 min



\*\*Large GLMM/GLM Loops\*\*

\- Run script 02\_models\_loop.R  

\- Output: model summaries for all GLMM results across response × distance × model × random structure, GLMs for mixed-ploidy species, and GLMs for mixed-ploidy genera; figures for model predictions; figures for effects

\- Expected run time: 5 min



\*\*Phylogenetic correction\*\*

\- Run script 03\_phylosig\_pglm.R

\- Output: Pagel's Lambda and phylogenetic GLM summary

\- Expected run time: 5 min





\## 4. Instructions

1\. Adjust filtering options, model variables, and predictor selections as needed.  

2\. Run scripts in R or RStudio.  



Optional:

\- Subset analyses by light preference, dispersal distance, or taxon

\- Change distance and ploidy versions

\- Skip GLMMs or pSEM for faster runs



\## 5. Contact

For questions regarding scripts or datasets, contact the author TZ.  



