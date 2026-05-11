# Phylogenetic GLMM:
# Does elevation predict ploidy?
# added to revised manuscript
# v 2026.4.14

library(ape)
library(dplyr)
library(here)
library(rio)
library(phyr)
library(phylolm)
library(phytools)

# 1. Load data
df      <- import(here("data", "pro", "Source_Data.xlsx"), which = 2)
log_raw <- import(here("data", "pg", "log.csv"))

phy_orig =  read.tree(here("data", "pg", "231108_pruned_phyloalps_tree.tre"))


# 2. Prepare species-to-tip lookup (direct matches only)
log_table <- log_raw |>
  filter(action == "direct match") |>
  select(spid = species, label = tip_label)

# 3. Prepare variables
df <- df |>
  mutate(
    elevation       = ifelse(belt == "N", NA, elevation),
    elevation_rel   = elevation - treeline,
    elevation       = scale(elevation_rel)[, 1],
    ploidy          = ifelse(ploidy_orig != 2, 1, 0),
    distance        = scale(dist_v0)[, 1],
    spid            = paste0(genus, "_", species)
  ) |>
  left_join(log_table, by = "spid")



# 4. Minimum sample filter (observations per species)
MIN_OBS <- 10

occ_counts   <- table(df$abbreviation)
keep_species <- names(occ_counts[occ_counts >= MIN_OBS])

df <- df |>
  filter(abbreviation %in% keep_species)

cat("Species retained after min-sample filter:", n_distinct(df$label), "\n")
cat("Observations retained:                   ", nrow(df), "\n")



# 5. Subset to columns needed for modelling & match tree tips
df <- df |>
  select(species = label, elevation, distance, ploidy) |>
  filter(species %in% phy_orig$tip.label)

cat("Species with tree match:                 ", n_distinct(df$species), "\n")
cat("Observations with tree match:            ", nrow(df), "\n")



# 6. Drop NAs (required for all models below)
df <- df |> na.omit()

cat("Species after dropping NAs:              ", n_distinct(df$species), "\n")
cat("Observations after dropping NAs:         ", nrow(df), "\n")



# 7. Prune tree to modelling species
model_species <- unique(df$species)
phy_pruned = drop.tip(phy_orig, setdiff(phy_orig$tip.label, model_species))
phy_pruned$edge.length

cat("Tips in pruned tree:                     ", length(phy_pruned$tip.label), "\n")
cat("Tips missing from data (should be 0):    ",
    length(setdiff(phy_pruned$tip.label, model_species)), "\n")

# ANALYSIS 1 — Phylogenetic GLM (phyloglm)
# using species-level summaries
# Aggregate to one row per species (proportion polyploid)
df_sp <- df |>
  group_by(species) |>
  summarise(
    ploidy    = mean(ploidy),        # proportion polyploid
    elevation = mean(elevation, na.rm = TRUE),
    distance  = max(distance,  na.rm = TRUE),
    .groups   = "drop"
  ) |>
# filter(ploidy == 0 |ploidy == 1) |> # unifrom
  filter(species %in% phy_pruned$tip.label)
#plot(phy_sp)

phy_sp <- drop.tip(phy_pruned, setdiff(phy_pruned$tip.label, df_sp$species))

cat("Species in phyloglm dataset:             ", nrow(df_sp), "\n")

# Note: phyloglm requires ploidy as numeric 0/1
# Using proportion here; to binary if strictly required
df_sp$ploidy_bin <- ifelse(df_sp$ploidy > 0.5, 1, 0)


rownames(df_sp) <- df_sp$species

# model 
mod_phyloglm <- phyloglm(
  ploidy_bin ~ elevation,
  data             = df_sp,
  phy              = phy_sp,
  method           = "logistic_MPLE",
  btol             = 50,
  log.alpha.bound  = 15
  
)
mod_phyloglm$alpha
summary(mod_phyloglm)

# ANALYSIS 2 — Pagel's lambda (phylosig)
# Tests whether ploidy shows phylogenetic signal across species

# Requires a single numeric value per species
ploidy_vec <- setNames(df_sp$ploidy, df_sp$species)

# Pagel's Lambda
lambda_pagel <- phylosig(
  phy_sp,
  ploidy_vec,
  method = "lambda",
  test   = TRUE
)

print(lambda_pagel)