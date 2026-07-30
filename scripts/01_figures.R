## descriptive stat
## tests
## except for large GLMMs/GLM loops
## figures for MS
## TZ last update 10/04/2026

# SETUP ----
## packages
library(tidyverse)
library(rio)
library(terra)
library(tidyterra)
library(janitor)
library(patchwork)
library(rstatix) 
library(ggpubr)
library(rstatix)
library(lme4)
library(piecewiseSEM)
library(sjPlot)
library(ggeffects)
library(DiagrammeR) # paths diagram
library(extrafont)
library(openxlsx)

source_wb <- createWorkbook()
loadfonts(device = "win")


# DIRECTORIES ----
project_dir <- here::here()

## paths
data_path <- file.path(project_dir, "data")
gis_path  <- file.path(project_dir, "GIS_data")
output_path <- file.path(project_dir, "output")

## options
# ploidy_ version one of "ploidy_orig" or "ploidy_alt"
ploidy_version = "ploidy_orig"

# subset: one of "all", "far", "close"
subset = "all"

# distance version "dist_v0", "dist_vA", "dist_v_B"
distance_version = "dist_v0"

## output directory (subset and distance version)
output_dir = paste0(output_path, "/", subset,"_",distance_version)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

## color palettes
pal_belt = c("#006147", "#80A552", "#E7D47A", "#D09D4D", "#B56721")
pal_region = c(
  "SLA" = "#721f81",
  "NLA" = "#1f78b4",
  "CEA" = "#e31a1c"
)
effect_cols = c("#8654AF",'#Cfa0fb') 

## theme ggplot
theme <- theme_classic(base_family = "Arial", base_size = 10) +
  theme(panel.border = element_rect(color = "black", 
                                    fill = NA),
        legend.background = element_blank(),
        legend.box.background = element_blank())

# IMPORT DATA ----
## dataset ----
data <- import(paste0(data_path, "/pro/Source_Data_rev.xlsx"), which = 2)

### optional layers for map ----
alp = vect(paste0(gis_path,"/countries/borders_Alps.shp"))
kref = vect(paste0(gis_path,"/refugia/k_refs_complete.gpkg"))
sref = vect(paste0(gis_path,"/refugia/s_refs_complete.gpkg"))
subunits = vect(paste0(gis_path,"/refugia/EA/EA.shp"))
points = vect(paste0(gis_path,"/mt_geo_84.shp"))

# PREPARE DATA ----
data <- data |> clean_names()

## relative elevation to timberline ----
data <- data |> 
  mutate(loc = paste0(mountain_code, "-", belt),
         elevation = ifelse(belt == "N", NA, elevation),
         elevation_orig = elevation,
         elevation_relative = elevation_orig - treeline,
         elevation = elevation_relative,
        ploidy_inferred = as.logical(as.numeric(ploidy_inferred)),
         subunit = as.factor(subunit),
         refugium = as.factor(refugium),
        dispersal_ability = as.factor(dispersal_ability),
        family = as.factor(family),
        polyploid_orig = ifelse(ploidy_orig != 2, 1, 0),
        polyploid_alt = ifelse(ploidy_alt != 2, 1, 0)
  )
# data |> write.csv("temp.csv")
levels(data$subunit)  = c("CEA", "NLA", "SLA")
levels(data$refugium) = c("g", "r")

data$subunit = fct_relevel(data$subunit, c('SLA', 'NLA', "CEA"))

## species richness ----
mt<- data |> 
  mutate(grouping = paste0(abbreviation, "_", mountain_code)) |>
  distinct(grouping, .keep_all = TRUE) |> ungroup() |> group_by(mountain_code) |> tally()

b <- data |> mutate(grouping = paste0(abbreviation, "_", mountain_code, "-", belt)) |>
  distinct(grouping, .keep_all = TRUE) |> ungroup() |>
  group_by(mountain_code, belt) |> distinct() |> tally()


su_belt <- data |> 
  mutate(grouping = paste0(abbreviation, "_", mountain_code, subunit, belt)) |>
  distinct(grouping, .keep_all = TRUE) |> ungroup() |> 
  group_by(subunit, mountain_code, belt) |> tally() |> 
  group_by(subunit, belt) |> summarise(mean = mean(n))


# avg number of species per mountain
mt$n |> mean()

# avg number of species per belt
b$n |> mean()



## scale distances and elevation for binomial GLM(M)s ----
data_model <- data |>
  mutate(across(c(dist_v0, dist_v_a, dist_v_b, elevation), ~ scale(.)[, 1], .names =
                  "{.col}"))


if (subset == "far") {
  data = data|>
    dplyr::filter(dispersal_ability == "far") |> droplevels()
  data_model = data_model|>
    dplyr::filter(dispersal_ability == "far") |> droplevels()
  n_samples_far = data$sample_id |> unique() |> length()
  n_species_far = data$abbreviation |> unique() |> length()
}

if (subset == "close") {
  data = data|>
    dplyr::filter(dispersal_ability == "close") |> droplevels()
  data_model = data_model|>
    dplyr::filter(dispersal_ability == "close") |> droplevels()
  n_samples_close = data$sample_id |> unique() |> length()
  n_species_close = data$abbreviation |> unique() |> length()
}

## rename ploidy helper ----
ploidy = data[ploidy_version]
colnames(ploidy) = "ploidy"
data$ploidy = as.factor(ploidy$ploidy)

## rename dist helper ----
distance = data[distance_version]
colnames(distance) = "distance"
data$distance = distance$distance

## calculate counts per belt ----
# keep one sample per belt for each species, number ofsamples per ploidy level
counts_belt <- data |>
  mutate(grouping = paste0(abbreviation, "_", mountain_code, "-", belt)) |>
  distinct(grouping, .keep_all = TRUE) |>
  dplyr::select(
    loc,
    sample_id,
    ploidy,
    mountain_code,
    belt,
    family,
    x,
    y,
    distance, subunit, refugium
  ) |>
  group_by(loc) |>
  mutate(dist_m = mean(distance, na.rm = T)) |> # avg distance for the belt
  group_by(mountain_code, belt, ploidy, dist_m, subunit, refugium) |>
  summarise(nums = n(), .groups = "drop") |> 
  filter(belt != "N") |> # exclude belt N
  pivot_wider(names_from = ploidy,
              values_from = nums,
              names_prefix = "cyt_")  # reshape

# Substitute all the NA values with 0
counts_belt[is.na(counts_belt)] <- 0

## summaries
# total number of species, number and frequency of diploids and higher ploidy levels
counts_belt <- counts_belt %>%
  mutate(total_spec = rowSums(dplyr::select(., contains("cyt_")), na.rm = TRUE),
    n_px = total_spec - cyt_2,
    f_2x = cyt_2 / total_spec,
    f_px = n_px / total_spec,
    belt = as.factor(belt),
    mountain_code = as.factor(mountain_code)
  ) 

## dataframe for plotting by subunit, belt level
counts_geo_belt <- counts_belt |>
  group_by(subunit) |> 
  mutate(subunit = fct_relevel(subunit, c('SLA', 'NLA', "CEA")))

## subset only complete transects ---- 
## frequency of polyploids with complete transects
counts_belt_complete <- counts_belt |>
  group_by(mountain_code) |>
  filter(n_distinct(belt) == 5) |>  # keep only mountains with all 5 belts
  ungroup() |>
  droplevels() |>
  ungroup()

## calculate mean distance at mountain level ---
counts_per_mountain <- data |>
  mutate(grouping = paste0(abbreviation, "_", mountain_code, "-", ploidy)) |>
  distinct(grouping, .keep_all = TRUE) |>
  dplyr::select(
    mountain_code,
    distance, refugium, subunit, ploidy) |>
  group_by(mountain_code, subunit, refugium) |>
  mutate(dist_m = mean(distance, na.rm = T)) |>
  group_by(mountain_code, ploidy, dist_m, subunit, refugium) |>
  summarise(nums = n(), .groups = "drop") |> 
  pivot_wider(names_from = ploidy,
              values_from = nums,
              names_prefix = "cyt_")


# substitute all the NA values with 0
counts_per_mountain[is.na(counts_per_mountain)] <- 0

# summaries
counts_per_mountain <- counts_per_mountain %>%
  mutate(
    total_spec = rowSums(dplyr::select(., contains("cyt_")), na.rm = TRUE),
    n_px = total_spec - cyt_2,
    f_2x = cyt_2 / total_spec,
    f_px = n_px / total_spec
  ) |> 
  mutate(dst = sqrt(dist_m))

## dst sqrt() and Z transformed
counts_per_mountain$dst = scale(counts_per_mountain$dst)[, 1]

## dataframe for plotting by subunit, mountain level
counts_geo_mountain <- counts_per_mountain |>
  group_by(subunit) |>
  ungroup()

## data for glm and glmms ----
# correct polyploid version
data_model <- data_model |> 
  mutate(polyploid = if (ploidy_version == "ploidy_orig") {
    polyploid_orig
  } else {
    polyploid_alt
  })

# correct ploidy version
ploidy = data_model[ploidy_version]
colnames(ploidy) = "ploidy"
data_model$ploidy = as.factor(ploidy$ploidy)

# correct distance version
distance = data_model[distance_version]
colnames(distance) = "distance"
data_model$distance = distance$distance


# ANALYSES and FIGURES ----
## belt and frequency ----
### Friedman test, frequency and belt ----
test_all = counts_belt_complete |>
    ungroup() |>
    friedman_test(f_px ~ belt | mountain_code)
counts_belt_complete$f_px |> max()
test_all$statistic




### Fig 2a: boxplot belt with Friedman test ----

# Minimal source for Fig 2a (belt frequency)
fig2a_source <- counts_belt_complete |>
  dplyr::select(mountain_code, belt, f_px, total_spec)

addWorksheet(source_wb, "Fig2a")
writeData(source_wb, "Fig2a", fig2a_source)

(bxp_belt_freq <- counts_belt_complete |>
    ggplot(aes(x = belt, y = f_px)) +
    # Jittered points with size mapped to total_spec
    geom_jitter(
      aes(size = total_spec, color = belt),
      shape = 1,
      width = 0.4,
      height = 0,
      alpha = .6
    ) +
    geom_boxplot(
      outliers = F,
      outlier.shape = 4,
      outlier.alpha = 1,
      outlier.size = 1,
      fill = NA
    ) +
    ylab("Frequency of polyploids") +
    xlab("Belt") +
    theme +
    theme(
      legend.position = c(.9, .2),
      legend.title = element_text(size = 8),
      legend.text  = element_text(size = 7),
      legend.key.size = unit(.8, "lines")
    ) +
    scale_colour_manual(values = pal_belt, guide = "none") +
        scale_radius(range = c(0, 3), # scale radius of points
                 name = "Number of \nspecies") +
    scale_y_continuous(limits = c(0, .51), expand = expansion(mult = c(0, .1))) +
    coord_flip() +
    stat_friedman_test(
      size = 2,
      aes(wid = mountain_code),
      label.x = 1,
      label.y = .1,
      label = "p = {p.format},\nn = {n}"
    )
)

### pairwise wilcoxon signed-rank test with Bonferroni Holm correction####
(  pairwise_test <- wilcox_test(
    data = counts_belt_complete,
    f_px ~ belt,
    paired = T,
    p.adjust.method = "holm" ))

# annotate plot with pairwise comparisons result
# annotation hard coded based on pairwise_test res
(
   bxp_belt_freq_letters <- bxp_belt_freq +
    stat_summary(
      geom = 'text',
      #label = c("a", "ab", "b", "c", "c"), # alt
      #label = c("a", "a", "a", "a", "a"), # far
      #label = c("a", "b", "c", "d", "e"), # close
      label = c("a", "a", "b", "c", "c"), # all
      fun.y = max,
      vjust = 0,
      hjust = 0,
      mapping = aes(y = 0.5)
    )
)

## save 
ggsave(
  create.dir = T,
  plot = bxp_belt_freq_letters,
  filename = paste0(ploidy_version, "_2a.pdf"),
  dpi = 600,
  device = cairo_pdf,
  height = 75,
  width = 75,
  units = "mm",
  path = paste0(output_dir, "/", ploidy_version)
)

## save 
# ggsave(
#   create.dir = T,
#   plot = bxp_belt_freq_letters,
#   filename = paste0(ploidy_version, "_2a.png"),
#   dpi = 600,
#   units = "mm",
#   height = 75,
#   width = 75,
#   path = paste0(output_dir, "/", ploidy_version)
# )

## belt and frequency refugia vs glaciated ---- 
# frequency of polyploids with complete transects
counts_geo_complete <- counts_geo_belt |>
  group_by(mountain_code) |>
  filter(n_distinct(belt) == 5) |>  # keep only mountains with all 5 belts
  ungroup() |>
  droplevels()

### Fig 2b: boxplot frequency elevation by glaciated vs refugia ----
counts_geo_belt_filtered <-  counts_geo_belt |>
  group_by(mountain_code) |>
  filter(n_distinct(belt) == 5) |>  # keep only mountains with all 5 belts
  ungroup() |>
  droplevels() |>
  mutate(refugium = factor(
    refugium,
    levels = c("g", "r"),
    labels = c("Glaciated", "Refugium")
  ))

# Minimal source for Fig 2b (belt x refugium frequency)
fig2b_source <- counts_geo_belt_filtered |>
  dplyr::select(mountain_code, belt, refugium, f_px, total_spec)

addWorksheet(source_wb, "Fig2b")
writeData(source_wb, "Fig2b", fig2b_source)

pos <- position_dodge(width = .9)

region_belt_plot_ref <- counts_geo_belt_filtered |>
  ggplot(aes(x = refugium, y = f_px)) +
  geom_boxplot(
    aes(group = interaction(belt, refugium), color = belt),
    outliers = F,
    outlier.shape = 4,
    outlier.alpha = 1,
    position = pos,
    fill = NA
  ) +
  ylab("Frequency of polyploids") +
  xlab("") +
  theme +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 7),
    legend.key.size = unit(.8, "lines")
  ) +
  scale_colour_manual(values = pal_belt,
                      name = "Belt",
                      guide = "none") +
  scale_y_continuous(limits = c(0, .56), expand = expansion(mult = c(0, .1))) +
  stat_friedman_test(
    aes(wid = mountain_code, group = belt),
    label = "{p.signif}\nn = {n}",
    size = 2,
    label.y = .45
  ) +
  theme(legend.position = "none")

counts_geo_belt |>
  group_by(mountain_code) |>
  filter(n_distinct(belt) == 5) |>  # keep only mountains with all 5 belts
  ungroup() |>
  droplevels() |>
  mutate(refugium = factor(
    refugium,
    levels = c("g", "r"))) |>
  filter(refugium=="g") |>
  droplevels() |>
friedman_test(f_px ~ belt | mountain_code)

## save
ggsave(
  create.dir = T,
  plot = region_belt_plot_ref,
  filename = paste0(ploidy_version, "_2b.pdf"),
  device = cairo_pdf,
  dpi = 600,
  units = "mm",
  width = 75,
  height = 75,
  path = paste0(output_dir, "/", ploidy_version)
)

## save
# ggsave(
#   create.dir = T,
#   plot = region_belt_plot_ref,
#   filename = paste0(ploidy_version, "_2b.png"),
#   dpi = 600,
#   units = "mm",
#   width = 75,
#   height = 75,
#   path = paste0(output_dir, "/", ploidy_version)
# )

## flip coords
region_belt_plot_ref = region_belt_plot_ref +
  coord_flip()

## distance and polypoidy ----
### model polyploid and avg distance ---- 
model_mt_dist = glm(cbind(n_px, total_spec - n_px) ~ dst,
                    data = counts_geo_mountain,
                    family = binomial)

summary(model_mt_dist)

### Fig 3a: frequency and distance with model ----
(dst_fpx_mountain = plot_model(
    show.values = F,
    model_mt_dist,
    type = "pred",
    terms = "dst",
    title = "",
    colors = "black",
    transform = "exp",
    show.p = TRUE)     +
    geom_point(
      data = counts_geo_mountain,
      aes(x = dst, y = f_px),
      size = 1,
      stroke = NA,
      alpha = 0.6) +
    scale_color_manual(
    values = c("gray70", "black"),
  name = "",
  labels = c("Glaciated", "Refugium")) +
    ylab("Frequency of polyploids") +
    xlab("Cost-distance from refugia") +
    theme+
    theme(legend.position = "none")
)

# extract OR and significance from model
coefs <- summary(model_mt_dist)$coefficients
or_df <- data.frame(term = rownames(coefs),
                    OR = exp(coefs[, 1]),
                    p = coefs[, 4]) |>
  filter(term != "(Intercept)") |>
  mutate(sig = case_when(p < 0.001 ~ "***", 
                         p < 0.01  ~ "**", 
                         p < 0.05  ~ "*", 
                         TRUE      ~ ""))

# add annotation manually
(dst_fpx_mountain_annot = dst_fpx_mountain +
    geom_text(
      data = or_df,
      aes(
        x = -1,
        y = .4,
        label = paste0("OR = ", round(OR, 2), sig)),
      inherit.aes = FALSE,
      hjust = 0
    ))



## save
ggsave(
  create.dir = T,
  plot = dst_fpx_mountain_annot,
  filename = paste0(ploidy_version, "_3a.pdf"),
  dpi = 600,
  units = "mm",
  device = cairo_pdf,
  width = 75,
  height = 75,
  path = paste0(output_dir, "/", ploidy_version)
)

## save
# ggsave(
#   create.dir = T,
#   plot = dst_fpx_mountain_annot,
#   filename = paste0(ploidy_version, "_3a.png"),
#   dpi = 600,
#   units = "mm",
#   width = 75,
#   height = 75,
#   path = paste0(output_dir, "/", ploidy_version)
# )

## polyploidy and elevation ----
### glmm global (all transects) ----
model <- data_model |>
  glmer(formula = polyploid ~ elevation + (1 | mountain_code),
        family = binomial)

summary(model)
# extract OR
fixef(model)[2] |> exp()  # or for MS

### binomial glm for each transect----
results <- data_model |>
  group_by(mountain_code) |>
  summarize(slope = coef(glm(polyploid ~ elevation, family = binomial))[2]) |> 
  mutate(odds_ratio = exp(slope))

## add to spatial object to plot
# check correct data type
points$mountain_code <- as.numeric(points$mountain)
results$mountain_code <- as.numeric(results$mountain_code)

# merge
points_vect <- left_join(points, results)


# crop raster to extent of points
alp_crop = crop(alp, ext(9, 16.7, 45.2, 48.2))

### Fig 2c: map with ORs ----
(map = ggplot() +
   geom_spatvector(data = alp_crop, fill = "grey95") +
   geom_spatvector(
     data = subunits,
     color = NA,
     alpha = .2,
     aes(fill = geo)
   ) +
   scale_fill_manual(values = pal_region, guide = "none") +
   ggnewscale::new_scale_fill() +
   ggnewscale::new_scale_color() +
   geom_spatvector(
     data = kref,
     fill = "#1f78b499",
     color = "#165680ff"
   ) +
   geom_spatvector(
     data = sref,
     fill = "#e31a1c99",
     color = "#991114ff"
   ) +
   geom_spatvector(
     data = points_vect,
     aes(fill = odds_ratio),
     size = 3,
     shape = 21,
     stroke = .1
   ) +
   scale_fill_gradientn(
     name = "OR relative\nelevation",
     colors = c("#006147ff", "white", "#B56721"),
     values = scales::rescale(c(
       min(results$odds_ratio, na.rm = TRUE),
       1,
       max(results$odds_ratio, na.rm = TRUE)
     )),
     # position white at 1
     limits = c(
       min(results$odds_ratio, na.rm = TRUE),
       max(results$odds_ratio, na.rm = TRUE)
     )
   ) +
   theme +
   theme(
     legend.position = c(.9, .3),
     legend.title = element_text(size = 8),
     legend.text  = element_text(size = 7),
     legend.key.size = unit(.8, "lines"),
     axis.title = element_blank()
   ) +
   guides(
     shape = guide_legend(override.aes = list(size = 2)),
     color = guide_legend(override.aes = list(size = 2))
   ) +
   scale_x_continuous(labels = ~ .x, breaks = c(10, 12, 14, 16)) +
   scale_y_continuous(labels = ~ .x, breaks = c(46, 48)) +
   coord_sf(expand = F,
            datum = sf::st_crs(4326))
)

# Get attributes + XY
points_df <- terra::as.data.frame(points) |> 
  dplyr::select(mountain, x, y) |> 
  mutate(mountain_code = as.numeric(mountain))

# Minimal source for Fig. 2c (with coordinates)
fig2c_source <- results |>
  dplyr::select(mountain_code, odds_ratio) |>
  dplyr::left_join(points_df, by = "mountain_code") |>
  dplyr::select(mountain_code, x, y, odds_ratio)



addWorksheet(source_wb, "Fig2c")
writeData(source_wb, "Fig2c", fig2c_source)

## save
ggsave(
  create.dir = T,
  filename = paste0(ploidy_version, "_2c.pdf"),
  path = paste0(output_dir, "/", ploidy_version),
  device = cairo_pdf,
  dpi = 600,
  height = 75,
  units = "mm",
  width = 75
)

## final Fig 2 ----
Fig2_a <- bxp_belt_freq_letters +
  theme(
  legend.title = element_text(size = 8),
  legend.text  = element_text(size = 7),
  legend.key.size = unit(.8, "lines")
) +
  guides(shape = guide_legend(override.aes = list(size = 2))) +
  annotate(
    "text",
    x = 5.2,
    y = .02,
    label = "bold(a)",
    parse = T
  )

Fig2_a_x <- Fig2_a +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank())

Fig2_c <- map +
  theme(axis.title.x = element_blank()) +
  annotate(
    "text",
    x = 9.2,
    y = 48,
    label = "bold(c)",
    parse = T
  )
(
  Fig2_b <- region_belt_plot_ref +
    scale_x_discrete(position = "top") +
    theme(
    axis.text.y = element_text(angle = 90))+
    annotate(
      "text",
      x = 2.4,
      y = .01,
      label = "bold(b)",
      parse = T
    )
)

# combine panels
(Fig2 <- ((Fig2_a + Fig2_b) / (Fig2_c)) + plot_layout(widths = c(1, 1)))


## save
ggsave(
  plot = Fig2,
  create.dir = T,
  filename = paste0("fig_2_current_f_", ploidy_version, ".pdf"),
  path = paste0(output_dir, "/", ploidy_version),
  device = cairo_pdf,
  height = 200,
  units = "mm",
  width = 175
)
temph = Fig2_a + Fig2_b + plot_layout(widths = c(1, 1))
Fig2h <- temph + Fig2_c + plot_layout(ncol = 3, nrow = 1, widths = c(1, 1, 2.5))
ggsave(
  plot = Fig2h,
  create.dir = T,
  filename = paste0("fig_2_current_h_", ploidy_version, ".pdf"),
  path = paste0(output_dir, "/", ploidy_version),
  device = cairo_pdf,
  height = 60,
  units = "mm",
  width = 180
)

## distance and frequency ----
### Fig 3b: boxplot by subunit ----
(Fig3_b <- counts_geo_mountain |>
    ggplot(aes(x = subunit, y = f_px)) +
    geom_jitter(
      aes(x = subunit,
        y = f_px,
        size = total_spec,
        color = refugium),
      shape = 1,
      width = 0.4,
      height = 0,
      alpha = .6) +
    geom_boxplot(
      outliers = F,
      outlier.shape = 4,
      outlier.alpha = 1,
      outlier.size = 1,
      fill = NA) +
    scale_color_manual(
      values = c("gray70", "black"),
      name = "",
      labels = c("Glaciated", "Refugium")) +
    ylab("Frequency of polyploids") +
    xlab("Subunit") +
    theme+
    theme(panel.border = element_rect(color = "black", fill = NA)) +
    scale_radius(range = c(0, 3), # point size
                 name = "Number of \nspecies") +
    scale_y_continuous(limits = c(0, .7), expand = expansion(mult = c(0, 0)))
)

### anova for subunit and frequency ----
aov(formula = f_px ~ subunit, data = counts_geo_mountain) |> summary()

### pairwise t-test ----
pwc <-pairwise_t_test(data = counts_geo_mountain,
                      f_px ~ subunit, 
                      p.adjust.method = "holm") |>
  add_xy_position(x = "subunit")

## annotate 
(Fig3_b_ann <-  Fig3_b +
    stat_pvalue_manual(pwc, hide.ns = TRUE, step.increase = .2))

## fix scale
Fig3_a = dst_fpx_mountain_annot +
  scale_y_continuous(limits = c(0, .5), expand = expansion(mult = c(0, 0)))

## annotate and fix scale
Fig3_b_ann <- Fig3_b_ann + 
  scale_y_continuous(limits = c(0, .5), 
                     expand = expansion(mult =c(0, 0))) +
  theme(legend.title = element_text(size = 8),
    legend.text  = element_text(size = 7),
    legend.key.size = unit(.8, "lines"),
    axis.title.y = element_blank(),
    axis.text.y = element_blank())

## final Fig 3 ----
(fig_3 = (Fig3_a+
   annotate(
     "text",
     x = -1,
     y = .45,
     label = "bold(a)",
     parse = T
   ))+ (Fig3_b_ann + annotate(
     "text",
     x = 1,
     y = .45,
     label = "bold(b)",
     parse = T
   )))

# Minimal source for Fig 3a (distance vs frequency with model)
fig3_source <- counts_geo_mountain |>
  dplyr::select(mountain_code, subunit, refugium, dst, f_px, total_spec)

addWorksheet(source_wb, "Fig3")
writeData(source_wb, "Fig3", fig3_source)


## save
ggsave(plot = fig_3,
  create.dir = T,
  filename = paste0("fig_3_current_", ploidy_version, ".pdf"),
  path = paste0(output_dir, "/", ploidy_version),
  device = cairo_pdf,
  height = 75,
  units = "mm",
  width = 150)

## save
# ggsave(
#   plot = fig_3,
#   create.dir = T,
#   filename = paste0("fig_3_current_", ploidy_version, ".png"),
#   path = paste0(output_dir, "/", ploidy_version),
#   height = 75,
#   units = "mm",
#   width = 150,
#   dpi = 300
# )

## richness ----
(
  richness_subunit_ref <- counts_geo_belt |>
    ggplot(aes(
      color = refugium, x = belt, y = total_spec
    )) +
    geom_point(
      size = .5,
      alpha = .6,
      position = position_jitterdodge(jitter.width = .3)
    ) +
    geom_boxplot(
      outliers = F,
      outlier.shape = 4,
      outlier.alpha = 1,
      outlier.size = 1,
      fill = NA
    ) +    # Labels and theme
    theme +
    theme(legend.position = "none") +
    scale_color_manual(
      values = c("gray50", "black"),
      name = "",
      labels = c("Glaciated", "Refugium")
    ) +
    ylab("Number of \nspecies") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .1))) +
    facet_wrap( ~ subunit)
)
(
  diplo_subunit_ref <- counts_geo_belt |>
    ggplot(aes(
      color = refugium, x = belt, y = total_spec - n_px
    )) +
    geom_point(
      size = .5,
      alpha = .6,
      position = position_jitterdodge(jitter.width = .3)
    ) +
    geom_boxplot(
      outliers = F,
      outlier.shape = 4,
      outlier.alpha = 1,
      outlier.size = 1,
      fill = NA
    ) +    # Labels and theme
    theme +
    theme(legend.position = "none") +
    scale_color_manual(
      values = c("gray50", "black"),
      name = "",
      labels = c("Glaciated", "Refugium")
    ) +
    ylab("Number of diploids") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .1))) +
    facet_wrap( ~ subunit)
)


(
  poly_subunit_ref <- counts_geo_belt |>
    ggplot(aes(
      color = refugium, x = belt, y = n_px
    )) +
    geom_point(
      size = .5,
      alpha = .6,
      position = position_jitterdodge(jitter.width = .3)
    ) +
    geom_boxplot(
      outliers = F,
      outlier.shape = 4,
      outlier.alpha = 1,
      outlier.size = 1,
      fill = NA
    ) +    # Labels and theme
    theme +
    theme(legend.position = c(.9, .9)) +
    scale_color_manual(
      values = c("gray50", "black"),
      name = "",
      labels = c("Glaciated", "Refugium")
    ) +
    ylab("Number of polyploids") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .1))) +
    facet_wrap( ~ subunit)
)


# additional sep subunit
(
  richness_subunit <- counts_geo_belt |>
    ggplot(aes(x = belt, y = total_spec)) +
    geom_jitter(size = .5, alpha = .6) +
    geom_boxplot(
      outliers = F,
      outlier.shape = 4,
      outlier.alpha = 1,
      outlier.size = 1,
      fill = NA
    ) +    # Labels and theme
    theme +
    ylab("Number of \nspecies") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .1))) +
    facet_wrap( ~ subunit)
)
(
  diplo_subunit <- counts_geo_belt |>
    ggplot(aes(x = belt, y = total_spec - n_px)) +
    geom_jitter(size = .5, alpha = .6) +
    geom_boxplot(
      outliers = F,
      outlier.shape = 4,
      outlier.alpha = 1,
      outlier.size = 1,
      fill = NA
    ) +    # Labels and theme
    theme +
    ylab("Number of diploids") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .1))) +
    facet_wrap( ~ subunit)
)


(
  poly_subunit <- counts_geo_belt |>
    ggplot(aes(x = belt, y = n_px)) +
    geom_jitter(size = .5, alpha = .6) +
    geom_boxplot(
      outliers = F,
      outlier.shape = 4,
      outlier.alpha = 1,
      outlier.size = 1,
      fill = NA
    ) +    # Labels and theme
    theme +
    ylab("Number of polyploids") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .1))) +
    facet_wrap( ~ subunit)
)


# additional
(
  richness <- counts_geo_belt |>
    ggplot(aes(x = belt, y = total_spec)) +
    geom_jitter(size = .5, alpha = .6) +
    geom_boxplot(
      outliers = F,
      outlier.shape = 4,
      outlier.alpha = 1,
      outlier.size = 1,
      fill = NA
    ) +    # Labels and theme
    theme +
    ylab("Number of \nspecies") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .1)))
)

(
  diplo <- counts_geo_belt |>
    ggplot(aes(x = belt, y = total_spec - n_px)) +
    geom_jitter(size = .5, alpha = .6) +
    geom_boxplot(
      outliers = F,
      outlier.shape = 4,
      outlier.alpha = 1,
      outlier.size = 1,
      fill = NA
    ) +    # Labels and theme
    theme +
    ylab("Number of diploids") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .1)))
)


(
  poly <- counts_geo_belt |>
    ggplot(aes(x = belt, y = n_px)) +
    geom_jitter(alpha = .6, size = .5) +
    geom_boxplot(
      outliers = F,
      outlier.shape = 4,
      outlier.alpha = 1,
      outlier.size = 1,
      fill = NA
    ) +    # Labels and theme
    theme +
    ylab("Number of polyploids") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .1)))
)

richness + diplo + poly

## save
ggsave(
  create.dir = T,
  filename = paste0("div_", ploidy_version, ".pdf"),
  path = paste0(output_dir, "/", ploidy_version),
  device = cairo_pdf,
  height = 100,
  units = "mm",
  width = 210
)

richness_subunit / diplo_subunit / poly_subunit

## save
ggsave(
  create.dir = T,
  filename = paste0("div_subunits_", ploidy_version, ".pdf"),
  device = cairo_pdf,
  path = paste0(output_dir, "/", ploidy_version),
  height = 297,
  units = "mm",
  width = 210
)

richness_subunit_ref / diplo_subunit_ref / poly_subunit_ref

## save
ggsave(
  create.dir = T,
  filename = paste0("div_subunits_refs_", ploidy_version, ".pdf"),
  path = paste0(output_dir, "/", ploidy_version),
  device = cairo_pdf,
  height = 297,
  units = "mm",
  width = 210
)

## supplementary figures 2-3 ----
### a: belts ----
supplement_a <- Fig2_a

### b: distance ----
supplement_b <- Fig3_a +
  annotate("text", x = -1.5, y = .48, label = "bold(b)", parse = T)

### c: effects ----
# formulas
f_glmer <- as.formula("polyploid ~ elevation * distance +  (1|mountain_code)")
f_glmer_fam <- as.formula("polyploid ~ elevation * distance +  (1|mountain_code) + (1|family)")

# models
glmm_pred <- glmer(f_glmer, data = data_model, family = binomial)
glmm_fam_pred <- glmer(f_glmer_fam, data = data_model, family = binomial)



# plot
(supplement_c <- plot_models(
  glmm_pred, glmm_fam_pred,
  show.values = TRUE,
  show.p = TRUE,
  vline.color = "gray",
  colors = effect_cols,
  legend.title = "Random intercepts",
  m.labels = c("Mountain", "Mountain, family"),
  axis.labels = c("Relative elevation: distance", "Distance", "Relative elevation"),
  dot.size = 2, wrap.labels = 4, wrap.legend.title = 4, p.shape = F) +
  ylim(.85, 1.38)+
  theme +
  theme(legend.position = c(.8,.2),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7.5),
        legend.key.size = unit(.8, "lines")) +
  annotate("text", x = 3.45, y = .85, label = "bold(c)", parse = T))

### final fig ----
supplement_a + supplement_b + supplement_c

## save
ggsave(path = paste0(output_dir,"/", ploidy_version),
  filename = paste0("fig_subset_supp_", subset, ".pdf"),  dpi = 600,
  device = cairo_pdf,units = "mm",
  width = 180, height = 67.5
)

## model comparison ----
# for comparisons LRT between models without and with interaction
f_glmer_add <- as.formula("polyploid ~ elevation + distance +  (1|mountain_code)")
f_glmer_fam_add <- as.formula("polyploid ~ elevation + distance +  (1|mountain_code) + (1|family)")
glmm_pred_add <- glmer(f_glmer_add, data = data_model, family = binomial)
glmm_fam_pred_add <- glmer(f_glmer_fam_add, data = data_model, family = binomial)

# comparison sig = interaction model better
(anova(glmm_pred_add, glmm_pred)) 

# comparison family
(anova(glmm_fam_pred_add, glmm_fam_pred))

## piecewise SEM ----
## data for psem
str(data_model)
df <- data_model |> 
  dplyr::select(polyploid, distance, elevation, dispersal_ability, mountain_code, belt) |> 
  mutate(
         dispersal = as.factor(dispersal_ability))
df = na.omit(df)

# convert categorical traits to numeric 0/1 for SEM
df$disp_far <- ifelse(df$dispersal == "far", 1, 0)

### models ----
# GLMM for polyploid frequency with interaction between elevation and distance (change accordingly)
mod_poly <- glmer(polyploid ~ elevation * distance + disp_far +  (1|mountain_code), 
                  family = binomial, data = df)

# GLMs for indirect effects (traits as responses)
mod_disp <- glmer(disp_far ~ elevation + (1|mountain_code), family = binomial, data = df)

# source data
addWorksheet(source_wb, "Fig5")
writeData(source_wb, "Fig5", df)


# Include correlated error between traits
sem_model <- psem(
  mod_poly,
  mod_disp
#  mod_forest
)

summary(mod_poly
        )
# summary
sem_summary <- summary(sem_model, conserve = TRUE)
rsquared(sem_model)
sem_summary

# print(sem_summary)

# direct effects: from mod_poly summary
direct_effects <- coef(summary(mod_poly))
print(direct_effects)

### Fig 5: plotting pSEM ----
(b = sem_summary)
b$coefficients
b$coefficients |> write.csv("edges_all_old.csv")


# define edges for paths
edges <- b$coefficients |> 
  select(Response, Predictor, P.Value, Std.Estimate) |> 
  mutate(
    is_cov = str_starts(Response, "~~"),  # covariance row starts ~~
    from = if_else(is_cov,
      str_remove(Response, "~~"),  # remove ~~
      Predictor), # start from predictor
    to = if_else(is_cov,
      str_remove(Predictor, "~~"),
      Response
    ),  # go to response
    est  = Std.Estimate,
    pval = P.Value,
    color = case_when(pval >= 0.05 ~ "grey75",
      est > 0      ~ "black",
      est < 0      ~ "red"
    ), # define color based on p and neg/pos effect
    lty = if_else(is_cov, "dashed", "solid"))  |> # linetype different for cov.
  mutate(color = ifelse(is_cov,  "black",color)) |> # color different for cov.
  select(from, to, est, pval, color, lty, is_cov)

# nodes with numerical id and labels (fror mapping names to ids)
nodes <- data.frame(
  id = seq_along(unique(c(edges$from, edges$to))),
  label = unique(c(edges$from, edges$to))
)

# Lookup table for mapping names to integer IDs
id_lookup <- setNames(nodes$id, nodes$label)

# mapping
edges_int <- edges %>%
  mutate(
    from_id = id_lookup[from],
    to_id   = id_lookup[to]
  )

# nodes dataframe (DiagrammeR)
nodes_df <- create_node_df(
  n = nrow(nodes),
  label = nodes$label
)

# edges dataframe (DiagrammeR)
edges_df <- create_edge_df(
  from  = edges_int$from_id,
  to    = edges_int$to_id,
  rel   = ifelse(edges_int$is_cov, "cov", "reg"),
  color = edges_int$color,
  style = edges_int$lty,
  label = ifelse(edges_int$pval < 0.05, round(edges_int$est, 2), "")
)

# make graph based on nodes and edges dataframe
graph <- create_graph(
  nodes_df = nodes_df,
  edges_df = edges_df,
  directed = TRUE) |> 
  set_edge_attrs(
    edge_attr = "arrowhead",
    values = ifelse(edges_int$is_cov, "none", "normal")) |> 
    set_edge_attrs(
    edge_attr = "penwidth",
    values = ifelse(edges_int$is_cov, 1, 1 + abs(edges_int$est)*5) )  |> # scale 
  add_global_graph_attrs("layout", "dot", "graph") |> 
  add_global_graph_attrs("rankdir", "LR", "graph") |> 
  set_node_attrs(
    node_attr = "shape",
    values = "box") |> 
  set_node_attrs(
    node_attr = "width",
    values = 1  # small nodes
  ) |> 
  set_node_attrs(
    node_attr = "height",
    values = 0.5
  ) |> 
  set_node_attrs(
    node_attr = "fillcolor",
    values = "white"
  ) |>  # fill boxes
  set_node_attrs(
    node_attr = "color",
    values = "black"
  ) #|> # color boxes
  # set_node_attrs(
  #   node_attr = "style",
  #   values = "filled"
  # )
graph <- create_graph(
  nodes_df = nodes_df,
  edges_df = edges_df,
  directed = TRUE
) |>
  set_edge_attrs(edge_attr = "arrowhead",
                 values = ifelse(edges_int$is_cov, "none", "normal")) |>
  set_edge_attrs(edge_attr = "penwidth",
                 values = ifelse(edges_int$is_cov, 1, 1 + abs(edges_int$est) * 5)) |>

  add_global_graph_attrs("layout",  "dot", "graph") |>
  add_global_graph_attrs("rankdir", "TB",  "graph") |>
  add_global_graph_attrs("ranksep", "0.2", "graph") |>
  add_global_graph_attrs("nodesep", "0.1", "graph") |>
  set_node_attrs(node_attr = "shape",     values = "box") |>
  set_node_attrs(node_attr = "fixedsize", values = "true") |>
  set_node_attrs(node_attr = "width",     values = 1) |>
  set_node_attrs(node_attr = "height",    values = 0.7) |>
  set_node_attrs(node_attr = "style",     values = "filled") |>
  set_node_attrs(node_attr = "fillcolor", values = "white") |>
  set_node_attrs(node_attr = "color",     values = "black")
# render (plot)
render_graph(graph)



## save pdf
export_graph(graph = graph, "paths_d.pdf", file_type = "pdf")

## comparison mixed ploidy species ----
## now matches version in MS with aggregates
## data
df <- data_model |> 
  select(aggregate_grp, ploidy, mountain_code, elevation_relative) |> 
  mutate(ploidy = as.numeric(as.character(ploidy))) 


colnames(df) <- c("species", "ploidy", "mountain", "elevation")
df = na.omit(df)

# Keep only species with both diploids and polyploids
species_counts <- df |>
  group_by(species) |>
  summarise(has_2x = any(ploidy == 2),
            has_poly = any(ploidy >= 4),
            is_mixed = any(length(unique(ploidy)) > 1), 
            .groups="drop")

valid_species <- species_counts |>
  filter(is_mixed) |>
  pull(species)

df_filtered <- df |> filter(species %in% valid_species)

# Remove cytotypes with <10 occurrences
df_filtered_counts <- df_filtered |>
  ungroup() |> 
  group_by(species, ploidy) |>
  mutate(n_occ = n()) |>
  filter(n_occ >= 10) |>
  ungroup()

# Keep only species that still have ≥2 cytotypes after filtering
species_valid <- df_filtered_counts |>
  group_by(species) |>
  summarise(n_cytotypes = n_distinct(ploidy), .groups = "drop") |>
  filter(n_cytotypes > 1) |>
  pull(species)

df_final <- df_filtered_counts |>
  filter(species %in% species_valid) |> 
  select(species, ploidy, elevation)

## Keep two most frequent ploidy levels
df_final_A <- df_final %>%
  group_by(species, ploidy) %>%
  summarise(n_occ = n(), .groups = "drop") %>%
  arrange(species, desc(n_occ), ploidy) %>%
  group_by(species) %>%
  slice_head(n = 2) %>%                # top 2 ploidies
  inner_join(df_final, by = c("species", "ploidy")) %>%
  ungroup()

df_final = df_final_A
df_final$species |> unique()
### mean elevation ----
mean_elev <- df_final |>
  group_by(species, ploidy) |>
  summarise(mean_elev = mean(elevation), .groups="drop")

# Differences (corrected & original)
diff_elev <- mean_elev |>
  group_by(species) |>
  summarise(diff = mean_elev[ploidy == max(ploidy)] - 
              mean_elev[ploidy == min(ploidy)])

### Monte Carlo permutatation----
set.seed(123)
n_perm <- 10000
perm_diffs <- replicate(n_perm, {
  df_final |>
    group_by(species) |>
    mutate(ploidy_perm = sample(ploidy)) |>
    group_by(species, ploidy_perm) |>
    summarise(mean_elev = mean(elevation), .groups="drop") |>
    group_by(species) |>
    summarise(diff = mean_elev[ploidy_perm == max(ploidy_perm)] - 
                mean_elev[ploidy_perm == min(ploidy_perm)]) |>
    pull(diff) |>
    mean()
})

### test obs vs permutation ----
obs_diff <- mean(diff_elev$diff)
p_val <- mean(perm_diffs >= obs_diff)

obs_diff
p_val

# Minimal source for Fig 4a (GLMM effects) — model input rows used
fig4_source <- data_model |>
  dplyr::select(polyploid, elevation, distance, mountain_code, family) |>
  tidyr::drop_na()

addWorksheet(source_wb, "Fig4")
writeData(source_wb, "Fig4", fig4_source)

### Fig 4a: effects ----
# same as supplementc c but change label
(p_left <- plot_models(
      glmm_pred, glmm_fam_pred,
      show.values = TRUE,
      show.p = TRUE,
      p.shape = FALSE,
      vline.color = "gray",
      colors = effect_cols,
      legend.title = "Random intercepts",
      m.labels = c("Mountain", "Mountain, family"),
      axis.labels = c("Relative elevation: distance", "Distance", "Relative elevation"),
      std.response = FALSE,
      dot.size = 2, wrap.labels = 4, wrap.legend.title = 4
    ) +
      ylim(.85, 1.2) +
      theme +
      theme(legend.position = c(.8,.2),
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 7.5),
            legend.key.size = unit(.8, "lines")) +
      annotate("text", x = 3.45, y = .85, label = "bold(a)", parse = T))
    
    
### Fig 4b ----
# relative elevation form model predictions
mu_elev <- mean(data_model$elevation_relative, na.rm = TRUE)
sd_elev <- sd(data_model$elevation_relative, na.rm = TRUE)
elev_raw <- c(-500, -250, 0, 250, 500)
elev_z <- (elev_raw - mu_elev) / sd_elev


fig4b <- plot_model(glmm_pred,
                 type = "pred",
                 terms = c("distance", paste0("elevation[", paste(elev_z, collapse = ","), "]")),
                 show.values = TRUE,
                 title = "",
                 axis.title = c("Cost-distance [scaled and centered]",
                                "Polyploidy probability (predicted)")) +
  aes(color = group, fill = group) +
  scale_fill_manual(values = pal_belt) +
  scale_color_manual(values = pal_belt, name = "Relative elevation", labels = elev_raw) +
  theme +
  ylim(.1, .5) +
  scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(.1, .1)))+
  theme(
    legend.position = c(0.3, 0.7),
    # axis.title.y = element_blank(),
    #axis.text.y = element_blank(),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7.5),
    legend.key.size = unit(.8, "lines")
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 1), reverse = TRUE),
    color = guide_legend(override.aes = list(size = 1), reverse = TRUE)
  ) +
  annotate("text", x = -1, y = .4,   label = "bold(b)", parse = T)

## final Fig 4 ----
(final_plot <- p_left|fig4b)

## save
ggsave(
  filename = paste0("Fig_4_glmm.pdf"),
  device = cairo_pdf,
  plot = final_plot,
  dpi = 600, width = 150, height = 80, units = "mm",
  path = paste0(output_dir, "/", ploidy_version)
)
(final_plot <- p_left/fig4b)
ggsave(
  filename = paste0("Fig_4_glmm_v.pdf"),
  device = cairo_pdf,
  plot = final_plot,
  dpi = 600, width = 80, height = 140, units = "mm",
  path = paste0(output_dir, "/", ploidy_version)
)



## check
# n_samples_close
# n_samples_far
# 
# n_species_close
# n_species_far

saveWorkbook(
  source_wb,
  file.path(output_dir, "Source_Data_final.xlsx"),
  overwrite = TRUE
)

