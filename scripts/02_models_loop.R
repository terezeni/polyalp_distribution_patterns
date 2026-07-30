## large GLMMs/GLM loops all
## TZ last update 01/2026

# SETUP ----
library(tidyverse)
library(rio)
library(janitor)
library(lme4)
library(patchwork)
library(ggpubr)
library(car)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(ggeffects)
library(dotwhisker)
library(visreg)
library(performance)
library(MuMIn)
library(broom)
library(broom.mixed)
library(purrr)
library(extrafont)
loadfonts(device = "win")

## THEME, COLORS ----
theme <- theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA))

pal_belt <- c("#006147", "#80A552","#E7D47A", "#D09D4D", "#B56721")
new_cols <- c("#D55E00", "#7570b3")  
mycols <- c("#1b9e77","#d95f02","#7570b3")


## PATHS ----
project_dir <- here::here()

## paths
data_path <- file.path(project_dir, "data")
output_path <- file.path(project_dir, "model_output")
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
effects_dir <- paste0(output_path, "/effects/")
plots_dir <- paste0(output_path, "/plots/")
if (!dir.exists(effects_dir)) dir.create(effects_dir, recursive = TRUE)
mxsl <- paste0(plots_dir, "/glm_taxa/mxsl/")
if (!dir.exists(mxsl)) dir.create(mxsl, recursive = TRUE)
gen <- paste0(plots_dir, "/glm_taxa/genus/")
if (!dir.exists(gen)) dir.create(gen, recursive = TRUE)

# DATA ----
data <- import(paste0(data_path,
                      "/pro/Source_Data_rev.xlsx"), which = 2) |> 
  clean_names() |> 
  mutate(ref = paste0(mountain_code, "-", belt),
         mountain = as.factor(mountain_code),
         elevation = ifelse(belt == "N", NA, elevation)) # set elevation to NA fro "N" belt samples

# elevation relative to timberline
data <- data |> 
  mutate(elevation_orig = elevation,
         elevation_relative = elevation_orig - treeline,
         elevation = elevation_relative)

# Scale numeric vars
data <- data |>
  mutate(across(c(dist_v0, dist_v_a, dist_v_b, elevation), 
                ~scale(.)[,1], .names="{.col}")) |>
  mutate(
    family = as.factor(family),
    polyploid_orig = ifelse(ploidy_orig != 2, 1, 0),
    polyploid_alt = ifelse(ploidy_alt != 2, 1, 0)
  )


## filter subsets ----
# # L_complete
# data = data |>
#   filter(light_preference == "light")
#filter(light_preference == "shade")

# # dispersal
# data = data |>
#   # filter(dispersal == "far")
#   filter(dispersal == "close")

# MODELLING ----

# variables for models ----
# ALL
response_vars <- c("polyploid_orig", "polyploid_alt")
dist_vars <- c("dist_v0", "dist_v_a", "dist_v_b")

# # quicker for subsets only one ploidy and one distance
#response_vars <- c("polyploid_orig")
#dist_vars <- c("dist_v0")

## model formulas ----
# all
model_formulas <- list(
  elevation_only  = function(dist) paste0("~ elevation"),
  distance_only   = function(dist) paste0("~ ", dist),
  additive = function(dist) paste0("~ elevation + ", dist),
  interaction     = function(dist) paste0("~ elevation * ", dist)
)

# taxa (genus or sl mx ploidy)
model_formulas_taxa <- list(
  elevation_only  = "~ elevation",
  distance_only = "~ distance"
)

## initialize lists for loop for models ----
# all
results_list <- list()
glmer_plot_list <- list()

# mixed ploidy genera
results_list_g <- list()
glm_plot_list_g <- list()

# mixed ploidy species/aggregates
results_list_mxsl <- list()
glm_plot_list_mxsl <- list()

## random structure options ----
random_structures <- list(
  family_mountain = "(1|family) + (1|mountain)",
  mountain_only   = "(1|mountain)"
)

## loop fitting models ----
for (resp in response_vars) {
  for (dist in dist_vars) {
    for (model_name in names(model_formulas)) {
      for (rand_name in names(random_structures)) {
        
        print(paste("Fitting:", resp, dist, model_name, rand_name))
        
        # formula string
        f_glmer <- as.formula(
          paste0(resp, " ", model_formulas[[model_name]](dist), " + ", random_structures[[rand_name]])
        )
        
        # subset to complete cases
        subset_m <- data[complete.cases(data[, all.vars(f_glmer)]), ]
        
        # GLMM
        glmer_model <- try(glmer(f_glmer, data = subset_m, family = binomial), silent = TRUE)
        
        # store results (some checks)
        if (!inherits(glmer_model, "try-error")) {
          model_id <- paste(resp, dist, model_name, rand_name, sep = "_")
          
          glmer_plot_list[[model_id]] <- glmer_model
          r2_vals <- performance::r2(glmer_model)
          
          results_list[[model_id]] <-
            tidy(glmer_model, effects = "fixed", conf.int = TRUE) |>
            mutate(
              response = resp,
              model_type = "glmer",
              dist_var = dist,
              model_name = model_name,
              random_structure = rand_name,
              AIC = AIC(glmer_model),
              R2_marginal = r2_vals$R2_marginal,
              R2_conditional = r2_vals$R2_conditional,
              sample_size = nobs(glmer_model)
            )
        } else {
          warning(paste("Model failed for:", resp, dist, model_name, rand_name))
        }
      }
    }
  }
}

## export results table ----
supp_results_all <- bind_rows(results_list)
export(supp_results_all, 
            paste0(output_path, 
                   "/tables/supplement_glmer_results_combined.xlsx"))

# only one ploidy definition and one distance
response_vars <- c("polyploid_orig")
dist_vars <- c("dist_v0")

# MODELLING MX TAXA ----
# grouping factor for aggregates (if no aggregate abbreviation kept)
mxsl_groups <- unique(data$aggregate_grp) |> sort()
data$distance <- data$dist_v0 # not other for mxpl


## loop ----
for (taxon in mxsl_groups) {
  subset_t <- data[data$aggregate_grp == taxon,]
  
  for (resp in response_vars) {
    
    for (model_name in names(model_formulas_taxa)) {
      f_fixed <- as.formula(paste0(resp, " ", model_formulas_taxa[[model_name]]))
      subset_m <- subset_t[complete.cases(subset_t[, all.vars(f_fixed)]), ]
     ### check data ----     
      # response variable is not uniform (mixed ploidy)
      if (length(unique(subset_m[[resp]])) == 1) {
        print(paste("Skipping iteration for", resp, "because the response is uniform."))
        next  # Skip this iteration and go to the next one
      }
      
      # both 0 and 1 have at least 10 samples
      response_counts <- table(subset_m[[resp]])
      if (any(response_counts < 10)) {
        print(paste("Skipping iteration for", resp, "because one of the response classes has fewer than 10 samples."))
        next  
      }
      
      
      # separation
      predictor <- subset_m[[all.vars(f_fixed)[2]]]  # predictor directly from the formula
      
      # sort the predictor and response to detect separation
      sorted_data <- data.frame(predictor = predictor, response = subset_m[[resp]])
      sorted_data <- sorted_data[order(sorted_data$predictor), ]
      
      # any point where all 0s come before all 1s (or opposite)
      response_sorted <- sorted_data$response
      separation_detected <- all(cummax(response_sorted) == response_sorted) || all(cummin(response_sorted) == response_sorted)
      
      if (separation_detected) {
        print(paste("Skipping iteration for", resp, "because separation detected in predictor:", all.vars(f_fixed)[2]))
        next  
      }
      
      ## GLM ----
      glm_model <- try(glm(f_fixed, data = subset_m, family = binomial), silent = TRUE)
      if (!inherits(glm_model, "try-error")) {
        glm_plot_list_mxsl[[paste(resp, taxon,  model_name, sep = "_")]] <- glm_model
        results_list_mxsl[[paste0(resp, taxon, "_glm_", "_", model_name)]] <-
          tidy(glm_model, effects = "fixed", conf.int = TRUE) |>
          mutate(response = resp,
                 taxon = taxon,
                 model_type = "glm",
                 model_name = model_name,
                 AIC = AIC(glm_model),
                 R2_marginal = performance::r2(glm_model)$R2,
                 R2_conditional = NA,
                 sample_size = nobs(glm_model))
      }
      
      
    } 
    
  }
  
}

## export results ----
supp_results_mxsl <- bind_rows(results_list_mxsl)
export(supp_results_mxsl, 
            paste0(output_path, 
                   "/tables/supplement_glmer_mxsl_results_combined.xlsx"))

 
## filter models based on response, taxon ----
filter_models <- function(model_list, resp, taxon, pred = NULL) {
  # base pattern: resp + taxon
  pattern <- paste0("^", resp, "_", taxon, "_")
  
  # in case append predictor to pattern
  if (!is.null(pred)) {
    pattern <- paste0(pattern, pred)
  }
  
  model_list[grepl(pattern, names(model_list))]
}

## loop over predictors for plotting ----
response_vars <- c("polyploid_orig")
predictors <- c("elevation","distance")  # loop over predictors

for (pred in predictors) {
  # one multpage pdf per predictor
  pdf(file = paste0(mxsl,"/mxsl_glm_pred_", pred, ".pdf"), 
      width = 8.27, height = 11.69) 
  
  plot_list <- list() 
  for (taxon in mxsl_groups) {
    for (resp in response_vars) {
      
      # filter models for this response, taxon, and predictor
      glm_sub <- filter_models(glm_plot_list_mxsl, resp, taxon, pred)
      
      if (length(glm_sub) > 0) {
        # subset data for taxon
        subset_t <- data[data$aggregate_grp == taxon,]
        
        # extract model summary
        glm_summary <- summary(glm_sub[[1]])
        coefficients <- glm_summary$coefficients
        odds_ratios <- exp(coefficients[, 1])
        p_values <- coefficients[, 4]
        
        # predictor
        model_pred <- setdiff(rownames(coefficients), "(Intercept)")[1]
        p_value <- p_values[model_pred]
        odds_ratio <- odds_ratios[model_pred]
        line_color <- ifelse(p_value < 0.05,
                             ifelse(odds_ratio > 1, 
                                    "#B56721", "#006147"), "grey50")
        
        annotation_text <- paste("OR = ", round(odds_ratio, 2), 
                                 "\np-value = ", round(p_value, 3))
        
        p_glm <- plot_model(glm_sub[[1]], 
                            colors = line_color,
                            type = "pred",        
                            terms = model_pred,         
                            axis.title = c(paste(model_pred, "[centered and scaled]"), 
                                           "probability of being polyploid")) +
          geom_point(data = subset_t, aes_string(x = model_pred, y = resp), alpha = 0.6, size = 0.5) +
          theme+
          ggtitle(paste("taxon:", taxon, "\nresponse:", resp)) +  # Title with taxon and response
          theme(plot.title = element_text(size = 10), axis.title = element_text(size = 10)) +
          ggtitle(paste("taxon:", taxon, "\nresponse:", resp)) +
          ylim(0,1) +
          annotate("text", y = 0.7, 
                   x = quantile(subset_t[[model_pred]], probs = 0.2, na.rm = TRUE), 
                   label = annotation_text, size = 3, color = line_color, hjust = 0, vjust = 0)
        
        plot_list[[length(plot_list) + 1]] <- p_glm
      } else {
        message(paste("Skipping plotting for taxon:", taxon, "Response:", resp, "Model is not valid."))
      }
    }
  }
  
  # Print plots in pages
  plots_per_page <- 12
  for (i in seq(1, length(plot_list), by = plots_per_page)) {
    current_plots <- plot_list[i:min(i + plots_per_page - 1, length(plot_list))]
    combined_plot <- wrap_plots(plotlist = current_plots, ncol = 3, nrow = 4, byrow = FALSE)
    print(combined_plot)
  }
  
  dev.off()  # Close PDF for this predictor
}

# MODELLING GENERA ----
# get species list for loop, min 20 observations
genus_valid <- 
  data |>
  group_by(genus) |> 
  summarise(n = n()) |> 
  filter(n>20)

## loop ----
for (taxon in genus_valid$genus) {
  subset_t <- data[data$genus == taxon,]
  
  for (resp in response_vars) {
    
    for (model_name in names(model_formulas_taxa)) {
      print("started")
      f_fixed <- as.formula(paste0(resp, " ", model_formulas_taxa[[model_name]]))
      subset_m <- subset_t[complete.cases(subset_t[, all.vars(f_fixed)]), ]
### check data ----
      # Check if response variable is uniform (i.e., all values are the same)
      if (length(unique(subset_m[[resp]])) == 1) {
        print(paste("Skipping iteration for", resp, "because the response is uniform."))
        next  # Skip this iteration and go to the next one
      }
      # Check if both 0 and 1 have at least 10 samples
      response_counts <- table(subset_m[[resp]])
      if (any(response_counts < 10)) {
        print(paste("Skipping iteration for", resp, "because one of the response classes has fewer than 10 samples."))
        next  # Skip this iteration if the condition is met
      }
      # separation
      predictor <- subset_m[[all.vars(f_fixed)[2]]]  # Extract the predictor directly from the formula
      
      # sort the predictor and response to detect separation
      sorted_data <- data.frame(predictor = predictor, response = subset_m[[resp]])
      sorted_data <- sorted_data[order(sorted_data$predictor), ]
      
      # Check if there is any point where all 0s come before all 1s, or vice versa
      response_sorted <- sorted_data$response
      separation_detected <- all(cummax(response_sorted) == response_sorted) || all(cummin(response_sorted) == response_sorted)
      
      if (separation_detected) {
        print(paste("Skipping iteration for", resp, "because separation detected in predictor:", all.vars(f_fixed)[2]))
        next  # Skip this iteration if separation is detected
      }
      
      ## GLM ----
      glm_model <- try(glm(f_fixed, data = subset_m, family = binomial), silent = TRUE)
      if (!inherits(glm_model, "try-error")) {
        glm_plot_list_g[[paste(resp, taxon,  model_name, sep = "_")]] <- glm_model
        results_list_g[[paste0(resp, taxon, "_glm_", "_", model_name)]] <-
          tidy(glm_model, effects = "fixed", conf.int = TRUE) |>
          mutate(response = resp,
                 taxon = taxon,
                 model_type = "glm",
                 model_name = model_name,
                 AIC = AIC(glm_model),
                 R2_marginal = performance::r2(glm_model)$R2,
                 R2_conditional = NA,
                 sample_size = nobs(glm_model))
      }
      
      
    } # formulas
    
  }
  
}


## export results ----
supp_results_g <- bind_rows(results_list_g)
export(supp_results_g, 
            paste0(output_path, 
                   "/tables/supplement_glm_genus_results_combined.xlsx"))

# response and predictors
response_vars <- c("polyploid_orig")
predictors <- c("elevation","distance")  # loop over predictors

## loop over predictors for plotting ---- 
for (pred in predictors) {
  
  # separate PDF for each predictor
  pdf(file = paste0(gen,"/mxgenus_glm_prediction_", pred, ".pdf"), 
      width = 8.27, height = 11.69)  # A4 size
  
  plot_list <- list()  # initialize plot list for current predictor
  
  for (taxon in genus_valid$genus) {
    for (resp in response_vars) {
      
      # filter models for this response, taxon, and predictor
      glm_sub <- filter_models(glm_plot_list_g, resp, taxon, pred)
      
      if (length(glm_sub) > 0 && inherits(glm_sub[[1]], "glm")) {
        subset_t <- data[data$genus == taxon,]
        
        # Extract model summary
        glm_summary <- summary(glm_sub[[1]])
        coefficients <- glm_summary$coefficients
        odds_ratios <- exp(coefficients[, 1])
        p_values <- coefficients[, 4]
        
        # Use the predictor 
        model_pred <- setdiff(rownames(coefficients), "(Intercept)")[1]
        p_value <- p_values[model_pred]
        odds_ratio <- odds_ratios[model_pred]
        line_color <- ifelse(p_value < 0.05, ifelse(odds_ratio > 1, "#B56721", "#006147"), "grey50")
        
        annotation_text <- paste("OR = ", round(odds_ratio, 2), "\np-value = ", round(p_value, 3))
        
        p_glm <- plot_model(glm_sub[[1]], 
                            colors = line_color,
                            type = "pred",        
                            terms = model_pred,         
                            axis.title = c(paste(model_pred, "[centered and scaled]"), 
                                           "probability of being polyploid")) +
          theme +  
          geom_point(data = subset_t, aes_string(x = model_pred, y = resp), alpha = 0.6, size = 0.5) +theme+
          ggtitle(paste("taxon:", taxon, "\nresponse:", resp)) +  # Title with taxon and response
          theme(plot.title = element_text(size = 10), axis.title = element_text(size = 10)) +
          ggtitle(paste("taxon:", taxon, "\nresponse:", resp)) +
          ylim(0,1) +
          annotate("text", y = 0.7, 
                   x = quantile(subset_t[[model_pred]], probs = 0.2, na.rm = TRUE), 
                   label = annotation_text, size = 3, color = line_color, hjust = 0, vjust = 0)
        
        plot_list[[length(plot_list) + 1]] <- p_glm
      } else {
        message(paste("Skipping plotting for taxon:", taxon, "Response:", resp, "Model is not valid."))
      }
    }
  }
  
  # Print plots in pages
  plots_per_page <- 12
  for (i in seq(1, length(plot_list), by = plots_per_page)) {
    current_plots <- plot_list[i:min(i + plots_per_page - 1, length(plot_list))]
    combined_plot <- wrap_plots(plotlist = current_plots, ncol = 3, nrow = 4, byrow = FALSE)
    print(combined_plot)
  }
  
  dev.off()  # Close PDF for this predictor
}

# add taxa information
## use info from aggregate_grp
mxsl_info = data |> 
  select(aggregate_grp,dispersal_ability, family) |> 
  unique()

# elevation table
supp_results_mxslt_e <- left_join(supp_results_mxsl, mxsl_info, by = c("taxon" = "aggregate_grp")) |> 
  filter(term == "elevation")

# distance table
supp_results_mxslt_d <- left_join(supp_results_mxsl, mxsl_info, by = c("taxon" = "aggregate_grp")) |> 
  filter(term == "distance")

# line color based on p-value and odds ratio
supp_results_mxslt_e$color = ifelse(supp_results_mxslt_e$p.value < 0.05, 
                                    ifelse(exp(supp_results_mxslt_e$estimate) > 1, 
                                           "#B56721", 
                                           "#006147"),
                                    "grey50")
supp_results_mxslt_d$color = ifelse(supp_results_mxslt_d$p.value < 0.05,
                                    ifelse(exp(supp_results_mxslt_d$estimate) > 1, 
                                           "#B56721", 
                                           "#006147"),
                                    "grey50")

# plot by group slopes or odd ratios
or_sl_e <- ggplot(data = supp_results_mxslt_e, aes(x = dispersal_ability, y = exp(estimate)))+
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = .4, color = "gray50")+
  geom_violin(trim = T, scale= "count", width = 1, linewidth = .1, fill = NA)+
  geom_jitter(aes(color = color),  alpha = .4, stroke = NA, height = 0, width = .1)+
  geom_boxplot(outliers = F, width = .1, linewidth =.2, fill = NA)+
  scale_colour_identity()+#
  ylab("OR elevation")+
  xlab("dispersal")+
  theme


# plot by group slopes or odd ratios
or_sl_d <- ggplot(data = supp_results_mxslt_d, aes(x = dispersal_ability, y = exp(estimate)))+
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = .4, color = "gray50")+
  geom_violin(trim = T, scale= "count", width = 1, linewidth = .1, fill = NA)+
  geom_jitter(aes(color = color),  alpha = .4, stroke = NA, height = 0, width = .1)+
  geom_boxplot(outliers = F, width = .1, linewidth =.2, fill = NA)+
  scale_colour_identity()+#
  ylab("OR distance")+
  xlab("dispersal")+
  theme

or_sl_e + or_sl_d + plot_annotation(tag_levels = 'a')

ggsave(filename = "effects_mxpl_groups.pdf", width = 6, height = 3, path = effects_dir)
ggsave(filename = "effects_mxpl_groups.png", width = 6, height = 3, path = effects_dir)


## test differences from no effect
supp_results_mxslt_e %>%
  group_by(dispersal_ability) %>%
  summarise(
    n = n(),
    median_or = median(exp(estimate), na.rm = TRUE),
    p_wilcox = wilcox.test(exp(estimate), mu = 1, alternative = "two.sided")$p.value
  )

## test differences from no effect
supp_results_mxslt_d %>%
  group_by(dispersal_ability) %>%
  summarise(
    n = n(),
    median_or = median(exp(estimate), na.rm = TRUE),
    p_wilcox = wilcox.test(exp(estimate), mu = 1, alternative = "two.sided")$p.value
  )


supp_results_mxslt_e %>%
  group_by(dispersal_ability) |> write_clip()
supp_results_mxslt_d %>%
  group_by(dispersal_ability) |> write_clip()