# This script aggregates the fatty acid data and performs
# several univariate and multivariate analyses aimed at
# relating fatty acid community compositions with sewage
# indicators in Lake Baikal.

library(tidyverse)
library(viridis)
library(viridisLite)
library(vegan)
library(ggpubr)
library(ggrepel)
library(ggtext)
library(fs)
library(broom)
library(MixSIAR)

# Check if figures directory exists 
# If not, create the figures directory
sub_dir <- "figures"
output_dir <- file.path(here::here(), sub_dir)

if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Dir 'figures' already exists!")
}

# Here, we will define a function to perform a permutational analysis of the 
# of the data. In general, this function permutes the response variable 5,000 times, 
# creates a linear model with the permuted data, and extracts the p-value and R-squared 
# value from the model. Once the interation is finished, the 5,000 p-values and R-squared
# values are plotted on a histogram and compared to the model generated with non-permuted
# data.

permute_data_analytics <- function(data, metric, full_model, metric_plot_title, predictor){
  for(i in 1:5000){
    
    # First permute the response variable. The variable is supplied by the user.
    permuted_data <- data %>%
      select(paste(metric)) %>%
      rename("permuted_metric" = paste(metric)) %>%
      sample_frac(size = 1) %>%
      cbind(., data)
    
    # This step takes the permuted data and builds the appropriate model with the 
    # desired predictor. 
    if(predictor == "idw_pop"){
      permuted_model <- lm((permuted_metric) ~ log10(distance_weighted_population),
                           data = permuted_data)
    }else if(predictor == "ppcp"){
      permuted_model <- lm((permuted_metric) ~ log10(ppcp_sum),
                           data = permuted_data)
    }else{
      return(cat("You entered an invalid predictor"))
    }
    
    # If this iteration is the first, then function creates two repos for the 
    # p-value and r-squared values. Note that this step requires the broom package
    # be installed. 
    if(i == 1){
      glance_repo <- glance(permuted_model)
      tidy_repo <- tidy(permuted_model)
    } else {
      glance_repo <- rbind(glance_repo, glance(permuted_model))
      tidy_repo <- rbind(tidy_repo, tidy(permuted_model)) %>%
        filter(term != "(Intercept)")
    }
  }
  
  # This step removes all intercept coefficients from the repo. 
  tidy_full_model <- tidy(full_model) %>%
    filter(term != "(Intercept)")
  
  # This step plots the p-value histogram figure. 
  permuted_plot_pval <- ggplot() +
    geom_histogram(data = tidy_repo, 
                   mapping = aes(p.value), 
                   fill = "white",
                   color = "black") +
    geom_vline(data = tidy_full_model,
               mapping = aes(xintercept = p.value),
               linetype = "dashed", size = 2, color = viridis(30)[15]) +
    ggtitle(paste(metric_plot_title, " (", 
                  round(length(tidy_repo$p.value[tidy_repo$p.value <= tidy_full_model$p.value])/5000, 3)*100, 
                  "% of models had lower p-value)", sep = "")) +
    ylab("Frequency") +
    xlab("P-Value") +
    theme_minimal() +
    theme(plot.title = ggtext::element_markdown())
  
  # This code chunk plots the histogram of the R-squared values.
  permuted_plot_rsquared <- ggplot() +
    geom_histogram(data = glance_repo, 
                   mapping = aes(r.squared), 
                   fill = "white",
                   color = "black") +
    geom_vline(data = glance(full_model),
               mapping = aes(xintercept = r.squared),
               linetype = "dashed", size = 2, color = viridis(30)[15]) +
    ggtitle(paste(metric_plot_title, " (", 
                  round(length(glance_repo$r.squared[glance_repo$r.squared >= glance(full_model)$r.squared])/5000, 3)*100, 
                  "% of models had higher R<sup>2</sup>)", sep = "")) +
    ylab("") +
    xlab("R Squared") +
    theme_minimal() +
    theme(plot.title = ggtext::element_markdown())
  
  # The two resulting figures are returned as a list. 
  return(list(permuted_plot_pval, permuted_plot_rsquared))
}


# 1. Load the data --------------------------------------------------------

fatty_acid <- read.csv(file = "../cleaned_data/fatty_acid.csv",
                       header = TRUE, stringsAsFactors = FALSE)

site_information <- read.csv(file = "../cleaned_data/site_information.csv",
                     header = TRUE, stringsAsFactors = FALSE)

distance <- read.csv(file = "../cleaned_data/distance_weighted_population_metrics.csv",
                     header = TRUE, stringsAsFactors = FALSE)

nutrients <- read.csv(file = "../cleaned_data/nutrients.csv",
                      header = TRUE, stringsAsFactors = FALSE)

site_info_dist <- full_join(x = site_information, y = distance)

ppcp <- read.csv(file = "../cleaned_data/ppcp.csv",
                 header = TRUE, stringsAsFactors = FALSE)

ppcp_site_info_dist <- full_join(x = ppcp, y = site_info_dist, by = c("site"))

fatty_acid_ppcp_site_info_dist <- inner_join(x = fatty_acid, y = ppcp_site_info_dist,
                                        by = c("site")) %>%
  unite(taxon, c("Genus", "Species"), remove = FALSE)

# 1.5 Create folder for tables --------------------------------------------

sub_dir <- "tables"
output_dir <- file.path(here::here(), sub_dir)

if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Dir 'tables' already exists!")
}

# 2. Overall Fatty Acid Analysis ------------------------------------------


# 2.1 Whole fatty acid spectrum -------------------------------------------

# Create a dataframe of the entire fatty acid spectrum
fatty_acid_whole_wide <- fatty_acid %>%
  filter(Genus != "Hyalella") %>%
  select(site:c24_0) %>%
  gather(key = fatty_acid, value = concentration, c12_0:c24_0) %>%
  group_by(site, Genus, Species) %>%
  mutate(total_fatty_acid = sum(concentration),
         prop_fatty_acid = concentration / total_fatty_acid) %>%
  select(-concentration, - total_fatty_acid) %>%
  spread(key = fatty_acid, value = prop_fatty_acid) %>%
  unite(taxon, c("Genus", "Species"), sep = "_", remove = FALSE)

# Define amphipod taxa in a vector to italicize
amphipods <- c("Eulimnogammarus verrucosus", "Eulimnogammarus vittatus",
               "Pallasea cancellus", "Eulimnogammarus cyaneus")

# Identify mean, variance, and coefficient of variation across all
# sites for each taxonomic grouping
mean_var <- fatty_acid_whole_wide %>%
  pivot_longer(cols = a_15_0:i_17_0, names_to = "fa") %>%
  group_by(fa) %>%
  summarize(
    Mean = mean(value),
    Variance = var(value)
  ) %>%
  mutate(Var_Mean_Ratio = Variance / Mean) %>%
  arrange(desc(Var_Mean_Ratio))

write.csv(x = mean_var, file = "../tables/cv_all_fa.csv", row.names = FALSE)

# Perform NMDS
whole_fatty_acid_metaMDS <- metaMDS(comm = fatty_acid_whole_wide[, 5:63],
                                    distance = "bray", k = 2, try = 100)
whole_fatty_acid_metaMDS

# Extract scores for each site and add additional data for analysis and plotting
data_scores <- as.data.frame(scores(whole_fatty_acid_metaMDS)) %>%
  mutate(site = fatty_acid_whole_wide$site,
         taxon = fatty_acid_whole_wide$taxon,
         taxon = gsub(pattern = "_", replacement = " ", x = taxon),
         taxon = gsub(pattern = "NA", replacement = "", x = taxon),
         taxon = gsub(pattern = "Draparnaldia", replacement = "Draparnaldia spp.", x = taxon),
         taxon = trimws(taxon))

# Pull species scores from NMDS
species_scores <- as.data.frame(scores(x = whole_fatty_acid_metaMDS, display = "species"))
species_scores$species <- rownames(species_scores)

# Create plot
# This figure is associated with figure S1 in the associated manuscript
nmds <- ggplot() +
  geom_point(data = data_scores %>%
             mutate(taxon = ifelse(taxon == "Draparnaldia spp.",
                              yes = "*Draparnaldia* spp.",
                              no = taxon), 
                    taxon = ifelse(taxon %in% amphipods,
                            yes = paste0("*", taxon, "*"),
                            no = taxon)),
             mapping = aes(x = NMDS1, y = NMDS2, fill = taxon),
             size = 10, alpha = .6, shape = 21, stroke = 2, color = "grey60") +
  scale_fill_viridis_d(option = "plasma") +
  geom_text_repel(data = species_scores %>%
                    filter(species %in% c("c18_3w3", "c18_1w9", "c18_2w6", "c16_0", 
                                          "c14_0", "c20_5w3", "c16_1w7")) %>%
                    mutate(NMDS2 = ifelse(test = species == "c14_0", yes = NMDS2+0.04, 
                                          no = NMDS2),
                           NMDS2 = ifelse(test = species == "c16_1w7", yes = NMDS2-0.02, 
                                          no = NMDS2),
                           species = gsub(pattern = "c", replacement = "", 
                                          x = species),
                           species = gsub(pattern = "_", replacement = ":", 
                                          x = species),
                           species = gsub(pattern = "w", replacement = "\U03C9", 
                                          x = species)), 
                  aes(x = NMDS1, y = NMDS2, label = species),
                  size = 10) +
  annotate("label", x = 0.4, y = 0.45,
           label = paste("Stress: ", round(whole_fatty_acid_metaMDS$stress, digits = 3)),
           size = 10) +
  theme_minimal() +
  theme(legend.position = "bottom",
        title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_markdown(size = 12))
nmds

ggsave(filename = "all_species_all_FA.png", plot = nmds, device = "png",
       path = "../figures/", width = 12, height = 10, units = "in", dpi = 300)


# 2.2 Redo analysis with only essential fatty acids -----------------------

# Create a dataframe of only the essential fatty acids
fatty_acid_essential_wide <- fatty_acid %>%
  filter(Genus != "Hyalella") %>%
  select(site, Genus, Species, c18_3w3, c18_4w3, c20_5w3, c22_5w3, c22_6w3,
         c18_2w6, c18_2w6t, c20_4w6) %>%
  gather(key = fatty_acid, value = concentration, c18_3w3:c20_4w6) %>%
  group_by(site, Genus, Species) %>%
  mutate(total_fatty_acid = sum(concentration),
         prop_fatty_acid = concentration / total_fatty_acid) %>%
  select(-concentration, - total_fatty_acid) %>%
  spread(key = fatty_acid, value = prop_fatty_acid) %>%
  unite(col = taxon, c("Genus", "Species"), sep = "_", remove = FALSE)

# Identify mean, variance, and coefficient of variation across all
# sites for each taxonomic grouping
mean_var <- fatty_acid_essential_wide %>%
  pivot_longer(cols = c18_2w6:c22_6w3, names_to = "fa") %>%
  group_by(fa) %>%
  summarize(
    Mean = mean(value),
    Variance = var(value)
  ) %>%
  mutate(Var_Mean_Ratio = Variance / Mean) %>%
  arrange(desc(Var_Mean_Ratio))

write.csv(x = mean_var[order(mean_var$Var_Mean_Ratio, decreasing = TRUE), ], 
          file = "../tables/cv_efa.csv", row.names = FALSE)


# Perform NMDS
essential_fatty_acid_metaMDS <- metaMDS(comm = fatty_acid_essential_wide[, 5:12],
                                        distance = "bray", k = 2, try = 100)
essential_fatty_acid_metaMDS

# Extract data scores for plotting and analysis of the NMDS
data_scores <- as.data.frame(scores(essential_fatty_acid_metaMDS)) %>%
  mutate(site = fatty_acid_essential_wide$site,
         taxon = fatty_acid_essential_wide$taxon,
         taxon = gsub(pattern = "_", replacement = " ", x = taxon),
         taxon = gsub(pattern = "NA", replacement = "", x = taxon),
         taxon = gsub(pattern = "Draparnaldia", replacement = "Draparnaldia spp.", x = taxon),
         taxon = trimws(taxon))

# Pull species scores from NMDS
species_scores <- as.data.frame(scores(x = essential_fatty_acid_metaMDS, display = "species"))
species_scores$species <- rownames(species_scores)

# Plot NMDS for all species but only Essential Fatty Acids
# This plot is figure S2 in the associated ms.
nmds <- ggplot() +
  geom_point(data = data_scores %>%
               mutate(taxon = ifelse(taxon == "Draparnaldia spp.",
                                     yes = "*Draparnaldia* spp.",
                                     no = taxon), 
                      taxon = ifelse(taxon %in% amphipods,
                                     yes = paste0("*", taxon, "*"),
                                     no = taxon)),
             mapping = aes(x = NMDS1, y = NMDS2, fill = taxon),
             size = 10, alpha = .6, shape = 21, stroke = 2, color = "grey70") +
  scale_fill_viridis_d(option = "plasma") +
  geom_text_repel(data = species_scores %>%
                    filter(species %in% c("c18_3w3", "c18_2w6", "c20_5w3")) %>%
                    mutate(species = gsub(pattern = "c", replacement = "", 
                                          x = species),
                           species = gsub(pattern = "_", replacement = ":", 
                                          x = species),
                           species = gsub(pattern = "w", replacement = "\U03C9", 
                                          x = species)), 
                  aes(x = NMDS1, y = NMDS2, label = species),
                  size = 10) +
  annotate("label", x = -0.5, y = 0.35,
           label = paste("Stress: ",
                         round(essential_fatty_acid_metaMDS$stress, digits = 3)),
           size = 10) +
  theme_minimal() +
  theme(legend.position = "bottom",
        title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_markdown(size = 12))
nmds

ggsave(filename = "all_species_essential_FA.png", plot = nmds, device = "png",
       path = "../figures/", width = 12, height = 10, units = "in",
       dpi = 300)


# 3. Correlating fatty acids with sewage ----------------------------------

# Create proportions of each fatty acid in the dataset
fatty_acid_prop_ppcp_site_info_dist <- fatty_acid_ppcp_site_info_dist %>%
  gather(key = fatty_acid, value = concentration, c12_0:c24_0) %>%
  group_by(site, Genus, Species) %>%
  mutate(total_fatty_acid = sum(concentration),
         prop_fatty_acid = concentration / total_fatty_acid) %>%
  select(-concentration, -total_fatty_acid) %>%
  spread(key = fatty_acid, value = prop_fatty_acid)

# 3.1 First do filamentous:diatom FA --------------------------------------

# First compare essential fatty acid ratios in periphyton with total PPCP concentration
peri_ppcp_lm <- lm((c18_3w3 + c18_1w9 + c18_2w6 + c16_0)/(c16_1w7 + c20_5w3 + c14_0 + c16_0)  ~ log10(ppcp_sum),
                   data = filter(fatty_acid_prop_ppcp_site_info_dist, Genus == "Periphyton"))

summary(peri_ppcp_lm)

# Also compare essential fatty acid ratios in periphyton with IDW population
peri_pop_lm <- lm((c18_3w3 + c18_1w9 + c18_2w6 + c16_0)/(c16_1w7 + c20_5w3 + c14_0 + c16_0)  ~ log10(distance_weighted_population),
                   data = filter(fatty_acid_prop_ppcp_site_info_dist, Genus == "Periphyton"))

summary(peri_pop_lm)

peri_ppcp_lm_permute <- permute_data_analytics(data = fatty_acid_prop_ppcp_site_info_dist %>%
                                                 ungroup() %>%
                                                 filter(Genus == "Periphyton") %>%
                                                 mutate(fil_diatom_fa_ratio = (c18_3w3 + c18_1w9 + c18_2w6 + c16_0)/(c16_1w7 + c20_5w3 + c14_0 + c16_0)), 
                                               metric = "fil_diatom_fa_ratio", 
                                               metric_plot_title = "Periphyton vs. &#91;Total PPCP&#93;",
                                               full_model = peri_ppcp_lm, 
                                               predictor = "ppcp")

peri_pop_lm_permute <- permute_data_analytics(data = fatty_acid_prop_ppcp_site_info_dist %>%
                                                 ungroup() %>%
                                                 filter(Genus == "Periphyton") %>%
                                                 mutate(fil_diatom_fa_ratio = (c18_3w3 + c18_1w9 + c18_2w6 + c16_0)/(c16_1w7 + c20_5w3 + c14_0 + c16_0)), 
                                               metric = "fil_diatom_fa_ratio", 
                                               metric_plot_title = "Periphyton vs. IDW Pop",
                                               full_model = peri_pop_lm, 
                                               predictor = "idw_pop")

# Second compare essential fatty acid ratios in E. verrucosus with total PPCP conentration
eulimnogammarus_verrucosus_ppcp_lm <- 
  lm((c18_3w3 + c18_1w9 + c18_2w6 + c16_0)/(c16_1w7 + c20_5w3 + c14_0 + c16_0) ~ log10(ppcp_sum),
     data = filter(fatty_acid_prop_ppcp_site_info_dist, 
                   taxon == "Eulimnogammarus_verrucosus"))

summary(eulimnogammarus_verrucosus_ppcp_lm)

# Next compare essential fatty acid ratios in E. verrucosus with IDW population
eulimnogammarus_verrucosus_pop_lm <- 
  lm((c18_3w3 + c18_1w9 + c18_2w6 + c16_0)/(c16_1w7 + c20_5w3 + c14_0 + c16_0) ~ log10(distance_weighted_population),
     data = filter(fatty_acid_prop_ppcp_site_info_dist, 
                   taxon == "Eulimnogammarus_verrucosus"))

summary(eulimnogammarus_verrucosus_pop_lm)

eulimnogammarus_verrucosus_ppcp_lm_permute <- permute_data_analytics(data = fatty_acid_prop_ppcp_site_info_dist %>%
                                                 ungroup() %>%
                                                 filter(taxon == "Eulimnogammarus_verrucosus") %>%
                                                 mutate(fil_diatom_fa_ratio = (c18_3w3 + c18_1w9 + c18_2w6 + c16_0)/(c16_1w7 + c20_5w3 + c14_0 + c16_0)), 
                                               metric = "fil_diatom_fa_ratio", 
                                               metric_plot_title = "*Eulimnogammarus verrucosus* vs. &#91;Total PPCP&#93;",
                                               full_model = eulimnogammarus_verrucosus_ppcp_lm, 
                                               predictor = "ppcp")

eulimnogammarus_verrucosus_pop_lm_permute <- permute_data_analytics(data = fatty_acid_prop_ppcp_site_info_dist %>%
                                                                       ungroup() %>%
                                                                       filter(taxon == "Eulimnogammarus_verrucosus") %>%
                                                                       mutate(fil_diatom_fa_ratio = (c18_3w3 + c18_1w9 + c18_2w6 + c16_0)/(c16_1w7 + c20_5w3 + c14_0 + c16_0)), 
                                                                     metric = "fil_diatom_fa_ratio", 
                                                                     metric_plot_title = "*Eulimnogammarus verrucosus* vs. IDW Pop",
                                                                     full_model = eulimnogammarus_verrucosus_pop_lm, 
                                                                     predictor = "idw_pop")

# Third compare essential fatty acid ratios in E. vittatus  with total PPCP conentration
eulimnogammarus_vittatus_ppcp_lm <- 
  lm((c18_3w3 + c18_1w9 + c18_2w6 + c16_0)/(c16_1w7 + c20_5w3 + c14_0 + c16_0) ~ log10(ppcp_sum),
     data = filter(fatty_acid_prop_ppcp_site_info_dist, 
                   taxon == "Eulimnogammarus_vittatus"))

summary(eulimnogammarus_vittatus_ppcp_lm)

eulimnogammarus_vittatus_pop_lm <- 
  lm((c18_3w3 + c18_1w9 + c18_2w6 + c16_0)/(c16_1w7 + c20_5w3 + c14_0 + c16_0) ~ log10(distance_weighted_population),
     data = filter(fatty_acid_prop_ppcp_site_info_dist, 
                   taxon == "Eulimnogammarus_vittatus"))

summary(eulimnogammarus_vittatus_pop_lm)

eulimnogammarus_vittatus_ppcp_lm_permute <- permute_data_analytics(data = fatty_acid_prop_ppcp_site_info_dist %>%
                                                                       ungroup() %>%
                                                                       filter(taxon == "Eulimnogammarus_vittatus") %>%
                                                                       mutate(fil_diatom_fa_ratio = (c18_3w3 + c18_1w9 + c18_2w6 + c16_0)/(c16_1w7 + c20_5w3 + c14_0 + c16_0)), 
                                                                     metric = "fil_diatom_fa_ratio", 
                                                                     metric_plot_title = "*Eulimnogammarus vittatus* vs.  &#91;Total PPCP&#93;",
                                                                     full_model = eulimnogammarus_vittatus_ppcp_lm, 
                                                                     predictor = "ppcp")

eulimnogammarus_vittatus_pop_lm_permute <- permute_data_analytics(data = fatty_acid_prop_ppcp_site_info_dist %>%
                                                                    ungroup() %>%
                                                                    filter(taxon == "Eulimnogammarus_vittatus") %>%
                                                                    mutate(fil_diatom_fa_ratio = (c18_3w3 + c18_1w9 + c18_2w6 + c16_0)/(c16_1w7 + c20_5w3 + c14_0 + c16_0)), 
                                                                  metric = "fil_diatom_fa_ratio", 
                                                                  metric_plot_title = "*Eulimnogammarus vittatus* vs.  IDW Pop",
                                                                  full_model = eulimnogammarus_vitatus_pop_lm, 
                                                                  predictor = "idw_pop")
# Now we will build the plots for FA ~ PPCPs

# Extract model R-squared values from each linear model
r_squared <- c(summary(peri_ppcp_lm)$r.squared,
               summary(eulimnogammarus_verrucosus_ppcp_lm)$r.squared,
               summary(eulimnogammarus_vittatus_ppcp_lm)$r.squared)

# Extract model p-values from each linear model
p_values <- c(summary(peri_ppcp_lm)$coefficients[2,4],
              summary(eulimnogammarus_verrucosus_ppcp_lm)$coefficients[2,4],
              summary(eulimnogammarus_vittatus_ppcp_lm)$coefficients[2,4])

n_samps <- c(length(peri_ppcp_lm$residuals),
             length(eulimnogammarus_verrucosus_ppcp_lm$residuals),
             length(eulimnogammarus_vittatus_ppcp_lm$residuals))

taxon <- c("Periphyton", "*Eulimnogammarus verrucosus*", "*Eulimnogammarus vittatus*")

labels <- data.frame(taxon, p_values, r_squared, n_samps) %>%
  mutate(label = paste0("p-value: ",
                        round(p_values, 3),
                        "<br>R<sup>2</sup>: ",
                        round(r_squared, 3),
                        "<br> N = ", n_samps))

# This figure is Figure 7 within the body of the associated manuscript.

ppcp_filamentous_diatom_fa_plot <- fatty_acid_prop_ppcp_site_info_dist %>%
    filter(taxon %in% c("Eulimnogammarus_verrucosus", "Eulimnogammarus_vittatus",
                        "Periphyton_NA")) %>%
    mutate(taxon = ifelse(test = taxon == "Eulimnogammarus_verrucosus",
                          yes = "*Eulimnogammarus verrucosus*", no = taxon),
           taxon = ifelse(test = taxon == "Eulimnogammarus_vittatus",
                          yes = "*Eulimnogammarus vittatus*", no = taxon),
           taxon = ifelse(test = taxon == "Periphyton_NA",
                          yes = "Periphyton", no = taxon)) %>%
    ggplot(aes(x = (ppcp_sum), 
               y = (c18_3w3 + c18_1w9 + c18_2w6 + c16_0)/(c16_1w7 + c20_5w3 + c14_0 + c16_0))) +
    geom_point(size = 3) +
    scale_x_log10() +
    facet_wrap(~ taxon) +
    geom_smooth(method = "lm") +
    geom_richtext(data = labels %>% filter(taxon != "Periphyton"), 
               mapping = aes(label = label, 
                             x = 0.01, y = 1.75), size = 7) +
    geom_richtext(data = labels %>% filter(taxon == "Periphyton"), 
                  mapping = aes(label = label, 
                                x = 0.01, y = 1.15), size = 7) +
    xlab(label = "Total PPCP (\u00b5g/L)") +
    ylab(label = expression(frac(18:3~omega~3 + 18:1~omega~9 + 18:2~omega~6 + 16:0, 
                                 16:1~omega~7 + 20:5~omega~3 + 16:0 + 14:0))) +
    theme_bw() +
    theme(text = element_text(size = 20),
          axis.title = element_text(size = 24),
          strip.text = element_markdown(size = 24))

ggsave(filename = "ppcp_filamentous_diatom_fa_plot.png", plot =  ppcp_filamentous_diatom_fa_plot, device = "png",
       path = "../figures/", width = 16, height = 8, units = "in")

# Now we will build the plots for FA ~ IDW population

# Extract model R-squared values from each linear model
r_squared <- c(summary(peri_pop_lm)$r.squared,
               summary(eulimnogammarus_verrucosus_pop_lm)$r.squared,
               summary(eulimnogammarus_vittatus_pop_lm)$r.squared)

# Extract model p-values from each linear model
p_values <- c(summary(peri_pop_lm)$coefficients[2,4],
              summary(eulimnogammarus_verrucosus_pop_lm)$coefficients[2,4],
              summary(eulimnogammarus_vittatus_pop_lm)$coefficients[2,4])

n_samps <- c(length(peri_pop_lm$residuals),
             length(eulimnogammarus_verrucosus_pop_lm$residuals),
             length(eulimnogammarus_vittatus_pop_lm$residuals))

taxon <- c("Periphyton", "*Eulimnogammarus verrucosus*", "*Eulimnogammarus vittatus*")

labels <- data.frame(taxon, p_values, r_squared, n_samps) %>%
  mutate(label = paste0("p-value: ",
                        round(p_values, 3),
                        "<br>R<sup>2</sup>: ",
                        round(r_squared, 3),
                        "<br> N = ", n_samps))

# This figure is Figure 7 within the body of the associated manuscript.

pop_filamentous_diatom_fa_plot <- fatty_acid_prop_ppcp_site_info_dist %>%
  filter(taxon %in% c("Eulimnogammarus_verrucosus", "Eulimnogammarus_vittatus",
                      "Periphyton_NA")) %>%
  mutate(taxon = ifelse(test = taxon == "Eulimnogammarus_verrucosus",
                        yes = "*Eulimnogammarus verrucosus*", no = taxon),
         taxon = ifelse(test = taxon == "Eulimnogammarus_vittatus",
                        yes = "*Eulimnogammarus vittatus*", no = taxon),
         taxon = ifelse(test = taxon == "Periphyton_NA",
                        yes = "Periphyton", no = taxon)) %>%
  ggplot(aes(x = (distance_weighted_population), 
             y = (c18_3w3 + c18_1w9 + c18_2w6 + c16_0)/(c16_1w7 + c20_5w3 + c14_0 + c16_0))) +
  geom_point(size = 3) +
  facet_wrap(~ taxon) +
  geom_smooth(method = "lm") +
  scale_x_log10() +
  geom_richtext(data = labels %>% filter(taxon != "Periphyton"), 
                mapping = aes(label = label, 
                              x = 1300, y = 1.75), size = 7) +
  geom_richtext(data = labels %>% filter(taxon == "Periphyton"), 
                mapping = aes(label = label, 
                              x = 1300, y = 1.15), size = 7) +
  xlab(label = "log10(IDW Population)") +
  ylab(label = expression(frac(18:3~omega~3 + 18:1~omega~9 + 18:2~omega~6 + 16:0, 
                               16:1~omega~7 + 20:5~omega~3 + 16:0 + 14:0))) +
  theme_bw() +
  theme(text = element_text(size = 20),
        strip.text = element_markdown())

ggsave(filename = "pop_filamentous_diatom_fa_plot.png", plot =  pop_filamentous_diatom_fa_plot, device = "png",
       path = "../figures/", width = 12, height = 6, units = "in")

# 3.2 Second do filamentous:diatom EFAs -----------------------------------


# Create linear models to explore which sewage indicators relate with fatty acid profiles
# Our analysis considers C18 fatty acids in comparison to C20,22 fatty acids.
# These two fatty acid groups roughly reflect Green:Diatom algal signature,
# which should inrease with an increasing sewage signal.

# First compare essential fatty acid ratios in periphyton with total PPCP concentration
peri_ppcp_lm <- lm(((c18_3w3 + c18_2w6) / (c20_5w3)) ~ log10(ppcp_sum),
               data = filter(fatty_acid_prop_ppcp_site_info_dist, Genus == "Periphyton"))

summary(peri_ppcp_lm)

# Next compare essential fatty acid ratios in periphyton with IDW population
peri_pop_lm <- lm(((c18_3w3 + c18_2w6) / (c20_5w3)) ~ log10(distance_weighted_population),
                   data = filter(fatty_acid_prop_ppcp_site_info_dist, Genus == "Periphyton"))

summary(peri_pop_lm)

peri_efa_ppcp_lm_permute <- permute_data_analytics(data = fatty_acid_prop_ppcp_site_info_dist %>%
                                                 ungroup() %>%
                                                 filter(Genus == "Periphyton") %>%
                                                 mutate(fil_diatom_fa_ratio = (c18_3w3 + c18_2w6) / (c20_5w3)), 
                                               metric = "fil_diatom_fa_ratio", 
                                               metric_plot_title = "Periphyton vs. &#91;Total PPCP&#93;",
                                               full_model = peri_ppcp_lm, 
                                               predictor = "ppcp")

peri_efa_pop_lm_permute <- permute_data_analytics(data = fatty_acid_prop_ppcp_site_info_dist %>%
                                                     ungroup() %>%
                                                     filter(Genus == "Periphyton") %>%
                                                     mutate(fil_diatom_fa_ratio = (c18_3w3 + c18_2w6) / (c20_5w3)), 
                                                   metric = "fil_diatom_fa_ratio", 
                                                   metric_plot_title = "Periphyton vs. IDW Pop",
                                                   full_model = peri_pop_lm, 
                                                   predictor = "idw_pop")

# Second compare essential fatty acid ratios in E. verrucosus with total PPCP conentration
eulimnogammarus_verrucosus_ppcp_lm <- lm(((c18_3w3 + c18_2w6) / (c20_5w3)) ~ log10(ppcp_sum),
                                         data = filter(fatty_acid_prop_ppcp_site_info_dist, 
                                                       taxon == "Eulimnogammarus_verrucosus"))

summary(eulimnogammarus_verrucosus_ppcp_lm)

# Second compare essential fatty acid ratios in E. verrucosus with IDW Population
eulimnogammarus_verrucosus_pop_lm <- lm(((c18_3w3 + c18_2w6) / (c20_5w3)) ~ log10(distance_weighted_population),
                                         data = filter(fatty_acid_prop_ppcp_site_info_dist, 
                                                       taxon == "Eulimnogammarus_verrucosus"))

summary(eulimnogammarus_verrucosus_pop_lm)

eulimnogammarus_verrucosus_efa_ppcp_lm_permute <- permute_data_analytics(data = fatty_acid_prop_ppcp_site_info_dist %>%
                                                                       ungroup() %>%
                                                                       filter(taxon == "Eulimnogammarus_verrucosus") %>%
                                                                       mutate(fil_diatom_fa_ratio = (c18_3w3 + c18_2w6) / (c20_5w3)), 
                                                                     metric = "fil_diatom_fa_ratio", 
                                                                     metric_plot_title = "*Eulimnogammarus verrucosus* vs. &#91;Total PPCP&#93;",
                                                                     full_model = eulimnogammarus_verrucosus_ppcp_lm, 
                                                                     predictor = "ppcp")

eulimnogammarus_verrucosus_efa_pop_lm_permute <- permute_data_analytics(data = fatty_acid_prop_ppcp_site_info_dist %>%
                                                                           ungroup() %>%
                                                                           filter(taxon == "Eulimnogammarus_verrucosus") %>%
                                                                           mutate(fil_diatom_fa_ratio = (c18_3w3 + c18_2w6) / (c20_5w3)), 
                                                                         metric = "fil_diatom_fa_ratio", 
                                                                         metric_plot_title = "*Eulimnogammarus verrucosus* vs. IDW pop",
                                                                         full_model = eulimnogammarus_verrucosus_pop_lm, 
                                                                         predictor = "idw_pop")

# Third compare essential fatty acid ratios in E. vittatus  with total PPCP conentration
eulimnogammarus_vittatus_ppcp_lm <- lm(((c18_3w3 + c18_2w6) / (c20_5w3)) ~ log10(ppcp_sum),
                                         data = filter(fatty_acid_prop_ppcp_site_info_dist, 
                                                       taxon == "Eulimnogammarus_vittatus"))

summary(eulimnogammarus_vittatus_ppcp_lm)

# Next compare essential fatty acid ratios in E. vittatus  with IDW population
eulimnogammarus_vittatus_pop_lm <- lm(((c18_3w3 + c18_2w6) / (c20_5w3)) ~ log10(distance_weighted_population),
                                       data = filter(fatty_acid_prop_ppcp_site_info_dist, 
                                                     taxon == "Eulimnogammarus_vittatus"))

summary(eulimnogammarus_vittatus_pop_lm)

eulimnogammarus_vittatus_efa_ppcp_lm_permute <- permute_data_analytics(data = fatty_acid_prop_ppcp_site_info_dist %>%
                                                                    ungroup() %>%
                                                                    filter(taxon == "Eulimnogammarus_vittatus") %>%
                                                                    mutate(fil_diatom_fa_ratio = (c18_3w3 + c18_2w6) / (c20_5w3)), 
                                                                  metric = "fil_diatom_fa_ratio", 
                                                                  metric_plot_title = "*Eulimnogammarus vittatus* vs. &#91;Total PPCP&#93;",
                                                                  full_model = eulimnogammarus_vitatus_ppcp_lm, 
                                                                  predictor = "ppcp")

eulimnogammarus_vittatus_efa_pop_lm_permute <- permute_data_analytics(data = fatty_acid_prop_ppcp_site_info_dist %>%
                                                                         ungroup() %>%
                                                                         filter(taxon == "Eulimnogammarus_vittatus") %>%
                                                                         mutate(fil_diatom_fa_ratio = (c18_3w3 + c18_2w6) / (c20_5w3)), 
                                                                       metric = "fil_diatom_fa_ratio", 
                                                                       metric_plot_title = "*Eulimnogammarus vittatus* vs. IDW Pop",
                                                                       full_model = eulimnogammarus_vitatus_pop_lm, 
                                                                       predictor = "idw_pop")

# Extract model R-squared values from each linear model
r_squared <- c(summary(peri_ppcp_lm)$r.squared,
               summary(eulimnogammarus_verrucosus_ppcp_lm)$r.squared,
               summary(eulimnogammarus_vittatus_ppcp_lm)$r.squared)

# Extract model p-values from each linear model
p_values <- c(summary(peri_ppcp_lm)$coefficients[2,4],
              summary(eulimnogammarus_verrucosus_ppcp_lm)$coefficients[2,4],
              summary(eulimnogammarus_vittatus_ppcp_lm)$coefficients[2,4])

n_samps <- c(length(peri_ppcp_lm$residuals),
             length(eulimnogammarus_verrucosus_ppcp_lm$residuals),
             length(eulimnogammarus_vittatus_ppcp_lm$residuals))

taxon <- c("Periphyton", "*Eulimnogammarus verrucosus*", "*Eulimnogammarus vittatus*")

labels <- data.frame(taxon, p_values, r_squared, n_samps) %>%
  mutate(label = paste0("p-value: ",
                        round(p_values, 3),
                        "<br>R<sup>2</sup>: ",
                        round(r_squared, 3),
                        "<br> N = ", n_samps))

# This figure is Figure 7 within the body of the associated manuscript.
ppcp_efa_plot <- fatty_acid_prop_ppcp_site_info_dist %>%
  filter(taxon %in% c("Eulimnogammarus_verrucosus", "Eulimnogammarus_vittatus",
                      "Periphyton_NA")) %>%
  mutate(taxon = ifelse(test = taxon == "Eulimnogammarus_verrucosus",
                        yes = "*Eulimnogammarus verrucosus*", no = taxon),
         taxon = ifelse(test = taxon == "Eulimnogammarus_vittatus",
                        yes = "*Eulimnogammarus vittatus*", no = taxon),
         taxon = ifelse(test = taxon == "Periphyton_NA",
                        yes = "Periphyton", no = taxon)) %>%
  ggplot(aes(x = (ppcp_sum), y = ((c18_3w3 + c18_2w6) / (c20_5w3)))) +
  scale_x_log10() +
  geom_point(size = 3) +
  facet_wrap(~ taxon) +
  geom_smooth(method = "lm") +
  geom_richtext(data = labels %>% filter(taxon != "Periphyton"), 
                mapping = aes(label = label, x = 0.01, y = 5), size = 7) +
  geom_richtext(data = labels %>% filter(taxon == "Periphyton"), 
             mapping = aes(label = label, x = 0.01, y = 2.5), size = 7) +
  xlab(label = "Total PPCP (\u00b5g/L)") +
  ylab(label = expression(frac(18:3~omega~3 + 18:2~omega~6, 20:5~omega~3 ))) +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.title = element_text(size = 24),
        strip.text = element_markdown(size = 24))

ggsave(filename = "ppcp_efa_plot.png", plot = ppcp_efa_plot, device = "png",
       path = "../figures/", width = 16, height = 8, units = "in")

# Now we will build the plots for FA ~ IDW population

# Extract model R-squared values from each linear model
r_squared <- c(summary(peri_pop_lm)$r.squared,
               summary(eulimnogammarus_verrucosus_pop_lm)$r.squared,
               summary(eulimnogammarus_vittatus_pop_lm)$r.squared)

# Extract model p-values from each linear model
p_values <- c(summary(peri_pop_lm)$coefficients[2,4],
              summary(eulimnogammarus_verrucosus_pop_lm)$coefficients[2,4],
              summary(eulimnogammarus_vittatus_pop_lm)$coefficients[2,4])

n_samps <- c(length(peri_pop_lm$residuals),
             length(eulimnogammarus_verrucosus_pop_lm$residuals),
             length(eulimnogammarus_vittatus_pop_lm$residuals))

taxon <- c("Periphyton", "*Eulimnogammarus verrucosus*", "*Eulimnogammarus vittatus*")

labels <- data.frame(taxon, p_values, r_squared, n_samps) %>%
  mutate(label = paste0("p-value: ",
                        round(p_values, 3),
                        "<br>R<sup>2</sup>: ",
                        round(r_squared, 3),
                        "<br> N = ", n_samps))

# This figure is Figure 7 within the body of the associated manuscript.

pop_efa_plot <- fatty_acid_prop_ppcp_site_info_dist %>%
  filter(taxon %in% c("Eulimnogammarus_verrucosus", "Eulimnogammarus_vittatus",
                      "Periphyton_NA")) %>%
  mutate(taxon = ifelse(test = taxon == "Eulimnogammarus_verrucosus",
                        yes = "*Eulimnogammarus verrucosus*", no = taxon),
         taxon = ifelse(test = taxon == "Eulimnogammarus_vittatus",
                        yes = "*Eulimnogammarus vittatus*", no = taxon),
         taxon = ifelse(test = taxon == "Periphyton_NA",
                        yes = "Periphyton", no = taxon)) %>%
  ggplot(aes(x = (distance_weighted_population), 
             y = (c18_3w3 + c18_2w6)/(c20_5w3))) +
  geom_point(size = 3) +
  facet_wrap(~ taxon) +
  geom_smooth(method = "lm") +
  geom_richtext(data = labels %>% filter(taxon != "Periphyton"), 
                mapping = aes(label = label, 
                              x = 1300, y = 2.5), size = 7) +
  geom_richtext(data = labels %>% filter(taxon == "Periphyton"), 
                mapping = aes(label = label, 
                              x = 1300, y = 1.15), size = 7) +
  xlab(label = "log10(IDW Population)") +
  ylab(label = expression(frac(18:3~omega~3 + 18:2~omega~6, 
                               20:5~omega~3))) +
  theme_bw() +
  theme(text = element_text(size = 20),
        strip.text = element_markdown())

ggsave(filename = "pop_efa_plot.png", plot =  pop_efa_plot, device = "png",
       path = "../figures/", width = 12, height = 6, units = "in")


# 3.3 Combine the plots into one ------------------------------------------

arranged_plots <- ggarrange(ppcp_filamentous_diatom_fa_plot, ppcp_efa_plot, ncol = 1, nrow = 2, 
                            labels = "AUTO", font.label = list(size = 20, face = "bold"))

ggsave(filename = "combined_ppcp_fattyy_acids.png", plot = arranged_plots, device = "png", path = "../figures/", 
       width = 16, height = 12, units = "in")

arranged_plots <- ggarrange(pop_filamentous_diatom_fa_plot, pop_efa_plot, ncol = 1, nrow = 2, 
                            labels = "AUTO", font.label = list(size = 20, face = "bold"))

ggsave(filename = "combined_pop_fatty_acids.png", plot = arranged_plots, device = "png", path = "../figures/", 
       width = 12, height = 12, units = "in")

permuted_fil_dia_plots <- c(peri_ppcp_lm_permute, 
                            eulimnogammarus_verrucosus_ppcp_lm_permute, 
                            eulimnogammarus_vitatus_ppcp_lm_permute)

arranged_plots <- ggarrange(plotlist = permuted_fil_dia_plots, ncol = 2, nrow = 3, 
                             font.label = list(size = 20, face = "bold"))

ggsave(filename = "permuted_fil_dia_fatty_acids.png", plot = arranged_plots, 
       device = "png", path = "../figures/", 
       width = 14, height = 12, units = "in")

permuted_fil_dia_plots <- c(peri_pop_lm_permute, 
                            eulimnogammarus_verrucosus_pop_lm_permute, 
                            eulimnogammarus_vitatus_pop_lm_permute)

arranged_plots <- ggarrange(plotlist = permuted_fil_dia_plots, ncol = 2, nrow = 3, 
                            font.label = list(size = 20, face = "bold"))

ggsave(filename = "permuted_pop_fil_dia_fatty_acids.png", plot = arranged_plots, 
       device = "png", path = "../figures/", 
       width = 14, height = 12, units = "in")

permuted_efa_plots <- c(peri_efa_ppcp_lm_permute, 
                        eulimnogammarus_verrucosus_efa_ppcp_lm_permute, 
                        eulimnogammarus_vittatus_efa_ppcp_lm_permute)

arranged_plots <- ggarrange(plotlist = permuted_efa_plots, ncol = 2, nrow = 3, 
                            font.label = list(size = 20, face = "bold"))

ggsave(filename = "permuted_efa_fatty_acids.png", plot = arranged_plots, 
       device = "png", path = "../figures/", 
       width = 14, height = 12, units = "in")

permuted_efa_plots <- c(peri_efa_pop_lm_permute, 
                        eulimnogammarus_verrucosus_efa_pop_lm_permute, 
                        eulimnogammarus_vittatus_efa_pop_lm_permute)

arranged_plots <- ggarrange(plotlist = permuted_efa_plots, ncol = 2, nrow = 3, 
                            font.label = list(size = 20, face = "bold"))

ggsave(filename = "permuted_pop_efa_fatty_acids.png", plot = arranged_plots, 
       device = "png", path = "../figures/", 
       width = 14, height = 12, units = "in")


# 4. Analyses of fatty acid groups ----------------------------------------


# 4.1 Define the fatty acids we will be considering -----------------------

safa <- c("c12_0", "c14_0", "a_15_0",  "c15_0", "i_15_0",
          "c16_0", "a_17_0", "c17_0", "i_17_0", "c18_0", "c20_0", "c22_0", "c24_0")

mufa <- c("c14_1n5", "c15_1w7", "c17_1w7",
          "c16_1w5", "c16_1w6", "c16_1w7", "c16_1w8", "c16_1w9",
          "c18_1w7", "c18_1w9", "c20_1w7", "c20_1w9", "c22_1w7", "c22_1w9")

scufa <- c("c16_2w4", "c16_2w6", "c16_2w7", "c16_3w3", "c16_3w4", "c16_3w6", "c16_4w1", "c16_4w3", 
           "c18_2w6", "c18_2w6t", "c18_3w3", "c18_3w6", "c18_4w3", "c18_4w4", "c18_5w3")

lcufa <- c("c20_2w5_11", "c20_2w5_13", "c20_2w6", "c20_3w3", "c20_3w6", "c20_4w3", "c20_4w6", "c20_5w3",
           "c22_2w6", "c22_3w3", "c22_4w3", "c22_4w6", "c22_5w3", "c22_5w6", "c22_6w3")

complete_fatty_acid_repo <- data.frame(rbind(c("SAFA", paste(safa, collapse = ", ")),
                                             c("MUFA", paste(mufa, collapse = ", ")),
                                             c("SCUFA", paste(scufa, collapse = ", ")),
                                             c("LCUFA", paste(lcufa, collapse = ", ")))) %>%
  rename("fatty_acid_group" = "X1",
         "fatty_acid_included" = "X2")

write.csv(x = complete_fatty_acid_repo, 
          file = "../tables/complete_fatty_acid_repo.csv", 
          row.names = FALSE)


# 4.2 Create table of mean fatty acid proportions -------------------------

sample_count <- fatty_acid_ppcp_site_info_dist %>% 
  group_by(Genus, Species) %>% 
  count()

fatty_acid_type_props_summary_table <- fatty_acid_prop_ppcp_site_info_dist %>%
  gather(fatty_acid, fa_prop, a_15_0:i_17_0) %>%
  mutate(fatty_acid_type = ifelse(fatty_acid %in% safa, "SAFA", "Branched"),
         fatty_acid_type = ifelse(fatty_acid %in% mufa, "MUFA", fatty_acid_type),
         fatty_acid_type = ifelse(fatty_acid %in% scufa, "SCUFA", fatty_acid_type),
         fatty_acid_type = ifelse(fatty_acid %in% lcufa, "LCUFA", fatty_acid_type)) %>%
  select(site, Genus, Species, fatty_acid_type, fa_prop) %>%
  group_by(site, Genus, Species, fatty_acid_type) %>%
  summarize(sum_fa_prop = sum(fa_prop)) %>%
  ungroup() %>%
  group_by(Genus, Species, fatty_acid_type) %>%
  summarize(mean_fa_prop = mean(sum_fa_prop),
            sd_fa_prop = sd(sum_fa_prop)) %>%
  inner_join(., sample_count) %>%
  filter(Genus != "Hyalella") %>%
  pivot_wider(names_from = "fatty_acid_type", values_from = c("mean_fa_prop", "sd_fa_prop"))

write.csv(x = fatty_acid_type_props_summary_table, 
          file = "../tables/fatty_acid_type_props_summary_table.csv", 
          row.names = FALSE)


# 4.3 plot of fatty acid proportions over PPCP concentration --------------

fatty_acid_type_props_plot <- fatty_acid_prop_ppcp_site_info_dist %>%
  gather(fatty_acid, fa_prop, a_15_0:i_17_0) %>%
  mutate(fatty_acid_type = ifelse(fatty_acid %in% safa, "SAFA", "Branched"),
         fatty_acid_type = ifelse(fatty_acid %in% mufa, "MUFA", fatty_acid_type),
         fatty_acid_type = ifelse(fatty_acid %in% scufa, "SCPUFA", fatty_acid_type),
         fatty_acid_type = ifelse(fatty_acid %in% lcufa, "LCPUFA", fatty_acid_type)) %>%
  mutate(taxon = gsub(pattern = "_", replacement = " ", x = taxon),
         taxon = gsub(pattern = "NA", replacement = "", x = taxon),
         taxon = ifelse(test = grepl(pattern = "Draparnaldia", x = taxon),
                        yes = "*Draparnaldia* spp.", no = taxon),
         taxon = ifelse(test = !taxon %in% c("Periphyton ", "Snail ", "*Draparnaldia* spp."),
                        yes = paste0("*", taxon, "*"), no = taxon)) %>%
  select(site, ppcp_sum, distance_weighted_population, taxon, 
         Genus, Species, fatty_acid_type, fa_prop) %>%
  group_by(site, ppcp_sum, distance_weighted_population, taxon,
           Genus, Species, fatty_acid_type, ) %>%
  summarize(sum_fa_prop = sum(fa_prop)) %>%
  ungroup() %>%
  filter(Genus != "Hyalella") %>%
  ggplot(aes(x = log10(ppcp_sum), y = sum_fa_prop, color = taxon)) +
  geom_point(size = 4) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_viridis_d(option = "viridis") +
  facet_grid(~fatty_acid_type) +
  ylab("Total proportion") +
  xlab("log10([Total PPCP])") +
  theme_bw() +
  theme(legend.position = "bottom",
        title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        legend.text = element_markdown(size = 20),
        strip.text = element_text(size = 16))

ggsave(filename = "fatty_acid_type_props_plot.png", plot = fatty_acid_type_props_plot, device = "png",
       path = "../figures/", height = 8, width = 15, units = "in")


# 5. Multivariate analysis of filamentous and diatom FA -------------------


# 5.1 Primary producer analysis -------------------------------------------

periphyton_fatty_acids <- fatty_acid_prop_ppcp_site_info_dist %>%
  filter(Genus %in% c("Periphyton", "Draparnaldia")) %>%
  select(site:distance_weighted_population, c18_3w3, c16_0, c18_1w9, c18_2w6, 
         c16_1w7, c20_5w3, c14_0)

peri_nmds <- metaMDS(comm = periphyton_fatty_acids[ , 19:25], distance = "bray", try = 100, k = 1)

# Extract scores for each site and add additional data for analysis and plotting
data_scores <- as.data.frame(scores(peri_nmds)) %>%
  mutate(site = periphyton_fatty_acids$site,
         taxon = periphyton_fatty_acids$taxon,
         taxon = gsub(pattern = "_", replacement = " ", x = taxon),
         taxon = gsub(pattern = "NA", replacement = "", x = taxon),
         taxon = gsub(pattern = "Draparnaldia", replacement = "*Draparnaldia* spp.", x = taxon),
         ppcp_sum = periphyton_fatty_acids$ppcp_sum)

species_scores <- as.data.frame(scores(peri_nmds, display = "species"))
species_scores$species <- rownames(species_scores)

# Create plot
# This figure is associated with figure S1 in the associated manuscript
nmds <- ggplot() +
  geom_point(data = data_scores, aes(x = NMDS1, y = 0.5,
                                     size = ppcp_sum, fill = taxon),
             alpha = .5, shape = 21, stroke = 5, color = "grey70") +
  scale_fill_manual(values = viridis(20)[c(4, 10)]) +
  geom_text_repel(data = species_scores %>%
                    mutate(species = gsub(pattern = "c", replacement = "", 
                                          x = species),
                           species = gsub(pattern = "_", replacement = ":", 
                                          x = species),
                           species = gsub(pattern = "w", replacement = "\U03C9", 
                                          x = species)), 
                  aes(x = NMDS1, y = 1, label = species),
                  size = 10, segment.size = NA) +
  scale_size_continuous(name = "[Total PPCP]", range = c(10, 30)) +
  guides(shape = guide_legend(override.aes = list(size=10))) + 
  ylab("") + 
  ylim(c(0,1.25)) +
  annotate("label", x = 0, y = 0.25,
           label = paste("Stress: ", round(peri_nmds$stress, digits = 3)),
           size = 10) +
  theme_minimal() +
  guides(fill = guide_legend(override.aes = list(size = 15))) +
  theme(legend.position = "right",
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_markdown(size = 20),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.key.height = unit(0.75, "in"))

ggsave(filename = "filamentous_diatom_nmds_peri_draparnaldia.png", plot = nmds, device = "png", 
       path = "../figures/", width = 16, height = 8, units = "in")

# 5.2 Macroinvertebrate analysis ------------------------------------------

invert_fatty_acids <- fatty_acid_prop_ppcp_site_info_dist %>%
  filter(!(Genus %in% c("Draparnaldia", "Periphyton", "Snail", "Hyalella"))) %>%
  select(site:distance_weighted_population, c18_3w3, c16_0, c18_1w9, c18_2w6,  
         c16_1w7, c20_5w3, c14_0)

invert_nmds <- metaMDS(comm = invert_fatty_acids[ , 19:25], distance = "bray", try = 100, k = 2)

# Extract scores for each site and add additional data for analysis and plotting
data_scores <- as.data.frame(scores(invert_nmds)) %>%
  mutate(site = invert_fatty_acids$site,
         taxon = invert_fatty_acids$taxon,
         taxon = gsub(pattern = "_", replacement = " ", x = taxon),
         taxon = gsub(pattern = "NA", replacement = "", x = taxon),
         taxon = paste0("*", taxon, "*"),
         ppcp_sum = invert_fatty_acids$ppcp_sum)

species_scores <- as.data.frame(scores(invert_nmds, display = "species"))
species_scores$species <- rownames(species_scores)

# Create plot
# This figure is associated with figure S1 in the associated manuscript
nmds <- ggplot() +
  geom_point(data = data_scores, aes(x = NMDS1, y = NMDS2, fill = taxon,
                                     size = ppcp_sum),
             alpha = .5, shape = 21, stroke = 5, color = "grey70") +
  scale_fill_manual(name = "taxon", values = viridis(30)[c(4, 9, 17, 24)]) +
  geom_text_repel(data = species_scores %>%
                    mutate(species = gsub(pattern = "c", replacement = "", 
                                          x = species),
                           species = gsub(pattern = "_", replacement = ":", 
                                          x = species),
                           species = gsub(pattern = "w", replacement = "\U03C9", 
                                          x = species)), 
                  aes(x = NMDS1, y = NMDS2, label = species),
                  size = 10, segment.size = NA) +
  scale_size_continuous(name = "[Total PPCP]", range = c(10,30)) +
  guides(shape = guide_legend(override.aes = list(size=10))) + 
  annotate("label", x = 0.2, y = 0.2,
           label = paste("Stress: ", round(invert_nmds$stress, digits = 3)),
           size = 10) +
  theme_minimal() +
  guides(fill = guide_legend(override.aes = list(size = 10))) +
  theme(legend.position = "right",
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_markdown(size = 18),
        legend.key.height = unit(0.5, "in"))

ggsave(filename = "filamentous_diatom_nmds_amphipods.png", plot = nmds, device = "png", 
       path = "../figures/", width = 16, height = 10, units = "in")
