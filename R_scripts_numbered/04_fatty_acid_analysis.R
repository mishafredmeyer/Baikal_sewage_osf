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


# 1. Load the data --------------------------------------------------------

fatty_acid <- read.csv(file = "../cleaned_data/fatty_acid.csv",
                       header = TRUE, stringsAsFactors = FALSE)

metadata <- read.csv(file = "../cleaned_data/metadata.csv",
                     header = TRUE, stringsAsFactors = FALSE)

distance <- read.csv(file = "../cleaned_data/distance_weighted_population_metrics.csv",
                     header = TRUE, stringsAsFactors = FALSE)

nutrients <- read.csv(file = "../cleaned_data/nutrients.csv",
                      header = TRUE, stringsAsFactors = FALSE)

metadata_dist <- full_join(x = metadata, y = distance)

ppcp <- read.csv(file = "../cleaned_data/ppcp.csv",
                 header = TRUE, stringsAsFactors = FALSE)

ppcp_meta_dist <- full_join(x = ppcp, y = metadata_dist, by = c("site"))

fatty_acid_ppcp_meta_dist <- inner_join(x = fatty_acid, y = ppcp_meta_dist,
                                        by = c("site")) %>%
  unite(taxon, c("Genus", "Species"), remove = FALSE)


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
  column_to_rownames(var = "fa") %>%
  arrange(desc(Var_Mean_Ratio))

write.csv(x = mean_var, file = "../tables/cv_all_fa.csv", row.names = TRUE)

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
         taxon = gsub(pattern = "Drapa", replacement = "Drapa spp.", x = taxon))

# Pull species scores from NMDS
species_scores <- as.data.frame(scores(x = whole_fatty_acid_metaMDS, display = "species"))
species_scores$species <- rownames(species_scores)

# Create plot
# This figure is associated with figure S1 in the associated manuscript
nmds <- ggplot() +
  geom_point(data = data_scores, aes(x = NMDS1, y = NMDS2, shape = taxon),
             size = 10, alpha = .75) +
  #scale_color_manual(values = viridis(69)[c(1, 13, 25, 33, 38, 50, 60, 69)]) +
  scale_shape_manual(values = c(1:7)) + 
  geom_text_repel(data = species_scores %>%
                    filter(species %in% c("c18_3w3", "c18_1w9", "c18_2w6", "c16_0", 
                           "c14_0", "c20_5w3", "c161w7")), 
                  aes(x = NMDS1, y = NMDS2, label = species),
                  size = 5) +
  ggtitle("NMDS with Entire FA Spectrum") +
  annotate("label", x = 0.4, y = 0.45,
           label = paste("Stress: ", round(whole_fatty_acid_metaMDS$stress, digits = 3)),
           size = 10) +
  theme_minimal() +
  theme(legend.position = "bottom",
        title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 12))
nmds

ggsave(filename = "all_species_all_FA_symbols.png", plot = nmds, device = "png",
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
  column_to_rownames(var = "fa") %>%
  arrange(desc(Var_Mean_Ratio))

write.csv(x = mean_var[order(mean_var$Var_Mean_Ratio, decreasing = TRUE), ], 
          file = "../tables/cv_efa.csv", row.names = TRUE)


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
         taxon = gsub(pattern = "Drapa", replacement = "Drapa spp.", x = taxon))

# Pull species scores from NMDS
species_scores <- as.data.frame(scores(x = essential_fatty_acid_metaMDS, display = "species"))
species_scores$species <- rownames(species_scores)

# Plot NMDS for all species but only Essential Fatty Acids
# This plot is figure S2 in the associated ms.
nmds <- ggplot() +
  geom_point(data = data_scores, aes(x = NMDS1, y = NMDS2, shape = taxon),
             size = 10, alpha = .75) +
  #scale_fill_manual(values = inferno(69)[c(1, 13, 18, 20, 27, 33, 38, 45, 55)]) +
  scale_shape_manual(values = c(1:7)) + 
  geom_text_repel(data = species_scores %>%
                    filter(species %in% c("c18_3w3", "c18_2w6", "c20_5w3")), 
                  aes(x = NMDS1, y = NMDS2, label = species),
                  size = 5) +
  ggtitle("NMDS with Essential FA Spectrum") +
  annotate("label", x = -0.5, y = 0.35,
           label = paste("Stress: ",
                         round(essential_fatty_acid_metaMDS$stress, digits = 3)),
           size = 10) +
  theme_minimal() +
  theme(legend.position = "bottom",
        title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 12))
nmds

ggsave(filename = "all_species_essential_FA_symbol.png", plot = nmds, device = "png",
       path = "../figures/", width = 12, height = 10, units = "in",
       dpi = 300)


# 3. Correlating fatty acids with sewage ----------------------------------

# Create proportions of each fatty acid in the dataset
fatty_acid_prop_ppcp_meta_dist <- fatty_acid_ppcp_meta_dist %>%
  gather(key = fatty_acid, value = concentration, c12_0:c24_0) %>%
  group_by(site, Genus, Species) %>%
  mutate(total_fatty_acid = sum(concentration),
         prop_fatty_acid = concentration / total_fatty_acid) %>%
  select(-concentration, -total_fatty_acid) %>%
  spread(key = fatty_acid, value = prop_fatty_acid)

# 3.1 First do filamentous:diatom FA --------------------------------------

# First compare essential fatty acid ratios in periphyton with total PPCP concentration
peri_ppcp_lm <- lm((c18_3w3 + c18_1w9 + c18_2w6 + c16_0)/(c16_1w7 + c20_5w3 + c14_0 + c16_0)  ~ log10(ppcp_sum),
                   data = filter(fatty_acid_prop_ppcp_meta_dist, Genus == "Periphyton"))

summary(peri_ppcp_lm)

# Second compare essential fatty acid ratios in E. verrucosus with total PPCP conentration
eulimnogammarus_verrucosus_ppcp_lm <- 
  lm((c18_3w3 + c18_1w9 + c18_2w6 + c16_0)/(c16_1w7 + c20_5w3 + c14_0 + c16_0) ~ log10(ppcp_sum),
     data = filter(fatty_acid_prop_ppcp_meta_dist, 
                   taxon == "Eulimnogammarus_verrucosus"))

summary(eulimnogammarus_verrucosus_ppcp_lm)

# Third compare essential fatty acid ratios in E. vittatus  with total PPCP conentration
eulimnogammarus_vitatus_ppcp_lm <- 
  lm((c18_3w3 + c18_1w9 + c18_2w6 + c16_0)/(c16_1w7 + c20_5w3 + c14_0 + c16_0) ~ log10(ppcp_sum),
     data = filter(fatty_acid_prop_ppcp_meta_dist, 
                   taxon == "Eulimnogammarus_vittatus"))

summary(eulimnogammarus_vitatus_ppcp_lm)

# Extract model R-squared values from each linear model
r_squared <- c(summary(peri_ppcp_lm)$r.squared,
               summary(eulimnogammarus_verrucosus_ppcp_lm)$r.squared,
               summary(eulimnogammarus_vitatus_ppcp_lm)$r.squared)

# Extract model p-values from each linear model
p_values <- c(summary(peri_ppcp_lm)$coefficients[2,4],
              summary(eulimnogammarus_verrucosus_ppcp_lm)$coefficients[2,4],
              summary(eulimnogammarus_vitatus_ppcp_lm)$coefficients[2,4])

taxon <- c("Periphyton", "Eulimnogammarus verrucosus", "Eulimnogammarus vittatus")

labels <- data.frame(taxon, p_values, r_squared) %>%
  mutate(label = paste0("p-value: ",
                        round(p_values, 3),
                        "\nR-squared: ",
                        round(r_squared, 3)))

# This figure is Figure 7 within the body of the associated manuscript.

ppcp_filamentous_diatom_fa_plot <- fatty_acid_prop_ppcp_meta_dist %>%
    filter(taxon %in% c("Eulimnogammarus_verrucosus", "Eulimnogammarus_vittatus",
                        "Periphyton_NA")) %>%
    mutate(taxon = ifelse(test = taxon == "Eulimnogammarus_verrucosus",
                          yes = "Eulimnogammarus verrucosus", no = taxon),
           taxon = ifelse(test = taxon == "Eulimnogammarus_vittatus",
                          yes = "Eulimnogammarus vittatus", no = taxon),
           taxon = ifelse(test = taxon == "Periphyton_NA",
                          yes = "Periphyton", no = taxon)) %>%
    ggplot(aes(x = log10(ppcp_sum), 
               y = (c18_3w3 + c18_1w9 + c18_2w6 + c16_0)/(c16_1w7 + c20_5w3 + c14_0 + c16_0))) +
    geom_point(size = 3) +
    facet_wrap(~ taxon) +
    geom_smooth(method = "lm") +
    geom_label(data = labels %>% filter(taxon != "Periphyton"), aes(label = label, x = -2.0, y = 1.75), size = 4) +
    geom_label(data = labels %>% filter(taxon == "Periphyton"), aes(label = label, x = -2.0, y = 1.15), size = 4) +
    xlab(label = "log10([Total PPCP])") +
    ylab(label = expression(frac(18:3~omega~3 + 18:1~omega~9 + 18:2~omega~6 + 16:0, 
                                 16:1~omega~7 + 20:5~omega~3 + 16:0 + 14:0))) +
    theme_bw() +
    theme(text = element_text(size = 20))

ggsave(filename = "ppcp_filamentous_diatom_fa_plot.png", plot = ppcp_filamentous_diatom_fa_plot, device = "png",
       path = "../figures/", width = 12, height = 6, units = "in")

# 3.2 Second do filamentous:diatom EFAs -----------------------------------


# Create linear models to explore which sewage indicators relate with fatty acid profiles
# Our analysis considers C18 fatty acids in comparison to C20,22 fatty acids.
# These two fatty acid groups roughly reflect Green:Diatom algal signature,
# which should inrease with an increasing sewage signal.

# First compare essential fatty acid ratios in periphyton with total PPCP concentration
peri_ppcp_lm <- lm(((c18_3w3 + c18_2w6) / (c20_5w3)) ~ log10(ppcp_sum),
               data = filter(fatty_acid_prop_ppcp_meta_dist, Genus == "Periphyton"))

summary(peri_ppcp_lm)

# Second compare essential fatty acid ratios in E. verrucosus with total PPCP conentration
eulimnogammarus_verrucosus_ppcp_lm <- lm(((c18_3w3 + c18_2w6) / (c20_5w3)) ~ log10(ppcp_sum),
                                         data = filter(fatty_acid_prop_ppcp_meta_dist, 
                                                       taxon == "Eulimnogammarus_verrucosus"))

summary(eulimnogammarus_verrucosus_ppcp_lm)

# Third compare essential fatty acid ratios in E. vittatus  with total PPCP conentration
eulimnogammarus_vitatus_ppcp_lm <- lm(((c18_3w3 + c18_2w6) / (c20_5w3)) ~ log10(ppcp_sum),
                                         data = filter(fatty_acid_prop_ppcp_meta_dist, 
                                                       taxon == "Eulimnogammarus_vittatus"))

summary(eulimnogammarus_vitatus_ppcp_lm)

# Extract model R-squared values from each linear model
r_squared <- c(summary(peri_ppcp_lm)$r.squared,
               summary(eulimnogammarus_verrucosus_ppcp_lm)$r.squared,
               summary(eulimnogammarus_vitatus_ppcp_lm)$r.squared)

# Extract model p-values from each linear model
p_values <- c(summary(peri_ppcp_lm)$coefficients[2,4],
              summary(eulimnogammarus_verrucosus_ppcp_lm)$coefficients[2,4],
              summary(eulimnogammarus_vitatus_ppcp_lm)$coefficients[2,4])

taxon <- c("Periphyton", "Eulimnogammarus verrucosus", "Eulimnogammarus vittatus")

labels <- data.frame(taxon, p_values, r_squared) %>%
  mutate(label = paste0("p-value: ",
                        round(p_values, 3),
                        "\nR-squared: ",
                        round(r_squared, 3)))

# This figure is Figure 7 within the body of the associated manuscript.
ppcp_efa_plot <- fatty_acid_prop_ppcp_meta_dist %>%
  filter(taxon %in% c("Eulimnogammarus_verrucosus", "Eulimnogammarus_vittatus",
                      "Periphyton_NA")) %>%
  mutate(taxon = ifelse(test = taxon == "Eulimnogammarus_verrucosus",
                        yes = "Eulimnogammarus verrucosus", no = taxon),
         taxon = ifelse(test = taxon == "Eulimnogammarus_vittatus",
                        yes = "Eulimnogammarus vittatus", no = taxon),
         taxon = ifelse(test = taxon == "Periphyton_NA",
                        yes = "Periphyton", no = taxon)) %>%
  ggplot(aes(x = log10(ppcp_sum), y = ((c18_3w3 + c18_2w6) / (c20_5w3)))) +
  geom_point(size = 3) +
  facet_wrap(~ taxon) +
  geom_smooth(method = "lm") +
  geom_label(data = labels %>% filter(taxon != "Periphyton"), aes(label = label, x = -2.0, y = 5), size = 4) +
  geom_label(data = labels %>% filter(taxon == "Periphyton"), aes(label = label, x = -2.0, y = 2.5), size = 4) +
  xlab(label = "log10([Total PPCP])") +
  ylab(label = expression(frac(18:3~omega~3 + 18:2~omega~6, 20:5~omega~3 ))) +
  theme_bw() +
  theme(text = element_text(size = 20))

ggsave(filename = "ppcp_efa_plot.png", plot = ppcp_efa_plot, device = "png",
       path = "../figures/", width = 12, height = 6, units = "in")


# 3.3 Combine the plots into one ------------------------------------------

arranged_plots <- ggarrange(ppcp_filamentous_diatom_fa_plot, ppcp_efa_plot, ncol = 1, nrow = 2, 
                            labels = "AUTO", font.label = list(size = 20, face = "bold"))

ggsave(filename = "combined_ppcp_fattyy_acids.png", plot = arranged_plots, device = "png", path = "../figures/", 
       width = 12, height = 12, units = "in")


# 4. Analyses of fatty acid groups ----------------------------------------


# 4.1 Define the fatty acids we will be considering -----------------------

safa <- c("c12_0", "c14_0", "a_15_0",  "c15_0", "i_15_0",
          "c16_0", "a_17_0", "c17_0", "i_17_0", "c18_0", "c20_0", "c22_0", "c24_0")

mufa <- c("c14_1n5", "c15_1w7", "c17_1n7",
          "c16_1w5", "c16_1w6", "c16_1w7", "c16_1w8", "c16_1w9",
          "c18_1w7", "c18_1w9", "c20_1w7", "c20_1w9", "c22_1w7", "c22_1w9")

scufa <- c("c16_2w4", "c16_2w6", "c16_2w7", "c16_3w3", "c16_3w4", "c16_3w6", "c16_4w1", "c16_4w3", 
           "c18_2w6", "c18_2w6t", "c18_3w3", "c18_3w6", "c18_4w3", "c18_4w4", "c18_5w3")

lcufa <- c("c20_2_5_11", "c20_2_5_13", "c20_2w6", "c20_3w3", "c20_3w6", "c20_4w3", "c20_4w6", "c20_5w3",
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

sample_count <- fatty_acid_ppcp_meta_dist %>% 
  group_by(Genus, Species) %>% 
  count()

fatty_acid_type_props_summary_table <- fatty_acid_prop_ppcp_meta_dist %>%
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
  summarize(mean_fa_prop = mean(sum_fa_prop)) %>%
  inner_join(., sample_count) %>%
  filter(Genus != "Hyalella") %>%
  spread(fatty_acid_type, mean_fa_prop)

write.csv(x = fatty_acid_type_props_summary_table, 
          file = "../tables/fatty_acid_type_props_summary_table.csv", 
          row.names = FALSE)


# 4.3 plot of fatty acid proportions over PPCP concentration --------------

fatty_acid_type_props_plot <- fatty_acid_prop_ppcp_meta_dist %>%
  gather(fatty_acid, fa_prop, a_15_0:i_17_0) %>%
  mutate(fatty_acid_type = ifelse(fatty_acid %in% safa, "SAFA", "Branched"),
         fatty_acid_type = ifelse(fatty_acid %in% mufa, "MUFA", fatty_acid_type),
         fatty_acid_type = ifelse(fatty_acid %in% scufa, "SCPUFA", fatty_acid_type),
         fatty_acid_type = ifelse(fatty_acid %in% lcufa, "LCPUFA", fatty_acid_type)) %>%
  mutate(taxon = gsub(pattern = "_", replacement = " ", x = taxon),
         taxon = gsub(pattern = "NA", replacement = "", x = taxon),
         taxon = ifelse(test = grepl(pattern = "Drapa", x = taxon),
                        yes = "Drapa spp.", no = taxon)) %>%
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
  theme_bw() +
  theme(legend.position = "bottom",
        title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 16))

ggsave(filename = "fatty_acid_type_props_plot.png", plot = fatty_acid_type_props_plot, device = "png",
       path = "../figures/", height = 8, width = 15, units = "in")


# 5. Multivariate analysis of filamentous and diatom FA -------------------


# 5.1 Primary producer analysis -------------------------------------------

periphyton_fatty_acids <- fatty_acid_prop_ppcp_meta_dist %>%
  filter(Genus %in% c("Periphyton", "Drapa")) %>%
  select(site:distance_weighted_population, c18_3w3, c16_0, c18_1w9, c18_2w6, 
         c16_1w7, c20_5w3, c14_0)

peri_nmds <- metaMDS(comm = periphyton_fatty_acids[ , 19:25], distance = "bray", try = 100, k = 1)

# Extract scores for each site and add additional data for analysis and plotting
data_scores <- as.data.frame(scores(peri_nmds)) %>%
  mutate(site = periphyton_fatty_acids$site,
         taxon = periphyton_fatty_acids$taxon,
         taxon = gsub(pattern = "_", replacement = " ", x = taxon),
         taxon = gsub(pattern = "NA", replacement = "", x = taxon),
         taxon = gsub(pattern = "Drapa", replacement = "Drapa spp.", x = taxon),
         ppcp_sum = periphyton_fatty_acids$ppcp_sum)

species_scores <- as.data.frame(scores(peri_nmds, display = "species"))

# Create plot
# This figure is associated with figure S1 in the associated manuscript
nmds <- ggplot() +
  geom_point(data = data_scores, aes(x = NMDS1, y = 0.5, shape = taxon,
                                     size = ppcp_sum),
             alpha = .5) +
  #geom_point(data = species_scores, aes(x = NMDS1, y = 1), size = 4) +
  geom_text_repel(data = species_scores, aes(x = NMDS1, y = 1,
                                             label = rownames(species_scores)),
                  size = 10, segment.size = NA) +
  scale_size_continuous(name = "[Total PPCP]", range = c(5,20)) +
  guides(shape = guide_legend(override.aes = list(size=10))) + 
  ggtitle("NMDS with Filamentous:Diatom Fatty Acids") +
  ylab("") + 
  ylim(c(0,1.25)) +
  annotate("label", x = 0, y = 0.25,
           label = paste("Stress: ", round(peri_nmds$stress, digits = 3)),
           size = 10) +
  theme_minimal() +
  theme(legend.position = "right",
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

ggsave(filename = "filamentous_diatom_nmds_peri_drapa.png", plot = nmds, device = "png", 
       path = "../figures/", width = 16, height = 8, units = "in")

# 5.2 Macroinvertebrate analysis ------------------------------------------

invert_fatty_acids <- fatty_acid_prop_ppcp_meta_dist %>%
  filter(!(Genus %in% c("Drapa", "Periphyton", "Snail", "Hyalella"))) %>%
  select(site:distance_weighted_population, c18_3w3, c16_0, c18_1w9, c18_2w6,  
         c16_1w7, c20_5w3, c14_0)

invert_nmds <- metaMDS(comm = invert_fatty_acids[ , 19:25], distance = "bray", try = 100, k = 2)

# Extract scores for each site and add additional data for analysis and plotting
data_scores <- as.data.frame(scores(invert_nmds)) %>%
  mutate(site = invert_fatty_acids$site,
         taxon = invert_fatty_acids$taxon,
         taxon = gsub(pattern = "_", replacement = " ", x = taxon),
         taxon = gsub(pattern = "NA", replacement = "", x = taxon),
         taxon = gsub(pattern = "Drapa", replacement = "Drapa spp.", x = taxon),
         ppcp_sum = invert_fatty_acids$ppcp_sum)

species_scores <- as.data.frame(scores(invert_nmds, display = "species"))

# Create plot
# This figure is associated with figure S1 in the associated manuscript
nmds <- ggplot() +
  geom_point(data = data_scores, aes(x = NMDS1, y = NMDS2, shape = taxon,
                                     size = ppcp_sum),
             alpha = .5) +
  scale_shape_manual(name = "taxon", values = c(15:18)) +
  #geom_point(data = species_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_text_repel(data = species_scores %>%
                    mutate(NMDS1 = ifelse(test = rownames(species_scores) == "c18_2w6",
                                          yes = NMDS1-0.04, no = NMDS1)), 
                  aes(x = NMDS1, y = NMDS2, label = rownames(species_scores)),
                  size = 10, segment.size = NA) +
  scale_size_continuous(name = "[Total PPCP]", range = c(5,20)) +
  guides(shape = guide_legend(override.aes = list(size=10))) + 
  ggtitle("NMDS with Filamentous:Diatom Fatty Acids") +
  annotate("label", x = 0.2, y = 0.2,
           label = paste("Stress: ", round(invert_nmds$stress, digits = 3)),
           size = 10) +
  theme_minimal() +
  theme(legend.position = "right",
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 14))

ggsave(filename = "filamentous_diatom_nmds_amphipods.png", plot = nmds, device = "png", 
       path = "../figures/", width = 16, height = 10, units = "in")

