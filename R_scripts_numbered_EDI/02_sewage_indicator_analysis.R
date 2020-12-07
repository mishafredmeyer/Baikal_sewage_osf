# This script is meant to perform univariate analysis of various sewage 
# indicators along Lake Baikal's southwestern shoreline. Each of the sewage
# indicators are regressed against log-transformed inverse distance weighted 
# population using a linear model. For each indicator, an individual plot is
# produced, and then all individual plots are merged into a single, 
# final figure.

# 1. Load packages --------------------------------------------------------

library(tidyverse)
library(viridis)
library(viridisLite)
library(ggpubr)
library(fs)

# Check if figures directory exists 
# If not, create the figures directory
sub_dir <- "figures"
output_dir <- file.path(here::here(), sub_dir)

if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Dir 'figures' already exists!")
}

# 2. Load data ------------------------------------------------------------

# PPCP Data
ppcp <- read.csv(file = "../cleaned_data/ppcp.csv",
                 header = TRUE, stringsAsFactors = FALSE)

# Nutrient data
nutrients <- read.csv(file = "../cleaned_data/nutrients.csv",
                      header = TRUE, stringsAsFactors = FALSE)

# Stable isotope data
stable_isotopes <- read.csv(file = "../cleaned_data/stable_isotopes.csv",
                            header = TRUE, stringsAsFactors = FALSE)

# Chlorophyll a data
chlorophylla <- read.csv(file = "../cleaned_data/chlorophylla.csv",
                         header = TRUE, stringsAsFactors = FALSE)

# Microplastics data
microplastics <- read.csv(file = "../cleaned_data/microplastics.csv",
                          header = TRUE, stringsAsFactors = FALSE)

# Site metadata
metadata <- read.csv(file = "../cleaned_data/metadata.csv",
                     header = TRUE, stringsAsFactors = FALSE)

# Site distance data
distance <- read.csv(file = "../cleaned_data/distance_weighted_population_metrics.csv",
                     header = TRUE, stringsAsFactors = FALSE)

# Join site metadata with distance data
metadata_dist <- full_join(x = metadata, y = distance, by = "site")

# 3. PPCP analysis --------------------------------------------------------

# Join PPCP data with metadata/distance and create two custom metrics
ppcp_meta_dist <- full_join(x = ppcp, y = metadata_dist, by = "site")

# Analyze total PPCPs as a function of population intensity
ppcp_PI_model <- lm(log10(ppcp_sum) ~ log10(distance_weighted_population),
                    data = ppcp_meta_dist)

# View model results
summary(ppcp_PI_model)

# Plot linear model
ppcp_PI_plot <- ggplot(data = ppcp_meta_dist,
                       aes(x = log10(distance_weighted_population),
                           y =  log10(ppcp_sum))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  ylab("log10([Total PPCP])") +
  xlab("log10(IDW Population)") +
  ggtitle("PPCP vs. IDW Population") +
  annotate(geom = "label", x = 3.5, y = -2.65,
           label = paste("p-value: ",
                         round(summary(ppcp_PI_model)$coefficients[2, 4], 3),
                         "\nR-squared: ",
                         round(summary(ppcp_PI_model)$r.squared, 3)),
           parse = FALSE) +
  theme_minimal()

# Export plot
ggsave(filename = "../figures/ppcp_PI_plot.png", plot = ppcp_PI_plot,
       device = "png", width = 18, height = 12, units = "in")


# 4. Nutrient analysis ----------------------------------------------------

# Join nutrient data with metadata/distance and create custom metric
nutrients_meta_dist <- full_join(x = nutrients, y = metadata_dist, by = "site")


# 4.1 Phosphorus ----------------------------------------------------------

# Analyze phosphorus as a function of population intensity
phosphorus_PI_model <- lm(log10(mean_tp_mg_dm3) ~ log10(distance_weighted_population),
                          data = nutrients_meta_dist)

# View model results
summary(phosphorus_PI_model)

# Plot linear model
phosphorus_PI_plot <- ggplot(data = nutrients_meta_dist,
                             aes(x = log10(distance_weighted_population),
                                 y = log10(mean_tp_mg_dm3))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  ylab("log10([Total Phosphorus])") +
  xlab("log10(IDW Population)") +
  ggtitle("Phosphorus vs. IDW Population") +
  annotate(geom = "label", x = 3.0, y = -1.45,
           label = paste0("p-value: ",
                          round(summary(phosphorus_PI_model)$coefficients[2, 4], 3),
                          "\nR-squared: ",
                          round(summary(phosphorus_PI_model)$r.squared, 3))) +
  theme_minimal()

# Export plot
ggsave("phosphorus_PI_plot.png", phosphorus_PI_plot, device = "png",
       path = "../figures", width = 18, height = 12, units = "in")


# 4.2 Nitrate -------------------------------------------------------------

# Analyze nitrate as a function of population intensity
nitrate_PI_model <- lm(log10(mean_no3_mg_dm3) ~ log10(distance_weighted_population),
                       data = nutrients_meta_dist)

# View model results
summary(nitrate_PI_model)

# Plot linear model
nitrate_PI_plot <- ggplot(data = nutrients_meta_dist,
                          aes(x = log10(distance_weighted_population),
                              y = log10(mean_no3_mg_dm3))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  ylab("log10([Total Nitrate])") +
  xlab("log10(IDW Population)") +
  ggtitle("Nitrate vs. IDW Population") +
  annotate("label", x = 3.00, y = -0.8,
           label = paste0("p-value: ",
                          round(summary(nitrate_PI_model)$coefficients[2, 4], 3),
                          "\nR-squared: ",
                          round(summary(nitrate_PI_model)$r.squared, 3))) +
  theme_minimal()

# Export plot
ggsave("nitrate_PI_plot.png", nitrate_PI_plot, device = "png",
       path = "../figures", width = 18, height = 12, units = "in")


# 4.3 Ammonium ------------------------------------------------------------

# Analyze ammonium as a function of population intensity
ammonium_PI_model <- lm(log10(mean_nh4_mg_dm3) ~ log10(distance_weighted_population),
                        data = nutrients_meta_dist)

# View model results
summary(ammonium_PI_model)

# Plot linear model
ammonium_PI_plot <- ggplot(nutrients_meta_dist,
                           aes(x = log10(distance_weighted_population),
                               y = log10(mean_nh4_mg_dm3))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  ylab("log10([Total Ammonium])") +
  xlab("log10(IDW Population)") +
  ggtitle("Ammonium vs. IDW Population") +
  annotate(geom = "label", x = 3.00, y = -1.25,
           label = paste0("p-value: ",
                          round(summary(ammonium_PI_model)$coefficients[2, 4], 3),
                          "\nR-squared: ",
                          round(summary(ammonium_PI_model)$r.squared, 3))) +
  theme_minimal()

# Export plot
ggsave(filename = "../figures/ammonium_PI_plot.png", plot = ammonium_PI_plot,
       device = "png", width = 18, height = 12, units = "in")


# 5. Stable isotopes analysis ---------------------------------------------

# Join stable isotope data with metadata/distance and create custom metric
stable_isotopes_meta_dist <- full_join(x = stable_isotopes, y = metadata_dist,
                                       by = "site")

# 5.1 N15 -----------------------------------------------------------------

# Analyze N15 as a function of population intensity
n15_PI_model <- lm(log10(N15) ~ log10(distance_weighted_population),
                   data = stable_isotopes_meta_dist[stable_isotopes$Genus != "Periphyton", ])

# View model results
summary(n15_PI_model)

# Plot linear model
n15_PI_plot <- ggplot(data = stable_isotopes_meta_dist[stable_isotopes$Genus != "Periphyton", ],
                      aes(x = log10(distance_weighted_population), y = log10(N15))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  ylab(expression(paste("log10(",delta^{15}, "N (\u2030))"))) +
  xlab("log10(IDW Population)") +
  ggtitle(expression(paste(delta^{15}, "N \u2030  vs. IDW Population"))) +
  xlim(c(2.75, 3.75)) +
  annotate(geom = "label", x = 3.00, y = 0.89,
           label = paste0("p-value: ",
                          round(summary(n15_PI_model)$coefficients[2, 4], 3),
                          "\nR-squared: ",
                          round(summary(n15_PI_model)$r.squared, 3))) +
  theme_minimal()

# Export plot
ggsave(filename = "../figures/n15_PI_plot.png", plot = n15_PI_plot,
       device = "png", width = 18, height = 12, units = "in")


# 5.2 C13 -----------------------------------------------------------------

# Analyze C13 as a function of population intensity
c13_PI_model <- lm((C13) ~ log10(distance_weighted_population),
                   data = stable_isotopes_meta_dist[stable_isotopes$Genus != "Periphyton", ])

# View model results
summary(c13_PI_model)

# Plot linear model
c13_PI_plot <- ggplot(stable_isotopes_meta_dist[stable_isotopes$Genus != "Periphyton", ],
                      aes(x = log10(distance_weighted_population), y = (C13))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  ylab(expression(paste("log10(",delta^{13}, "C (\u2030))"))) +
  xlab("log10(IDW Population)") +
  ggtitle(expression(paste(delta^{13}, "C \u2030 vs. IDW Population"))) +
  xlim(c(2.75, 3.75)) +
  annotate(geom = "label", x = 3.00, y = -17,
           label = paste0("p-value: ",
                          round(summary(c13_PI_model)$coefficients[2, 4], 3),
                          "\nR-squared: ",
                          round(summary(c13_PI_model)$r.squared, 3))) +
  theme_minimal()

# Export plot
ggsave(filename = "../figures/c13_PI_plot.png", plot = c13_PI_plot,
       device = "png", width = 18, height = 12, units = "in")


# 6. Chlorophyll a analysis -----------------------------------------------

# Join Chl a data with metadata/distance and create custom metric
chlorophylla_meta_dist <- full_join(x = chlorophylla, y = metadata_dist,
                                    by = "site")

# Analyze chl a as a function of population intensity
chlorophylla_PI_model <- lm((mean_chlorophylla) ~ log10(distance_weighted_population),
                            data = chlorophylla_meta_dist)

# View model results
summary(chlorophylla_PI_model)

# Plot linear model
chlorophylla_PI_plot <- ggplot(data = chlorophylla_meta_dist,
                               aes(x = log10(distance_weighted_population),
                                   y = (mean_chlorophylla))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  ylab("Chlorophyll a") +
  xlab("log10(IDW Population)") +
  ggtitle("Chlorophyll a vs. IDW Population") +
  annotate(geom = "label", x = 4.00, y = 1.3,
           label = paste0("p-value: ",
                          round(summary(chlorophylla_PI_model)$coefficients[2, 4], 3),
                          "\nR-squared: ",
                          round(summary(chlorophylla_PI_model)$r.squared, 3))) +
  theme_minimal()

# Export plot
ggsave(filename = "../figures/chlorophylla_PI_plot.png",
       plot = chlorophylla_PI_plot, device = "png", width = 18, height = 12,
       units = "in")


# 7. Microplastics analysis -----------------------------------------------

# Format microplastics data before join
microplastics <- microplastics %>%
  group_by(site) %>%
  summarize(mean_total = mean(x = total_microplastics, na.rm = TRUE),
            mean_density = mean(x = density, na.rm = TRUE),
            mean_fragment_density = mean(x = fragment_density, na.rm = TRUE),
            mean_fiber_density = mean(x = fiber_density, na.rm = TRUE),
            mean_bead_density = mean(x = bead_density, na.rm = TRUE))

# Join microplastics data with metadata/distance and create custom metric
microplastics_meta_dist <- full_join(x = microplastics, y = metadata_dist,
                                     by = "site")


# 7.1 Mean total microplastics --------------------------------------------

# Analyze mean total microplastics as a function of population intensity
microplastics_total_PI_model <- lm((mean_total) ~
                                     log10(distance_weighted_population),
                                   data = microplastics_meta_dist)

# View model results
summary(microplastics_total_PI_model)

# Plot linear model
microplastics_total_PI_plot <- ggplot(data = microplastics_meta_dist,
                                      aes(x = log10(distance_weighted_population),
                                          y = (mean_total))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  ylab("Mean Total Microplastics") +
  xlab("log10(IDW Population)") +
  ggtitle("Total Microplastics vs. IDW Population") +
  annotate(geom = "label", x = 3.00, y = 5,
           label = paste0("p-value: ",
                          round(summary(microplastics_total_PI_model)$coefficients[2, 4], 3),
                          "\nR-squared: ",
                          round(summary(microplastics_total_PI_model)$r.squared, 3))) +
  theme_minimal()

# Export plot
ggsave(filename = "../figures/microplastics_total_PI_plot.png",
       plot = microplastics_total_PI_plot, device = "png", width = 18,
       height = 12, units = "in")


# 7.2 Mean microplastic density -------------------------------------------

# Analyze mean microplastic density as a function of population intensity
microplastics_density_PI_model <- lm((mean_density) ~
                                       log10(distance_weighted_population),
                                     data = microplastics_meta_dist)

# View model results
summary(microplastics_density_PI_model)

# Plot linear model
microplastics_density_PI_plot <- ggplot(data = microplastics_meta_dist,
                                        aes(x = log10(distance_weighted_population),
                                            y = (mean_density))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  ylab("Mean Microplastic Density") +
  xlab("log10(IDW Population)") +
  ggtitle("Microplastics Density vs. IDW Population") +
  annotate(geom = "label", x = 3.65, y = 0.004,
           label = paste0("p-value: ",
                          round(summary(microplastics_density_PI_model)$coefficients[2, 4], 3),
                          "\nR-squared: ",
                          round(summary(microplastics_density_PI_model)$r.squared, 3))) +
  theme_minimal()

# Export plot
ggsave(filename = "../figures/microplastics_density_PI_plot.png",
       plot = microplastics_density_PI_plot, device = "png", width = 18,
       height = 12, units = "in")


# 8. Combine plots --------------------------------------------------------

ggarrange(ppcp_PI_plot, n15_PI_plot, phosphorus_PI_plot, chlorophylla_PI_plot,
          nitrate_PI_plot, microplastics_total_PI_plot, ammonium_PI_plot,
          microplastics_density_PI_plot, ncol = 2, nrow = 4, labels = "AUTO") %>%
  ggexport(filename = "../figures/combined_plot.png",
           height = 1900, width = 1200, res = 120)
