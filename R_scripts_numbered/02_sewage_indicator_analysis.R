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
library(broom)
library(ggtext)

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

permute_data_analytics <- function(data, metric, full_model, metric_plot_title, log_transform_response = TRUE){
  for(i in 1:5000){
    
    # First permute the response variable. The variable is supplied by the user.
    permuted_data <- data %>%
      select(paste(metric)) %>%
      rename("permuted_metric" = paste(metric)) %>%
      sample_frac(size = 1) %>%
      cbind(., data)
    
    # If the user has specified to log-transform the variable, this step will 
    # actually perform the log-transform. 
    if(log_transform_response){
      permuted_model <- lm(log10(permuted_metric) ~ log10(distance_weighted_population),
                           data = permuted_data)
    }else{
      permuted_model <- lm((permuted_metric) ~ log10(distance_weighted_population),
                           data = permuted_data)
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
    theme_minimal()
  
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
    theme(plot.title = element_markdown())
  
  # The two resulting figures are returned as a list. 
  return(list(permuted_plot_pval, permuted_plot_rsquared))
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

# Site site_information
site_information <- read.csv(file = "../cleaned_data/site_information.csv",
                     header = TRUE, stringsAsFactors = FALSE)

# Site distance data
distance <- read.csv(file = "../cleaned_data/distance_weighted_population_metrics.csv",
                     header = TRUE, stringsAsFactors = FALSE)

# Join site site_information with distance data
site_information_dist <- full_join(x = site_information, y = distance, by = "site")

# 3. PPCP analysis --------------------------------------------------------

# Join PPCP data with site_information/distance and create two custom metrics
ppcp_site_info_dist <- full_join(x = ppcp, y = site_information_dist, by = "site")

# Analyze total PPCPs as a function of population intensity
ppcp_PI_model <- lm(log10(ppcp_sum) ~ log10(distance_weighted_population),
                    data = ppcp_site_info_dist)

# View model results
summary(ppcp_PI_model)

# Plot linear model
ppcp_PI_plot <- ggplot(data = ppcp_site_info_dist,
                       aes(x = (distance_weighted_population),
                           y =  (ppcp_sum))) +
  geom_point(size = 6) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_y_log10()+
  scale_x_log10() +
  ylab("Total PPCP (\u00b5g/L)") +
  xlab("IDW Population (people)") +
  ggtitle("PPCP vs. IDW Population") +
  geom_richtext(data = glance(ppcp_PI_model), 
                mapping = aes(label = paste0("p-value: ",
                                             round(p.value, 3),
                                             "<br>R<sup>2</sup>: ",
                                             round(r.squared, 3),
                                             "<br> N = ", nobs), x = 6000, y = 0.1), size = 10) +
  theme_minimal() +
  theme(strip.text = element_markdown(),
        axis.text = element_text(size = 28),
        axis.title = element_text(size = 32),
        plot.title = element_text(size = 34))

# Export plot
ggsave(filename = "../figures/ppcp_PI_plot.png", plot = ppcp_PI_plot,
       device = "png", width = 18, height = 12, units = "in")

ppcp_permute_plots <- permute_data_analytics(data = ppcp_site_info_dist, 
                                             metric = "ppcp_sum", 
                                             metric_plot_title = "PPCP vs. IDW Population", 
                                             full_model = ppcp_PI_model)


# 4. Nutrient analysis ----------------------------------------------------

# Join nutrient data with site_information/distance and create custom metric
nutrients_site_info_dist <- full_join(x = nutrients, y = site_information_dist, by = "site")


# 4.1 Phosphorus ----------------------------------------------------------

# Analyze phosphorus as a function of population intensity
phosphorus_PI_model <- lm(log10(mean_tp_mg_dm3) ~ log10(distance_weighted_population),
                          data = nutrients_site_info_dist)

# View model results
summary(phosphorus_PI_model)

# Plot linear model
phosphorus_PI_plot <- ggplot(data = nutrients_site_info_dist,
                             aes(x = (distance_weighted_population),
                                 y = (mean_tp_mg_dm3))) +
  geom_point(size = 6) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_y_log10()+
  scale_x_log10() +
  ylab("Phosphorus (mg/L)") +
  xlab("IDW Population (people)") +
  ggtitle("Total Phosphorus vs. IDW Population") +
  geom_richtext(data = glance(phosphorus_PI_model), 
                mapping = aes(label = paste0("p-value: ",
                                             round(p.value, 3),
                                             "<br>R<sup>2</sup>: ",
                                             round(r.squared, 3),
                                             "<br> N = ", nobs), x = 6000, y = 0.04), size = 10) +
  theme_minimal() +
  theme(strip.text = element_markdown(),
        axis.text = element_text(size = 28),
        axis.title = element_text(size = 32),
        plot.title = element_text(size = 34))

# Export plot
ggsave("phosphorus_PI_plot.png", phosphorus_PI_plot, device = "png",
       path = "../figures", width = 18, height = 12, units = "in")

phosphorus_permute_plots <- permute_data_analytics(data = nutrients_site_info_dist, 
                                                   metric = "mean_tp_mg_dm3", 
                                                   metric_plot_title = "Total Phosphorus vs. IDW Population", 
                                                   full_model = phosphorus_PI_model)


# 4.2 Nitrate -------------------------------------------------------------

# Analyze nitrate as a function of population intensity
nitrate_PI_model <- lm(log10(mean_no3_mg_dm3) ~ log10(distance_weighted_population),
                       data = nutrients_site_info_dist)

# View model results
summary(nitrate_PI_model)

# Plot linear model
nitrate_PI_plot <- ggplot(data = nutrients_site_info_dist,
                          aes(x = (distance_weighted_population),
                              y = (mean_no3_mg_dm3))) +
  geom_point(size = 6) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_y_log10()+
  scale_x_log10() +
  ylab("Nitrate (mg/L)") +
  xlab("IDW Population (people)") +
  ggtitle("Nitrate vs. IDW Population") +
  geom_richtext(data = glance(nitrate_PI_model), 
                mapping = aes(label = paste0("p-value: ",
                                             round(p.value, 3),
                                             "<br>R<sup>2</sup>: ",
                                             round(r.squared, 3),
                                             "<br> N = ", nobs), x = 6000, y = 0.15), size = 10) +
  theme_minimal() +
  theme(strip.text = element_markdown(),
        axis.text = element_text(size = 28),
        axis.title = element_text(size = 32),
        plot.title = element_text(size = 34))

# Export plot
ggsave("nitrate_PI_plot.png", nitrate_PI_plot, device = "png",
       path = "../figures", width = 18, height = 12, units = "in")

nitrate_permute_plots <- permute_data_analytics(data = nutrients_site_info_dist, 
                                                   metric = "mean_no3_mg_dm3", 
                                                   metric_plot_title = "Nitrate vs. IDW Population", 
                                                   full_model = nitrate_PI_model)


# 4.3 Ammonium ------------------------------------------------------------

# Analyze ammonium as a function of population intensity
ammonium_PI_model <- lm(log10(mean_nh4_mg_dm3) ~ log10(distance_weighted_population),
                        data = nutrients_site_info_dist)

# View model results
summary(ammonium_PI_model)

# Plot linear model
ammonium_PI_plot <- ggplot(nutrients_site_info_dist,
                           aes(x = (distance_weighted_population),
                               y = (mean_nh4_mg_dm3))) +
  geom_point(size = 6) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_y_log10()+
  scale_x_log10() +
  ylab("Ammonium (mg/L)") +
  xlab("IDW Population (people)") +
  ggtitle("Ammonium vs. IDW Population") +
  geom_richtext(data = glance(ammonium_PI_model), 
                mapping = aes(label = paste0("p-value: ",
                                             round(p.value, 3),
                                             "<br>R<sup>2</sup>: ",
                                             round(r.squared, 3),
                                             "<br> N = ", nobs), x = 6000, y = 0.03), size = 10) +
  theme_minimal() +
  theme(strip.text = element_markdown(),
        axis.text = element_text(size = 28),
        axis.title = element_text(size = 32),
        plot.title = element_text(size = 34))

# Export plot
ggsave(filename = "../figures/ammonium_PI_plot.png", plot = ammonium_PI_plot,
       device = "png", width = 18, height = 12, units = "in")

ammonium_permute_plots <- permute_data_analytics(data = nutrients_site_info_dist, 
                                                metric = "mean_nh4_mg_dm3", 
                                                metric_plot_title = "Ammonium vs. IDW Population", 
                                                full_model = ammonium_PI_model)

# 5. Stable isotopes analysis ---------------------------------------------

# Join stable isotope data with site_information/distance and create custom metric
stable_isotopes_site_info_dist <- inner_join(x = stable_isotopes, y = site_information_dist,
                                       by = "site")

# 5.1 N15 -----------------------------------------------------------------

# Analyze N15 as a function of population intensity
n15_PI_model <- lm((N15) ~ log10(distance_weighted_population),
                   data = stable_isotopes_site_info_dist[stable_isotopes_site_info_dist$Genus != "Periphyton", ])

# View model results
summary(n15_PI_model)

# Plot linear model
n15_PI_plot <- ggplot(data = stable_isotopes_site_info_dist[stable_isotopes_site_info_dist$Genus != "Periphyton", ],
                      aes(x = (distance_weighted_population), y = (N15))) +
  geom_point(size = 6) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_x_log10(limits = c(min(site_information_dist[,c("distance_weighted_population")]),
                           max(site_information_dist[,c("distance_weighted_population")]))) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab("IDW Population (people)") +
  ggtitle(expression(paste(delta^{15}, "N \u2030  vs. IDW Population"))) +
  geom_richtext(data = glance(n15_PI_model), 
                mapping = aes(label = paste0("p-value: ",
                                             round(p.value, 3),
                                             "<br>R<sup>2</sup>: ",
                                             round(r.squared, 3),
                                             "<br> N = ", nobs), x = 2000, y = 7.0), size = 10) +
  theme_minimal() +
  theme(strip.text = element_markdown(),
        axis.text = element_text(size = 28),
        axis.title = element_text(size = 32),
        plot.title = element_text(size = 34))
  

# Export plot
ggsave(filename = "../figures/n15_PI_plot.png", plot = n15_PI_plot,
       device = "png", width = 18, height = 12, units = "in")

n15_permute_plots <- permute_data_analytics(data = stable_isotopes_site_info_dist, 
                                                 metric = "N15", 
                                                 metric_plot_title = "N \u2030 vs. IDW Population", 
                                                 full_model = n15_PI_model, 
                                                 log_transform_response = FALSE)


# 5.2 C13 -----------------------------------------------------------------

# Analyze C13 as a function of population intensity
c13_PI_model <- lm((C13) ~ log10(distance_weighted_population),
                   data = stable_isotopes_site_info_dist[stable_isotopes$Genus != "Periphyton", ])

# View model results
summary(c13_PI_model)

# Plot linear model
c13_PI_plot <- ggplot(stable_isotopes_site_info_dist[stable_isotopes$Genus != "Periphyton", ],
                      aes(x = (distance_weighted_population), y = (C13))) +
  geom_point(size = 6) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_x_log10(limits = c(min(site_information_dist[,c("distance_weighted_population")]),
                           max(site_information_dist[,c("distance_weighted_population")]))) +
  ylab(expression(paste(delta^{13}, "C (\u2030)"))) +
  xlab("IDW Population (people)") +
  ggtitle(expression(paste(delta^{13}, "C \u2030  vs. IDW Population"))) +
  geom_richtext(data = glance(c13_PI_model), 
                mapping = aes(label = paste0("p-value: ",
                                             round(p.value, 3),
                                             "<br>R<sup>2</sup>: ",
                                             round(r.squared, 3),
                                             "<br> N = ", nobs), x = 2000, y = -17), size = 10) +
  theme_minimal() +
  theme(strip.text = element_markdown(),
        axis.text = element_text(size = 28),
        axis.title = element_text(size = 32),
        plot.title = element_text(size = 34))

# Export plot
ggsave(filename = "../figures/c13_PI_plot.png", plot = c13_PI_plot,
       device = "png", width = 18, height = 12, units = "in")

c13_permute_plots <- permute_data_analytics(data = stable_isotopes_site_info_dist, 
                                            metric = "C13", 
                                            metric_plot_title = "C \u2030 vs. IDW Population", 
                                            full_model = c13_PI_model, 
                                            log_transform_response = FALSE)

# 6. Chlorophyll a analysis -----------------------------------------------

# Join Chl a data with site_information/distance and create custom metric
chlorophylla_site_info_dist <- full_join(x = chlorophylla, y = site_information_dist,
                                    by = "site")

# Analyze chl a as a function of population intensity
chlorophylla_PI_model <- lm((mean_chlorophylla) ~ log10(distance_weighted_population),
                            data = chlorophylla_site_info_dist)

# View model results
summary(chlorophylla_PI_model)

# Plot linear model
chlorophylla_PI_plot <- ggplot(data = chlorophylla_site_info_dist,
                               aes(x = (distance_weighted_population),
                                   y = (mean_chlorophylla))) +
  geom_point(size = 6) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_y_log10()+
  scale_x_log10() +
  ylab("Chlorophyll a (mg/L)") +
  xlab("IDW Population (people)") +
  ggtitle("Chlorophyll a vs. IDW Population") +
  geom_richtext(data = glance(chlorophylla_PI_model), 
                mapping = aes(label = paste0("p-value: ",
                                             round(p.value, 3),
                                             "<br>R<sup>2</sup>: ",
                                             round(r.squared, 3),
                                             "<br> N = ", nobs), x = 6000, y = 1), size = 10) +
  theme_minimal() +
  theme(strip.text = element_markdown(),
        axis.text = element_text(size = 28),
        axis.title = element_text(size = 32),
        plot.title = element_text(size = 34))

# Export plot
ggsave(filename = "../figures/chlorophylla_PI_plot.png",
       plot = chlorophylla_PI_plot, device = "png", width = 18, height = 12,
       units = "in")

chlorophylla_permute_plots <- permute_data_analytics(data = chlorophylla_site_info_dist, 
                                            metric = "mean_chlorophylla", 
                                            metric_plot_title = "Chlorophyll a vs. IDW Population", 
                                            full_model = chlorophylla_PI_model, 
                                            log_transform_response = FALSE)

# 7. Microplastics analysis -----------------------------------------------

# Format microplastics data before join
microplastics <- microplastics %>%
  group_by(site) %>%
  summarize(mean_total = mean(x = total_microplastics, na.rm = TRUE),
            mean_density = mean(x = density, na.rm = TRUE),
            mean_fragment_density = mean(x = fragment_density, na.rm = TRUE),
            mean_fiber_density = mean(x = fiber_density, na.rm = TRUE),
            mean_bead_density = mean(x = bead_density, na.rm = TRUE))

# Join microplastics data with site_information/distance and create custom metric
microplastics_site_info_dist <- full_join(x = microplastics, y = site_information_dist,
                                     by = "site")


# 7.1 Mean total microplastics --------------------------------------------

# Analyze mean total microplastics as a function of population intensity
microplastics_total_PI_model <- lm((mean_total) ~
                                     log10(distance_weighted_population),
                                   data = microplastics_site_info_dist)

# View model results
summary(microplastics_total_PI_model)

## During review, a reviewer suggested we trying a poisson regression for MP

summary(glm((mean_total) ~ log10(distance_weighted_population),
   family = "poisson",
   data = microplastics_site_info_dist))

# Plot linear model
microplastics_total_PI_plot <- ggplot(data = microplastics_site_info_dist,
                                      aes(x = (distance_weighted_population),
                                          y = (mean_total))) +
  geom_point(size = 6) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_y_log10()+
  scale_x_log10() +
  ylab("Microplastics (MPs)") +
  xlab("IDW Population (people)") +
  ggtitle("Total Microplastics vs. IDW Population") +
  geom_richtext(data = glance(microplastics_total_PI_model), 
                mapping = aes(label = paste0("p-value: ",
                                             round(p.value, 3),
                                             "<br>R<sup>2</sup>: ",
                                             round(r.squared, 3),
                                             "<br> N = ", nobs), x = 6000, y = 1), size = 10) +
  theme_minimal() +
  theme(strip.text = element_markdown(),
        axis.text = element_text(size = 28),
        axis.title = element_text(size = 32),
        plot.title = element_text(size = 34))

# Export plot
ggsave(filename = "../figures/microplastics_total_PI_plot.png",
       plot = microplastics_total_PI_plot, device = "png", width = 18,
       height = 12, units = "in")

microplastics_total_permute_plots <- permute_data_analytics(data = microplastics_site_info_dist, 
                                                    metric = "mean_total", 
                                                    metric_plot_title = "Total Microplastics vs. IDW Population", 
                                                    full_model = microplastics_total_PI_model, 
                                                    log_transform_response = FALSE)

# 7.2 Mean microplastic density -------------------------------------------

# Analyze mean microplastic density as a function of population intensity
microplastics_density_PI_model <- lm((mean_density) ~
                                       log10(distance_weighted_population),
                                     data = microplastics_site_info_dist)

# View model results
summary(microplastics_density_PI_model)

## During review, a reviewer suggested we trying a poisson regression 
## for MP Density

summary(glm((mean_density) ~ log10(distance_weighted_population),
            family = "poisson",
            data = microplastics_site_info_dist))

# Plot linear model
microplastics_density_PI_plot <- ggplot(data = microplastics_site_info_dist,
                                        aes(x = (distance_weighted_population),
                                            y = (mean_density))) +
  geom_point(size = 6) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_y_log10()+
  scale_x_log10() +
  ylab("Microplastics Density (MPs/L)") +
  xlab("IDW Population (people)") +
  ggtitle("Microplastics Density vs. IDW Population") +
  geom_richtext(data = glance(microplastics_density_PI_model), 
                mapping = aes(label = paste0("p-value: ",
                                             round(p.value, 3),
                                             "<br>R<sup>2</sup>: ",
                                             round(r.squared, 3),
                                             "<br> N = ", nobs), x = 6000, y = 0.001), size = 10) +
  theme_minimal() +
  theme(strip.text = element_markdown(),
        axis.text = element_text(size = 28),
        axis.title = element_text(size = 32),
        plot.title = element_text(size = 34))

# Export plot
ggsave(filename = "../figures/microplastics_density_PI_plot.png",
       plot = microplastics_density_PI_plot, device = "png", width = 18,
       height = 12, units = "in")

microplastics_density_permute_plots <- permute_data_analytics(data = microplastics_site_info_dist, 
                                                            metric = "mean_density", 
                                                            metric_plot_title = "Microplastics Density vs. IDW Population", 
                                                            full_model = microplastics_density_PI_model, 
                                                            log_transform_response = FALSE)


# 8. Combine plots --------------------------------------------------------

ggarrange(nitrate_PI_plot, ammonium_PI_plot, 
          phosphorus_PI_plot, chlorophylla_PI_plot,
          ppcp_PI_plot, n15_PI_plot, 
          microplastics_total_PI_plot, microplastics_density_PI_plot,
          ncol = 2, nrow = 4, labels = "AUTO", 
          font.label = list(size = 36)) %>%
  ggexport(filename = "../figures/combined_plot.png",
           height = 4000, width = 2400, res = 120)

combined_permuted_plots <- c(ppcp_permute_plots, n15_permute_plots, phosphorus_permute_plots, chlorophylla_permute_plots,
                             nitrate_permute_plots, ammonium_permute_plots, microplastics_total_permute_plots, 
                             microplastics_density_permute_plots)

ggarrange(plotlist =  combined_permuted_plots, ncol = 2, nrow = 8) %>%
  ggexport(filename = "../figures/permuted_combined_plot.png",
           height = 1900, width = 1600, res = 120)
