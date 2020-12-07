# This script aggregates stable isotopes with sewage
# indicator data in Lake Baikal to produce a stable isotope
# biplot for the associated manuscript. 

library(tidyverse)
library(viridis)
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

# 1. Load the data --------------------------------------------------------

stable_isotopes <- read.csv("../cleaned_data/stable_isotopes.csv",
                            header = TRUE)


# 2. Define IDW population groupings --------------------------------------

low <- c("BGO-1", "BGO-2", "BGO-3", "KD-1", "KD-2", "MS-1")
mod <- c("BK-2", "BK-3", "SM-1")
high <- c("BK-1", "EM-1", "LI-3", "LI-2", "LI-1")


# 3. Build and export stable isotope biplot -------------------------------

# To aggregate the data necessary to creat the stable isotopes 
# biplot, we first find the mean and then the standard deviaton
# of each C13 and N15 value for each taxon-IDW population grouping.

foodweb_data <- stable_isotopes %>%
  mutate(idw_group = ifelse(site %in% c(low, mod), "Low/Mod", NA),
         idw_group = ifelse(site %in% high, "High", idw_group),
         Genus = as.character(Genus),
         Genus = ifelse(test = Genus != "Periphyton",
                        yes = paste0("*", as.character(Genus), "*"),
                        no = Genus)) %>%
  unite(taxon_idw_pop, c("Genus", "idw_group"), sep = " + ") %>%
  filter(!grepl("NA", taxon_idw_pop)) %>%
  select(taxon_idw_pop, C13, N15) %>%
  group_by(taxon_idw_pop) %>%
  summarize(mean_N15 = mean(N15),
            mean_C13 = mean(C13),
            sd_N15 = sd(N15),
            sd_C13 = sd(C13))


# Step 4: Build the stable isotopes biplot -------------------------------


foodweb_plot <- ggplot(data = foodweb_data, 
                       aes(mean_C13, mean_N15, 
                           color = taxon_idw_pop)) +
  geom_point(size = 6) +
  geom_errorbar(aes(ymin = mean_N15-sd_N15, 
                    ymax = mean_N15+sd_N15), size = 1) +
  geom_errorbarh(aes(xmin = mean_C13-sd_C13, 
                     xmax = mean_C13+sd_C13), size = 1) +
  labs(color = "Taxon + IDW Population Group") +
  scale_color_viridis_d(begin = 0,
                        end = 0.9,
                        option = "magma") +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  theme_minimal() +
  theme(text = element_text(size = 16),
        legend.text = element_markdown())

ggsave("../figures/stable_isotopes_biplot.png", foodweb_plot,
       device = "png", width = 9, height = 6, units = "in")
