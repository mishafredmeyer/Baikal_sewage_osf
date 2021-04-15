# This script is meant to perform multivariate and some univariate analysis of
# macroinvertebrates and periphyton community composition along Lake Baikal's
# southwestern shoreline. These community compositions are then related to
# sewage indicator concentrations and inverse distance weight population
# along the shoreline.

# 1. Load packages --------------------------------------------------------

library(tidyverse)
library(ggplot2)
library(viridisLite)
library(vegan)
library(factoextra)
library(cluster)
library(ggpubr)
library(ggrepel)
library(pvclust)

# Pull in function for plotting correlations
source("panel_cor_function.R")

# Check if figures directory exists 
# If not, create the figures directory
sub_dir <- "figures"
output_dir <- file.path(here::here(), sub_dir)

if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Dir 'figures' already exists!")
}

# 2. Load the data --------------------------------------------------------

# Invertebrate data
invertebrates <- read.csv("../cleaned_data/invertebrates.csv",
                          header = TRUE, stringsAsFactors = FALSE)

# Periphyton data
periphyton <- read.csv("../cleaned_data/periphyton.csv",
                       header = TRUE, stringsAsFactors = FALSE)

# Site metadata
metadata <- read.csv("../cleaned_data/metadata.csv",
                     header = TRUE, stringsAsFactors = FALSE)

# Site distance data
distance <- read.csv("../cleaned_data/distance_weighted_population_metrics.csv",
                     header = TRUE, stringsAsFactors = FALSE)

# Join site metadata with distance data
metadata_dist <- full_join(x = metadata, y = distance, by = "site")

# PPCP data
ppcp <- read.csv("../cleaned_data/ppcp.csv",
                 header = TRUE, stringsAsFactors = FALSE)

# Join PPCP data with metadata/distance and select sites of interest
ppcp_meta_dist <- full_join(x = ppcp, y = metadata_dist, by = "site") %>%
  filter(!(site %in% c("OS-1", "OS-2", "OS-3")))


# Define low/moderate/high sites
# These groupings are derived from the inverse distance weighted population metric,
# where sites with "low" IDW population are recorded in the "low" vector, moderate in "mod",
# and high in the "high" vector. These same delineations are used in both the periphyton
# and invertebrate analysis.

low <- c("BGO-1", "BGO-2", "BGO-3", "KD-1", "KD-2", "MS-1")
mod <- c("BK-2", "BK-3", "SM-1")
high <- c("BK-1", "EM-1", "LI-3", "LI-2", "LI-1")

# Small Function to add some plot labels at the end
add_label <- function(xfrac, yfrac, label, pos = 4, ...) {
  u <- par("usr")
  x <- u[1] + xfrac * (u[2] - u[1])
  y <- u[4] - yfrac * (u[4] - u[3])
  text(x, y, label, pos = pos, ...)
}

# 3. Periphyton analysis ---------------------------------------------------

# Join periphyton data with metadata/distance and create custom metric
periphyton_meta_dist <- full_join(x = periphyton, y = ppcp_meta_dist,
                                  by = "site")


# 3.1 Univariate analysis -------------------------------------------------

# Convert periphyton data to long format and add total count
periphyton_meta_dist_long <- periphyton_meta_dist %>%
  gather(key = taxon, value = count, desmidales:ulothrix) %>%
  filter(!(taxon %in% c("desmidales", "pediastrum"))) %>%
  group_by(site) %>%
  mutate(total_count = sum(count))

# Rework Site column as a factor
# Remove data with Site = NA
periphyton_long_clean <- periphyton_meta_dist_long %>% 
  filter(!is.na(site)) %>%
  mutate(site = factor(x = site,
                       levels = c("BGO-2",  "BGO-3", "KD-1",
                                  "KD-2", "BGO-1", "MS-1",
                                  "BK-3", "BK-2", "SM-1",
                                  "BK-1", "EM-1", "LI-1",
                                  "LI-2", "LI-3")))

# Plot periphyton counts as a function of site and taxa
periphyton_meta_dist_plot <- ggplot(data = periphyton_long_clean) +
  geom_bar(aes(x = site, y = total_count), fill = "grey80", stat = "identity") +
  geom_bar(aes(x = site, y = count, fill = taxon), stat = "identity") +
  scale_fill_viridis_d(option = "inferno") +
  facet_wrap(~ taxon) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Number of cells") +
  xlab("Site (Arranged by increasing inverse distance-weighted population)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(size = 20),
        strip.text.x = element_text(size = 27),
        strip.background = element_rect(fill = "white"),
        panel.background = element_rect(color = "black"),
        axis.title = element_text(size = 30),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(margin = margin(0, 20, 0, 0)),
        axis.title.x = element_text(size = 24, margin = margin(20, 0, 0, 0)),
        legend.text = element_text(size = 16),
        axis.ticks.x = element_blank())

# This plot is not included in the associated manuscript but is included to
# offer insight for multivariate analyeses.

ggsave(filename = "periphyton_univariate.png", plot = periphyton_meta_dist_plot,
       device = "png", path = "../figures/", width = 11, height = 8.5,
       units = "in")


# 3.2 Multivariate analysis -----------------------------------------------

# Clean dataset and add IDW_pop_group
periphyton_meta_dist_wide <- periphyton_meta_dist %>%
  filter(!(site %in% c("OS-1", "OS-2", "OS-3"))) %>%
  mutate(IDW_pop_group = ifelse(test = site %in% high,
                                yes = "High", no = "NULL"),
         IDW_pop_group = ifelse(test = (site %in% c(low,mod)), 
                                yes = "Low/Mod", no = IDW_pop_group))

# Define community for NMDS
peri_community <- periphyton_meta_dist_wide %>%
  select(diatom:ulothrix)

# Run NMDS
periphyton_nmds <- metaMDS(comm = peri_community,  try = 100, distance = "bray")
periphyton_nmds

# Pull scores from NMDS and add site data
data_scores <- as.data.frame(scores(x = periphyton_nmds))
data_scores$site <- periphyton_meta_dist_wide %>%
  pull(site)

# Join scores with PPCP data and code into PI groups
data_scores <- inner_join(x = data_scores, y = ppcp_meta_dist, 
                          by = "site") %>%
  mutate(IDW_pop_group = ifelse(test = site %in% high,
                                yes = "High", no = "NULL"),
         IDW_pop_group = ifelse(test = (site %in% c(low,mod)), 
                                yes = "Low/Mod", no = IDW_pop_group))

# Rework IDW_pop_group column as a factor
data_scores$IDW_pop_group <- factor(x = data_scores$IDW_pop_group,
                                    levels = c("High", "Low/Mod"))

# Pull species scores from NMDS
species_scores <- as.data.frame(scores(x = periphyton_nmds, display = "species"))
species_scores$species <- rownames(species_scores)

# Plot the NMDS
periphyton_IDW_pop_group_plot <- ggplot() +
  geom_point(data = data_scores,
             aes(x = NMDS1, y = NMDS2, 
                 size = log10(distance_weighted_population),
                 color = IDW_pop_group)) +
  scale_size_continuous(range = c(12, 28), guide = FALSE) +
  scale_color_manual(values = inferno(15)[c(3, 8, 11)],
                     name = "IDW Population Grouping") +
  geom_text_repel(data = species_scores %>% 
                    filter(species %in% c("spirogyra", "ulothrix", "diatom")) %>%
                    mutate(species = ifelse(species == "spirogyra", 
                                            yes = "Spirogyra",
                                            no = species),
                           species = ifelse(species == "ulothrix", 
                                            yes = "Ulothrix",
                                            no = species),
                           species = ifelse(species %in% c("Ulothrix", "Spirogyra"), 
                                            yes = paste0("italic(", species, ")"),
                                            no = species)), 
                  aes(x = NMDS1, y = NMDS2, label = species), 
            size = 9, parse = TRUE) + 
  guides(colour = guide_legend(override.aes = list(size = 10))) +
  annotate("label", x = 0, y = -0.35, size = 10,
           label = paste("Stress: ", round(periphyton_nmds$stress, digits = 3))) +
  theme_minimal() +
  theme(legend.position = "right",
        legend.key=element_blank(),
        text = element_text(size = 24))

periphyton_IDW_pop_group_plot

# Export plot
ggsave(filename = "../figures/periphyton_IDW_pop_group_plot.png",
       plot = periphyton_IDW_pop_group_plot, device = "png", height = 10,
       width = 20, dpi = 300)

# Vizualize optimum cluster number with wss method
peri_cluster_wss <- fviz_nbclust(x = as.matrix(vegdist(x = peri_community, method = "bray",
                                                       diag = TRUE, upper = TRUE)),
                                 FUNcluster = cluster::pam, method = "wss")
peri_cluster_wss

peri_cluster_sil <- fviz_nbclust(x = as.matrix(vegdist(x = peri_community, method = "bray",
                                                       diag = TRUE, upper = TRUE)),
                                 FUNcluster = cluster::pam, method = "silhouette")
peri_cluster_sil

peri_hclust <- pvclust(as.matrix(vegdist(x = peri_community, method = "bray",
                                           diag = TRUE, upper = TRUE)), method.hclust = "median", 
                         method.dist = "correlation")

plot(peri_hclust)
pvrect(peri_hclust, alpha = 0.90)

# Export plots
# These plots correspond with figures S1a & S1c in the associated manuscript.
ggsave(filename = "../figures/periphyton_pam_analysis_wss.png",
       plot = peri_cluster_wss, device = "png", height = 10, width = 20,
       dpi = 300)

ggsave(filename = "../figures/periphyton_pam_analysis_sil.png",
       plot = peri_cluster_sil, device = "png", height = 10, width = 20,
       dpi = 300)

ggsave(filename = "../figures/periphyton_hclust_analysis.png",
       plot = plot(peri_hclust), device = "png", height = 10, width = 20,
       dpi = 300)

# Run PERMANOVA
adonis2(formula = periphyton_meta_dist_wide[, 2:7]
       ~ log10(periphyton_meta_dist_wide[, 21]),
       data = periphyton_meta_dist_wide,
       method = "bray", permutations = 999)

adonis2(formula = periphyton_meta_dist_wide[, 2:7]
        ~ periphyton_meta_dist_wide[, 22],
        data = periphyton_meta_dist_wide,
        method = "bray", permutations = 999)

simper_results <- simper(comm = periphyton_meta_dist_wide[, 2:7],
                         group = periphyton_meta_dist_wide[, 22],
                         permutations = 999)

summary(simper_results)


# 4. Invertebrate analysis ---------------------------------------------------


# 4.1 Univariate analysis -------------------------------------------------

# Convert invertebrate data to long format and add total count
invertebrate_meta_dist <- full_join(x = invertebrates, y = ppcp_meta_dist,
                                    by = "site")

# Define amphipod genera
amphipods <- c("Eulimnogammarus", "Poekilogammarus", "Pallasea",
               "Cryptoropus", "Brandtia")
# Define mollusc genera
molluscs <- c("Acroloxidae", "Baicaliidae",  "Benedictidae", "Planorbidae",
              "Valvatidae")

# Make invertebrate data long format, sum counts, and group genera
invertebrates_long <- invertebrate_meta_dist %>%
  select(site:Valvatidae) %>%
  gather(key = taxon, value = count, Acroloxidae:Valvatidae) %>%
  filter(!grepl("juvenile", taxon, ignore.case = TRUE)) %>%
  separate(col = taxon, into = c("Genus", "Species", "Subspecies"), sep = "_") %>%
  group_by(site, Genus) %>%
  summarize(total_Genus = sum(count)) %>%
  ungroup() %>%
  group_by(site) %>%
  mutate(total_site = sum(total_Genus)) %>%
  unique() %>%
  mutate(grouping = ifelse(test = Genus %in% amphipods,
                           yes = "Amphipod", no = NA),
         grouping = ifelse(test = Genus %in% molluscs,
                           yes = "Mollusc", no = grouping),
         grouping = ifelse(test = Genus == "Asellidae",
                           yes = "Isopod", no = grouping),
         grouping = ifelse(test = Genus == "Caddisflies",
                           yes = "Caddisflies", no = grouping),
         grouping = ifelse(test = Genus == "flatworms",
                           yes = "Planaria", no = grouping),
         grouping = ifelse(test = Genus == "Leeches",
                           yes = "Hirudinea", no = grouping)) %>%
  filter(!is.na(grouping))

# Rework Site column as a factor
invertebrates_long$site <- factor(x = invertebrates_long$site,
                                  levels = c("BGO-2",  "BGO-3", "KD-1",
                                             "KD-2", "BGO-1", "MS-1",
                                             "BK-3", "BK-2", "SM-1",
                                             "BK-1", "EM-1", "LI-1",
                                             "LI-2", "LI-3"))

# Plot invert counts by site and genus
ggplot(invertebrates_long) +
  geom_bar(aes(x = site, y = total_site), alpha = 0.5, stat = "identity") +
  geom_bar(aes(x = site, y = total_Genus, fill = Genus), stat = "identity") +
  facet_wrap(~ Genus) +
  xlab("Site (Arranged by increasing inverse distance-weighted population)") +
  ylab("Number of Individuals") +
  theme_classic() +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        strip.background = element_rect(fill = "white"),
        panel.background = element_rect(color = "black"),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(margin = margin(0, 20, 0, 0)),
        axis.title.x = element_text(margin = margin(20, 0, 0, 0)),
        legend.text = element_text(size = 16))

# Create a cleaned up wide format invertebrate genus dataset
invertebrates_condensed_wide <- invertebrate_meta_dist %>%
  select(site:Valvatidae) %>%
  gather(key = taxon, value = count, Acroloxidae:Valvatidae) %>%
  filter(!grepl("juvenile", taxon, ignore.case = TRUE)) %>%
  separate(col = taxon, into = c("Genus", "Species", "Subspecies")) %>%
  group_by(site, Genus) %>%
  summarize(total_Genus = sum(count)) %>%
  ungroup() %>%
  group_by(site) %>%
  mutate(total_site = sum(total_Genus)) %>%
  unique() %>%
  spread(key = Genus, value = total_Genus)

# Check how data are formatted
head(invertebrates_condensed_wide)

# Check for correlations
pairs(invertebrates_condensed_wide[, 3:17], upper.panel = panel.cor)


# Remove poorly preserved genera from dataset
invertebrates_well_preserved_long <- invertebrate_meta_dist %>%
  select(site:Valvatidae) %>%
  gather(key = taxon, value = count, Acroloxidae:Valvatidae) %>%
  filter(!grepl("juvenile", taxon, ignore.case = TRUE)) %>%
  separate(col = taxon, into = c("Genus", "Species", "Subspecies")) %>%
  group_by(site, Genus) %>%
  summarize(total_Genus = sum(count)) %>%
  ungroup() %>%
  group_by(site) %>%
  mutate(total_site = sum(total_Genus)) %>%
  unique() %>%
  ungroup()

# Rework Genus and Site columns as factors
invertebrates_well_preserved_long_ordered <- invertebrates_well_preserved_long %>%
  mutate(Genus = factor(x = Genus,
                        levels = c("Brandtia", "Cryptoropus", "Eulimnogammarus",
                                   "Pallasea", "Poekilogammarus",
                                   "Acroloxidae", "Baicaliidae", "Benedictidae",
                                   "Maackia", "Planorbidae", "Valvatidae",
                                   "Asellidae", "Caddisflies", "Flatworms",
                                   "Leeches")),
         site = factor(x = site,
                       levels = c("BGO-2",  "BGO-3", "KD-1",
                                  "KD-2", "BGO-1", "MS-1",
                                  "BK-3", "BK-2", "SM-1",
                                  "BK-1", "EM-1", "LI-1",
                                  "LI-2", "LI-3")))

# Prep invert data for plotting and then do so
invertebrate_well_preserved_plot <- invertebrates_well_preserved_long_ordered %>%
  mutate(Group = ifelse(test = Genus %in% amphipods,
                        yes = "Amphipod", no = Genus),
         Group = ifelse(test = Genus %in% molluscs,
                        yes = "Mollusc", no = Group),
         Genus = factor(x = Genus,
                        levels = c("Brandtia", "Cryptoropus", "Eulimnogammarus",
                                   "Pallasea", "Poekilogammarus",
                                   "Acroloxidae", "Baicaliidae", "Benedictidae",
                                   "Maackia", "Planorbidae", "Valvatidae",
                                   "Asellidae", "Caddisflies", "Flatworms",
                                   "Leeches"))) %>%
  ggplot() +
  geom_bar(aes(x = site, y = total_site), alpha = 0.5, stat = "identity") +
  geom_bar(aes(x = site, y = total_Genus, fill = as.factor(Group)),
           stat = "identity") +
  scale_fill_viridis_d(option = "inferno") +
  facet_wrap(~ Genus) +
  xlab("Site (Arranged by increasing IDW population)") +
  ylab("Number of Individuals") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        strip.background = element_rect(fill = "white"),
        panel.background = element_rect(color = "black"),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(margin = margin(0, 20, 0, 0)),
        axis.title.x = element_text(margin = margin(20, 0, 0, 0)),
        legend.text = element_text(size = 16))

# Export plot
ggsave(filename = "../figures/invertebrate_well_preserved_plot.png",
       plot = invertebrate_well_preserved_plot, device = "png",
       width = 18, height = 12, dpi = 300)


# 4.2 Multivariate analysis -----------------------------------------------

# Create a cleaned up wide format invertebrate genus dataset with well preserved
# genera. In this step we also remove rare species (i.e., taxa representing less
# than 1% of the intercommunity abundance).
invertebrates_well_preserved_wide <- invertebrate_meta_dist %>%
  select(site:Valvatidae) %>%
  gather(key = taxon, value = count, Acroloxidae:Valvatidae) %>%
  filter(!grepl("juvenile", taxon)) %>%
  separate(col = taxon, into = c("Genus", "Species", "Subspecies")) %>%
  group_by(site, Genus) %>%
  summarize(total_Genus = sum(count)) %>%
  ungroup() %>%
  group_by(site) %>%
  mutate(total_site = sum(total_Genus),
         prop_Genus = total_Genus / total_site) %>%
  ungroup() %>%
  group_by(Genus) %>%
  mutate(mean_prop = mean(prop_Genus, na.rm = TRUE)) %>%
  filter(mean_prop >= 0.01) %>%
  select(-mean_prop, -prop_Genus) %>%
  spread(key = Genus, value = total_Genus)

# Square-root transform
inver_comm <- sqrt(invertebrates_well_preserved_wide[, 3:14])

# Vizualize optimum cluster number for invert community
# Sometimes this can produce an error that reads:
# Error in .Call.graphics(C_palette2, .Call(C_palette2, NULL)) : 
# invalid graphics state.
# This can be fixed by running the command: dev.off()
invert_cluster_wss <- fviz_nbclust(x = as.matrix(vegdist(x = inver_comm, method = "bray",
                                                         diag = TRUE, upper = TRUE)),
                                   FUNcluster = cluster::pam, method = "wss")
invert_cluster_wss

invert_cluster_sil <- fviz_nbclust(x = as.matrix(vegdist(x = inver_comm, method = "bray",
                                                         diag = TRUE, upper = TRUE)),
                                   FUNcluster = cluster::pam, method = "silhouette", print.summary = FALSE)
invert_cluster_sil

invert_hclust <- pvclust(as.matrix(vegdist(x = inver_comm, method = "bray",
                                   diag = TRUE, upper = TRUE)), method.hclust = "median", 
                         method.dist = "correlation")

plot(invert_hclust)
pvrect(invert_hclust, alpha = 0.90)

# Export plots
ggsave(filename = "../figures/invertebrate_pam_analysis_wss.png",
       plot = invert_cluster_wss, device = "png",
       width = 18, height = 12, dpi = 300)

ggsave(filename = "../figures/invertebrate_pam_analysis_sil.png",
       plot = invert_cluster_sil, device = "png",
       width = 18, height = 12, dpi = 300)

ggsave(filename = "../figures/invertebrate_hclust_analysis.png",
       plot = plot(invert_hclust), device = "png",
       width = 18, height = 12, dpi = 300)


# Run invert NMDS
invertebrates_metaMDS <- metaMDS(comm = inver_comm, distance = "bray", try = 100)
invertebrates_metaMDS

# Pull scores from NMDS and add site data
data_scores <- as.data.frame(scores(x = invertebrates_metaMDS, display = "sites"))
data_scores$site <- invertebrates_well_preserved_wide$site

# Join scores with PPCP data and code into IDW_pop_group
data_scores <- full_join(x = data_scores, y = ppcp_meta_dist, by = "site") %>%
  mutate(IDW_pop_group = ifelse(test = site %in% c(low, mod),
                                yes = "Low/Mod", no = NA),
         IDW_pop_group = ifelse(test = site %in% high,
                                yes = "High", no = IDW_pop_group),
         IDW_pop_group = factor(x = IDW_pop_group,
                                levels = c("High", "Low/Mod")))

species_scores <- as.data.frame(scores(x = invertebrates_metaMDS, display = "species"))
species_scores$species <- rownames(species_scores)

# Plot NMDS
inverts_well_preserved_nmds <- ggplot() +
  geom_point(data = data_scores,
             aes(x = NMDS1, y = NMDS2, size = log10(distance_weighted_population),
                 color = IDW_pop_group)) +
  scale_size_continuous(range = c(12, 28), guide = FALSE) +
  scale_color_manual(values = inferno(15)[c(3, 8, 11, 14)],
                     name = "IDW Population Grouping") +
  guides(colour = guide_legend(override.aes = list(size = 10))) +
  geom_text_repel(data = species_scores %>% 
                    filter(species %in% c("Eulimnogammarus", "Poekilogammarus", "Valvatidae",
                                          "Caddisflies", "Brandtia", "Planorbidae", "Baicaliidae",
                                          "Cryptoropus", "Flatworms")) %>%
                    mutate(NMDS1 = ifelse(test = species == "Poekilogammarus", 
                                                 yes = NMDS1+0.015, no = NMDS1),
                           NMDS2 = ifelse(test = species == "Poekilogammarus", 
                                          yes = NMDS2+0.075, no = NMDS2),
                           NMDS1 = ifelse(test = species == "Eulimnogammarus", 
                                                 yes = NMDS1-0.08, no = NMDS1),
                           NMDS1 = ifelse(test = species == "Eulimnogammarus", 
                                          yes = NMDS1-0.06, no = NMDS1),
                           NMDS2 = ifelse(test = species == "Flatworms", 
                                          yes = NMDS2-0.15, no = NMDS2),
                           species = ifelse(test = species %in% amphipods, 
                                            yes = paste0("italic(", species, ")"),
                                            no = species)),
                  aes(x = NMDS1, y = NMDS2, label = species), 
                  size = 8, parse = TRUE) + 
  coord_equal() +
  annotate("label", x = -0.15, y = 0.4, size = 10,
           label = paste("Stress: ",
                         round(invertebrates_metaMDS$stress, digits = 3))) +
  theme(legend.position = "right",
        legend.key=element_blank(),
        strip.text.x = element_text(size = 20, color = "grey80"),
        text = element_text(size = 24),
        axis.title.y = element_text(margin = margin(0, 20, 0, 0)),
        axis.title.x = element_text(margin = margin(20, 0, 0, 0)),
        panel.background = element_rect("white"),
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor = element_line(colour = "grey80"),
        axis.ticks = element_line(color = "grey80"))

inverts_well_preserved_nmds

# Export plot. This plot is Figure 4 in the associated manuscript.
ggsave(filename = "../figures/inverts_well_preserved_nmds.png",
       plot = inverts_well_preserved_nmds, device = "png",
       height = 10, width = 20, dpi = 300)

# Re-join the metadata and distance data back in with invert data
inverts_well_preserved_meta_dist_wide <- full_join(x = invertebrates_well_preserved_wide,
                                                   y = ppcp_meta_dist,
                                                   by = "site") %>%
  mutate(IDW_pop_group = ifelse(test = site %in% c(low, mod),
                                yes = "Low/Mod", no = NA),
         IDW_pop_group = ifelse(test = site %in% high,
                                yes = "High", no = IDW_pop_group),
         distance_weighted_population = log10(distance_weighted_population)) %>%
  data.frame()

# Run PERMANOVA
adonis2(formula = sqrt(inverts_well_preserved_meta_dist_wide[, 3:14]) ~
          inverts_well_preserved_meta_dist_wide[, 28],
        data = inverts_well_preserved_meta_dist_wide,
        method = "bray")

adonis2(formula = sqrt(inverts_well_preserved_meta_dist_wide[, 3:14]) ~
         inverts_well_preserved_meta_dist_wide[, 29],
       data = inverts_well_preserved_meta_dist_wide,
       method = "bray")

# Run SIMPER
simper_results <- simper(comm = sqrt(inverts_well_preserved_meta_dist_wide[, 3:14]),
                         group = inverts_well_preserved_meta_dist_wide[, 29],
                         permutations = 999)

summary(simper_results)

# Combine pam wss and silhouette plots into one.
# This plot is associated with Figure S1 in the associated manuscript.

ggarrange(peri_cluster_wss, invert_cluster_wss,
          ncol = 2, nrow = 1, labels = "AUTO") %>%
  ggexport(filename = "../figures/keams_wss_combined_plot.png",
           hieght = 1800, width = 1200, res = 120)

ggarrange(peri_cluster_sil, invert_cluster_sil, 
          ncol = 2, nrow = 1, labels = "AUTO") %>%
  ggexport(filename = "../figures/keams_sil_combined_plot.png",
           hieght = 1800, width = 1200, res = 120)

# pvclust plots are not easily converted to a grob object, 
# so this is method for outputting the combined plot is a workaround 
png(filename="../figures/hclust_combined_plot.png", 
    width = 10, height = 8, units = "in", res = 120)
par(mfrow = c(1,2))
plot(peri_hclust, main = "Periphyton Cluster Dendrogram")
pvrect(peri_hclust, alpha = 0.9)
add_label(0.02, 0.07, "A")
plot(invert_hclust, main = "Macroinvertebrate Cluster Dendrogram")
pvrect(invert_hclust, alpha = 0.9, )
add_label(0.02, 0.07, "B")
dev.off()

ggarrange(periphyton_IDW_pop_group_plot, inverts_well_preserved_nmds,
          ncol = 1, nrow = 2, labels = "AUTO", 
          font.label = list(size = 36, face = "bold")) %>%
  ggexport(filename = "../figures/nmds_combined_plot.png", height = 3000, width = 3000, res = 120)
