# This script performs two food web inferences studies using 
# stable isotopes and fatty acids. First the script aggregates 
# stable isotopes with sewage indicator data in Lake Baikal to 
# produce a stable isotope biplot for the associated manuscript. 
# Second, the script performs a Bayesian mixing model using fatty acid 
# data collected at Lake Baikal, already published fatty acid data for 
# Lake Baikal primary producers, and antarctic amphipod trophic discrimination
# factors (the most closely related taxa with publicly available data at the time
# of this analysis). 

library(tidyverse)
library(viridis)
library(ggtext)
library(car)
library(MixSIAR)
library(ggpubr)

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
            sd_C13 = sd(C13),
            n = length(taxon_idw_pop)) %>%
  mutate(taxon_idw_pop = paste(taxon_idw_pop, " (N=", n, ")", sep = ""))


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
  ggtitle("Stable Isotope Isospace") +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  theme_minimal() +
  theme(legend.position = "bottom",
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_markdown(size = 16),
        legend.key.height = unit(0.5, "in"),
        legend.key.width = unit(0.25, "in")) +
  guides(color = guide_legend(nrow = 3))


ggsave("../figures/stable_isotopes_biplot.png", foodweb_plot,
       device = "png", width = 9, height = 6, units = "in")


# Step 5: Perform signifance tests ----------------------------------------

foodweb_analysis_data <- stable_isotopes %>%
  mutate(idw_group = ifelse(site %in% c(low, mod), "Low/Mod", NA),
         idw_group = ifelse(site %in% high, "High", idw_group),
         Genus = as.character(Genus),
         Genus = ifelse(test = Genus != "Periphyton",
                        yes = paste0("*", as.character(Genus), "*"),
                        no = Genus)) %>%
  unite(taxon_idw_pop, c("Genus", "idw_group"), sep = " + ") %>%
  filter(!grepl("NA", taxon_idw_pop)) %>%
  select(taxon_idw_pop, C13, N15)

Anova(lm(formula = N15 ~ taxon_idw_pop, 
         data = foodweb_analysis_data %>%
           filter(grepl(pattern = "Periphyton", x = taxon_idw_pop, ignore.case = TRUE))), 
      type = "II")

Anova(lm(formula = N15 ~ taxon_idw_pop, 
         data = foodweb_analysis_data %>%
           filter(grepl(pattern = "Eulimnogammarus", x = taxon_idw_pop, ignore.case = TRUE))), 
      type = "II")

# 6. Building a Bayesian Mixing Model for trophic interactions ------------


# 6.1 Aggregate data into required formats --------------------------------

# First create the producer file
# These data come from two main sources:
# Ulothrix Fatty Acid data come from Osipova et al 2009
# Diatom Fatty Acid data from Shishlyannikov et al 2019. These data were for 
# S. acus and A. baicalensis, so we average between these two taxa to create 
# a characteristic Baikalian diatom signature. 
# For the producers and consumers, we also remove certain fatty acids, so that 
# the final files only contain fatty acids present in both producers and consumers.

producer_file <- read.csv("../original_data/fatty_acid_producer_literature.csv") %>%
  filter(grepl(pattern = "prop", x = taxon)) %>%
  mutate(taxon = as.character(taxon),
         taxon = ifelse(grepl(pattern = "ulothrix", x = taxon), yes = "Ulothrix", no = taxon),
         taxon = ifelse(grepl(pattern = "acus", x = taxon), yes = "Diatom", no = taxon),
         taxon = ifelse(grepl(pattern = "baicalensis", x = taxon), yes = "Diatom", no = taxon)) %>%
  pivot_longer(cols = c(c14_0:std22_6w3), names_to = "fatty_acid", values_to = "prop") %>%
  filter(!grepl(pattern = "22_5w3*|*16_2w4", x = fatty_acid)) %>%
  group_by(taxon, fatty_acid) %>%
  summarize(prop_mean = mean(prop, na.rm = TRUE)) %>%
  pivot_wider(names_from = fatty_acid, values_from = prop_mean) 

unique_producer_fa <- colnames(producer_file)[grepl(pattern = "c", x = colnames(producer_file))]

drapa_data <- read.csv("../cleaned_data/fatty_acid.csv") %>%
  unite(col = "Genus_Species", c("Genus", "Species"), sep = "_") %>%
  filter(Genus_Species %in% c("Draparnaldia_NA")) %>%
  mutate(Genus_Species = ifelse(Genus_Species == "Draparnaldia_NA", "Draparnaldia", Genus_Species)) %>%
  pivot_longer(cols = c(c14_0:c22_6w3), names_to = "fatty_acid", values_to = "concentration") %>%
  filter(fatty_acid %in% unique_producer_fa) %>%
  group_by(site, Genus_Species) %>%
  mutate(sum_concentration = sum(concentration),
         proportion = concentration/sum_concentration) %>%
  select(-concentration, -sum_concentration) %>%
  ungroup() %>%
  group_by(Genus_Species, fatty_acid) %>%
  summarize(mean_proportion = mean(proportion)*100,
            sd_proportion = sd(proportion)*100) %>%
  ungroup() %>%
  pivot_wider(names_from = c(fatty_acid) , 
              values_from = c(mean_proportion, sd_proportion)) %>%
  ungroup() %>%
  rename("taxon" = "Genus_Species")

names(drapa_data) <- names(drapa_data) %>%
  gsub(pattern = "mean_proportion_", replacement = "", .) %>%
  gsub(pattern = "sd_proportion_c", replacement = "std", .)

producers_whole <- rbind(data.frame(producer_file), drapa_data)

producers_mixsiar <- producers_whole %>%
  ## adding this step for MixSIAR
  pivot_longer(cols = c14_0:std22_6w3, names_to = "fatty_acid", values_to = "prop") %>%
  mutate(fatty_acid = gsub(pattern = "std", x = fatty_acid, replacement = "SDc"),
         fatty_acid = ifelse(!grepl(pattern = "SDc", x = fatty_acid), paste0("Mean", fatty_acid), fatty_acid),
         prop = ifelse(grepl(pattern = "SDc", x = fatty_acid) & prop == 0, 0.000001, prop)) %>%
  pivot_wider(names_from = fatty_acid, values_from = prop) %>%
  mutate(n = ifelse(taxon == "Diatom", 2, NA),
         n = ifelse(taxon == "Ulothrix", 2, n),
         n = ifelse(taxon == "Draparnaldia", 7, n))

write.csv(x = producers_mixsiar, 
          file = "../cleaned_data/fatty_acid_producer_mixsiar.csv", 
          row.names = FALSE)

unique_producer_fa <- colnames(producer_file)[grepl(pattern = "c", x = colnames(producer_file))]

# Second create the consumer file
# These data come from our own samples, which we collected as part of this 
# samping campaign. Because our main analyses were restricted to 
# E. vittatus and E. verrucosus, we restrict this analysis to those 
# species. 

consumer_file <- read.csv("../cleaned_data/fatty_acid.csv") %>%
  select(site, Genus, Species, paste(unique_producer_fa)) %>%
  unite(col = "Genus_Species", c("Genus", "Species"), sep = "_") %>%
  filter(Genus_Species %in% c("Eulimnogammarus_verrucosus	", "Eulimnogammarus_vittatus")) %>%
  pivot_longer(cols = c(c14_0:c22_6w3), names_to = "fatty_acid", values_to = "concentration") %>%
  group_by(site, Genus_Species) %>%
  mutate(sum_concentration = sum(concentration),
         proportion = (concentration/sum_concentration)*100) %>%
  select(-concentration, -sum_concentration) %>%
  pivot_wider(names_from = fatty_acid, values_from = proportion)

write.csv(x = consumer_file, 
          file = "../cleaned_data/fatty_acid_consumer_formatted.csv",
          row.names = FALSE)

# Third create the trophic discrimination factor (TDF) file
# These data come from Schram et al (2019), which used marine 
# antartic amphipods in grazing experiments. Because Schram et al (2019)
# considered Diatoms, Red, and Brown algae as potential resources, we average
# among all three, so that the final TDFs for our analysis are identical for 
# all three potential resources defined in the producer file.

unique_discrim_fa <- gsub(pattern = "c", replacement = "C", x = unique_producer_fa, ignore.case = FALSE)
unique_discrim_fa <- gsub(pattern = "_", replacement = ".", x = unique_discrim_fa, ignore.case = FALSE)
unique_discrim_fa <- unique_discrim_fa[!(unique_discrim_fa %in% c("C16.2w4", "C22.5w3"))]

amphipod_discrim <- read.csv("../original_data/B022_2018FA.csv") %>%
  select(diet_type, tissue, time, meas_type, paste(unique_discrim_fa)) %>%
  filter(tissue == "amphipod" & time == "final" & meas_type == "prop") %>%
  filter(diet_type != "Palmeria") %>%
  group_by(diet_type) %>%
  summarize(across(C14.0:C22.6w3, list(Mean = mean, SD = sd), .names = "{.fn}{.col}")) %>%
  pivot_longer(cols = c(MeanC14.0:SDC22.6w3), names_to = "fatty_acid", values_to = "prop") %>%
  mutate(fatty_acid = gsub(pattern = "C", replacement = "c", x = fatty_acid, ignore.case = FALSE),
         fatty_acid = gsub(pattern = "\\.", replacement = "_", x = fatty_acid, ignore.case = FALSE)) %>%
  pivot_wider(names_from = fatty_acid, values_from = prop)

alga_discrim <- read.csv("../original_data/B022_2018FA.csv") %>%
  select(diet_type, tissue, time, meas_type, paste(unique_discrim_fa)) %>%
  filter(tissue == "alga" & time == "final" & meas_type == "prop") %>%
  group_by(diet_type) %>%
  summarize(across(C14.0:C22.6w3, list(Mean = mean, SD = sd), .names = "{.fn}{.col}")) %>%
  pivot_longer(cols = c(MeanC14.0:SDC22.6w3), names_to = "fatty_acid", values_to = "prop") %>%
  mutate(fatty_acid = gsub(pattern = "C", replacement = "c", x = fatty_acid, ignore.case = FALSE),
         fatty_acid = gsub(pattern = "\\.", replacement = "_", x = fatty_acid, ignore.case = FALSE)) %>%
  pivot_wider(names_from = fatty_acid, values_from = prop)

tdf_raw <- amphipod_discrim[ , -1] - alga_discrim[ , -1]

tdf <- tdf_raw %>%
  cbind(alga_discrim [, 1], .) %>%
  pivot_longer(cols = c(Meanc14_0:SDc22_6w3), names_to = "fatty_acid", values_to = "prop") %>%
  separate(col = "fatty_acid", into = c("stat", "fatty_acid"), sep = "c", remove = TRUE) %>%
  group_by(stat, fatty_acid) %>%
  mutate(mean_fa = mean(prop)) %>%
  select(-prop) %>%
  unite(col = "fatty_acid", stat:fatty_acid, remove = TRUE, sep = "c") %>%
  pivot_wider(names_from = fatty_acid, values_from = mean_fa) %>%
  mutate(diet_type = as.character(diet_type),
         diet_type = ifelse(diet_type == "Desmarestia", "Ulothrix", diet_type),
         diet_type = ifelse(diet_type == "Himantothallus", "Draparnaldia", diet_type),
         diet_type = ifelse(diet_type == "diatom", "Diatom", diet_type))

write.csv(x = tdf, 
          file = "../cleaned_data/tdf_formatted.csv",
          row.names = FALSE)


# 6.2 Define specifications for Bayesian Mixing Model ---------------------

modelname = "MixSIR"                        # underlying model can be "mixSIR", "SIAR" 
mcmc.chainLength <- as.integer(100000)      # post-burn 
mcmc.burn <- as.integer(50000) 
mcmc.thin = 50
mcmc.chains = 3 


# 6.3 Load data for model building ----------------------------------------

# First load consumer data

mix <- load_mix_data(filename = "../cleaned_data/fatty_acid_consumer_formatted.csv",
                     iso_names = unique_producer_fa,
                     factors=NULL,
                     fac_random=NULL,
                     fac_nested=NULL,
                     cont_effects=NULL)

# Second load resource data

source <- load_source_data(filename = "../cleaned_data/fatty_acid_producer_mixsiar.csv",
                           source_factors=NULL,
                           conc_dep=FALSE,
                           data_type="means",
                           mix)

# Third load TDF data

discr <- load_discr_data(filename = "../cleaned_data/tdf_formatted.csv", mix)

# Our model will use an uniformed prior. Just be sure, we can plot our prior isospace

plot_prior(alpha.prior=1,source)


# 6.4 Building the model --------------------------------------------------

# Write the JAGS model file
model_filename <- "MixSIAR_model.txt"
resid_err <- FALSE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

jags.1 <- run_model(run="test", mix, source, discr, model_filename)

jags.1 <- run_model(run="normal", mix, source, discr, model_filename)

# Now plot the prior distributions 

p_outputs <- data.frame(jags.1$BUGSoutput$sims.list$p.global) %>%
  rename("Diatom" = X1,
         "Draparnaldia" = X2,
         "Ulothrix" = X3) %>%
  pivot_longer(cols = c(Diatom:Ulothrix), names_to = "resources", values_to = "posteriors") %>%
  mutate(resources = gsub(pattern = "Draparnaldia", replacement = "*Draparnaldia* spp.", x = resources),
         resources = gsub(pattern = "Ulothrix", replacement = "*Ulothrix* spp.", x = resources))


posterior_plot <- ggplot(p_outputs) +
  geom_density(aes(x = posteriors, y = ..scaled.., fill = resources, color = resources), alpha = 0.65) +
  scale_fill_manual(values = viridis(30)[c(15, 22, 28)], name = "Resources") +
  scale_color_manual(values = viridis(30)[c(15, 22, 28)], name = "Resources") +
  ggtitle("Posterior Density Plots of Diet Proportions ") +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Proportion of diet") +
  ylab("Scaled frequency of posterior") +
  theme_minimal() +
  theme(legend.position = "bottom",
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_markdown(size = 16),
        legend.key.height = unit(0.5, "in"),
        legend.key.width = unit(0.25, "in"))

posterior_boxplot <- ggplot(p_outputs) +
  geom_violin(aes(x = resources, y = posteriors), alpha = 0.65, fill = "grey70") +
  geom_boxplot(aes(x = resources, y = posteriors), 
               alpha = 0.65, width = 0.2, outlier.alpha = 0) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  ylab("Proportion of diet") +
  theme_minimal() +
  theme(legend.position = "bottom",
        title = element_text(size = 20),
        axis.text.x = element_markdown(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.title.x = element_blank(),
        legend.text = element_markdown(size = 16))

ggsave(filename = "mixsiar_posterior_results.png", plot = posterior_plot, device = "png", 
       path = "../figures/", width = 8, height = 6, units = "in")


arranged_plots <- ggarrange(foodweb_plot, posterior_plot,
                            font.label = list(size = 28, face = "bold"), ncol = 2, nrow = 1,
                            labels = "AUTO")

ggsave(filename = "sia_mixsiar_results.png", plot = arranged_plots, device = "png", 
       path = "../figures/", width = 24, height = 12, units = "in")


# 6.5: Sensitivty Analysis ------------------------------------------------

# Because our TDFs are based on Antarctic marine taxa, we performed a sensitivity
# analysis by increasing the SD for each fatty acid by 5%, 10%, and 25% to 
# qualitatively and quantitatively understand how results may change. 


# 6.5.1: 5% increase ------------------------------------------------------

# Test case for 50% increase
tdf_fifty_percent <- read.csv(file = "../cleaned_data/tdf_formatted.csv") %>%
  pivot_longer(cols = c(Meanc14_0:SDc22_6w3), names_to = "fatty_acid", values_to = "prop") %>%
  mutate(prop = ifelse(grepl(pattern = "SD", x = fatty_acid, ), prop*2, prop)) %>%
  pivot_wider(names_from = fatty_acid, values_from = prop)

write.csv(x = tdf_fifty_percent, 
          file = "../cleaned_data/tdf_fifty_percent.csv",
          row.names = FALSE)


# 6.2 Define specifications for Bayesian Mixing Model ---------------------

modelname = "MixSIR"                        # underlying model can be "mixSIR", "SIAR" 
mcmc.chainLength <- as.integer(100000)      # post-burn 
mcmc.burn <- as.integer(50000) 
mcmc.thin = 50
mcmc.chains = 3 


# Load data for model building ----------------------------------------

# First load consumer data

mix <- load_mix_data(filename = "../cleaned_data/fatty_acid_consumer_formatted.csv",
                     iso_names = unique_producer_fa,
                     factors=NULL,
                     fac_random=NULL,
                     fac_nested=NULL,
                     cont_effects=NULL)

# Second load resource data

source <- load_source_data(filename = "../cleaned_data/fatty_acid_producer_mixsiar.csv",
                           source_factors=NULL,
                           conc_dep=FALSE,
                           data_type="means",
                           mix)

# Third load TDF data

discr <- load_discr_data(filename = "../cleaned_data/tdf_fifty_percent.csv", mix)

# Our model will use an uniformed prior. Just be sure, we can plot our prior isospace

plot_prior(alpha.prior=1,source)


# 6.4 Building the model --------------------------------------------------

# Write the JAGS model file
model_filename <- "MixSIAR_model.txt"
resid_err <- FALSE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

jags.1 <- run_model(run="test", mix, source, discr, model_filename)

jags.1 <- run_model(run="normal", mix, source, discr, model_filename)

# Now plot the prior distributions 

p_outputs <- data.frame(jags.1$BUGSoutput$sims.list$p.global) %>%
  rename("Diatom" = X1,
         "Draparnaldia" = X2,
         "Ulothrix" = X3) %>%
  pivot_longer(cols = c(Diatom:Ulothrix), names_to = "resources", values_to = "posteriors") %>%
  mutate(resources = gsub(pattern = "Draparnaldia", replacement = "*Draparnaldia* spp.", x = resources),
         resources = gsub(pattern = "Ulothrix", replacement = "*Ulothrix* spp.", x = resources))


poster_plot_five_percent <- ggplot(p_outputs) +
  geom_density(aes(x = posteriors, y = ..scaled.., fill = resources, color = resources), alpha = 0.65) +
  scale_fill_manual(values = viridis(30)[c(15, 22, 28)], name = "Resources") +
  scale_color_manual(values = viridis(30)[c(15, 22, 28)], name = "Resources") +
  ggtitle("Posterior Density Plots of Diet Proportions: Doubling TDF standard deviation") +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Proportion of diet") +
  ylab("Scaled frequency of posterior") +
  theme_minimal() +
  theme(legend.position = "bottom",
        title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_markdown(size = 16),
        legend.key.height = unit(0.5, "in"),
        legend.key.width = unit(0.25, "in"))
