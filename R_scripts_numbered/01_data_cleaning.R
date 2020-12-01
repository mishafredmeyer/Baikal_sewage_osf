## This script aggregates data from the cleaned but disaggregated forms from
## script 00_disaggregated_data_cleaning.R. This script is meant to prepare 
## the data for analysis in scripts 02-06. 

# 1. Load packages --------------------------------------------------------

library(tidyverse)
library(lubridate)
library(stringr)

sub_dir <- "cleaned_data"
output_dir <- file.path(here::here(), sub_dir)

if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Dir 'cleaned_data' already exists!")
}

# 2. Load and clean PPCP data ---------------------------------------------

# Raw data
ppcp_orig <- read.csv(file = "../clean_disaggregated_data/ppcp.csv",
                      header = TRUE)

# Select variables and sites of interest, sum PPCPs
ppcp <- ppcp_orig %>%
  select(site, caffeine, acetaminophen, paraxanthine, cotinine) %>%
  group_by(site) %>%
  mutate(ppcp_sum = caffeine + acetaminophen + paraxanthine + cotinine)
  
# Take a look
head(ppcp)

# Export new version of the data
write.csv(x = ppcp, file = "../cleaned_data/ppcp.csv",
          row.names = FALSE)


# 2. Load and clean nutrient data -----------------------------------------

nutrients_orig <- read.csv(file = "../clean_disaggregated_data/nutrients.csv",
                           header = TRUE)

# Take nutrients averages by site
nutrients <- nutrients_orig %>%
  group_by(site) %>%
  summarize(across(-replicate, mean, .names = "mean_{col}"))

head(nutrients)

write.csv(x = nutrients, file = "../cleaned_data/nutrients.csv",
          row.names = FALSE)


# 3. Load and clean chlorophyll a data ------------------------------------

chla_orig <- read.csv(file = "../clean_disaggregated_data/chlorophylla.csv",
                      header = TRUE)

# Select sites of interest and take average chl a by site
chlorophylla <- chla_orig %>%
  select(site, chl_conc) %>%
  group_by(site) %>%
  summarize(mean_chlorophylla = mean(chl_conc))

head(chlorophylla)

write.csv(x = chlorophylla, file = "../cleaned_data/chlorophylla.csv",
          row.names = FALSE)


# 4. Load and clean lat/long metadata -------------------------------------

metadata_orig <- read.csv(file = "../clean_disaggregated_data/metadata.csv",
                          header = TRUE)

# Select columns of interest
metadata <- metadata_orig %>%
  select(site, lat, long, depth_m, distance_to_shore_m, air_temp_celsius, 
         surface_temp_celsius, mid_temp_celsius, bottom_temp_celsius) 

head(metadata)

write.csv(x = metadata, file = "../cleaned_data/metadata.csv",
          row.names = FALSE)


# 5. Load and clean macroinvertebrate data --------------------------------

inverts_orig <- read.csv(file = "../clean_disaggregated_data/invertebrates.csv",
                         header = TRUE)

# Take mean counts by taxon
inverts_wide <- inverts_orig %>%
  group_by(site) %>%
  summarize(across(-replicate, mean))

head(inverts_wide)

write.csv(x = inverts_wide, file = "../cleaned_data/invertebrates.csv",
          row.names = FALSE)


# 6. Load and clean periphyton data ---------------------------------------

periphyton_orig <- read.csv(file = "../clean_disaggregated_data/periphyton.csv",
                            header = TRUE)

# Take mean counts by taxon
periphyton_summarized <- periphyton_orig %>%
  select(-contains("filament")) %>%
  filter(!is.na(diatom)) %>%
  group_by(site) %>%
  summarize(across(-one_of(c("replicate", "comments", "subsamples_counted")), mean)) %>%
  select(site, desmidales, diatom, pediastrum, spirogyra, tetrasporales, ulothrix)

head(periphyton_summarized)

write.csv(x = periphyton_summarized, file = "../cleaned_data/periphyton.csv",
          row.names = FALSE)

# 7. Load and clean stable isotope data -----------------------------------

stable_isotopes_orig <- read.csv(file = "../clean_disaggregated_data/stable_isotopes.csv",
                                 header = TRUE)

# Remove comments column
stable_isotopes <- stable_isotopes_orig %>%
  select(-comments)

head(stable_isotopes)

write.csv(x = stable_isotopes, file = "../cleaned_data/stable_isotopes.csv",
          row.names = FALSE)


# 8. Load and clean fatty acid data ---------------------------------------

fatty_acid_orig <- read.csv(file = "../clean_disaggregated_data/fatty_acid.csv",
                            header = TRUE)

# Remove comments column because our analysis focuses on proportions 
fatty_acids <- fatty_acid_orig %>%
  select(-comments)

head(fatty_acids)

write.csv(x = fatty_acids, file = "../cleaned_data/fatty_acid.csv",
          row.names = FALSE)


# 9. Load and clean microplastics data -----------------------------------

microplastics_orig <- read.csv(file = "../clean_disaggregated_data/microplastics.csv",
                               header = TRUE)

# Isolate samples that are not controls
microplastics_uncorrected <- microplastics_orig %>%
  select(-comments) %>%
  filter(replicate != "C")

# Separate out the controls so that they can be removed from
# the experimental counts
microplastics_controls <- microplastics_orig %>%
  select(-comments) %>%
  filter(replicate == "C") %>%
  group_by(site) %>%
  summarize(fiber_controls = mean(fibers),
            fragment_controls = mean(fragments),
            beads_controls = mean(beads))

microplastics_corrected <- left_join(x = microplastics_uncorrected, microplastics_controls,
                                     by = "site") %>%
  mutate(fragments_corrected = fragments - fragment_controls,
         fragments_corrected = ifelse(test = fragments_corrected < 0,
                                      yes = 0, no = fragments_corrected),
         fibers_corrected = fibers - fiber_controls,
         fibers_corrected = ifelse(test = fibers_corrected < 0,
                                   yes = 0, no = fibers_corrected),
         beads_corrected = beads - beads_controls,
         beads_corrected = ifelse(test = beads_corrected < 0,
                                  yes = 0, no = beads_corrected),
         total_microplastics = fragments_corrected + fibers_corrected + beads_corrected,
         density = total_microplastics / volume_filtered_ml,
         fragment_density = fragments_corrected / volume_filtered_ml,
         fiber_density = fibers_corrected / volume_filtered_ml,
         bead_density = beads_corrected / volume_filtered_ml) %>%
  select(site, replicate, total_microplastics, density,
         fragment_density, fiber_density, bead_density)

head(microplastics_corrected)

write.csv(x = microplastics_corrected, file = "../cleaned_data/microplastics.csv",
          row.names = FALSE)


# 11. Load and clean calculated distance metrics --------------------------

distance_orig <- read.csv(file = "../clean_disaggregated_data/distance_weighted_population_metrics.csv",
                          header = TRUE)

# All looks good from previous aggregation.

head(distance_orig)

write.csv(x = distance_orig, file = "../cleaned_data/distance_weighted_population_metrics.csv",
          row.names = FALSE)
