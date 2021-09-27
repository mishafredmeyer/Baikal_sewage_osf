# This script aggregates data to produce tables used within the
# associated manuscript.

library(tidyverse)
library(stringr)

# Check if tables directory exists 
# If not, create the tables directory
sub_dir <- "tables"
output_dir <- file.path(here::here(), sub_dir)

if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Dir 'tables' already exists!")
}


# 1. Load the data --------------------------------------------------------

site_information <- read.csv(file = "../cleaned_data/site_information.csv",
                     header = TRUE, stringsAsFactors = FALSE)

ppcp <- read.csv(file = "../cleaned_data/ppcp.csv",
                 header = TRUE, stringsAsFactors = FALSE)

nutrients <- read.csv(file = "../clean_disaggregated_data/nutrients.csv",
                      header = TRUE, stringsAsFactors = FALSE)

microplastics <- read.csv(file = "../clean_disaggregated_data/microplastics.csv",
                          header = TRUE, stringsAsFactors = FALSE)

distance <- read.csv(file = "../cleaned_data/distance_weighted_population_metrics.csv",
                     header = TRUE, stringsAsFactors = FALSE)


# 2. Site information table --------------------------------------------------

site_information_formatted <- site_information %>%
  select(-mid_temp_celsius, -bottom_temp_celsius) %>%
  rename(Latitude = lat,
         Longitude = long,
         Depth_m = depth_m,
         Distance_to_shore_m = distance_to_shore_m,
         Air_Temperature_C = air_temp_celsius,
         Surface_Temperature_C = surface_temp_celsius,
         Site = site,
         Time_of_Sampling_hh_mm = time) %>%
  mutate(Adjacent_Population = ifelse(test = Site %in% c("BGO-2", "BGO-3"),
                                      yes = 600, no = 0),
         Adjacent_Population = ifelse(test = Site %in% c("BK-1", "BK-2", "BK-3"),
                                      yes = 80, no = Adjacent_Population),
         Adjacent_Population = ifelse(test = Site %in% c("LI-1", "LI-2", "LI-3"),
                                      yes = 2000, no = Adjacent_Population),
         Adjacent_Population = ifelse(test = Site %in% c("OS-1", "OS-2", "OS-3"),
                                      yes = NA, no = Adjacent_Population),
         Sampling_Date = paste(year, sprintf(fmt = "%02d", month), day, sep = "")) %>%
  select(Site:Surface_Temperature_C, Adjacent_Population, Sampling_Date, Time_of_Sampling_hh_mm)

# This table is meant to be associated with Table 1 in the main
# body of the manuscript.

write.csv(x = site_information_formatted, file = "../tables/site_information_table1.csv",
          row.names = FALSE)


# 3. PPCP table -----------------------------------------------------------

ppcp_formatted <- ppcp %>%
  rename(Total_PPCP_Concentration = ppcp_sum)

write.csv(x = ppcp_formatted, file = "../tables/ppcp.csv",
          row.names = FALSE)


# 4. Nutrients table ------------------------------------------------------

nutrients_formatted <- nutrients %>%
  select(-tpo43_mg_dm3) %>%
  group_by(site) %>%
  summarize(NH4_mg_dm3_mean = mean(nh4_mg_dm3),
            NH4_mg_dm3_sd = sd(nh4_mg_dm3),
            NO3_mg_dm3_mean = mean(no3_mg_dm3),
            NO3_mg_dm3_sd = sd(no3_mg_dm3),
            TP_mg_dm3_mean = mean(tp_mg_dm3),
            TP_mg_dm3_sd = sd(tp_mg_dm3))

write.csv(x = nutrients_formatted, file = "../tables/nutrients.csv",
          row.names = FALSE)


# 5. Microplastics table --------------------------------------------------

# Isolate samples that are not controls
microplastics_uncorrected <- microplastics %>%
  select(-comments) %>%
  filter(replicate != "C")

# Separate out the controls so that they can be removed from
# the experimental counts
microplastics_controls <- microplastics %>%
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

microplastics_formatted <- microplastics_corrected %>%
  group_by(site) %>%
  summarize(Microplastic_mean_density_microplastics_per_L = mean(density),
            Microplastic_sd_density_microplastics_per_L = sd(density),
            Fragment_mean_density_microplastics_per_L = mean(fragment_density),
            Fragment_sd_density_microplastics_per_L = sd(fragment_density),
            Fiber_mean_density_microplastics_per_L = mean(fiber_density),
            Fiber_sd_density_microplastics_per_L = sd(fiber_density),
            Bead_mean_density_microplastics_per_L = mean(bead_density),
            Bead_sd_density_microplastics_per_L = sd(bead_density)) 

write.csv(x = microplastics_formatted, file = "../tables/microplastics.csv",
          row.names = FALSE)


# 6. Combined table -------------------------------------------------------

low <- c("BGO-1", "BGO-2", "BGO-3", "KD-1", "KD-2", "MS-1", "OS-1")
mod <- c("BK-2", "BK-3", "SM-1", "OS-3")
high <- c("BK-1", "EM-1", "LI-3", "LI-2", "LI-1", "OS-2")

# This table is associated with Table 2 in the manuscript

meta_nut_ppcp_mp <- full_join(x = site_information_formatted,
                              y = nutrients_formatted,
                              by = c("Site" = "site")) %>%
  full_join(x = ., y = ppcp_formatted, by = c("Site" = "site")) %>%
  full_join(x = ., y = microplastics_formatted, by = c("Site" = "site")) %>%
  full_join(x = ., y = distance, by = c("Site" = "site")) %>%
  mutate(Categorical_distance_weighted_population = ifelse(test = Site %in% low,
                                                           yes = "Low", no = NA),
         Categorical_distance_weighted_population = ifelse(test = Site %in% mod,
                                                           yes = "Mod",
                                                           no = Categorical_distance_weighted_population),
         Categorical_distance_weighted_population = ifelse(test = Site %in% high,
                                                           yes = "High",
                                                           no = Categorical_distance_weighted_population)) %>%
  select(Site, NH4_mg_dm3_mean:cotinine, 
         Microplastic_mean_density_microplastics_per_L:Bead_sd_density_microplastics_per_L,
         distance_weighted_population,
         Categorical_distance_weighted_population) 

write_csv(x = meta_nut_ppcp_mp, path = "../tables/combined_table2.csv", col_names = TRUE, )
