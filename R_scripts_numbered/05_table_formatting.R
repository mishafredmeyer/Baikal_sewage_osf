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

metadata <- read.csv(file = "../cleaned_data/metadata.csv",
                     header = TRUE, stringsAsFactors = FALSE)

ppcp <- read.csv(file = "../cleaned_data/ppcp.csv",
                 header = TRUE, stringsAsFactors = FALSE)

nutrients <- read.csv(file = "../cleaned_data/nutrients.csv",
                      header = TRUE, stringsAsFactors = FALSE)

microplastics <- read.csv(file = "../cleaned_data/microplastics.csv",
                          header = TRUE, stringsAsFactors = FALSE)

distance <- read.csv(file = "../cleaned_data/distance_weighted_population_metrics.csv",
                     header = TRUE, stringsAsFactors = FALSE)


# 2. Site metadata table --------------------------------------------------

metadata_formatted <- metadata %>%
  select(-mid_temp_celsius, -bottom_temp_celsius) %>%
  rename(Latitude = lat,
         Longitude = long,
         Depth_m = depth_m,
         Distance_to_shore_m = distance_to_shore_m,
         Air_Temperature_C = air_temp_celsius,
         Surface_Temperature_C = surface_temp_celsius,
         Site = site) %>%
  mutate(Adjacent_Population = ifelse(test = Site %in% c("BGO-2", "BGO-3"),
                                      yes = 600, no = 0),
         Adjacent_Population = ifelse(test = Site %in% c("BK-1", "BK-2", "BK-3"),
                                      yes = 80, no = Adjacent_Population),
         Adjacent_Population = ifelse(test = Site %in% c("LI-1", "LI-2", "LI-3"),
                                      yes = 2000, no = Adjacent_Population),
         Adjacent_Population = ifelse(test = Site %in% c("OS-1", "OS-2", "OS-3"),
                                      yes = NA, no = Adjacent_Population))

# This table is meant to be associated with Table 1 in the main
# body of the manuscript.

write.csv(x = metadata_formatted, file = "../tables/metadata_table1.csv",
          row.names = FALSE)


# 3. PPCP table -----------------------------------------------------------

ppcp_formatted <- ppcp %>%
  rename(Total_PPCP_Concentration = ppcp_sum)

write.csv(x = ppcp_formatted, file = "../tables/ppcp.csv",
          row.names = FALSE)


# 4. Nutrients table ------------------------------------------------------

nutrients_formatted <- nutrients %>%
  select(-mean_tpo43_mg_dm3) %>%
  rename(NH4_mg_dm3 = mean_nh4_mg_dm3,
         NO3_mg_dm3 = mean_no3_mg_dm3,
         TP_mg_dm3 = mean_tp_mg_dm3)

write.csv(x = nutrients_formatted, file = "../tables/nutrients.csv",
          row.names = FALSE)


# 5. Microplastics table --------------------------------------------------

microplastics_formatted <- microplastics %>%
  group_by(site) %>%
  summarize(Microplastic_density_microplastics_per_L = mean(density),
            Fragment_density_microplastics_per_L = mean(fragment_density),
            Fiber_density_microplastics_per_L = mean(fiber_density),
            Bead_density_microplastics_per_L = mean(bead_density)) 

write.csv(x = microplastics_formatted, file = "../tables/microplastics.csv",
          row.names = FALSE)


# 6. Combined table -------------------------------------------------------

low <- c("BGO-1", "BGO-2", "BGO-3", "KD-1", "KD-2", "MS-1", "OS-1")
mod <- c("BK-2", "BK-3", "SM-1", "OS-3")
high <- c("BK-1", "EM-1", "LI-3", "LI-2", "LI-1", "OS-2")

# This table is associated with Table 2 in the manuscript

meta_nut_ppcp_mp <- full_join(x = metadata_formatted,
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
  select(Site, NH4_mg_dm3:cotinine, 
         Microplastic_density_microplastics_per_L:Bead_density_microplastics_per_L,
         distance_weighted_population,
         Categorical_distance_weighted_population)

write.csv(x = meta_nut_ppcp_mp, file = "../tables/combined_table2.csv",
          row.names = FALSE)
