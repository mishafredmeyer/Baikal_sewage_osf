## This script reads in the raw data and formats the data such that misspelled 
## taxa are correctly spelled, column headers are uniform, and data can be 
## shared easily with collaborators. This script is meant to prime the data 
## for script 01_data_cleaning.R, where data are averaged and prepped for 
## the appropriate analysis in scripts 02-07. 


# 1. Load packages --------------------------------------------------------

library(tidyverse)
library(lubridate)
library(stringr)
library(janitor)
library(sf)
library(spdplyr)


# 2. Load and clean PPCP data ---------------------------------------------

# Station IDs
st_ids <- c("BK-1", "BGO-3", "BK-2", "BK-3", "BGO-1", "BGO-2", "KD-1",
            "KD-2", "MS-1", "SM-1", "OS-1", "OS-2", "EM-1", "LI-1", "LI-2",
            "LI-3", "OS-3")


# Raw data
ppcp_orig <- read.csv(file = "../original_data/PPCP_Baikal_orig_20180524.csv",
                      header = TRUE)

# Select variables and sites of interest, sum PPCPs
ppcp <- ppcp_orig %>%
  clean_names() %>%
  filter(sample_id %in% st_ids) %>%
  group_by(sample_id) %>%
  rename(paraxanthine = x1_7_dimethylxanthine, 
         site = sample_id) %>%
  mutate(collection_year = year(mdy(collection_date)),
         collection_month = month(mdy(collection_date)),
         collection_day = day(mdy(collection_date)),
         analysis_year = year(mdy(analysis_date)),
         analysis_month = month(mdy(analysis_date)),
         analysis_day = day(mdy(analysis_date))) %>%
  select(-collection_date, -analysis_date, -batch)

# Take a look
head(ppcp)

# Export new version of the data
write.csv(x = ppcp, file = "../clean_disaggregated_data/ppcp.csv",
          row.names = FALSE)


# 2. Load and clean nutrient data -----------------------------------------

nutrients_orig <- read.csv(file = "../original_data/baikal_nearshore_nutrient_data_20151009.csv",
                           header = TRUE)

# Take nutrients, clean column names, and standardize relational column
nutrients <- nutrients_orig %>%
  clean_names(case = "snake") %>%
  rename(site = sample)

head(nutrients)

write.csv(x = nutrients, file = "../clean_disaggregated_data/nutrients.csv",
          row.names = FALSE)


# 3. Load and clean chlorophyll a data ------------------------------------

chla_orig <- read.csv(file = "../original_data/chlorophyll_20170117.csv",
                      header = TRUE)

# Select sites of interest and take average chl a by site
chlorophylla <- chla_orig %>%
  select(-Chlorophyll) %>%
  filter(Station %in% st_ids) %>%
  clean_names(case = "snake") %>%
  # Correct unintended formatting as dates
  mutate(replicate = as.character(replicate),
         replicate = ifelse(test = grepl(pattern = "Jan", x = replicate), 
                            yes = 1, no = replicate),
         replicate = ifelse(test = grepl(pattern = "Feb", x = replicate), 
                            yes = 2, no = replicate),
         replicate = ifelse(test = grepl(pattern = "Mar", x = replicate), 
                            yes = 3, no = replicate)) %>%
  rename(site = station)

head(chlorophylla)

write.csv(x = chlorophylla, file = "../clean_disaggregated_data/chlorophylla.csv",
          row.names = FALSE)


# 4. Load and clean lat/long metadata -------------------------------------

metadata_orig <- read.csv(file = "../original_data/baikal_nearshore_metadata_201508.csv",
                          header = TRUE)

# Select columns of interest
metadata <- metadata_orig %>%
  rename(common_locaiton = location,
         local_site = site,
         site = loc_site,
         site_description = site_desc,
         depth_m = depth,
         distance_to_shore_m = dist_to_shore,
         air_temp_celsius = air_temp,
         surface_temp_celsius = surface_temp,
         mid_temp_celsius = mid_temp,
         bottom_temp_celsius = bottom_temp) %>%
  select(year, month, day, time, site, lat, long, site_description, distance_to_shore_m, depth_m,
         air_temp_celsius, surface_temp_celsius, mid_temp_celsius, bottom_temp_celsius, 
         comments, shore_photo, substrate_photo, sponges, brandtia)

head(metadata)

write.csv(x = metadata, file = "../clean_disaggregated_data/metadata.csv",
          row.names = FALSE)


# 5. Load and clean macroinvertebrate data --------------------------------

inverts_orig <- read.csv(file = "../original_data/macroinvert_community_QAQC_mfm_20171108.csv",
                         header = TRUE)

# Clean up naming in data and add up counts by taxon
# During enumeration, choronomids and and oligochaetes were poorly preserved, 
# so we resolved after a few samples were counted to remove them from the dataset. 
# Additionally, hyallela_cziarnianski was misidentified in some samples. Rather than
# potentially erroneously include those counts, we remove this taxon. 
# LI-3 was partitioned into two samples due to the high number of individuals. 
# This step combines that replicate into one row. 
inverts_wide <- inverts_orig %>%
  select(-X) %>%
  gather(key = site, value = count, MS1.3:BK1.3) %>%
  rename(taxon = Invertebrate) %>%
  mutate(site = gsub(pattern = ".", replacement = "_", x = site, fixed = TRUE),
         taxon = trimws(taxon),
         taxon = gsub(pattern = " ", replacement = "_", x = taxon),
         count = ifelse(test = is.na(count), yes = 0, no = count)) %>%
  separate(col = site, into = c("location", "replicate", "duplicate"),
           remove = FALSE) %>%
  filter(!(taxon %in% c("Propapaidae", "choronomids", "hyallela_cziarnianski"))) %>%
  mutate(replicate = ifelse(test = replicate %in% c("1B1", "1B2"),
                            yes = "1", no = replicate),
         taxon = ifelse(test = grepl(pattern = "Brandtia_latissima", x = taxon),
                        yes = "Brandtia_latissima", no = taxon),
         taxon = ifelse(test = grepl(pattern = "Benedictidate", x = taxon),
                        yes = "Benedictidae", no = taxon)) %>%
  group_by(location, replicate, taxon) %>%
  summarize(sum_count = sum(count)) %>%
  ungroup() %>%
  separate(col = taxon, into = c("Genus", "Species", "Subspecies")) %>%
  mutate(Genus = ifelse(test = Genus == "E",
                        yes = "Eulimnogammarus", no = Genus),
         Genus = ifelse(test = Genus == "Eulimno",
                        yes = "Eulimnogammarus", no = Genus),
         Genus = ifelse(test = Genus == "flatworms",
                        yes = "Flatworms", no = Genus),
         Genus = ifelse(test = Genus == "caddisflies",
                        yes = "Caddisflies", no = Genus),
         Genus = ifelse(test = Genus == "pallasea",
                        yes = "Pallasea", no = Genus),
         Genus = ifelse(test = Genus == "hyallela",
                        yes = "Hyallela", no = Genus),
         Genus = ifelse(test = Genus == "poekilo",
                        yes = "Poekilogammarus", no = Genus),
         Genus = ifelse(test = Genus == "Poekilo",
                        yes = "Poekilogammarus", no = Genus),
         Genus = ifelse(test = Genus == "valvatidae",
                        yes = "Valvatidae", no = Genus),
         Species = ifelse(test = Species == "spp",
                          yes = NA, no = Species),
         Species = ifelse(test = Species == "Juveniles",
                          yes = "juveniles", no = Species),
         Species = ifelse(test = Species == "maacki",
                          yes = "maackii", no = Species),
         Subspecies = ifelse(test = Subspecies == "virids",
                          yes = "viridis", no = Subspecies)) %>%
  unite(col = "taxon", Genus, Species, Subspecies) %>%
  mutate(taxon = gsub(pattern = "_NA_NA", replacement = "", x = taxon),
         taxon = gsub(pattern = "_NA", replacement = "", x = taxon),
         taxon = ifelse(test = taxon == "flatworms_",
                        yes = "Flatworms", no = taxon)) %>%
  spread(key = taxon, value = sum_count) %>%
  select(-Total) %>%
  separate(col = location, into = c("place", "number"), sep = -1) %>%
  unite(col = "site", place, number, sep = "-")
  
head(inverts_wide)

write.csv(x = inverts_wide, file = "../clean_disaggregated_data/invertebrates.csv",
          row.names = FALSE)


# 6. Load and clean periphyton data ---------------------------------------

periphyton_orig <- read.csv(file = "../original_data/periphyton_20180917.csv",
                            header = TRUE)

# Make long format with relational columns
periphyton_wide <- periphyton_orig %>%
  rename("replicate" = "rep",
         "subsamples_counted" = "counts",
         "tetrasporales" = "tetraporales") %>%
  separate(col = site, into = c("Location", "Number"), sep = -1) %>%
  unite(col = "site", Location, Number, sep = "-") %>%
  select(-date, -Lyngbya) %>%
  clean_names()

head(periphyton_wide)

write.csv(x = periphyton_wide, file = "../clean_disaggregated_data/periphyton.csv",
          row.names = FALSE)


# 7. Load and clean stable isotope data -----------------------------------

stable_isotopes_orig <- read.csv(file = "../original_data/sia_results_mfm_20170509.csv",
                                 header = TRUE)

# Parse identifier data into site and taxonomic data
# In this analysis, "Splash zone" was used as the shorthand 
# for periphyton. We relabel "Splash zone" as periphyton as such. 
stable_isotopes <- stable_isotopes_orig %>%
  separate(col = Identifier, into = c("site", "Genus", "Species"), sep = " ") %>%
  mutate(Genus = ifelse(test = Genus == "E.",
                        yes = "Eulimnogammarus", Genus),
         Genus = ifelse(test = Genus == "P.",
                        yes = "Pallasea", no = Genus),
         Genus = ifelse(test = Genus == "Sp.",
                        yes = "Periphyton", no = Genus),
         Species = ifelse(test = Species == "can",
                          yes = "cancellus", no = Species),
         Species = ifelse(test = Species == "ver",
                          yes = "verrucosus", no = Species),
         Species = ifelse(test = Species == "vitetus",
                          yes = "vittatus", no = Species),
         Species = ifelse(test = Species == "veruossus",
                          yes = "verrucosus", no = Species),
         Species = ifelse(test = Species == "cyan",
                          yes = "cyaneus", no = Species),
         Species = ifelse(test = Species == "zone",
                          yes = NA, no = Species))

head(stable_isotopes)

write.csv(x = stable_isotopes, file = "../clean_disaggregated_data/stable_isotopes.csv",
          row.names = FALSE)


# 8. Load and clean fatty acid data ---------------------------------------

fatty_acid_orig <- read.csv(file = "../original_data/BaikalFAs_wt_20180322.csv",
                            header = TRUE)

# Parse spp column into taxonomic data
fatty_acid <- fatty_acid_orig %>%
  clean_names() %>%
  select(-gc_id, -sample) %>%
  separate(col = spp, into = c("Genus", "Species")) %>%
  mutate(Genus = ifelse(test = Genus == "E" & Species == "ver",
                        yes = "Eulimnogammarus", no = Genus),
         Genus = ifelse(test = Genus == "E" & Species == "vitatus",
                        yes = "Eulimnogammarus", no = Genus),
         Genus = ifelse(test = Genus == "E" & Species == "cyan",
                        yes = "Eulimnogammarus", no = Genus),
         Genus = ifelse(test = Genus == "P" & Species == "can",
                        yes = "Pallasea", no = Genus),
         Genus = ifelse(test = Genus == "Spl",
                        yes = "Periphyton", no = Genus),
         Species = ifelse(test = Genus == "Eulimnogammarus" & Species == "ver",
                          yes = "verrucosus", no = Species),
         Species = ifelse(test = Genus == "Eulimnogammarus" & Species == "vitatus",
                          yes = "vittatus", no = Species),
         Species = ifelse(test = Genus == "Eulimnogammarus" & Species == "cyan",
                          yes = "cyaneus", no = Species),
         Species = ifelse(test = Genus == "Pallasea" & Species == "can",
                          yes = "cancellus", no = Species),
         Species = ifelse(test = Species == "zone",
                          yes = NA, no = Species)) %>%
  rename(site = location) %>%
  select(site:c24_0, comments)

head(fatty_acid)

write.csv(x = fatty_acid, file = "../clean_disaggregated_data/fatty_acid.csv",
          row.names = FALSE)


# 9. Load and clean total lipid data --------------------------------------

total_lipid_orig <- read.csv(file = "../original_data/Baikal.total.lipid.mfm.20180322.csv",
                             header = TRUE)

# Parse sample.id column into site and taxonomic data
total_lipid <- total_lipid_orig %>%
  select(-sample.num) %>%
  separate(col = sample.id, into = c("site", "SPP"), sep = "\\,") %>%
  separate(col = SPP, into = c("Genus", "Species"), sep = "[.]") %>%
  mutate(Genus = ifelse(test = Genus == " E",
                        yes = "Eulimnogammarus", no = Genus),
         Genus = ifelse(test = Genus == "E",
                        yes = "Eulimnogammarus", no = Genus),
         Genus = ifelse(test = Genus == "E" & Species == "ver",
                        yes = "Eulimnogammarus", no = Genus),
         Genus = ifelse(test = Genus == "E" & Species == "vitatus",
                        yes = "Eulimnogammarus", no = Genus),
         Genus = ifelse(test = Genus == "E" & Species == "cyan",
                        yes = "Eulimnogammarus", no = Genus),
         Genus = ifelse(test = Genus == " P" & Species == "can",
                        yes = "Pallasea", no = Genus),
         Genus = ifelse(test = Genus == " Spl",
                        yes = "Periphyton", no = Genus),
         Genus = ifelse(test = grepl(pattern = "Drapa", x = Genus),
                        yes = "Drapa", no = Genus),
         Genus = ifelse(test = grepl(pattern = "Hyalella", x = Genus),
                        yes = "Hyalella", no = Genus),
         Genus = ifelse(test = grepl(pattern = "Snails", x = Genus),
                        yes = "Snails", no = Genus),
         Species = ifelse(test = Genus == "Eulimnogammarus" & Species == "ver",
                          yes = "verrucosus", no = Species),
         Species = ifelse(test = Genus == "Eulimnogammarus" & Species == "ever",
                          yes = "verrucosus", no = Species),
         Species = ifelse(test = Genus == "Eulimnogammarus" & grepl(pattern = "verucossus", x = Species),
                          yes = "verrucosus", no = Species),
         Species = ifelse(test = Genus == "Eulimnogammarus" & grepl(pattern = "vitatus", x = Species),
                          yes = "vitatus", no = Species),
         Species = ifelse(test = Genus == "Eulimnogammarus" & Species == "cyan",
                          yes = "cyaneus", no = Species),
         Species = ifelse(test = Genus == "Eulimnogammarus" & grepl(pattern = "cyan", x = Species),
                          yes = "cyaneus", no = Species),
         Species = ifelse(test = Genus == "Pallasea" & Species == "can",
                          yes = "cancellus", no = Species),
         Species = ifelse(test = Species == "zone",
                          yes = NA, no = Species)) %>%
  rename(total_lipid_mg_per_g = total.lipid.mg.g)

head(total_lipid)

write.csv(x = total_lipid, file = "../clean_disaggregated_data/total_lipid.csv",
          row.names = FALSE)


# 10. Load and clean microplastics data -----------------------------------

microplastics_orig <- read.csv(file = "../original_data/microplastics_mfm_20171010.csv",
                               header = TRUE)

# Here we convert the volume per filtration and the number of times
# that volume was filtered for a total volume filtered. 

microplastics <- microplastics_orig %>%
  select(-date) %>%
  unite(col = "site", location, site, sep = "-") %>%
  mutate(volume_filtered_ml = volume * volume_rep) %>%
  rename("replicate" = "rep") %>%
  select(-volume, -volume_rep)

head(microplastics)

write.csv(x = microplastics, file = "../clean_disaggregated_data/microplastics.csv",
          row.names = FALSE)


# 11. Load and clean calculated distance metrics --------------------------

# Load the shapefiles and metadata containing information of the development
# polygons and the lake shorelines.

baikal_shapefile <- sf::st_read(dsn = "../clean_disaggregated_data/Baikal_shapefile.kml")


metadata <- read.csv(file = "../clean_disaggregated_data/metadata.csv",
                     header = TRUE, stringsAsFactors = FALSE)


# Extract and format development area data 
loc_areas <- baikal_shapefile %>%
  filter(!grepl("shoreline", Name)) %>%
  st_area() %>%
  enframe(name = NULL) %>%
  rename(development_area_m2 = value) %>%
  mutate(development_area_m2 = as.numeric(development_area_m2),
         development_area_km2 = development_area_m2 / 1000000) %>%
  cbind(baikal_shapefile$Name) %>%
  rename(site = `baikal_shapefile$Name`) %>%
  mutate(site = as.character(site),
         site = ifelse(test = site == "Bolshoe Goloustnoe",
                       yes = "BGO", no = site),
         site = ifelse(test = site == "Bolshie Koty",
                       yes = "BK", no = site),
         site = ifelse(test = site == "Listvyanka",
                       yes = "LI", no = site)) %>%
  filter(site %in% c("BK", "LI", "BGO"))


# Extract Development names 
shoreline_names <- baikal_shapefile %>%
  filter(grepl(pattern = "shoreline", x = Name)) %>%
  as_tibble()


# Extract shoreline length info and join with area
loc_shoreline_area_length <- baikal_shapefile %>%
  filter(grepl("shoreline", Name)) %>%
  st_length() %>%
  as_tibble() %>%
  rename(development_shoreline_length_m = value) %>%
  mutate(development_shoreline_length_m = as.numeric(development_shoreline_length_m),
         development_shoreline_length_km = development_shoreline_length_m / 1000) %>%
  cbind(shoreline_names$Name) %>%
  rename(site = `shoreline_names$Name`) %>%
  mutate(site = as.character(site),
         site = ifelse(grepl("Bolshoe Goloustnoe", site), "BGO", site),
         site = ifelse(grepl("Bolshie Koty", site), "BK", site),
         site = ifelse(grepl("Listvyanka", site), "LI", site)) %>%
  filter(site %in% c("BK", "LI", "BGO")) %>%
  full_join(x = ., y = loc_areas, by = c("site"))


# Calculate centroids for each developed site 
baikal_shapefile$centroids <- st_centroid(x = baikal_shapefile)


# Convert metadata into a spatial object 
site_loc <- metadata %>%
  dplyr::select(site, lat, long)

site_loc_pts <- st_as_sf(site_loc[, 2:3], coords = c("long", "lat"),
                         crs = 4326)


# Join all data together and format into a dataframe
locs_centroids <- st_distance(x = site_loc_pts,
                              y = baikal_shapefile$centroids$geometry[1:3]) %>%
  as_tibble() %>%
  rename(BGO = V1,
         BK = V2,
         LI = V3) %>%
  mutate(BGO = as.numeric(BGO),
         BK = as.numeric(BK),
         LI = as.numeric(LI)) %>%
  cbind(., site_loc) %>%
  gather(key = nearest_neighbor, value = distance, BGO:LI) %>%
  group_by(site, lat, long) %>%
  mutate(distance_km = distance / 1000,
         population = ifelse(test = nearest_neighbor == "BGO",
                             yes = 600, no = NA),
         population = ifelse(test = nearest_neighbor == "BK",
                             yes = 80, no = population),
         population = ifelse(test = nearest_neighbor == "LI",
                             yes = 5000, no = population)) %>%
  left_join(x = ., y = loc_shoreline_area_length,
            by = c("nearest_neighbor" = "site")) %>%
  mutate(distance_weighted_population = ((population * development_shoreline_length_km) /
                                           development_area_km2) / distance_km) %>%
  group_by(site) %>%
  summarize(distance_weighted_population = sum(distance_weighted_population)) %>%
  arrange(distance_weighted_population)

write.csv(x = locs_centroids,
          file = "../clean_disaggregated_data/distance_weighted_population_metrics.csv",
          row.names = FALSE)
