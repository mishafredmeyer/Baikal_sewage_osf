# This script creates the map of sampling locations for the associated manuscript.
# The script also merges inverse distance weighted population metrics with sampling
# location to reveal the IDW population gradient over space.

library(tidyverse)
library(OpenStreetMap)
library(ggpubr)
library(cowplot)
library(ggrepel)
library(viridis)
library(ggsn)
library(sp)
library(ggspatial)
library(cowplot)
library(sf)


# Check if figures directory exists 
# If not, create the figures directory
sub_dir <- "figures"
output_dir <- file.path(here::here(), sub_dir)

if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Dir 'figures' already exists!")
}


# Step 1: Import data and join datasets -----------------------------------

site_information <- read.csv(file = "../cleaned_data/site_information.csv",
                          stringsAsFactors = FALSE)

distance <- read.csv(file = "../cleaned_data/distance_weighted_population_metrics.csv",
                     header = TRUE, stringsAsFactors = FALSE)

sample_points <- full_join(x = site_information,
                           y = distance,
                           by = c("site"))

# Make the locations Mercator
sample_points_merc <- projectMercator(lat = sample_points$lat,
                                      long = sample_points$long) %>%
  as.data.frame() %>%
  bind_cols(., sample_points)


# Step 2: Start building the map ------------------------------------------

# Get a basemap
# Options: https://www.r-bloggers.com/the-openstreetmap-package-opens-up/
base_map <- openmap(upperLeft = c(55.915113, 102.2324553),
                    lowerRight = c(51.1800703, 110.8),
                    type = "esri") %>%
  openproj()

# Get a zoomed in map
base_map_zoom <- openmap(upperLeft = c(52.15, 104.75),
                         lowerRight = c(51.75, 105.55),
                         type = "bing", zoom = 11) %>%
  openproj() 


# Build the inset map with the basemap
inset_map <- autoplot.OpenStreetMap(base_map) +
  geom_point(data = sample_points_merc,
             aes(x = long, y = lat,
                 fill = log10(distance_weighted_population)),
             alpha = 1,  color = "grey70", size = 3,
             shape = 21) +
  scale_fill_viridis(option = "plasma") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.background = element_rect(fill = "snow1")) +
  annotate(geom = "text", label = "Irkutsk Oblast",
           x = 105, y = 54, color = "black", size = 5) +
  annotate(geom = "text", label = "Republic of\nBuryatiya",
           x = 109, y = 52.35, color = "black", size = 5) +
  coord_quickmap() + 
  xlab("") +
  ylab("") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# Build a close up map with satellite imagery & zoomed in map
zoom_map <- autoplot.OpenStreetMap(base_map_zoom) +
  geom_point(data = sample_points_merc, 
             aes(x = long, y = lat,
                 size = (distance_weighted_population),
                 fill = (distance_weighted_population)),
             alpha = 0.8,  color = "grey70", shape = 21,
             stroke = 2.5) +
  #geom_sf(data = st_as_sf(sample_points_merc, coords = c("long", "lat"), remove = FALSE), 
  #        aes(geometry = geometry, x = long, y = lat)) +
  scale_fill_viridis(option = "plasma", name = "IDW Population\n(Persons)", trans ="log10") +
  scale_size_continuous(range = c(1, 20), guide = "none", trans ="log10") +
  scalebar(data = sf::st_as_sf(sample_points_merc[ , -c(1:2)], coords = c("long", "lat")), 
           anchor = c(x = 105.43, y = 51.80), dist = 10,
           dist_unit = "km", transform = TRUE, model = "WGS84",
           st.bottom = FALSE, st.size = 8, st.dist = 0.1, 
           st.color = "white", height = 0.05, 
           box.color = "white") + 
  xlab("Longitude") +
  ylab("Latitude") +
  annotate(geom = "text", label = "Bolshoe Goloustnoe",
           x = 105.4, y = 52.06, color = "white", size = 10) +
  annotate(geom = "text", label = "Bolshie Koty",
           x = 105.075, y = 51.935, color = "white", size = 10) +
  annotate(geom = "text", label = "Listvyanka",
           x = 104.85, y = 51.9, color = "white", size = 10) +
  coord_quickmap() + 
  theme(legend.key.height = unit(1, "in"),
        legend.key.width = unit(0.65, "in"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        panel.background = element_blank(),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 20)) 

# Combine both maps. Resource used for this:
# https://stackoverflow.com/questions/5219671/it-is-possible-to-create-inset-graphs
baikal_combine <- ggdraw() +
  draw_plot(zoom_map) +
  draw_plot(inset_map, x = 0.01, y = 0.57, width = .5, height = .33, scale = 1) 

ggarrange(baikal_combine, arranged_maps, ncol = 2, nrow = 1,
          widths = c(2.95,1), labels = c("A"), font.label = list(size = 20))

# This plot is associated with Figure 1 in the accompanying manuscript.
ggsave(filename = "../figures/baikal_map.png", 
       plot = ggarrange(baikal_combine, arranged_maps, 
                        ncol = 2, nrow = 1,
                        widths = c(2.95,1), labels = c("A"), 
                        font.label = list(size = 20)), 
       device = "png", width = 18, height = 9, units = "in")


# Step 3: Create individual maps for developments -------------------------

site_information <- read.csv(file = "../cleaned_data/site_information.csv",
                             stringsAsFactors = FALSE)

distance <- read.csv(file = "../cleaned_data/distance_weighted_population_metrics.csv",
                     header = TRUE, stringsAsFactors = FALSE)

baikal_shapefile <- sf::st_read(dsn = "../original_data/Baikal_shapefile.kml") %>%
  filter(!grepl(pattern = "shoreline", x = Name))

## Combine site informaiton and distance dfs
## Also make object an sf object

sample_points <- full_join(x = site_information,
                           y = distance,
                           by = c("site")) %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326)


# Step 2: Make a map of all of Baikal -------------------------------------

## Listvyanka map

listvyanka <- ggplot() +
  annotation_map_tile(
    type = 'osm',
    zoom = 15,
    cachedir = "data/map_tiles/") +
  geom_sf(data = baikal_shapefile %>%
            filter(Name == "Listvyanka"), 
          color = "black", alpha = 0.4, fill = viridis(10)[1]) +
  annotation_scale(location = "bl") +
  xlab("") +
  ylab("") +
  coord_sf(xlim = c(104.84, 104.89), ylim = c(51.84, 51.87), crs = 4326) +
  scale_x_continuous(breaks = seq(104.84, 104.89, by = 0.02)) +
  scale_y_continuous(breaks = seq(51.84, 51.87, by = 0.015)) +
  theme_bw() +
  theme(panel.background = element_blank(),
        axis.text = element_text(size = 12)) 


bolshie_koty <- ggplot() +
  annotation_map_tile(
    type = 'osm',
    zoom = 15,
    cachedir = "data/map_tiles/") +
  geom_sf(data = baikal_shapefile %>%
            filter(Name == "Bolshie Koty"), 
          color = "black", alpha = 0.4, fill = viridis(10)[1]) +
  annotation_scale(location = "br") +
  xlab("") +
  ylab("") +
  coord_sf(xlim = c(105.06, 105.086), ylim = c(51.9, 51.91), crs = 4326) +
  scale_x_continuous(breaks = seq(105.06, 105.086, by = 0.015)) +
  scale_y_continuous(breaks = seq(51.9, 51.91, by = 0.005)) +
  theme_bw() +
  theme(panel.background = element_blank(),
        axis.text = element_text(size = 12))

bolshoe_goloustnoe <- ggplot() +
  annotation_map_tile(
    type = 'osm',
    zoom = 15,
    cachedir = "data/map_tiles/") +
  geom_sf(data = baikal_shapefile %>%
            filter(Name == "Bolshoe Goloustnoe"), 
          color = "black", alpha = 0.4, fill = viridis(10)[1]) +
  annotation_scale(location = "br") +
  xlab("") +
  ylab("") +
  coord_sf(xlim = c(105.392, 105.4175), ylim = c(52.025, 52.052), crs = 4326) +
  scale_x_continuous(breaks = seq(105.392, 105.4175, by = 0.02)) +
  scale_y_continuous(breaks = seq(52.025, 52.052, by = 0.01)) +
  theme_bw() +
  theme(panel.background = element_blank(),
        axis.text = element_text(size = 12))


# Step 4: Join the maps and save the file ---------------------------------

arranged_maps <- ggarrange(listvyanka, bolshie_koty, bolshoe_goloustnoe, 
                           ncol = 1, nrow = 3, labels = c("B", "C", "D"), 
                           font.label = list(size = 20))

baikal_map <- ggdraw() +
  draw_plot(zoom_map) +
  draw_plot(plot = inset_map,
            x = -0.09, y = 0.54, width = .7, height = .4, scale = 1) 

ggsave(filename = "../figures/baikal_map.png", plot = baikal_map, device = "png",
       width = 18, height = 9, units = "in")

