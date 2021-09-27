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
library(sf)
library(ggspatial)
library(cowplot)


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

baikal_shapefile <- sf::st_read(dsn = "../original_data/Baikal_shapefile.kml") %>%
  filter(!grepl(pattern = "shoreline", x = Name))

## Combine site informaiton and distance dfs
## Also make object an sf object

sample_points <- full_join(x = site_information,
                           y = distance,
                           by = c("site")) %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326)


# Step 2: Make a map of all of Baikal -------------------------------------


full_baikal_map <- ggplot() +
  annotation_map_tile(
    type = "https://tile.openstreetmap.org/${z}/${x}/${y}.png",
    zoom = 6,
    cachedir = "data/map_tiles/") +
  # geom_sf(data = developments, color = "black", alpha = 0.5, fill = viridis(10)[2]) +
  geom_sf(data = sample_points,
          mapping = aes(fill = distance_weighted_population),
          color = "black", 
          pch = 21, size = 3, alpha = 0.7) +
  scale_fill_viridis(option = "plasma", name = "IDW Population", trans = "log10",
                     # breaks = c(10, 50, 100, 200),
                     # labels = c(10, 50, 100, 200),
                     guide = guide_colorbar(barwidth = 4,
                                            barheight = 20),
                     na.value = "grey50") +
  #scale_x_continuous(breaks = c(-114.45, -114.25, -114.05, -113.85)) +
  #scale_y_continuous(breaks = c(47.70, 47.82, 47.94, 48.06)) +
  coord_sf(xlim = c(102.2324553, 110.8), ylim = c(51.1800703, 55.915113), crs = 4326) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.background = element_rect(fill = "snow1")) +
  annotate(geom = "text", label = "Irkutsk Oblast",
           x = 105, y = 54, color = "black", size = 7) +
  annotate(geom = "text", label = "Republic of\nBuryatiya",
           x = 109, y = 52.4, color = "black", size = 7) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.background = element_rect(color = "black", fill="snow2", size=3))


# Step 3: Create the zoomed in map of the sampling area -------------------


zoomed_in_baikal_map <- ggplot() +
  annotation_map_tile(
    type = 'https://stamen-tiles.a.ssl.fastly.net/terrain-background/${z}/${x}/${y}.png',
    zoom = 11,
    cachedir = "data/map_tiles/") +
  geom_sf(data = baikal_shapefile, color = "black", alpha = 0.95, fill = viridis(10)[1]) +
  geom_sf(data = sample_points, 
          aes(fill = distance_weighted_population,
              size = distance_weighted_population), 
          alpha = 0.6,  color = "grey70", shape = 21,
          stroke = 1) +
  scale_fill_viridis(option = "plasma", name = "IDW Population", trans = "log10",
                     # breaks = c(10, 50, 100, 200),
                     # labels = c(10, 50, 100, 200),
                     guide = guide_colorbar(barwidth = 4,
                                            barheight = 20), 
                     na.value = "grey50") +
  scale_size_continuous(range = c(1, 15), guide = "none", trans = "log10") +
  annotate(geom = "text", label = "Bolshoe Goloustnoe",
           x = 105.35, y = 52.06, color = "white", size = 10) +
  annotate(geom = "text", label = "Bolshie Koty",
           x = 105.075, y = 51.935, color = "white", size = 10) +
  annotate(geom = "text", label = "Listvyanka",
           x = 104.85, y = 51.9, color = "white", size = 10) +
  annotation_scale(location = "bl") +
  xlab("") +
  ylab("") +
  coord_sf(xlim = c(104.75, 105.55), ylim = c(51.75, 52.15), crs = 4326) +
  theme_bw() +
  theme(legend.key.height = unit(2, "in"),
        legend.key.width = unit(0.85, "in"),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 36),
        panel.background = element_blank(),
        axis.title = element_text(size = 36),
        axis.text = element_text(size = 30)) 


# Step 4: Join the maps and save the file ---------------------------------


baikal_map <- ggdraw() +
  draw_plot(zoomed_in_baikal_map) +
  draw_plot(plot = full_baikal_map,
            x = 0, y = 0.62, width = .6, height = .35, scale = 1) 

ggsave(filename = "../figures/baikal_map.png", plot = baikal_map, device = "png",
       width = 18, height = 9, units = "in")
