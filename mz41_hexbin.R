
# Notes ----
# make hexbin map of BpSs for LF Map Zone 41 (ne Minnesota, nw Wisconsin)
# Randy Swaty
# May 28, 2025


# Set up ----


# load packages
library(exactextractr)
library(foreign)
library(htmltools)
library(htmlwidgets)
library(leaflet)
library(rlandfire)
library(sf)
library(terra)
library(tidyverse)



shp <- st_read("inputs/mz41_us.shp", quiet = TRUE) %>% 
  st_transform(crs = 5070) %>%
  st_union() %>%
  st_sf()


# read in .csv of CONUS-wide attributes

bps_conus_atts <- read_csv("inputs/LF20_BPS_220.csv")

# create an accessible pallette (from https://thenode.biologists.com/data-visualization-with-flying-colors/research/)
custom_palette <- c("#ffbf00", "#67d8ff", "#00be8a", "#ffff4f", "#0089d6", "#ff7100", "#f591c8", "#ededed", "#ffcd02", "#615641")




# Download BpS data ----

aoi <- getAOI(shp)

products <-  c("200BPS")
projection <- 5070  
resolution <- 30
email <- "rswaty@tnc.org" # Replace with your email address. 

# R specific arguments
save_file <- tempfile(fileext = ".zip")

# call API
ncal <- landfireAPIv2(
  products, 
  aoi, 
  projection, 
  resolution, 
  path = save_file,
  email = email)


# define the destination path
dest_file <- file.path("inputs", "landfire_data.zip")

# move and rename the file
file.rename(save_file, dest_file)

# create a temporary directory for unzipping
temp_dir <- tempfile()
dir.create(temp_dir)

# unzip the file into the temporary directory
unzip(dest_file, exdir = temp_dir)

# get the list of unzipped files
unzipped_files <- list.files(temp_dir, full.names = TRUE)

# rename each unzipped file to "landfire_data" with its full original extension
for (file in unzipped_files) {
  file_name <- basename(file)
  file_extension <- sub("^[^.]*", "", file_name)  # Extract the full extension
  new_file_path <- file.path("inputs", paste0("landfire_data", file_extension))
  file.rename(file, new_file_path)
}

# clean up the temporary directory
unlink(temp_dir, recursive = TRUE)

# GIS processing

# load downloaded stacked raster
US_200BPS <- rast("inputs/landfire_data.tif")


# crop and mask the us_200bps raster using the shapefile
bps_aoi <- US_200BPS %>%
  crop(shp) %>%
  mask(shp)

# set the levels of the raster to bps_conus_atts
levels(bps_aoi)[[1]] <- bps_conus_atts
# set the active category of the raster to "value"
activeCat(bps_aoi) <- "VALUE"

# extract values from the raster, remove na values, and convert to dataframe
bps_aoi_atts <- values(bps_aoi, dataframe = T, na.rm = T) %>%
  # create a frequency table of the values
  table(dnn = "VALUE") %>%
  # convert the table to a dataframe
  as.data.frame() %>%
  # convert all columns to character type
  mutate_all(as.character) %>%
  # convert all columns to integer type
  mutate_all(as.integer) %>%
  # join the dataframe with the raster categories
  left_join(cats(bps_aoi)[[1]], by = "VALUE") %>%
  # filter out rows with a frequency of 0
  filter(Freq != 0) %>%
  # calculate acres and relative percentage
  mutate(ACRES = round((Freq * 900 / 4046.86), 0),
         REL_PERCENT = round((Freq / sum(Freq)), 3) * 100) %>%
  # arrange the dataframe by relative percentage in descending order
  arrange(desc(REL_PERCENT))

# write the raster to a file with specified options
# writeRaster(bps_aoi, "outputs/bps_aoi.tif",
#             gdal = c("COMPRESS=NONE", "TFW=YES"),
#             datatype = "INT2S",
#             overwrite = T)
# 
# # write the attributes dataframe to a dbf file
# write.dbf(bps_aoi_atts, "outputs/bps_aoi.tif.vat.dbf")
# 
# # write the attributes dataframe to a csv file
# write.csv(bps_aoi_atts, "outputs/bps_aoi_attributes.csv")


# Hexbin processing ----

# Calculate the side length for hexagons with 200 km² area
side_length <- sqrt((2 * 200000000) / (3 * sqrt(3)))

# Create a hexagon grid within the shapefile's bounding box
hex_grid <- st_make_grid(shp, cellsize = side_length, square = FALSE, what = 'polygons')

# Convert the grid to an sf object
hex_grid <- st_sf(geometry = hex_grid)

# Add an index column to the hexagon grid
hex_grid <- hex_grid %>%
  mutate(index = row_number())


# Filter hexagons by centroid using st_intersects
hex_grid <- hex_grid %>%
  filter(lengths(st_intersects(st_centroid(hex_grid), shp)) > 0)


# Ensure attributes are correctly handled
hex_grid <- hex_grid %>%
  st_set_geometry(st_geometry(hex_grid))



# Perform the exact_extract operation
bps_majority_hex <- exact_extract(bps_aoi, hex_grid, 'majority', append_cols = "index") %>%
  left_join(select(bps_aoi_atts,
                   VALUE,
                   BPS_MODEL,
                   BPS_NAME),
            by = c('majority' = 'VALUE'))



# Join both BpS attributes to hex shapefile
hexs_bps<- hex_grid %>%
  left_join(bps_majority_hex, by = 'index') 




# Make hexbin leaflet ----

# get top 10 BpSs (probably unnecesary)
bps_name_10 <- bps_aoi_atts %>%
  group_by(BPS_NAME) %>%
  summarize(ACRES = sum(ACRES),
            REL_PERCENT = sum(REL_PERCENT)) %>%
  arrange(desc(REL_PERCENT)) %>%
  top_n(n = 10, wt = REL_PERCENT)

# get bps names for palette
bps_names <- unique(bps_name_10$BPS_NAME)

# make palette
bps_colors <- setNames(custom_palette[1:length(bps_names)], bps_names)

pal_bps <- colorFactor(
  palette = bps_colors,
  domain = hexs_bps$BPS_NAME
)


# Transform the hex_grid to WGS84
hexs_bps_4326 <- st_transform(hexs_bps, crs = 4326)



map <- leaflet(hexs_bps_4326) %>%
  addProviderTiles("CartoDB.Positron") %>%
  setView(lng = -93, lat = 47, zoom = 7) %>%
  addPolygons(
    fillColor = ~pal_bps(BPS_NAME),
    color = "#BDBDC3",
    weight = 1,
    opacity = 1,
    fillOpacity = 1.0,
    highlightOptions = highlightOptions(
      weight = 2,
      color = "#666",
      fillOpacity = 1.0,
      bringToFront = TRUE
    ),
    label = ~BPS_NAME,
    labelOptions = labelOptions(
      style = list("font-weight" = "normal", padding = "3px 8px"),
      textsize = "15px",
      direction = "auto"
    )) %>%
  addScaleBar(position = "bottomleft", options = scaleBarOptions(metric = TRUE, imperial = FALSE)) %>%

  addControl("<h3>Top LANDFIRE Biophysical Settings (historical ecosystems) per hexagon for Map Zone 41 <br> 
             Hexagons are ~200km, more information on BpSs at <a href='https://landfire.gov/vegetation/bps' target='_blank'>https://landfire.gov/vegetation/bps</a></h3>", 
             position = "topright", className = "map-title") %>%
  addControl(tags$style(HTML("
  .map-title {
    background-color: white;
    border: 2px solid #666;
    padding: 25px;
    border-radius: 5px;
    font-size: 16px;
    font-weight: bold;
    text-align: right; /* Center the title */
  }
")), position = "topright")



map


# Save the map as a self-contained HTML file
saveWidget(map, "bps_map.html", selfcontained = TRUE)

# Calculate annual acres burned for top 10 BpSs


annualFire <- bps_aoi_atts%>%
  mutate(annual_fire_acres = ((1/FRI_ALLFIR)*ACRES)) %>%
  filter(BPS_NAME != 'Open Water') %>%
  group_by(BPS_NAME) %>%
  summarize(acres = sum(annual_fire_acres)) %>%
  arrange(desc(acres)) %>%
  top_n(n = 10, wt = acres)

print(sum(annualFire$acres))

1,956,152

