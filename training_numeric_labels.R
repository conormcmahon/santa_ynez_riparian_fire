
# Add numeric indicator for class to vector training data, to be used in Google Earth Engine

library(raster)
library(rgdal)
library(here)

training <- readOGR(here::here("riparian_gee_training_v3.shp"))

training$numeric_class <- as.numeric(training$class)

writeOGR(training, dsn=here::here("riparian_gee_training_v4.shp"), layer="riparian_gee_training_v4", driver="ESRI Shapefile")
