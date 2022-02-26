
# Generate random training points to be manually labelled for phenology classifier 

library(raster)
library(rgdal)
library(here)

# This is a raster image which covers the spatial exent we're interested in classifying
pheno_image <- raster(here::here("pheno_image_example.tif"))

# Coordinates for each raster cell in the image
image_coordinates <- coordinates(pheno_image)

# Randomly sample points in image
num_points <- 1000
point_indices <- sample(1:nrow(image_coordinates), num_points)
point_coords <- image_coordinates[point_indices,]
# Create SpatialPointsDataFrame
point_data = data.frame(class = as.character(rep("", num_points)))
point_df <- SpatialPointsDataFrame(point_coords, point_data)
crs(point_df) <- crs(pheno_image)

# Output Points
writeOGR(point_df, dsn=here::here("training_data_random.shp"), layer="training_data_random", driver="ESRI Shapefile", overwrite_layer=TRUE)