
library(raster)
library(tidyverse)
library(janitor)
library(here)

ripdel <- stack(here::here("riparian_masks/santa_ynez_front_range_classes_2006.tif"),
                here::here("riparian_masks/santa_ynez_front_range_classes_2007.tif"),
                here::here("riparian_masks/santa_ynez_front_range_classes_2008.tif"),
                here::here("riparian_masks/santa_ynez_front_range_classes_2009.tif"),
                here::here("riparian_masks/santa_ynez_front_range_classes_2014.tif"),
                here::here("riparian_masks/santa_ynez_front_range_classes_2015.tif"),
                here::here("riparian_masks/santa_ynez_front_range_classes_2016.tif"),
                here::here("riparian_masks/santa_ynez_front_range_classes_2017.tif"),
                here::here("riparian_masks/santa_ynez_front_range_classes_2018.tif"),
                here::here("riparian_masks/santa_ynez_front_range_classes_2019.tif") )

# Load SRTM elevation data, reproject to riparian mask crs
srtm <- raster(here::here("..","srtm_south_coast_depressionless.tif"))
srtm <- projectRaster(srtm, ripdel)

# Bin SRTM data - 400 m bins
srtm_binned <- srtm
srtm_binned[srtm_binned < 400] <- 0
srtm_binned[srtm_binned >= 400 & srtm_binned < 800] <- 1
srtm_binned[srtm_binned >= 800 & srtm_binned < 1200] <- 2
srtm_binned[srtm_binned >= 1200] <- 3

# Combine rasters, convert to data frame
rip_elev <- stack(ripdel, srtm_binned)
rip_df <- as.data.frame(rip_elev)
rip_longer_df <- rip_df %>% 
  pivot_longer(1:10, names_to="year", values_to="riparian")
names(rip_longer_df) <- c("elevation_bin", "year", "riparian")
rip_longer_df$year <- substr(rip_longer_df$year,32,35)
rip_summary <- rip_longer_df %>%
  group_by(year, riparian, elevation_bin) %>%
  tally()
srtm_bin_counts <- as.data.frame(srtm_binned) %>%
  group_by(srtm_south_coast_depressionless) %>%
  tally()
srtm_index <- unlist(lapply(rip_summary$elevation_bin, function(x){which(x == srtm_bin_counts$srtm_south_coast_depressionless)}))
rip_summary$elev_total <- srtm_bin_counts[srtm_index,]$n
rip_summary$elev_frac <- rip_summary$n / rip_summary$elev_total


# Plot Outputs
ggplot(rip_summary %>% filter(riparian==1)) + 
  geom_line(aes(x=year, y=elev_frac, group=elevation_bin, col=elevation_bin))


# Read NDVI timeseries
phenoseries <- read_csv(here::here("summer_greenness","riparian_ndvi_all_landsat.csv"))
# deal with annoying formatting
image_indices <- (1:(nrow(phenoseries)/2))*2 - 1
pheno_df <- data.frame(image = phenoseries[image_indices,]$data,
                       ndvi = as.numeric(phenoseries[image_indices+1,]$data)) %>%
  drop_na()
pheno_df$year = substr(pheno_df$image, 13,16)
pheno_df$month = substr(pheno_df$image, 17,18)
pheno_df$day = substr(pheno_df$image, 19,20)
pheno_df$doy = lubridate::yday(as.Date(paste(pheno_df$year,pheno_df$month,pheno_df$day,sep="-")))

ggplot(pheno_df) %>% 
  # geom_line(aes(x=doy, y=ndvi, group=year, col=year))

spring_df <- pheno_df %>% filter(doy <= 120)
jan_median <- median((spring_df %>% filter(doy < 32))$NDVI)

set.seed(1)

ransacLines <- function(ndvi_data, iterations, max_residual)
{
  inliers_best <- 0
  model_best <- NA
  for(i in 1:iterations)
  {
    # Choose two random seed points for line model
    indices <- c(0,0)
    while(indices[1] == indices[2])
    {
      indices <- sample(1:(nrow(ndvi_data)), 2)
    }
    # Fit line to seed points
    line_points <- data.frame(doy = c(ndvi_data[indices[1],]$doy,
                                      ndvi_data[indices[2],]$doy),
                              ndvi = c(ndvi_data[indices[1],]$ndvi,
                                      ndvi_data[indices[2],]$ndvi))
    line_model <- lm(data=line_points, ndvi ~ doy)
    # Check how many inliers fit the model
    residuals <- predict(line_model, ndvi_data) - ndvi_data$ndvi
    inliers_count <- sum(abs(residuals) < max_residual)
    # If this is a new best model, update 
    if(inliers_count > inliers_best)
    {
      inliers_best <- inliers_count
      model_best <- line_model
    }
  }
  return(model_best)
}

# 7-period approach
spring_1985 <- pheno_df %>% filter(doy <= 120, doy > 30, year %in% seq(1985,1989))
spring_1990 <- pheno_df %>% filter(doy <= 120, doy > 30, year %in% seq(1990,1994))
spring_1995 <- pheno_df %>% filter(doy <= 120, doy > 30, year %in% seq(1995,1999))
spring_2000 <- pheno_df %>% filter(doy <= 120, doy > 30, year %in% seq(2000,2004))
spring_2005 <- pheno_df %>% filter(doy <= 120, doy > 30, year %in% seq(2005,2009))
spring_2010 <- pheno_df %>% filter(doy <= 120, doy > 30, year %in% seq(2010,2014))
spring_2015 <- pheno_df %>% filter(doy <= 120, doy > 30, year %in% seq(2015,2019))

max_iterations <- 5000
max_residual <- 0.06
model_1985 <- ransacLines(spring_1985, max_iterations, max_residual)
model_1990 <- ransacLines(spring_1990, max_iterations, max_residual)
model_1995 <- ransacLines(spring_1995, max_iterations, max_residual)
model_2000 <- ransacLines(spring_2000, max_iterations, max_residual)
model_2005 <- ransacLines(spring_2005, max_iterations, max_residual)
model_2010 <- ransacLines(spring_2010, max_iterations, max_residual)
model_2015 <- ransacLines(spring_2015, max_iterations, max_residual)

plot(spring_1985$doy, spring_1985$ndvi, col=(abs((spring_1985$ndvi - predict(model_1985, spring_1985)))<max_residual))
points(spring_1985$doy, predict(model_1985, spring_1985), col="red")

spring_data <- data.frame(doy = spring_1985$doy, 
                          ndvi = spring_1985$ndvi,
                          ndvi_p = predict(model_1985, spring_1985),
                          period = 1985)
new_data <- data.frame(doy = spring_1990$doy, 
                       ndvi = spring_1990$ndvi,
                       ndvi_p = predict(model_1990, spring_1990),
                       period = 1990)
spring_data <- rbind(spring_data, new_data)
new_data <- data.frame(doy = spring_1995$doy, 
                       ndvi = spring_1995$ndvi,
                       ndvi_p = predict(model_1995, spring_1995),
                       period = 1995)
spring_data <- rbind(spring_data, new_data)
new_data <- data.frame(doy = spring_2000$doy, 
                       ndvi = spring_2000$ndvi,
                       ndvi_p = predict(model_2000, spring_2000),
                       period = 2000)
spring_data <- rbind(spring_data, new_data)
new_data <- data.frame(doy = spring_2005$doy, 
                       ndvi = spring_2005$ndvi,
                       ndvi_p = predict(model_2005, spring_2005),
                       period = 2005)
spring_data <- rbind(spring_data, new_data)
new_data <- data.frame(doy = spring_2010$doy, 
                       ndvi = spring_2010$ndvi,
                       ndvi_p = predict(model_2010, spring_2010),
                       period = 2010)
spring_data <- rbind(spring_data, new_data)
new_data <- data.frame(doy = spring_2015$doy, 
                       ndvi = spring_2015$ndvi,
                       ndvi_p = predict(model_2015, spring_2015),
                       period = 2015)
spring_data <- rbind(spring_data, new_data)

ggplot(spring_data) + geom_line(aes(x=doy, y=ndvi_p, group=period, col=period))








# 4-period approach
spring_1984 <- pheno_df %>% filter(doy <= 120, doy > 30, year %in% seq(1984,1992))
spring_1993 <- pheno_df %>% filter(doy <= 120, doy > 30, year %in% seq(1993,2001))
spring_2002 <- pheno_df %>% filter(doy <= 120, doy > 30, year %in% seq(2002,2010))
spring_2011 <- pheno_df %>% filter(doy <= 120, doy > 30, year %in% seq(2011,2019))

model_1984 <- ransacLines(spring_1984, max_iterations, max_residual)
model_1993 <- ransacLines(spring_1993, max_iterations, max_residual)
model_2002 <- ransacLines(spring_2002, max_iterations, max_residual)
model_2011 <- ransacLines(spring_2011, max_iterations, max_residual)

spring_data <- data.frame(doy = spring_1984$doy, 
                          ndvi = spring_1984$ndvi,
                          ndvi_p = predict(model_1984, spring_1984),
                          period = 1984)
new_data <- data.frame(doy = spring_1993$doy, 
                       ndvi = spring_1993$ndvi,
                       ndvi_p = predict(model_1993, spring_1993),
                       period = 1993)
spring_data <- rbind(spring_data, new_data)
new_data <- data.frame(doy = spring_2002$doy, 
                       ndvi = spring_2002$ndvi,
                       ndvi_p = predict(model_2002, spring_2002),
                       period = 2002)
spring_data <- rbind(spring_data, new_data)
new_data <- data.frame(doy = spring_2011$doy, 
                       ndvi = spring_2011$ndvi,
                       ndvi_p = predict(model_2011, spring_2011),
                       period = 2011)
spring_data <- rbind(spring_data, new_data)

ggplot(spring_data) + geom_line(aes(x=doy, y=ndvi_p, group=period, col=period))
