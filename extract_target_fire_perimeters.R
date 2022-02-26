
# Script to extract target fires from overall NIFC fire perimeter database

library(rgdal)
library(here)

# Get all fire perimeters
all_fire_perimeters <- readOGR(here::here("fire_perimeters","InteragencyFirePerimeterHistory.shp"))
# Drop perimeters with NA for incident name
all_valid_fire_perimeters <- all_fire_perimeters[!is.na(all_fire_perimeters$INCIDENT),]

# Get and save individual fire perimeters for:
#  2008 Tea Fire
#    NOTE for some reason the Tea polygon is in here twice - so only taking the first one
tea <- all_valid_fire_perimeters[all_valid_fire_perimeters$INCIDENT == "TEA",]
writeOGR(tea[1,], here::here("fire_perimeters","tea.shp"), layer="tea", driver="ESRI Shapefile", overwrite=TRUE)
#  2009 Jesusita Fire
jesusita <- all_valid_fire_perimeters[all_valid_fire_perimeters$INCIDENT == "JESUSITA",]
writeOGR(jesusita, here::here("fire_perimeters","jesusita.shp"), layer="jesusita", driver="ESRI Shapefile", overwrite=TRUE)
#  2016 Sherpa Fire
sherpa <- all_valid_fire_perimeters[all_valid_fire_perimeters$INCIDENT == "SHERPA",]
writeOGR(sherpa, here::here("fire_perimeters","sherpa.shp"), layer="sherpa", driver="ESRI Shapefile", overwrite=TRUE)
#  2017 Whittier Fire
whittier <- all_valid_fire_perimeters[all_valid_fire_perimeters$INCIDENT == "WHITTIER",]
writeOGR(whittier, here::here("fire_perimeters","whittier.shp"), layer="whittier", driver="ESRI Shapefile", overwrite=TRUE)
#  2017 Thomas Fire
thomas <- all_valid_fire_perimeters[all_valid_fire_perimeters$INCIDENT == "THOMAS",]
writeOGR(thomas, here::here("fire_perimeters","thomas.shp"), layer="thomas", driver="ESRI Shapefile", overwrite=TRUE)
#  2019 Cave Fire
#   AGAIN here there are two nonstandard things - first Cave is not all caps, and there are extra extraneous records, too. 
cave <- all_valid_fire_perimeters[all_valid_fire_perimeters$INCIDENT == "Cave",]
writeOGR(cave[4,], here::here("fire_perimeters","cave.shp"), layer="cave", driver="ESRI Shapefile", overwrite=TRUE)
