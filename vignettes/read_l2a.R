library(rGEDI)
library(sf)

bbox_4326 <- st_bbox(c(xmin = 139.40119, ymin = 36.86163, 
                       xmax = 139.52543, ymax = 36.94841), 
                     crs = 4326)

# Study area boundary box coordinates
ul_lat<- 36.94841
lr_lat<- 36.86163
ul_lon<- 139.40119
lr_lon<- 139.52543

# Specifying the date range
daterange=c("2022-07-15","2022-07-31")

# Get path to GEDI data
gLevel1B<-gedifinder(product="GEDI01_B",ul_lat, ul_lon, lr_lat, lr_lon,version="002",daterange=daterange)
gLevel2A<-gedifinder(product="GEDI02_A",ul_lat, ul_lon, lr_lat, lr_lon,version="002",daterange=daterange)
gLevel2B<-gedifinder(product="GEDI02_B",ul_lat, ul_lon, lr_lat, lr_lon,version="002",daterange=daterange)

# Set output dir for downloading the files
outdir=getwd()

# Downloading GEDI data
gediDownload(filepath=gLevel1B[1],outdir=outdir)
gediDownload(filepath=gLevel2A[1],outdir=outdir)
gediDownload(filepath=gLevel2B[1],outdir=outdir)

### Read GEDI Level 1B data
gedilevel1b<-readLevel1B(level1Bpath = file.path(outdir,"GEDI01_B_2022198155950_O20359_03_T02979_02_005_03_V002.h5"))
level1bGeo<-getLevel1BGeo(level1b=gedilevel1b,select=c("elevation_bin0"))
head(level1bGeo)
dim(level1bGeo)
level1bGeo$shot_number<-as.character(level1bGeo$shot_number)

### Read GEDI Level 2A data
gedilevel2a<-readLevel2A(level2Apath = file.path(outdir,"GEDI02_A_2022198155950_O20359_03_T02979_02_003_03_V002.h5"))
level2AM<-getLevel2AM(gedilevel2a)
head(level2AM[,c("beam","shot_number","elev_highestreturn","elev_lowestmode","rh100")])
# Converting shot_number as "integer64" to "character"
level2AM$shot_number<-as.character(level2AM$shot_number)
level2AM_spdf <- sf::st_as_sf(
  level2AM,
  coords = c("lon_lowestmode", "lat_lowestmode"),
  crs = "epsg:4326"
)










level2AM_spdf_3857 <- level2AM_spdf |> st_transform(3857)
sf::write_sf(level2AM_spdf,file.path(outdir,"GEDI02_A_2022198155950_O20359_03_T02979_02_003_03_V002.gpkg"))
sf::write_sf(level2AM_spdf_3857,file.path(outdir,"GEDI02_A_2022198155950_O20359_03_T02979_02_003_03_V002_3857.gpkg"))


library(leaflet)
library(leafsync)

leaflet() %>%
  addCircleMarkers(st_coordinates(level2AM_spdf)[,1],
                   st_coordinates(level2AM_spdf)[,2],
                   radius = 1,
                   opacity = 1,
                   color = "red")  %>%
  addScaleBar(options = list(imperial = FALSE)) %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%
  addLegend(colors = "red", labels= "Samples",title ="GEDI Level1B")

########################################################
# try salpa

dem <- elevatr::get_elev_point(level2AM_spdf, src = "aws",prj=4326)

library(salpa)
library(sf)
library(terra)
library(exactextractr)
library(GA)
library(dtw)
library(tmap)
library(tmaptools)
library(elevatr)
library(raster)
library(ggplot2)
library(dplyr)
library(tidyr)

aligned_points <- linear_alignment(level2AM_spdf, 3857)

plot(level2AM_spdf|> dplyr::select(rh100))


# plot intaractively
# Create the interactive map
interactive_plot_level2AM_spdf <- function(data = level2AM_spdf) {
  library(leaflet)
  
  # Create the interactive map
  map <- leaflet(data = data) %>%
    addProviderTiles(providers$Esri.WorldImagery) %>%
    addCircleMarkers(
      color = "red",
      radius = 3,
      stroke = FALSE,
      fillOpacity = 0.8,
      popup = ~paste("Shot Number:", shot_number, "<br>", "rh100:", rh100)
    ) %>%
    addScaleBar(options = list(imperial = FALSE))
  
  return(map)
}

# Call the function to display the interactive map
interactive_plot_level2AM_spdf()
