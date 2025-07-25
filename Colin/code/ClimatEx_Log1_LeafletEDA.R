
# compare climatological maps from various source using the BC PRISM DEM domain
# Colin Mahony colin.mahony@gov.bc.ca

library(terra)
library(climr)
library(data.table)
library(leaflet)
library(RColorBrewer)
library(ranger)
library(rworldmap)

# function for preparing data
prep <- function(x, studyarea, element){
  x <- crop(x, studyarea)
  studyarea <- project(studyarea, x)
  x <- mask(x, studyarea)
  if(element=="Pr") values(x) <- log2(values(x))
  values(x)[!is.finite(values(x))] <- NA
  return(x)
}

bound <- function(x, breaks){
  values(x)[values(x)>max(breaks)] <- max(breaks)
  values(x)[values(x)<min(breaks)] <- min(breaks)
  return(x)
  }

elements <- c("Tmin", "Tmax", "Pr")
monthcodes <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
monthdays <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 365.25)
month.abb.lowercase <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")

e <- 1
m <- 7

# load the source STATION data for the BC prism
dir <- "//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/"
stn.info <- fread(paste(dir, "Stations/",c("Tmin", "Tmax", "Pr")[e],"_uscdn_8110.csv", sep="")) #read in
for (i in which(names(stn.info)%in%c(month.abb, "Annual"))) stn.info[get(names(stn.info)[i])==c(-9999), (i):=NA, ] # replace -9999 with NA
stn.info <- stn.info[-which(El_Flag=="@"),]
stn <- vect(stn.info, geom = c("Long", "Lat"), crs = "EPSG:4326")
vals <- extract(prism.bc, stn)
stn.info <- stn.info[!is.na(vals[, 2])]
stn.data <- stn.info[,get(month.abb[m])]
stn.data <- if(elements[e]=="Pr") log2(stn.data) else stn.data/10
stn.info <- stn.info[is.finite(stn.data),]
stn.data <- stn.data[is.finite(stn.data)]

#PRISM DEM
dir <- paste("//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/PRISM_dem/", sep="")
dem.bc <- rast(paste(dir, "PRISM_dem.asc", sep=""))

# load the BC PRISM  data for the variable
dir <- paste("//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/", sep="")
file <- list.files(dir, pattern=paste(c("tmin", "tmax", "pr")[e],"_.*._",m, ".tif", sep=""))
prism.bc <- rast(paste(dir, file, sep=""))
if(elements[e]=="Pr") values(prism.bc) <- log2(values(prism.bc))
values(prism.bc)[!is.finite(values(prism.bc))] <- NA

# load the ClimatEx WRF data for the variable
dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/"
wrfclimatex <- rast(paste0(dir, paste("ClimatExWRF_climatology_", c("tmin", "tmax", "pr")[e], "_latlon.tif", sep="")))[[m]]
if(e != 3) wrfclimatex <- wrfclimatex - 273.15
wrfclimatex <- prep(wrfclimatex, studyarea=prism.bc, element=elements[e])


# load the CONUSII WRF data for the variable
dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_CONUSII/"
wrfconus2 <- rast(paste0(dir, paste("conus2_climatology_", c("tmin", "tmax", "pr")[e], "_latlon.tif", sep="")))[[m]]
if(e != 3) wrfconus2 <- wrfconus2 - 273.15
wrfconus2 <- prep(wrfconus2, studyarea=prism.bc, element=elements[e])

# load the USask WRF data for the variable
dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_CCRN/monthly_clim_regridded/"
dem.usask <- rast(paste(dir, "HGT/HGT_regrid.nc", sep=""))
wrfusask <- rast(paste0(dir, paste(c("tmin", "tmax", "prec")[e], monthcodes[m], "regrid.nc", sep="_")))
wrfusask <- wrfusask*monthdays[m]
wrfusask <- prep(wrfusask, studyarea=prism.bc, element=elements[e])

# color scheme
combined <- c(values(prism.bc), values(wrfclimatex), stn.data)
combined <- combined[is.finite(combined)]
inc=diff(range(combined))/500
breaks=seq(quantile(combined, 0.001)-inc, quantile(combined, .999)+inc, inc)
ColScheme <- colorRampPalette(if(elements[e]=="Pr") brewer.pal(9, "YlGnBu") else rev(brewer.pal(11, "RdYlBu")))(length(breaks)-1)
ColPal <- colorBin(ColScheme, bins=breaks, na.color = "white")
ColPal.raster <- colorBin(ColScheme, bins=breaks, na.color = "transparent")

wrfclimatex <- bound(wrfclimatex, breaks=breaks)
wrfconus2 <- bound(wrfconus2, breaks=breaks)
wrfusask <- bound(wrfusask, breaks=breaks)
stn.data[stn.data>max(breaks)] <- max(breaks)
stn.data[stn.data<min(breaks)] <- min(breaks)

# leaflet map
labels <- paste(stn.info$Name, "(El. ", stn.info$Elevation, "m)", sep="")
map <- leaflet(stn.info) %>%
  addTiles(group = "basemap") %>%
  addProviderTiles('Esri.WorldImagery', group = "sat photo") %>%
  addRasterImage(prism.bc, colors = ColPal.raster, opacity = 1, maxBytes = 7 * 1024 * 1024, group = "BC PRISM") %>%
  addRasterImage(wrfclimatex, colors = ColPal.raster, opacity = 1, maxBytes = 7 * 1024 * 1024, group = "WRF (ClimatEx)") %>%
  addRasterImage(wrfusask, colors = ColPal.raster, opacity = 1, maxBytes = 7 * 1024 * 1024, group = "WRF (USask)") %>%
  addRasterImage(wrfconus2, colors = ColPal.raster, opacity = 1, maxBytes = 7 * 1024 * 1024, group = "WRF (CONUS II)") %>%
  addCircleMarkers(lng = ~Long, lat = ~Lat, color="black", fillColor = ~ ColPal(stn.data), opacity = 1, fillOpacity = 1, popup = labels, radius=6, weight=2, group = "Stations") %>%
  addLayersControl(
    baseGroups = c("basemap", "sat photo"),
    overlayGroups = c("WRF (CONUS II)", "WRF (USask)", "WRF (ClimatEx)", "BC PRISM", "Stations"),
    options = layersControlOptions(collapsed = FALSE)
  )
map
