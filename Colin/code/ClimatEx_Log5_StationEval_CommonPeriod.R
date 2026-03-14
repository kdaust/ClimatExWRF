
# quantitative evaluation of WRF climatologies against PRISM stations
# New version using a common period (WY2001-2015) for all datasets
# Colin Mahony colin.mahony@gov.bc.ca


library(rnaturalearth)
library(terra)
library(data.table)
library(sf)
library(scales)
library(RColorBrewer)
library(bcmaps)
library(plotrix) # for Taylor plots

source("Colin/code/util.R")
element.names <- c("Tmin", "Tmax", "Precipitation")

#######################
## Common data
#######################


#PRISM DEM
dir <- paste("C:/Users/CMAHONY/OneDrive - Government of BC/Data/PRISM_BC/PRISM_dem/", sep="")
dem.bc <- rast(paste(dir, "PRISM_dem.asc", sep=""))

slope.bc <- terrain(dem.bc, v = "slope", unit = "radians")   # or unit = "degrees" if preferred
aspect.bc <- terrain(dem.bc, v = "aspect", unit = "radians")
hill.bc <- shade(slope.bc, aspect.bc, angle = 40, direction = 270)

# Create an ocean mask
maskarea <- ext(dem.bc)+c(4,4,4,4)
maskpoly <- as.polygons(maskarea) |> st_as_sf()
st_crs(maskpoly) <- crs(dem.bc)
land <- vect(ne_download(scale = 10, type = "land", category = "physical", returnclass = "sf"))
land <- crop(land, maskarea)
land_union <- st_union(st_as_sf(land))
oceanmask_geom <- st_difference(st_geometry(maskpoly), st_geometry(land_union))
oceanmask <- st_sf(geometry = oceanmask_geom, crs = st_crs(maskpoly))
oceanmask <- vect(oceanmask)

bc <- vect(bc_bound())
bc <- project(bc, dem.bc)
region1 <- ext(bc)

# define the regions
regions <- c("BC", "SW", "SE", "North")
lon1 <- -129.5; lat1 = 48.2
region2 <- ext(c(lon1, lon1+7.5, lat1, lat1+3))
lon1 <- -121.5; lat1 = 48.9
region3 <- ext(c(lon1, lon1+7.5, lat1, lat1+5))
lon1 <- -134; lat1 = 54.5
region4 <- ext(c(lon1, lon1+14, lat1, lat1+5.5))

datasets <- c("prism", "wrfclimatex", "wrfusask", "wrfconus2")
datasets.names <- c("PRISM", "ClimatEx", "USask", "CONUSII")

###########################################################
## plot a key map
#######################

# load the station normals for the WY2000_WY2015 common period
e=1 
m=1
dir <- paste("C:/Users/CMAHONY/OneDrive - Government of BC/Data/bc_stn_dt_final/", sep="")
stn.info <- fread(paste0(dir, paste("bc_obs_stns_climatology_WY2001_WY2015_", c("tmin", "tmax", "pr")[e], ".csv", sep=""))) #read in
stn <- vect(stn.info, geom = c("lon", "lat"), crs = "EPSG:4326")
vals <- extract(prism.bc, stn)
stn.info <- stn.info[!is.na(vals[, 2])]
stn.data <- stn.info[,get(month.abb[m])]
stn.info <- stn.info[is.finite(stn.data),]
stn.data <- stn.data[is.finite(stn.data)]


if(e==1 & m==1){
  png(filename=paste("Colin/results/ClimatExEval.Log5.KeyMap.png",sep="."), type="cairo", units="in", width=6.5, height=5.8, pointsize=10, res=600)
  par(mar=c(0,0,0,0), mfrow=c(1,1))
  legend.args=list(text='Elevation (m)', side=2, font=2, line=0.5, cex=0.8)
  X <- dem.bc
  lim <- quantile(values(X), 0.99)
  values(X)[which(values(X)>lim)] <- lim
  # plot(hill, col=alpha(grey(0:100/100), 1), maxpixels=ncell(hill), legend=F)
  plot(X, col=terrain.colors(99), xaxt="n", yaxt="n")
  plot(crop(hill.bc,X), add=T, col=alpha(grey(0:100/100), 0.5), legend=F, legend.mar=0)
  plot(oceanmask, add=T, col="white", border=F)
  plot(bc, add=T)
  plot(crop(stn, bc), add=T, pch=16, col="gray30", cex=0.8, lwd=0.5)
  # mtext(paste("(a)", sep=""), side=1, line=-1.5, adj=0.005, font=2, cex=0.8)
  
  for(i in 1:length(regions)){
    region <- get(paste("region", i, sep=""))
    plot(region, add=T, lty=2, lwd=2)
    text(ext(region)[1], ext(region)[4]-0.25, regions[i], font=2, pos=4, offset=0.1)
  }
  l <- ext(X)
  rect(l[1], l[3], l[2], l[4], border = "black", lwd = 1)
  
  dev.off()
}



#######################
## 2 months for 3 elements
#######################

months <- c(1,7)
regions.select <- c(1,4)

png(filename=paste("Colin/results/ClimatExEval.Log5.TaylorPlots", paste0(month.abb[months], collapse=""), "png",sep="."), type="cairo", units="in", width=9, height=6.25, pointsize=10, res=600)
par(mfrow=c(2,3), mar=c(0,0,0,0), mgp=c(2,0.25, 0), tck=-0.01)
m=1
for(m in months){
  monthcode = monthcodes[m]
  
  e=1
  for(e in 1:3){
    element = elements[e]
    
    # load the BC PRISM  data for the variable
    dir <- paste("C:/Users/CMAHONY/OneDrive - Government of BC/Data/PRISM_BC/", sep="")
    file <- list.files(dir, pattern=paste(c("tmin", "tmax", "pr")[e],"_.*._",m, ".tif", sep=""))
    prism.bc <- rast(paste(dir, file, sep=""))
    
    # load the source STATION data for the BC prism
    dir <- paste("C:/Users/CMAHONY/OneDrive - Government of BC/Data/PRISM_BC/", sep="")
    stn.prism.info <- fread(paste(dir, "Stations/",c("Tmin", "Tmax", "Pr")[e],"_uscdn_8110.csv", sep="")) #read in
    for (i in which(names(stn.prism.info)%in%c(month.abb, "Annual"))) stn.prism.info[get(names(stn.prism.info)[i])==c(-9999), (i):=NA, ] # replace -9999 with NA
    stn.prism.info <- stn.prism.info[-which(El_Flag=="@"),]
    stn <- vect(stn.prism.info, geom = c("Long", "Lat"), crs = "EPSG:4326")
    vals <- extract(prism.bc, stn)
    stn.prism.info <- stn.prism.info[!is.na(vals[, 2])]
    stn.prism.data <- stn.prism.info[,get(month.abb[m])]
    stn.prism.data <- if(e==3) log2(stn.prism.data) else stn.prism.data/10
    stn.prism.info <- stn.prism.info[is.finite(stn.prism.data),]
    stn.prism.data <- stn.prism.data[is.finite(stn.prism.data)]
    
    # load the station normals for the WY2000_WY2015 common period
    dir <- paste("C:/Users/CMAHONY/OneDrive - Government of BC/Data/bc_stn_dt_final/", sep="")
    stn.info <- fread(paste0(dir, paste("bc_obs_stns_climatology_WY2001_WY2015_", c("tmin", "tmax", "pr")[e], ".csv", sep=""))) #read in
    stn <- vect(stn.info, geom = c("lon", "lat"), crs = "EPSG:4326")
    vals <- extract(prism.bc, stn)
    stn.info <- stn.info[!is.na(vals[, 2])]
    stn.data <- stn.info[,get(month.abb[m])]
    if(e==3) stn.data <- log2(stn.data)
    stn.info <- stn.info[is.finite(stn.data),]
    stn.data <- stn.data[is.finite(stn.data)]
    
    # load the ClimatEx WRF data for the variable
    dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/"
    wrfclimatex.bc <- rast(paste0(dir, paste("ClimatExWRF_climatology_WY2001_WY2015_", c("tmin", "tmax", "pr")[e], "_latlon.tif", sep="")))[[m]]
    if(e != 3) wrfclimatex.bc <- wrfclimatex.bc - 273.15
    wrfclimatex.bc <- crop(wrfclimatex.bc, prism.bc)
    dem.wrfclimatex <- rast(paste0(dir, "HGT_latlon.nc"))
    dem.wrfclimatex <- crop(dem.wrfclimatex, prism.bc)
    
    # load the USask WRF data for the variable
    dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_CCRN/monthly_clim_regridded/"
    dem.usask <- rast(paste(dir, "HGT/HGT_regrid.nc", sep=""))
    wrfusask.bc <- rast(paste0(dir, paste(c("tmin", "tmax", "prec")[e], monthcodes[m], "regrid.nc", sep="_")))
    wrfusask.bc <- crop(wrfusask.bc, prism.bc)
    if(e == 3) wrfusask.bc <- wrfusask.bc*monthdays[m]
    
    # load the conus2 WRF data for the variable
    dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_CONUSII/"
    # dem.conus2 <- rast(paste(dir, "HGT/HGT_regrid.nc", sep=""))
    wrfconus2.bc <- rast(paste0(dir, paste("conus2_climatology_WY2001_WY2015_", element, "_latlon.tif", sep="")))[[m]]
    if(e != 3) wrfconus2.bc <- wrfconus2.bc - 273.15
    wrfconus2.bc <- crop(wrfconus2.bc, prism.bc)
    
    #################################
    ## Taylor plots
    #################################
    
    # Plot PRISM against observations
    taylor.diagram(ref = 0:5, model = 0:5, main = paste(month.name[m], element.names[e]),
                   col = "black", pch = 1, cex = 1, normalize = TRUE, xlab="",
                   sd.arcs = TRUE, grad.corr.lines = TRUE, pos.cor = TRUE)
    
    colScheme <- c("black", "dodgerblue", "yellow", "red")
    
    # Add a legend
    legend("topright",
           legend = datasets.names,
           title="Datasets",
           fill = colScheme,      # colored boxes
           border = NULL,
           bty = "n",
           inset = c(0,-0.07))
    
    # second column: black shapes
    legend("topright",
           legend = regions[regions.select],
           title="Regions",
           bg = "white",
           pch = c(22,21,24,25)[regions.select],  # shapes
           pt.cex = 1.5,
           bty = "n",
           inset = c(0.3,-0.07))   # shift to the left so they appear in second column
    
    i=4
    for(i in regions.select){
      studyarea <- get(paste("region", i, sep=""))  
      caseStudy <- regions[i]
      
      ## DEM
      
      dem <- crop(dem.bc, studyarea)
      prism <- crop(prism.bc, studyarea)
      wrfclimatex<- project(wrfclimatex.bc, dem)
      wrfusask<- project(wrfusask.bc, dem)
      wrfconus2<- project(wrfconus2.bc, dem)
      dem.wrf <- project(dem.wrfclimatex, dem)
      
      stn.prism.vect <- vect(stn.prism.info, geom = c("Long", "Lat"), crs = "EPSG:4326")
      stn.prism.crop <- crop(stn.prism.vect, studyarea)
      stn.prism.values <- as.data.frame(stn.prism.crop)[,which(names(stn.prism.crop)==month.abb[m])]
      stn.prism.values <- if(e==3) log2(stn.prism.values) else stn.prism.values/10
      
      stn.vect <- vect(stn.info, geom = c("lon", "lat"), crs = "EPSG:4326")
      stn.crop <- crop(stn.vect, studyarea)
      stn.values <- as.data.frame(stn.crop)[,which(names(stn.crop)==month.abb[m])]
      stn.values <- if(e==3) log2(stn.values) else stn.values
      
      for(dataset in datasets){
        d=which(datasets==dataset)
        temp <- get(dataset)
        # plot(temp)
        # plot(stn.crop, add=T)
        stn.temp <- as.vector(unlist(extract(temp, if(dataset=="prism") stn.prism.crop else stn.crop)[2]))
        stn.temp <- if(e==3) log2(stn.temp) else stn.temp
        # assign(paste("stn", dataset, elements[e], monthcodes[m], sep="."), stn.temp)
        # assign(paste("error", dataset, elements[e], monthcodes[m], sep="."), stn.temp - stn.values)
        # Add points on the Taylor diagram
        # taylor.diagram(ref = stn.values, model = stn.temp, add = TRUE,
        #                col = colScheme[d], pch = c(16,15,17)[i], cex = 1.2, normalize = TRUE, pcex=1.5)
        taylor.diagram.filled(ref = if(dataset=="prism") stn.prism.values else stn.values, model = stn.temp, add = TRUE,
                              bg = colScheme[d], pch = c(22,21,24,25)[i], cex = 2, normalize = TRUE)
      }
      
      print(paste("region", i))
    }
    print(element)
  }
  print(month.abb[m])
}
dev.off()



#######################
## 4 months for 3 elements
#######################

months <- c(1,4,7,10)

png(filename=paste("Colin/results/ClimatExEval.Log5.TaylorPlots", "png",sep="."), type="cairo", units="in", width=9, height=12, pointsize=10, res=600)
par(mfrow=c(4,3), mar=c(0,0,0,0), mgp=c(2,0.25, 0), tck=-0.01)
m=1
for(m in months){
  monthcode = monthcodes[m]
  
  e=1
  for(e in 1:3){
    element = elements[e]
    
    # load the BC PRISM  data for the variable
    dir <- paste("C:/Users/CMAHONY/OneDrive - Government of BC/Data/PRISM_BC/", sep="")
    file <- list.files(dir, pattern=paste(c("tmin", "tmax", "pr")[e],"_.*._",m, ".tif", sep=""))
    prism.bc <- rast(paste(dir, file, sep=""))
    
    # load the source STATION data for the BC prism
    dir <- paste("C:/Users/CMAHONY/OneDrive - Government of BC/Data/PRISM_BC/", sep="")
    stn.prism.info <- fread(paste(dir, "Stations/",c("Tmin", "Tmax", "Pr")[e],"_uscdn_8110.csv", sep="")) #read in
    for (i in which(names(stn.prism.info)%in%c(month.abb, "Annual"))) stn.prism.info[get(names(stn.prism.info)[i])==c(-9999), (i):=NA, ] # replace -9999 with NA
    stn.prism.info <- stn.prism.info[-which(El_Flag=="@"),]
    stn <- vect(stn.prism.info, geom = c("Long", "Lat"), crs = "EPSG:4326")
    vals <- extract(prism.bc, stn)
    stn.prism.info <- stn.prism.info[!is.na(vals[, 2])]
    stn.prism.data <- stn.prism.info[,get(month.abb[m])]
    stn.prism.data <- if(e==3) log2(stn.prism.data) else stn.prism.data/10
    stn.prism.info <- stn.prism.info[is.finite(stn.prism.data),]
    stn.prism.data <- stn.prism.data[is.finite(stn.prism.data)]
    
    # load the station normals for the WY2000_WY2015 common period
    dir <- paste("C:/Users/CMAHONY/OneDrive - Government of BC/Data/bc_stn_dt_final/", sep="")
    stn.info <- fread(paste0(dir, paste("bc_obs_stns_climatology_WY2001_WY2015_", c("tmin", "tmax", "pr")[e], ".csv", sep=""))) #read in
    stn <- vect(stn.info, geom = c("lon", "lat"), crs = "EPSG:4326")
    vals <- extract(prism.bc, stn)
    stn.info <- stn.info[!is.na(vals[, 2])]
    stn.data <- stn.info[,get(month.abb[m])]
    if(e==3) stn.data <- log2(stn.data)
    stn.info <- stn.info[is.finite(stn.data),]
    stn.data <- stn.data[is.finite(stn.data)]
    
    # load the ClimatEx WRF data for the variable
    dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/"
    wrfclimatex.bc <- rast(paste0(dir, paste("ClimatExWRF_climatology_WY2001_WY2015_", c("tmin", "tmax", "pr")[e], "_latlon.tif", sep="")))[[m]]
    if(e != 3) wrfclimatex.bc <- wrfclimatex.bc - 273.15
    wrfclimatex.bc <- crop(wrfclimatex.bc, prism.bc)
    dem.wrfclimatex <- rast(paste0(dir, "HGT_latlon.nc"))
    dem.wrfclimatex <- crop(dem.wrfclimatex, prism.bc)
    
    # load the USask WRF data for the variable
    dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_CCRN/monthly_clim_regridded/"
    dem.usask <- rast(paste(dir, "HGT/HGT_regrid.nc", sep=""))
    wrfusask.bc <- rast(paste0(dir, paste(c("tmin", "tmax", "prec")[e], monthcodes[m], "regrid.nc", sep="_")))
    wrfusask.bc <- crop(wrfusask.bc, prism.bc)
    if(e == 3) wrfusask.bc <- wrfusask.bc*monthdays[m]
    
    # load the conus2 WRF data for the variable
    dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_CONUSII/"
    # dem.conus2 <- rast(paste(dir, "HGT/HGT_regrid.nc", sep=""))
    wrfconus2.bc <- rast(paste0(dir, paste("conus2_climatology_WY2001_WY2015_", element, "_latlon.tif", sep="")))[[m]]
    if(e != 3) wrfconus2.bc <- wrfconus2.bc - 273.15
    wrfconus2.bc <- crop(wrfconus2.bc, prism.bc)
    
    
    #################################
    ## Taylor plots
    #################################
    
    # Plot PRISM against observations
    taylor.diagram(ref = 0:5, model = 0:5, main = paste(month.name[m], element.names[e]),
                   col = "black", pch = 1, cex = 1, normalize = TRUE, xlab="",
                   sd.arcs = TRUE, grad.corr.lines = TRUE, pos.cor = TRUE)
    
    colScheme <- c("black", "dodgerblue", "yellow", "red")
    
    # Add a legend
    legend("topright",
           legend = datasets.names,
           title="Datasets",
           fill = colScheme,      # colored boxes
           border = NULL,
           bty = "n",
           inset = c(0,-0.07))
    
    # second column: black shapes
    legend("topright",
           legend = regions[regions.select],
           title="Regions",
           bg = "white",
           pch = c(22,21,24,25)[regions.select],  # shapes
           pt.cex = 1.5,
           bty = "n",
           inset = c(0.3,-0.07))   # shift to the left so they appear in second column
    
    i=4
    for(i in regions.select){
      studyarea <- get(paste("region", i, sep=""))  
      caseStudy <- regions[i]
      
      ## DEM
      
      dem <- crop(dem.bc, studyarea)
      prism <- crop(prism.bc, studyarea)
      wrfclimatex<- project(wrfclimatex.bc, dem)
      wrfusask<- project(wrfusask.bc, dem)
      wrfconus2<- project(wrfconus2.bc, dem)
      dem.wrf <- project(dem.wrfclimatex, dem)
      
      stn.prism.vect <- vect(stn.prism.info, geom = c("Long", "Lat"), crs = "EPSG:4326")
      stn.prism.crop <- crop(stn.prism.vect, studyarea)
      stn.prism.values <- as.data.frame(stn.prism.crop)[,which(names(stn.prism.crop)==month.abb[m])]
      stn.prism.values <- if(e==3) log2(stn.prism.values) else stn.prism.values/10
      
      stn.vect <- vect(stn.info, geom = c("lon", "lat"), crs = "EPSG:4326")
      stn.crop <- crop(stn.vect, studyarea)
      stn.values <- as.data.frame(stn.crop)[,which(names(stn.crop)==month.abb[m])]
      stn.values <- if(e==3) log2(stn.values) else stn.values
      
      for(dataset in datasets){
        d=which(datasets==dataset)
        temp <- get(dataset)
        # plot(temp)
        # plot(stn.crop, add=T)
        stn.temp <- as.vector(unlist(extract(temp, if(dataset=="prism") stn.prism.crop else stn.crop)[2]))
        stn.temp <- if(e==3) log2(stn.temp) else stn.temp
        # assign(paste("stn", dataset, elements[e], monthcodes[m], sep="."), stn.temp)
        # assign(paste("error", dataset, elements[e], monthcodes[m], sep="."), stn.temp - stn.values)
        # Add points on the Taylor diagram
        # taylor.diagram(ref = stn.values, model = stn.temp, add = TRUE,
        #                col = colScheme[d], pch = c(16,15,17)[i], cex = 1.2, normalize = TRUE, pcex=1.5)
        taylor.diagram.filled(ref = if(dataset=="prism") stn.prism.values else stn.values, model = stn.temp, add = TRUE,
                              bg = colScheme[d], pch = c(22,21,24,25)[i], cex = 2, normalize = TRUE)
      }
      
      # print(paste("region", i))
    }
    print(element)
  }
  print(month.abb[m])
}
dev.off()



#######################
## 12 months for 1 element
#######################

e=1
for(e in 1:3){
  element = elements[e]
  
  png(filename=paste("Colin/results/ClimatExEval.Log5.TaylorPlots.12months", element, "png",sep="."), type="cairo", units="in", width=9, height=12, pointsize=10, res=600)
  par(mfrow=c(4,3), mar=c(0,0,0,0), mgp=c(2,0.25, 0), tck=-0.01)
  m=1
  for(m in 1:12){
    monthcode = monthcodes[m]
    
    # load the BC PRISM  data for the variable
    dir <- paste("C:/Users/CMAHONY/OneDrive - Government of BC/Data/PRISM_BC/", sep="")
    file <- list.files(dir, pattern=paste(c("tmin", "tmax", "pr")[e],"_.*._",m, ".tif", sep=""))
    prism.bc <- rast(paste(dir, file, sep=""))
    
    # load the source STATION data for the BC prism
    dir <- paste("C:/Users/CMAHONY/OneDrive - Government of BC/Data/PRISM_BC/", sep="")
    stn.prism.info <- fread(paste(dir, "Stations/",c("Tmin", "Tmax", "Pr")[e],"_uscdn_8110.csv", sep="")) #read in
    for (i in which(names(stn.prism.info)%in%c(month.abb, "Annual"))) stn.prism.info[get(names(stn.prism.info)[i])==c(-9999), (i):=NA, ] # replace -9999 with NA
    stn.prism.info <- stn.prism.info[-which(El_Flag=="@"),]
    stn <- vect(stn.prism.info, geom = c("Long", "Lat"), crs = "EPSG:4326")
    vals <- extract(prism.bc, stn)
    stn.prism.info <- stn.prism.info[!is.na(vals[, 2])]
    stn.prism.data <- stn.prism.info[,get(month.abb[m])]
    stn.prism.data <- if(e==3) log2(stn.prism.data) else stn.prism.data/10
    stn.prism.info <- stn.prism.info[is.finite(stn.prism.data),]
    stn.prism.data <- stn.prism.data[is.finite(stn.prism.data)]
    
    # load the station normals for the WY2000_WY2015 common period
    dir <- paste("C:/Users/CMAHONY/OneDrive - Government of BC/Data/bc_stn_dt_final/", sep="")
    stn.info <- fread(paste0(dir, paste("bc_obs_stns_climatology_WY2001_WY2015_", c("tmin", "tmax", "pr")[e], ".csv", sep=""))) #read in
    stn <- vect(stn.info, geom = c("lon", "lat"), crs = "EPSG:4326")
    vals <- extract(prism.bc, stn)
    stn.info <- stn.info[!is.na(vals[, 2])]
    stn.data <- stn.info[,get(month.abb[m])]
    if(e==3) stn.data <- log2(stn.data)
    stn.info <- stn.info[is.finite(stn.data),]
    stn.data <- stn.data[is.finite(stn.data)]
    
    # load the ClimatEx WRF data for the variable
    dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/"
    wrfclimatex.bc <- rast(paste0(dir, paste("ClimatExWRF_climatology_WY2001_WY2015_", c("tmin", "tmax", "pr")[e], "_latlon.tif", sep="")))[[m]]
    if(e != 3) wrfclimatex.bc <- wrfclimatex.bc - 273.15
    wrfclimatex.bc <- crop(wrfclimatex.bc, prism.bc)
    dem.wrfclimatex <- rast(paste0(dir, "HGT_latlon.nc"))
    dem.wrfclimatex <- crop(dem.wrfclimatex, prism.bc)
    
    # load the USask WRF data for the variable
    dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_CCRN/monthly_clim_regridded/"
    dem.usask <- rast(paste(dir, "HGT/HGT_regrid.nc", sep=""))
    wrfusask.bc <- rast(paste0(dir, paste(c("tmin", "tmax", "prec")[e], monthcodes[m], "regrid.nc", sep="_")))
    wrfusask.bc <- crop(wrfusask.bc, prism.bc)
    if(e == 3) wrfusask.bc <- wrfusask.bc*monthdays[m]
    
    # load the conus2 WRF data for the variable
    dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_CONUSII/"
    # dem.conus2 <- rast(paste(dir, "HGT/HGT_regrid.nc", sep=""))
    wrfconus2.bc <- rast(paste0(dir, paste("conus2_climatology_WY2001_WY2015_", element, "_latlon.tif", sep="")))[[m]]
    if(e != 3) wrfconus2.bc <- wrfconus2.bc - 273.15
    wrfconus2.bc <- crop(wrfconus2.bc, prism.bc)
    
    
    #################################
    ## Taylor plots
    #################################
    
    # Plot PRISM against observations
    taylor.diagram(ref = 0:5, model = 0:5, main = paste(month.name[m], element.names[e]),
                   col = "black", pch = 1, cex = 1, normalize = TRUE, xlab="",
                   sd.arcs = TRUE, grad.corr.lines = TRUE, pos.cor = TRUE)
    
    colScheme <- c("black", "dodgerblue", "yellow", "red")
    
    # Add a legend
    legend("topright",
           legend = datasets.names,
           title="Datasets",
           fill = colScheme,      # colored boxes
           border = NULL,
           bty = "n",
           inset = c(0,-0.07))
    
    # second column: black shapes
    legend("topright",
           legend = regions[regions.select],
           title="Regions",
           bg = "white",
           pch = c(22,21,24,25)[regions.select],  # shapes
           pt.cex = 1.5,
           bty = "n",
           inset = c(0.3,-0.07))   # shift to the left so they appear in second column
    
    i=4
    for(i in regions.select){
      studyarea <- get(paste("region", i, sep=""))  
      caseStudy <- regions[i]
      
      ## DEM
      
      dem <- crop(dem.bc, studyarea)
      prism <- crop(prism.bc, studyarea)
      wrfclimatex<- project(wrfclimatex.bc, dem)
      wrfusask<- project(wrfusask.bc, dem)
      wrfconus2<- project(wrfconus2.bc, dem)
      dem.wrf <- project(dem.wrfclimatex, dem)
      
      stn.prism.vect <- vect(stn.prism.info, geom = c("Long", "Lat"), crs = "EPSG:4326")
      stn.prism.crop <- crop(stn.prism.vect, studyarea)
      stn.prism.values <- as.data.frame(stn.prism.crop)[,which(names(stn.prism.crop)==month.abb[m])]
      stn.prism.values <- if(e==3) log2(stn.prism.values) else stn.prism.values/10
      
      stn.vect <- vect(stn.info, geom = c("lon", "lat"), crs = "EPSG:4326")
      stn.crop <- crop(stn.vect, studyarea)
      stn.values <- as.data.frame(stn.crop)[,which(names(stn.crop)==month.abb[m])]
      stn.values <- if(e==3) log2(stn.values) else stn.values
      
      for(dataset in datasets){
        d=which(datasets==dataset)
        temp <- get(dataset)
        # plot(temp)
        # plot(stn.crop, add=T)
        stn.temp <- as.vector(unlist(extract(temp, if(dataset=="prism") stn.prism.crop else stn.crop)[2]))
        stn.temp <- if(e==3) log2(stn.temp) else stn.temp
        # assign(paste("stn", dataset, elements[e], monthcodes[m], sep="."), stn.temp)
        # assign(paste("error", dataset, elements[e], monthcodes[m], sep="."), stn.temp - stn.values)
        # Add points on the Taylor diagram
        # taylor.diagram(ref = stn.values, model = stn.temp, add = TRUE,
        #                col = colScheme[d], pch = c(16,15,17)[i], cex = 1.2, normalize = TRUE, pcex=1.5)
        taylor.diagram.filled(ref = if(dataset=="prism") stn.prism.values else stn.values, model = stn.temp, add = TRUE,
                              bg = colScheme[d], pch = c(22,21,24,25)[i], cex = 2, normalize = TRUE)
      }
      
      # print(paste("region", i))
    }
    print(month.abb[m])
  }
  dev.off()
  print(element)
}

#######################
## 6 months for 1 element
#######################

months <- seq(1,12, 2)

e=1
for(e in 1:3){
  element = elements[e]
  
  png(filename=paste("Colin/results/ClimatExEval.Log5.TaylorPlots", element, "png",sep="."), type="cairo", units="in", width=9, height=6.25, pointsize=10, res=600)
  par(mfrow=c(2,3), mar=c(0,0,0,0), mgp=c(2,0.25, 0), tck=-0.01)
  m=1
  for(m in months){
    monthcode = monthcodes[m]
    
    # load the BC PRISM  data for the variable
    dir <- paste("C:/Users/CMAHONY/OneDrive - Government of BC/Data/PRISM_BC/", sep="")
    file <- list.files(dir, pattern=paste(c("tmin", "tmax", "pr")[e],"_.*._",m, ".tif", sep=""))
    prism.bc <- rast(paste(dir, file, sep=""))
    
    # load the source STATION data for the BC prism
    dir <- paste("C:/Users/CMAHONY/OneDrive - Government of BC/Data/PRISM_BC/", sep="")
    stn.prism.info <- fread(paste(dir, "Stations/",c("Tmin", "Tmax", "Pr")[e],"_uscdn_8110.csv", sep="")) #read in
    for (i in which(names(stn.prism.info)%in%c(month.abb, "Annual"))) stn.prism.info[get(names(stn.prism.info)[i])==c(-9999), (i):=NA, ] # replace -9999 with NA
    stn.prism.info <- stn.prism.info[-which(El_Flag=="@"),]
    stn <- vect(stn.prism.info, geom = c("Long", "Lat"), crs = "EPSG:4326")
    vals <- extract(prism.bc, stn)
    stn.prism.info <- stn.prism.info[!is.na(vals[, 2])]
    stn.prism.data <- stn.prism.info[,get(month.abb[m])]
    stn.prism.data <- if(e==3) log2(stn.prism.data) else stn.prism.data/10
    stn.prism.info <- stn.prism.info[is.finite(stn.prism.data),]
    stn.prism.data <- stn.prism.data[is.finite(stn.prism.data)]
    
    # load the station normals for the WY2000_WY2015 common period
    dir <- paste("C:/Users/CMAHONY/OneDrive - Government of BC/Data/bc_stn_dt_final/", sep="")
    stn.info <- fread(paste0(dir, paste("bc_obs_stns_climatology_WY2001_WY2015_", c("tmin", "tmax", "pr")[e], ".csv", sep=""))) #read in
    stn <- vect(stn.info, geom = c("lon", "lat"), crs = "EPSG:4326")
    vals <- extract(prism.bc, stn)
    stn.info <- stn.info[!is.na(vals[, 2])]
    stn.data <- stn.info[,get(month.abb[m])]
    if(e==3) stn.data <- log2(stn.data)
    stn.info <- stn.info[is.finite(stn.data),]
    stn.data <- stn.data[is.finite(stn.data)]
    
    # load the ClimatEx WRF data for the variable
    dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/"
    wrfclimatex.bc <- rast(paste0(dir, paste("ClimatExWRF_climatology_WY2001_WY2015_", c("tmin", "tmax", "pr")[e], "_latlon.tif", sep="")))[[m]]
    if(e != 3) wrfclimatex.bc <- wrfclimatex.bc - 273.15
    wrfclimatex.bc <- crop(wrfclimatex.bc, prism.bc)
    dem.wrfclimatex <- rast(paste0(dir, "HGT_latlon.nc"))
    dem.wrfclimatex <- crop(dem.wrfclimatex, prism.bc)
    
    # load the USask WRF data for the variable
    dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_CCRN/monthly_clim_regridded/"
    dem.usask <- rast(paste(dir, "HGT/HGT_regrid.nc", sep=""))
    wrfusask.bc <- rast(paste0(dir, paste(c("tmin", "tmax", "prec")[e], monthcodes[m], "regrid.nc", sep="_")))
    wrfusask.bc <- crop(wrfusask.bc, prism.bc)
    if(e == 3) wrfusask.bc <- wrfusask.bc*monthdays[m]
    
    # load the conus2 WRF data for the variable
    dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_CONUSII/"
    # dem.conus2 <- rast(paste(dir, "HGT/HGT_regrid.nc", sep=""))
    wrfconus2.bc <- rast(paste0(dir, paste("conus2_climatology_WY2001_WY2015_", element, "_latlon.tif", sep="")))[[m]]
    if(e != 3) wrfconus2.bc <- wrfconus2.bc - 273.15
    wrfconus2.bc <- crop(wrfconus2.bc, prism.bc)
    
    
    #################################
    ## Taylor plots
    #################################
    
    # Plot PRISM against observations
    taylor.diagram(ref = 0:5, model = 0:5, main = paste(month.name[m], element.names[e]),
                   col = "black", pch = 1, cex = 1, normalize = TRUE, xlab="",
                   sd.arcs = TRUE, grad.corr.lines = TRUE, pos.cor = TRUE)
    
    colScheme <- c("black", "dodgerblue", "yellow", "red")
    
    # Add a legend
    legend("topright",
           legend = datasets.names,
           title="Datasets",
           fill = colScheme,      # colored boxes
           border = NULL,
           bty = "n",
           inset = c(0,-0.07))
    
    # second column: black shapes
    legend("topright",
           legend = regions[regions.select],
           title="Regions",
           bg = "white",
           pch = c(22,21,24,25)[regions.select],  # shapes
           pt.cex = 1.5,
           bty = "n",
           inset = c(0.3,-0.07))   # shift to the left so they appear in second column
    
    i=4
    for(i in regions.select){
      studyarea <- get(paste("region", i, sep=""))  
      caseStudy <- regions[i]
      
      ## DEM
      
      dem <- crop(dem.bc, studyarea)
      prism <- crop(prism.bc, studyarea)
      wrfclimatex<- project(wrfclimatex.bc, dem)
      wrfusask<- project(wrfusask.bc, dem)
      wrfconus2<- project(wrfconus2.bc, dem)
      dem.wrf <- project(dem.wrfclimatex, dem)
      
      stn.prism.vect <- vect(stn.prism.info, geom = c("Long", "Lat"), crs = "EPSG:4326")
      stn.prism.crop <- crop(stn.prism.vect, studyarea)
      stn.prism.values <- as.data.frame(stn.prism.crop)[,which(names(stn.prism.crop)==month.abb[m])]
      stn.prism.values <- if(e==3) log2(stn.prism.values) else stn.prism.values/10
      
      stn.vect <- vect(stn.info, geom = c("lon", "lat"), crs = "EPSG:4326")
      stn.crop <- crop(stn.vect, studyarea)
      stn.values <- as.data.frame(stn.crop)[,which(names(stn.crop)==month.abb[m])]
      stn.values <- if(e==3) log2(stn.values) else stn.values
      
      for(dataset in datasets){
        d=which(datasets==dataset)
        temp <- get(dataset)
        # plot(temp)
        # plot(stn.crop, add=T)
        stn.temp <- as.vector(unlist(extract(temp, if(dataset=="prism") stn.prism.crop else stn.crop)[2]))
        stn.temp <- if(e==3) log2(stn.temp) else stn.temp
        # assign(paste("stn", dataset, elements[e], monthcodes[m], sep="."), stn.temp)
        # assign(paste("error", dataset, elements[e], monthcodes[m], sep="."), stn.temp - stn.values)
        # Add points on the Taylor diagram
        # taylor.diagram(ref = stn.values, model = stn.temp, add = TRUE,
        #                col = colScheme[d], pch = c(16,15,17)[i], cex = 1.2, normalize = TRUE, pcex=1.5)
        taylor.diagram.filled(ref = if(dataset=="prism") stn.prism.values else stn.values, model = stn.temp, add = TRUE,
                              bg = colScheme[d], pch = c(22,21,24,25)[i], cex = 2, normalize = TRUE)
      }
      
      # print(paste("region", i))
    }
    print(month.abb[m])
  }
  dev.off()
  
  print(element)
}


#-----------------------------------------
## boxplots of climatological mean bias in all months
#-----------------------------------------

i=1
for(i in 1:length(regions)){
  
  region <- regions[i]
  studyarea <- get(paste("region", i, sep=""))  
  caseStudy <- regions[i]
  
  
  png(filename=paste("Colin/results/ClimatExEval.Log5.Boxplots_bias_", caseStudy, ".png",sep=""), type="cairo", units="in", width=6.5, height=7.5, pointsize=12, res=300)
  
  par(mfrow=c(3,1))
  
  e=1
  for (e in 1:3) {
    
    # storage vectors
    values <- c()
    groups <- c()
    
    m = 1
    for (m in 1:12) {
      
      # load the BC PRISM  data for the variable
      dir <- paste("C:/Users/CMAHONY/OneDrive - Government of BC/Data/PRISM_BC/", sep="")
      file <- list.files(dir, pattern=paste(c("tmin", "tmax", "pr")[e],"_.*._",m, ".tif", sep=""))
      prism.bc <- rast(paste(dir, file, sep=""))
      
      # load the station normals for the WY2000_WY2015 common period
      dir <- paste("C:/Users/CMAHONY/OneDrive - Government of BC/Data/bc_stn_dt_final/", sep="")
      stn.info <- fread(paste0(dir, paste("bc_obs_stns_climatology_WY2001_WY2015_", c("tmin", "tmax", "pr")[e], ".csv", sep=""))) #read in
      stn <- vect(stn.info, geom = c("lon", "lat"), crs = "EPSG:4326")
      vals <- extract(prism.bc, stn)
      stn.info <- stn.info[!is.na(vals[, 2])]
      stn.data <- stn.info[,get(month.abb[m])]
      if(e==3) stn.data <- log2(stn.data)
      stn.info <- stn.info[is.finite(stn.data),]
      stn.data <- stn.data[is.finite(stn.data)]
      
      # load the ClimatEx WRF data for the variable
      dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/"
      wrfclimatex.bc <- rast(paste0(dir, paste("ClimatExWRF_climatology_WY2001_WY2015_", c("tmin", "tmax", "pr")[e], "_latlon.tif", sep="")))[[m]]
      if(e != 3) wrfclimatex.bc <- wrfclimatex.bc - 273.15
      wrfclimatex.bc <- crop(wrfclimatex.bc, prism.bc)
      dem.wrfclimatex <- rast(paste0(dir, "HGT_latlon.nc"))
      dem.wrfclimatex <- crop(dem.wrfclimatex, prism.bc)
      
      # load the USask WRF data for the variable
      dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_CCRN/monthly_clim_regridded/"
      dem.usask <- rast(paste(dir, "HGT/HGT_regrid.nc", sep=""))
      wrfusask.bc <- rast(paste0(dir, paste(c("tmin", "tmax", "prec")[e], monthcodes[m], "regrid.nc", sep="_")))
      wrfusask.bc <- crop(wrfusask.bc, prism.bc)
      if(e == 3) wrfusask.bc <- wrfusask.bc*monthdays[m]
      
      # load the conus2 WRF data for the variable
      dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_CONUSII/"
      # dem.conus2 <- rast(paste(dir, "HGT/HGT_regrid.nc", sep=""))
      wrfconus2.bc <- rast(paste0(dir, paste("conus2_climatology_WY2001_WY2015_", c("tmin", "tmax", "pr")[e], "_latlon.tif", sep="")))[[m]]
      if(e != 3) wrfconus2.bc <- wrfconus2.bc - 273.15
      wrfconus2.bc <- crop(wrfconus2.bc, prism.bc)
      
      #regrid
      dem <- crop(dem.bc, studyarea)
      prism <- crop(prism.bc, studyarea)
      wrfclimatex<- project(wrfclimatex.bc, dem)
      wrfusask<- project(wrfusask.bc, dem)
      wrfconus2<- project(wrfconus2.bc, dem)
      dem.wrf <- project(dem.wrfclimatex, dem)
      
      # stations for selected region
      stn.vect <- vect(stn.info, geom = c("lon", "lat"), crs = "EPSG:4326")
      stn.crop <- crop(stn.vect, studyarea)
      stn.values <- as.data.frame(stn.crop)[,which(names(stn.crop)==month.abb[m])]
      stn.values <- if(e==3) log2(stn.values) else stn.values
      
      # extract from spatial data
      wrfclimatex.ptData <- as.vector(extract(wrfclimatex, stn.crop)[,2])
      wrfusask.ptData <- as.vector(extract(wrfusask, stn.crop)[,2])
      wrfconus2.ptData <- as.vector(extract(wrfconus2, stn.crop)[,2])
      
      if(e==3){
        wrfclimatex.ptData <- log2(wrfclimatex.ptData)
        wrfusask.ptData <- log2(wrfusask.ptData)
        wrfconus2.ptData <- log2(wrfconus2.ptData)
      }
      
      wrfclimatex.error <- wrfclimatex.ptData - stn.values
      wrfusask.error <- wrfusask.ptData - stn.values
      wrfconus2.error <- wrfconus2.ptData - stn.values
      
      # accumulate values
      values <- c(values, wrfclimatex.error, wrfusask.error, wrfconus2.error)
      
      # accumulate groups
      groups <- c(
        groups,
        rep(paste0("wrfclimatex-", m), length(wrfclimatex.ptData)),
        rep(paste0("wrfusask-", m), length(wrfusask.ptData)),
        rep(paste0("wrfconus2-", m), length(wrfconus2.ptData))
      )
      
      print(m)  
    }
    
    # duplicate values for annual boxplot
    values <- c(values, values)
    groups <- c(groups, sapply(strsplit(groups, "-"), `[`, 1))
    
    
    # ensure correct order
    groups <- factor(groups, levels=unique(groups))
    
    # 3 datasets per month
    n_per_month <- 3
    
    # position calculations
    month_indices <- rep(1:13, each=n_per_month)
    spacing <- 1 # you can adjust this value
    at_positions <- (month_indices - 1) * (n_per_month + spacing) + rep(1:n_per_month, times=13)
    
    # plot
    par(mar=c(c(0.1,0.1,2)[e],4,c(2,0.1,0.1)[e],1), mgp=c(1.75, 0.25, 0), tck=-0.01)
    boxplot(values ~ groups, outline = FALSE, xlab="",
            at=at_positions, 
            las=2,
            col=rep(c("lightblue", "gold1", "forestgreen"), 12),
            xaxt="n",
            yaxt="n", 
            ylab=paste0(element.names[e], " error", c(" (\u00B0C)", " (\u00B0C)", "")[e]))
    axis(2, at=pretty(values), labels = if(e==3) round(2^(pretty(values)), 2) else pretty(values), tck=-0.01, las=2)
    par(xpd=TRUE)
    if(e==1) legend("topright", legend=c("WRF-ClimatEx", "WRF-USask", "WRF-CONUSII"), 
                    # inset = c(-0, -0.125), 
                    # bty="n", 
                    ncol=3, fill = c("lightblue", "gold1", "forestgreen"))
    par(xpd=FALSE)
    if(e==1) title(main=paste("WY2001-2015 climatology error relative to weather stations", if(i==1) "" else paste(" - ", caseStudy, " Region")))
    mtext(paste0("(",letters[e],")"), side=3, line=-1.5, font=2, adj=0.01)
    # mtext(c("Tmin", "Tmax", "PPT")[e], side=3, line=-1, font=2, adj=0.03)
    
    if(e==3){
      # add month labels under clusters
      month_centers <- tapply(at_positions, month_indices, mean)
      axis(1, at=month_centers, labels=c(month.abb, "All"), las=1, cex.axis=1)
    }
    
    abline(h=0, col="grey")    
    
    print(e)  
  }
  dev.off()
  
  print(paste("region", i))
}
