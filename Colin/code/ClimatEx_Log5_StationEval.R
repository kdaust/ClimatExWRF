
# quantitative evaluation of WRF climatologies against PRISM stations
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

# define the regions
regions <- c("SW", "SE", "NI")
lon1 <- -129.5; lat1 = 48.2
region1 <- ext(c(lon1, lon1+7.5, lat1, lat1+3))
lon1 <- -121.5; lat1 = 48.9
region2 <- ext(c(lon1, lon1+7.5, lat1, lat1+5))
lon1 <- -134; lat1 = 54.5
region3 <- ext(c(lon1, lon1+14, lat1, lat1+5.5))

bc <- vect(bc_bound())
bc <- project(bc, dem.bc)

datasets <- c("prism", "wrfclimatex", "wrfusask", "wrfconus2")
datasets.names <- c("PRISM", "ClimatEx", "USask", "CONUSII")



#######################
## 2 months for 3 elements
#######################

months <- c(1,7)

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
    stn.info <- fread(paste(dir, "Stations/",c("Tmin", "Tmax", "Pr")[e],"_uscdn_8110.csv", sep="")) #read in
    for (i in which(names(stn.info)%in%c(month.abb, "Annual"))) stn.info[get(names(stn.info)[i])==c(-9999), (i):=NA, ] # replace -9999 with NA
    stn.info <- stn.info[-which(El_Flag=="@"),]
    stn <- vect(stn.info, geom = c("Long", "Lat"), crs = "EPSG:4326")
    vals <- extract(prism.bc, stn)
    stn.info <- stn.info[!is.na(vals[, 2])]
    stn.data <- stn.info[,get(month.abb[m])]
    stn.data <- if(e==3) log2(stn.data) else stn.data/10
    stn.info <- stn.info[is.finite(stn.data),]
    stn.data <- stn.data[is.finite(stn.data)]
    
    # load the ClimatEx WRF data for the variable
    dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/"
    wrfclimatex.bc <- rast(paste0(dir, paste("ClimatExWRF_climatology_", c("tmin", "tmax", "pr")[e], "_latlon.tif", sep="")))[[m]]
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
    wrfconus2.bc <- rast(paste0(dir, paste("conus2_climatology_", element, "_latlon.tif", sep="")))[[m]]
    if(e != 3) wrfconus2.bc <- wrfconus2.bc - 273.15
    wrfconus2.bc <- crop(wrfconus2.bc, prism.bc)
    
    
    ###########################################################
    ## for Tmin, plot a key map
    #######################
    
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
           legend = regions,
           title="Regions",
           bg = "white",
           pch = c(21,22,24),  # shapes
           pt.cex = 1.5,
           bty = "n",
           inset = c(0.3,-0.07))   # shift to the left so they appear in second column
    
    for(i in 1:length(regions)){
      studyarea <- get(paste("region", i, sep=""))  
      caseStudy <- regions[i]
      
      ## DEM
      
      dem <- crop(dem.bc, studyarea)
      prism <- crop(prism.bc, studyarea)
      wrfclimatex<- project(wrfclimatex.bc, dem)
      wrfusask<- project(wrfusask.bc, dem)
      wrfconus2<- project(wrfconus2.bc, dem)
      dem.wrf <- project(dem.wrfclimatex, dem)
      
      stn.vect <- vect(stn.info, geom = c("Long", "Lat"), crs = "EPSG:4326")
      stn.crop <- crop(stn.vect, studyarea)
      stn.values <- as.data.frame(stn.crop)[,which(names(stn.crop)==month.abb[m])]
      stn.values <- if(e==3) log2(stn.values) else stn.values/10
      
      for(dataset in datasets){
        d=which(datasets==dataset)
        temp <- get(dataset)
        # plot(temp)
        # plot(stn.crop, add=T)
        stn.temp <- as.vector(unlist(extract(temp, stn.crop)[2]))
        stn.temp <- if(e==3) log2(stn.temp) else stn.temp
        assign(paste("stn", dataset, elements[e], monthcodes[m], sep="."), stn.temp)
        assign(paste("error", dataset, elements[e], monthcodes[m], sep="."), stn.temp - stn.values)
        # Add points on the Taylor diagram
        # taylor.diagram(ref = stn.values, model = stn.temp, add = TRUE,
        #                col = colScheme[d], pch = c(16,15,17)[i], cex = 1.2, normalize = TRUE, pcex=1.5)
        taylor.diagram.filled(ref = stn.values, model = stn.temp, add = TRUE,
                              bg = colScheme[d], pch = c(21,22,24)[i], cex = 1.5, normalize = TRUE)
      }
      
      # print(paste("region", i))
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
    stn.info <- fread(paste(dir, "Stations/",c("Tmin", "Tmax", "Pr")[e],"_uscdn_8110.csv", sep="")) #read in
    for (i in which(names(stn.info)%in%c(month.abb, "Annual"))) stn.info[get(names(stn.info)[i])==c(-9999), (i):=NA, ] # replace -9999 with NA
    stn.info <- stn.info[-which(El_Flag=="@"),]
    stn <- vect(stn.info, geom = c("Long", "Lat"), crs = "EPSG:4326")
    vals <- extract(prism.bc, stn)
    stn.info <- stn.info[!is.na(vals[, 2])]
    stn.data <- stn.info[,get(month.abb[m])]
    stn.data <- if(e==3) log2(stn.data) else stn.data/10
    stn.info <- stn.info[is.finite(stn.data),]
    stn.data <- stn.data[is.finite(stn.data)]
    
    # load the ClimatEx WRF data for the variable
    dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/"
    wrfclimatex.bc <- rast(paste0(dir, paste("ClimatExWRF_climatology_", c("tmin", "tmax", "pr")[e], "_latlon.tif", sep="")))[[m]]
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
    wrfconus2.bc <- rast(paste0(dir, paste("conus2_climatology_", element, "_latlon.tif", sep="")))[[m]]
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
           legend = regions,
           title="Regions",
           bg = "white",
           pch = c(21,22,24),  # shapes
           pt.cex = 1.5,
           bty = "n",
           inset = c(0.3,-0.07))   # shift to the left so they appear in second column
    
    for(i in 1:length(regions)){
      studyarea <- get(paste("region", i, sep=""))  
      caseStudy <- regions[i]
      
      ## DEM
      
      dem <- crop(dem.bc, studyarea)
      prism <- crop(prism.bc, studyarea)
      wrfclimatex<- project(wrfclimatex.bc, dem)
      wrfusask<- project(wrfusask.bc, dem)
      wrfconus2<- project(wrfconus2.bc, dem)
      dem.wrf <- project(dem.wrfclimatex, dem)
      
      stn.vect <- vect(stn.info, geom = c("Long", "Lat"), crs = "EPSG:4326")
      stn.crop <- crop(stn.vect, studyarea)
      stn.values <- as.data.frame(stn.crop)[,which(names(stn.crop)==month.abb[m])]
      stn.values <- if(e==3) log2(stn.values) else stn.values/10
      
      for(dataset in datasets){
        d=which(datasets==dataset)
        temp <- get(dataset)
        # plot(temp)
        # plot(stn.crop, add=T)
        stn.temp <- as.vector(unlist(extract(temp, stn.crop)[2]))
        stn.temp <- if(e==3) log2(stn.temp) else stn.temp
        assign(paste("stn", dataset, elements[e], monthcodes[m], sep="."), stn.temp)
        assign(paste("error", dataset, elements[e], monthcodes[m], sep="."), stn.temp - stn.values)
        # Add points on the Taylor diagram
        # taylor.diagram(ref = stn.values, model = stn.temp, add = TRUE,
        #                col = colScheme[d], pch = c(16,15,17)[i], cex = 1.2, normalize = TRUE, pcex=1.5)
        taylor.diagram.filled(ref = stn.values, model = stn.temp, add = TRUE,
                              bg = colScheme[d], pch = c(21,22,24)[i], cex = 1.5, normalize = TRUE)
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
    stn.info <- fread(paste(dir, "Stations/",c("Tmin", "Tmax", "Pr")[e],"_uscdn_8110.csv", sep="")) #read in
    for (i in which(names(stn.info)%in%c(month.abb, "Annual"))) stn.info[get(names(stn.info)[i])==c(-9999), (i):=NA, ] # replace -9999 with NA
    stn.info <- stn.info[-which(El_Flag=="@"),]
    stn <- vect(stn.info, geom = c("Long", "Lat"), crs = "EPSG:4326")
    vals <- extract(prism.bc, stn)
    stn.info <- stn.info[!is.na(vals[, 2])]
    stn.data <- stn.info[,get(month.abb[m])]
    stn.data <- if(e==3) log2(stn.data) else stn.data/10
    stn.info <- stn.info[is.finite(stn.data),]
    stn.data <- stn.data[is.finite(stn.data)]
    
    # load the ClimatEx WRF data for the variable
    dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/"
    wrfclimatex.bc <- rast(paste0(dir, paste("ClimatExWRF_climatology_", c("tmin", "tmax", "pr")[e], "_latlon.tif", sep="")))[[m]]
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
    wrfconus2.bc <- rast(paste0(dir, paste("conus2_climatology_", element, "_latlon.tif", sep="")))[[m]]
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
           legend = regions,
           title="Regions",
           bg = "white",
           pch = c(21,22,24),  # shapes
           pt.cex = 1.5,
           bty = "n",
           inset = c(0.3,-0.07))   # shift to the left so they appear in second column
    
    for(i in 1:length(regions)){
      studyarea <- get(paste("region", i, sep=""))  
      caseStudy <- regions[i]
      
      ## DEM
      
      dem <- crop(dem.bc, studyarea)
      prism <- crop(prism.bc, studyarea)
      wrfclimatex<- project(wrfclimatex.bc, dem)
      wrfusask<- project(wrfusask.bc, dem)
      wrfconus2<- project(wrfconus2.bc, dem)
      dem.wrf <- project(dem.wrfclimatex, dem)
      
      stn.vect <- vect(stn.info, geom = c("Long", "Lat"), crs = "EPSG:4326")
      stn.crop <- crop(stn.vect, studyarea)
      stn.values <- as.data.frame(stn.crop)[,which(names(stn.crop)==month.abb[m])]
      stn.values <- if(e==3) log2(stn.values) else stn.values/10
      
      for(dataset in datasets){
        d=which(datasets==dataset)
        temp <- get(dataset)
        # plot(temp)
        # plot(stn.crop, add=T)
        stn.temp <- as.vector(unlist(extract(temp, stn.crop)[2]))
        stn.temp <- if(e==3) log2(stn.temp) else stn.temp
        assign(paste("stn", dataset, elements[e], monthcodes[m], sep="."), stn.temp)
        assign(paste("error", dataset, elements[e], monthcodes[m], sep="."), stn.temp - stn.values)
        # Add points on the Taylor diagram
        # taylor.diagram(ref = stn.values, model = stn.temp, add = TRUE,
        #                col = colScheme[d], pch = c(16,15,17)[i], cex = 1.2, normalize = TRUE, pcex=1.5)
        taylor.diagram.filled(ref = stn.values, model = stn.temp, add = TRUE,
                              bg = colScheme[d], pch = c(21,22,24)[i], cex = 1.5, normalize = TRUE)
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
    stn.info <- fread(paste(dir, "Stations/",c("Tmin", "Tmax", "Pr")[e],"_uscdn_8110.csv", sep="")) #read in
    for (i in which(names(stn.info)%in%c(month.abb, "Annual"))) stn.info[get(names(stn.info)[i])==c(-9999), (i):=NA, ] # replace -9999 with NA
    stn.info <- stn.info[-which(El_Flag=="@"),]
    stn <- vect(stn.info, geom = c("Long", "Lat"), crs = "EPSG:4326")
    vals <- extract(prism.bc, stn)
    stn.info <- stn.info[!is.na(vals[, 2])]
    stn.data <- stn.info[,get(month.abb[m])]
    stn.data <- if(e==3) log2(stn.data) else stn.data/10
    stn.info <- stn.info[is.finite(stn.data),]
    stn.data <- stn.data[is.finite(stn.data)]
    
    # load the ClimatEx WRF data for the variable
    dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/"
    wrfclimatex.bc <- rast(paste0(dir, paste("ClimatExWRF_climatology_", c("tmin", "tmax", "pr")[e], "_latlon.tif", sep="")))[[m]]
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
    wrfconus2.bc <- rast(paste0(dir, paste("conus2_climatology_", element, "_latlon.tif", sep="")))[[m]]
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
           legend = regions,
           title="Regions",
           bg = "white",
           pch = c(21,22,24),  # shapes
           pt.cex = 1.5,
           bty = "n",
           inset = c(0.3,-0.07))   # shift to the left so they appear in second column
    
    for(i in 1:length(regions)){
      studyarea <- get(paste("region", i, sep=""))  
      caseStudy <- regions[i]
      
      ## DEM
      
      dem <- crop(dem.bc, studyarea)
      prism <- crop(prism.bc, studyarea)
      wrfclimatex<- project(wrfclimatex.bc, dem)
      wrfusask<- project(wrfusask.bc, dem)
      wrfconus2<- project(wrfconus2.bc, dem)
      dem.wrf <- project(dem.wrfclimatex, dem)
      
      stn.vect <- vect(stn.info, geom = c("Long", "Lat"), crs = "EPSG:4326")
      stn.crop <- crop(stn.vect, studyarea)
      stn.values <- as.data.frame(stn.crop)[,which(names(stn.crop)==month.abb[m])]
      stn.values <- if(e==3) log2(stn.values) else stn.values/10
      
      for(dataset in datasets){
        d=which(datasets==dataset)
        temp <- get(dataset)
        # plot(temp)
        # plot(stn.crop, add=T)
        stn.temp <- as.vector(unlist(extract(temp, stn.crop)[2]))
        stn.temp <- if(e==3) log2(stn.temp) else stn.temp
        assign(paste("stn", dataset, elements[e], monthcodes[m], sep="."), stn.temp)
        assign(paste("error", dataset, elements[e], monthcodes[m], sep="."), stn.temp - stn.values)
        # Add points on the Taylor diagram
        # taylor.diagram(ref = stn.values, model = stn.temp, add = TRUE,
        #                col = colScheme[d], pch = c(16,15,17)[i], cex = 1.2, normalize = TRUE, pcex=1.5)
        taylor.diagram.filled(ref = stn.values, model = stn.temp, add = TRUE,
                              bg = colScheme[d], pch = c(21,22,24)[i], cex = 1.5, normalize = TRUE)
      }
      
      # print(paste("region", i))
    }
    print(month.abb[m])
  }
  dev.off()
  
  print(element)
}



