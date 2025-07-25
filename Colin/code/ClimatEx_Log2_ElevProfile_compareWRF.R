# Transects of BC showing maps and elevation profiles of PRISM and WRF climatologies
# Colin Mahony colin.mahony@gov.bc.ca

library(rnaturalearth)
library(terra)
library(data.table)
library(sf)
library(scales)
library(RColorBrewer)

monthcodes <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
month.abb.lowercase <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
elements <- c("tmin", "tmax", "pr")
element.names <- c("Tmin (\u00B0C)", "Tmax (\u00B0C)", "precip. (mm)")

e <- 1
m <- 7


#######################
## Common data
#######################

#PRISM DEM
dir <- paste("//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/PRISM_dem/", sep="")
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

# define the transect regions
regions <- c("SCBC", "SEBC", "CCBC", "NCBC", "NEBC")
lon1 <- -127
region1 <- ext(c(lon1, lon1+7.125, 49.25, 50.25))
lon1 <- -123.25
region2 <- ext(c(lon1, lon1+7.125, 50.35, 51.35))
lon1 <- -133
region3 <- ext(c(lon1, lon1+7.125, 52.5, 53.5))
lon1 <- -138
region4 <- ext(c(lon1, lon1+7.125, 58.5, 59.5))
lon1 <- -129.5
region5 <- ext(c(lon1, lon1+7.125, 57, 58))

#######################
## Key Map
#######################

png(filename=paste("results\\ClimatExEval.Log2.KeyMap.png",sep="."), type="cairo", units="in", width=6.5, height=5.8, pointsize=10, res=600)
par(mar=c(0,0,0,0))
legend.args=list(text='Elevation (m)', side=2, font=2, line=0.5, cex=0.8)
X <- dem.bc
lim <- quantile(values(X), 0.99)
values(X)[which(values(X)>lim)] <- lim
# plot(hill, col=alpha(grey(0:100/100), 1), maxpixels=ncell(hill), legend=F)
plot(X, col=terrain.colors(99), xaxt="n", yaxt="n")
plot(crop(hill.bc,X), add=T, col=alpha(grey(0:100/100), 0.5), legend=F, legend.mar=0)
plot(oceanmask, add=T, col="white", border=F)
# mtext(paste("(a)", sep=""), side=1, line=-1.5, adj=0.005, font=2, cex=0.8)

for(i in 1:5){
  region <- get(paste("region", i, sep=""))
  plot(region, add=T)
  text(ext(region)[1], ext(region)[4]-0.25, regions[i], font=2, pos=4, offset=0.1)
}
l <- ext(X)
rect(l[1], l[3], l[2], l[4], border = "black", lwd = 1)

dev.off()


#######################
## Elevation Profile
#######################

e=3
for(e in 1:2){
  element = elements[e]
  
  m=1
  for(m in c(1,7)){
    monthcode = monthcodes[m]
    
    # load the source STATION data for the BC prism
    dir <- "//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/"
    stn.info <- fread(paste(dir, "Stations/",c("Tmin", "Tmax", "Pr")[e],"_uscdn_8110.csv", sep="")) #read in
    for (i in which(names(stn.info)%in%c(month.abb, "Annual"))) stn.info[get(names(stn.info)[i])==c(-9999), (i):=NA, ] # replace -9999 with NA
    stn.info <- stn.info[-which(El_Flag=="@"),]
    stn.data <- stn.info[,get(month.abb[m])]
    stn.data <- if(e==3) log2(stn.data) else stn.data/10
    stn.info <- stn.info[is.finite(stn.data),]
    stn.data <- stn.data[is.finite(stn.data)]
    
    # load the BC PRISM  data for the variable
    dir <- paste("//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/", sep="")
    file <- list.files(dir, pattern=paste(c("tmin", "tmax", "pr")[e],"_.*._",m, ".tif", sep=""))
    prism.bc <- rast(paste(dir, file, sep=""))

    # load the ClimatEx WRF data for the variable
    dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/"
    wrfclimatex.bc <- rast(paste0(dir, paste("ClimatExWRF_climatology_", c("tmin", "tmax", "pr")[e], "_latlon.tif", sep="")))[[m]]
    if(e != 3) wrfclimatex.bc <- wrfclimatex.bc - 273.15
    wrfclimatex.bc <- crop(wrfclimatex.bc, prism.bc)
    dem.wrfclimatex <- rast(paste0(dir, "HGT_latlon.nc"))
    dem.wrfclimatex <- crop(dem.wrfclimatex, prism.bc)
    
    # load the USask WRF data for the variable
    dir <- "//objectstore2.nrs.bcgov/ffec/Climatologies/USask_WRF/monthly_clim_regridded/"
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

    
    #======================================
    # climate transects
    #======================================

    i=1
    for(i in 1:5){
      slopeData <- "slope.Agg"
      # slopeData <- "slope.wrfna"
      studyarea <- get(paste("region", i, sep=""))  
      caseStudy <- regions[i]
      
      ## DEM
      
      dem <- crop(dem.bc, studyarea)
      slope =  crop(slope.bc, studyarea)
      aspect =  crop(aspect.bc, studyarea)
      hill =  crop(hill.bc, studyarea)
      prism <- crop(prism.bc, studyarea)
      wrfclimatex<- project(wrfclimatex.bc, dem)
      wrfusask<- project(wrfusask.bc, dem)
      wrfconus2<- project(wrfconus2.bc, dem)
      dem.wrf <- project(dem.wrfclimatex, dem)
      
      
      stn_vect <- vect(stn.info, geom = c("Long", "Lat"), crs = "EPSG:4326")
      stn_crop <- crop(stn_vect, studyarea)
      stn_values <- as.data.frame(stn_crop)[,which(names(stn_crop)==month.abb[m])]
      stn_values <- if(e==3) log2(stn_values) else stn_values/10
      
      #################################
      ## plot
      #################################
      
      
      row <- c(60,60,70,70,60)[which(regions==caseStudy)]
      rowlat <- yFromRow(dem, row)
      y.dem <- values(dem, row = row, nrows = 1)
      y.dem[is.na(y.dem)] <- 0
      xvals <- 1:length(y.dem)
      
      datasets <- c("prism", "wrfclimatex", "wrfusask", "wrfconus2")
      datasets.names <- c("PRISM", "WRF-ClimatEx", "WRF-USask", "WRF-CONUSII")
      
      x1 <- get(datasets[1])
      x2 <- get(datasets[2])
      x3 <- get(datasets[3])
      x4 <- get(datasets[4])
      X <- dem
      
      x1 <- mask(x1, land)
      x2 <- mask(x2, land)
      x3 <- mask(x3, land)
      x4 <- mask(x4, land)
      X <- mask(X, land)
      hill <- mask(hill, land)
      
      if(e==3){
        values(x1) <- log2(values(x1))
        values(x1)[!is.finite(values(x1))] <- NA
        values(x2) <- log2(values(x2))
        values(x2)[!is.finite(values(x2))] <- NA
        values(x3) <- log2(values(x3))
        values(x3)[!is.finite(values(x3))] <- NA
        values(x4) <- log2(values(x4))
        values(x4)[!is.finite(values(x4))] <- NA
      }
      
      values(x1)[1:2] <- range(c(values(x1),values(x2),values(x3),values(x4)), na.rm = T)
      values(x2)[1:2] <- range(c(values(x1),values(x2),values(x3),values(x4)), na.rm = T)
      values(x3)[1:2] <- range(c(values(x1),values(x2),values(x3),values(x4)), na.rm = T)
      values(x4)[1:2] <- range(c(values(x1),values(x2),values(x3),values(x4)), na.rm = T)
      
      inc=0.05
      breaks=seq(min(c(values(x1), values(x2), values(x3), values(x4)), na.rm = T)-inc, max(c(values(x1), values(x2), values(x3),values(x4)), na.rm = T)+inc, inc)
      ColScheme <- colorRampPalette(if(e==3) brewer.pal(9, "YlGnBu") else rev(brewer.pal(11, "RdYlBu")))(length(breaks)-1)
      
      
      png(filename=paste("results\\ClimatExEval.Log2.ElevProfile", caseStudy, elements[e], monthcodes[m],"png",sep="."), type="cairo", units="in", width=c(6.27, 6.1, 5.9, 5.05, 5.25)[i], height=5.3*5/4, pointsize=10, res=600)
      mat <- matrix(1:10,5, byrow=T)   #define the plotting order
      layout(mat, widths=c(1,0.125), heights=c(1, 1,1,1))   #set up the multipanel plot
      par(mar=c(0.5,0.5,0.5,0))
      
      ## elevation plot
      
      lim <- quantile(values(X), 0.99, na.rm=T)
      values(X)[which(values(X)>lim)] <- lim
      plot(X, col=terrain.colors(99), mar = NA, axes = FALSE, frame = TRUE,  legend=F, xaxt="n", yaxt="n")
      plot(hill, add=T, col=alpha(grey(0:100/100), 0.5), legend=F)
      lines(c(-180, ext(X)[2]), c(rowlat,rowlat), col="black", lty=2)
      text(xFromCol(dem, 1), rowlat+0.05, "Cross-section", pos=4, cex=1.)
      mtext(paste("(",letters[1],")", sep=""), line=-1.5, adj=0.005, side=3, cex=1)
      box()
      
      ## Terrain legend
      xl <- 0.3; yb <- 0.1; xr <- 0.5; yt <- 0.9
      plot(1, type="n", xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
      rect(xl,  head(seq(yb,yt,(yt-yb)/length(terrain.colors(99))),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(terrain.colors(99))),-1),  border=NA, col=terrain.colors(99))
      rect(xl,  yb,  xr,  yt)
      labels=round(quantile(values(X), c(0, 0.25,0.5,0.75, 1), na.rm=T),0)
      text(rep(xr,length(labels)),seq(yb,yt,(yt-yb)/(length(labels)-1)),labels,pos=4,cex=0.9,font=1, offset=0.1)
      text(xl-0.15, mean(c(yb,yt))-0.04, paste("Elevation (m)"), srt=90, pos=3, cex=1.2, font=2)
      
      
      ## climate plots
      
      for(dataset in datasets[1:3]){
        
        plot(get(paste("x", which(datasets==dataset), sep="")), col=ColScheme, zlim = c(min(breaks), max(breaks)), mar = NA, axes = FALSE, frame = TRUE, legend=F, main="")
        lines(c(-180, 0), c(rowlat,rowlat), col="black", lty=2)
        mtext(datasets.names[which(datasets==dataset)], side=1, line=-1.5, adj=0.01, font=2, cex=0.8)
        # for(i in 1:length(col)){
        #   box.lon1 <- xFromCol(dem,col[i]-(oddRound(ngb/latFactor)-1)/2)
        #   box.lon2 <- xFromCol(dem,col[i]+(oddRound(ngb/latFactor)-1)/2)
        #   box.lat1 <- yFromRow(dem,row+(ngb-1)/2)
        #   box.lat2 <- yFromRow(dem,row-(ngb-1)/2)
        #   rect(box.lon1, box.lat1, box.lon2, box.lat2)
        #   text(xFromCol(dem,col[i]),  yFromRow(dem,row), letters[6:26][i], pos=3, cex=1, font=2, offset=0.8)
        # }
        mtext(paste("(",letters[2:26][which(datasets==dataset)],")", sep=""), line=-1.5, adj=0.005, side=3, cex=1)
        if(dataset=="prism"){ 
          plot(stn_crop, add=T, pch=21, bg=ColScheme[cut(stn_values, breaks=breaks)], cex=1.25, lwd=0.5)
        }
        box()
        
        ## legend
        xl <- 0.3; yb <- 0.1; xr <- 0.5; yt <- 0.9
        plot(1, type="n", xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
        rect(xl,  head(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  xr,  tail(seq(yb,yt,(yt-yb)/length(ColScheme)),-1),  border=NA, col=ColScheme)
        rect(xl,  yb,  xr,  yt)
        labels=round(quantile(breaks, c(0, 0.25,0.5,0.75, 1)),0)
        text(rep(xr,length(labels)),seq(yb,yt,(yt-yb)/(length(labels)-1)),labels,pos=4,cex=0.9,font=1, offset=0.1)
        text(xl-0.15, mean(c(yb,yt))-0.04, paste0(month.abb[m], ". ", element.names[e]), srt=90, pos=3, cex=1.2, font=2)
        
      }
      
      ## latitudinal cross section
      transect_wrfconus2 <- if(e==3) log2(values(wrfconus2, row = row, nrows = 1)) else values(wrfconus2, row = row, nrows = 1)
      transect_wrfclimatex <- if(e==3) log2(values(wrfclimatex, row = row, nrows = 1)) else values(wrfclimatex, row = row, nrows = 1)
      transect_wrfusask <- if(e==3) log2(values(wrfusask, row = row, nrows = 1)) else values(wrfusask, row = row, nrows = 1)
      transect_prism <- if(e==3) log2(values(prism, row = row, nrows = 1)) else values(prism, row = row, nrows = 1)
      transect_combined <- c(transect_wrfconus2, transect_wrfclimatex, transect_wrfusask, transect_prism)
      
      ylim.lower=c(0,0,0,0,0)[which(regions==caseStudy)]
      ylim.upper=c(5000, 5000, 5000, 5000, 5000)[which(regions==caseStudy)]
      shiftFactor <- min(transect_combined, na.rm=T)
      scaleFactor <- ylim.upper/diff(range(transect_combined, na.rm=T))*0.95

      par(mar=c(0.5,0.5,0.5,0), mgp=c(2, 1, 0))
      plot(xvals, y.dem, col="white", xaxs="i", yaxs="i", ylim=c(ylim.lower,ylim.upper), ylab="", xlab="Degrees longitude", yaxt="n", xaxt="n")
      latseq <- seq(ceiling(ext(dem)[1]), floor(ext(dem)[2]),2)
      # axis(1, at=colFromX(dem, latseq), labels=latseq, tck=0)
      par(mgp=c(2, 0.1, 0))
      
      labels <- pretty(transect_combined)
      axis(4, at=(labels-shiftFactor)*scaleFactor, labels=if(e==3) round(2^labels) else labels, tck=0, las=2)
      mtext(side = 4, line = 2, paste0(month.abb[m], ". ", element.names[e]), cex=0.8, font=2)
      polygon(c(xvals, length(y.dem), 1), c(y.dem, 0, 0), col="gray", border=F)
      lines(xvals, values(dem.wrf, row = row, nrows = 1), col="grey50", lwd=1, lty=2)
      lines(xvals, (transect_wrfusask-shiftFactor)*scaleFactor, col="grey50", lwd=4)
      lines(xvals, (transect_wrfconus2-shiftFactor)*scaleFactor, col="red", lwd=1.2)
      lines(xvals, (transect_wrfclimatex-shiftFactor)*scaleFactor, col="blue", lwd=2)
      lines(xvals, (transect_prism-shiftFactor)*scaleFactor, col="black", lwd=2)
      # mtext("Cross-section", side=1, line=-1.5, adj=0.01, font=2, cex=0.8)
      mtext(paste("(",letters[5],")", sep=""), line=-1.5, adj=0.005, side=3, cex=1)
      # for(i in 1:length(col)){
      #   left <- col[i]-(oddRound(ngb/latFactor)-1)/2
      #   right <- col[i]+(oddRound(ngb/latFactor)-1)/2
      #   arrows(left,300, right, 300, length=0.01, angle=90, code=3)
      #   text(right,  300, letters[6:26][i], pos=4, cex=1, font=2)
      # }
      legend(c("topright", "topright", "topright", "topright", "topright")[which(regions==caseStudy)], ylim.upper+100, legend=datasets.names, col=c("black", "blue", "grey50", "red"), lty=c(1,1,1,1), lwd=c(2,2,4,1.2) , bty="n", y.intersp = 0.8, cex=1)
      box()
      
      dev.off()
      
      print(paste("transect", i))
    }
    print(month.abb[m])
  }
  print(element)
}



