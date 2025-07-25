
# Evaluation of multiscale precipitation-elevation relationships in the ClimatEx WRF simulations
# Colin Mahony colin.mahony@gov.bc.ca

library(terra)
library(data.table)
library(sf)
library(scales)
library(RColorBrewer)
library(rnaturalearth)

monthcodes <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
month.abb.lowercase <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
monthdays <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
elements <- c("tmin", "tmax", "pr")
element.names <- c("Tmin (\u00B0C)", "Tmax (\u00B0C)", "precip. (mm/day)")
element.units <- c("\u00B0C", "\u00B0C", "mm/day")

e <- 1
m <- 7

#PRISM DEM
dir <- paste("//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/PRISM_dem/", sep="")
dem.prism <- rast(paste(dir, "PRISM_dem.asc", sep=""))
studyarea <- ext(dem.prism)

slope <- terrain(dem.prism, v = "slope", unit = "radians")   # or unit = "degrees" if preferred
aspect <- terrain(dem.prism, v = "aspect", unit = "radians")
hill <- shade(slope, aspect, angle = 40, direction = 270)

# Create an ocean mask
maskarea <- ext(dem.prism)+c(4,4,4,4)
maskpoly <- as.polygons(maskarea) |> st_as_sf()
st_crs(maskpoly) <- crs(dem.prism)
land <- vect(ne_download(scale = 10, type = "land", category = "physical", returnclass = "sf"))
land <- crop(land, maskarea)
land_union <- st_union(st_as_sf(land))
oceanmask_geom <- st_difference(st_geometry(maskpoly), st_geometry(land_union))
oceanmask <- st_sf(geometry = oceanmask_geom, crs = st_crs(maskpoly))
oceanmask <- vect(oceanmask)


dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/"
dem.wrfclimatex <- rast(paste0(dir, "HGT_latlon.nc"))
dem.wrfclimatex <- project(dem.wrfclimatex, dem.prism) # project reduced resolution dem back to prism resolution

dir <- "//objectstore2.nrs.bcgov/ffec/Climatologies/USask_WRF/monthly_clim_regridded/"
dem.wrfusask <- rast(paste(dir, "HGT/HGT_regrid.nc", sep=""))[[1]]
dem.wrfusask <- project(dem.wrfusask, dem.prism) # project reduced resolution dem back to prism resolution

dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_CONUSII/"
dem.wrfconus2 <- rast(paste0(dir, "conus2_HGT_latlon.nc"))
dem.wrfconus2 <- project(dem.wrfconus2, dem.prism) # project reduced resolution dem back to prism resolution


datasets <- c("prism", "wrfclimatex", "wrfusask", "wrfconus2")

dataset <- "wrfclimatex"
for(dataset in datasets){
  #============================================
  ## smoothed and residual elevation
  # latFactor <- cos(mean(extent(dem.wrfhim)[3:4])*(pi/180)) # latitudinal correction for longitudinal length of cell
  w <- round(27/1)
  # studyarea.reduced <- c(-128.64583, -115.09583, 49.22917, 50.77917) # this is what you would get with a w=27. i define it so that the study area is the same for smaller Ws)
  
  dem.d <- focal(get(paste("dem", dataset, sep=".")), w=matrix(1,2*w+1,2*w+1), fun=mean, na.rm=T)
  dem.d <- trim(dem.d)
  dem.o <- crop(get(paste("dem", dataset, sep=".")), dem.d)
  dem.r <- dem.o-dem.d
  
  # par(mfrow=c(1,1))
  # plot(dem.r)
  
  stn.el.o <- stn.info$Elevation
  
  
  e=3
  # for(e in 1:2){
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
      prism <- rast(paste(dir, file, sep=""))
      if(e ==3) prism <- prism/monthdays[m]
      
      # load the ClimatEx WRF data for the variable
      dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/"
      wrfclimatex <- rast(paste0(dir, paste("ClimatExWRF_climatology_", elements[e], "_latlon.tif", sep="")))[[m]]
      wrfclimatex <- if(e ==3) wrfclimatex/monthdays[m] else wrfclimatex - 273.15
      wrfclimatex <- project(wrfclimatex, prism)
 
      # load the USask WRF data for the variable
      if(dataset=="wrfusask"){
        dir <- "//objectstore2.nrs.bcgov/ffec/Climatologies/USask_WRF/monthly_clim_regridded/"
        wrfusask <- rast(paste0(dir, paste(c("tmin", "tmax", "prec")[e], monthcodes[m], "regrid.nc", sep="_")))
        wrfusask <- project(wrfusask, prism)
      }
      
      # load the CONUS II WRF data for the variable
      if(dataset=="wrfconus2"){
        dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_CONUSII/"
      wrfconus2 <- rast(paste0(dir, paste("conus2_climatology_", elements[e], "_latlon.tif", sep="")))[[m]]
      wrfconus2 <- if(e ==3) wrfconus2/monthdays[m] else wrfconus2 - 273.15
      wrfconus2 <- project(wrfconus2, prism)
      }
      
      stn_vect <- vect(stn.info, geom = c("Long", "Lat"), crs = "EPSG:4326")
      stn_crop <- crop(stn_vect, studyarea)
      stn_values <- as.data.frame(stn_crop)[,which(names(stn_crop)==month.abb[m])]
      stn_values <- if(e==3) log2(stn_values) else stn_values/10
      
      #============================================
      ## smoothed and residual climate element
      
      clim <- if(e ==3) log2(get(dataset)) else get(dataset)
      clim.d <- focal(clim, w=matrix(1,2*w+1,2*w+1), fun=mean, na.rm=T)
      clim.d <- trim(clim.d)
      clim.o <- crop(clim, clim.d)
      clim.r <- clim.o-clim.d
      # clim.r <- 100*(2^(clim.o-clim.d)-1)
      # plot(clim.d)
      
      stn.clim.o <- if(e ==3) log2(stn_values) else stn_values
      
      
      #============================================
      ## Scatterplots for distinct regions
      
      
      png(filename=paste("results\\ClimatExEval.LapseRate", element, monthcode,dataset, "png",sep="."), type="cairo", units="in", width=6, height=8, pointsize=10, res=600)
      mat <- rbind(c(rep(1,3),2), c(rep(1,3),3), c(rep(1,3),4), matrix(c(5:16),3, byrow=F))  #define the plotting order
      layout(mat, widths=c(1,1,1,1), heights=c(1,1,1,1,1,1))   #set up the multipanel plot
      
      par(mar=c(3,3,0.1,0.1))
      X <- dem.o
      lim <- quantile(values(X), 0.99, na.rm=T)
      values(X)[which(values(X)>lim)] <- lim
      plot(X, col=terrain.colors(99), mar = NA, axes = FALSE, frame = TRUE,  legend=F, xaxt="n", yaxt="n")
      plot(hill, add=T, col=alpha(grey(0:100/100), 0.5), legend=F)
      mtext(paste("(a)", sep=""), side=1, line=-1.5, adj=0.005, font=2, cex=0.8)
      
      regions <- c("SC", "CM",  "CC", "NC", "NR")
      lon=-125; lat=50; latFactor <- cos(lat*(pi/180))
      region1 <- ext(c(lon, lon+3, lat, lat+1.5))
      
      lon=-119; lat=51; latFactor <- cos(lat*(pi/180))
      region2 <- ext(c(lon, lon+3, lat, lat+1.5))
      
      lon=-130; lat=54; latFactor <- cos(lat*(pi/180))
      region3 <- ext(c(lon, lon+3, lat, lat+1.5))
      
      lon=-135; lat=58; latFactor <- cos(lat*(pi/180))
      region4 <- ext(c(lon, lon+3, lat, lat+1.5))
      
      lon=-128; lat=58; latFactor <- cos(lat*(pi/180))
      region5 <- ext(c(lon, lon+3, lat, lat+1.5))
      
      for(i in 1:5){
        region <- get(paste("region", i, sep=""))
        plot(region, add=T)
        text(mean(ext(region)[1:2]), ext(region)[4], regions[i], font=2, pos=3, offset=0.1)
        # if(stations==T) {
        #   stn.i <- crop(stn.combo, region)
        #   points(stn.i, cex=1, col="blue", pch=1)
        #   points(stn.i[!is.na(apply(data.frame(stn.i)[,4:9], 1,sum)),], cex=1, col="blue", pch=16)
        # }
      }
      
      for(i in 1:5){
        region <- get(paste("region", i, sep=""))
        dem.d.i <- crop(dem.d, region)
        dem.o.i <- crop(dem.o, region)
        dem.r.i <- crop(dem.r, region)
        clim.d.i <- crop(clim.d, region)
        clim.o.i <- crop(clim.o, region)
        clim.r.i <- crop(clim.r, region)
        
        par(mar=c(3,3,0.1,0.1), mgp=c(1.2, 0.2, 0), tck=-0.01)
        # s <- sample(1:length(values(clim.o.i)), length(values(clim.o.i))*0.3)
        s <-1:length(values(clim.o.i))
        y <- values(clim.o.i)[s]
        plot(values(dem.o.i)[s], y, yaxt="n", xaxs="i", pch=16, cex=0.2, xlab="Elevation (m)", ylab=paste0(month.abb[m], ". ", element, " (", element.units[e], ")"))
        if(e==3) axis(2, at=seq(-99,99,1), labels=2^(seq(-99,99,1)), las=2) else axis(2, at=pretty(y), labels=pretty(y)) 
        mtext(paste("(", letters[c(2,5,8,11,14)[i]], ") ", regions[i], sep=""), side=3, line=-1.5, adj=0.05, font=2, cex=0.8)
        
        if(e != 3){
          model1 <- lm(y ~ values(dem.o.i))
        slope1 <- coefficients(model1)[2]
        rSq1 <- summary(model1)$r.squared
        slope1 <- slope1*1000 # convert lapse rate in /m to /km
        abline(model1, col="blue")
        # mtext(paste("slope = lapse rate = ", round(slope,1), "%/100m", sep=""), line=-1.5, adj=0.025, side=3)
        # mtext(paste("R-square =", round(rSq, 2)), line=-2.5, adj=0.025, side=3)
        mtext(paste("lapse = ", round(slope1,1), "%/km", sep=""), line=-1.5, adj=0.95, side=3, cex=0.65, col="blue")
        }
        
        y <- values(clim.d.i)[s]
        plot(values(dem.d.i)[s], y, xlim=range(values(dem.o.i)[s], na.rm=T), yaxt="n", ylim=range(values(clim.o.i)[s], na.rm=T) , xaxs="i", pch=16, cex=0.2, xlab="Smoothed elevation (m)", ylab=paste0("Smoothed ", element))
        if(e==3) axis(2, at=seq(-99,99,1), labels=2^(seq(-99,99,1)), las=2) else axis(2, at=pretty(y), labels=pretty(y)) 
        mtext(paste("(", letters[c(3,6,9,12,15)[i]], ") ", regions[i], sep=""), side=3, line=-1.5, adj=0.05, font=2, cex=0.8)
        
        if(e != 3){
        model1 <- lm(y ~ values(dem.d.i))
        slope1 <- coefficients(model1)[2]
        rSq1 <- summary(model1)$r.squared
        slope1 <- slope1*1000 # convert lapse rate in /m to /km
        abline(model1, col="blue")
        # mtext(paste("slope = lapse rate = ", round(slope,1), "%/100m", sep=""), line=-1.5, adj=0.025, side=3)
        # mtext(paste("R-square =", round(rSq, 2)), line=-2.5, adj=0.025, side=3)
        mtext(paste("lapse = ", round(slope1,1), "%/km", sep=""), line=-1.5, adj=0.95, side=3, cex=0.65, col="blue")
        }
        
        y=values(clim.r.i)[s]
        plot(values(dem.r.i)[s], y, xaxs="i", yaxt="n", xlim=c(-1100, 1100), ylim=c(-2, 1.5), pch=16, cex=0.2, xlab="Elevation residual (m)", ylab=paste0(element, " residual"))
        if(e==3) axis(2, at=seq(-3,3), labels=paste(round(100*(2^(seq(-3,3)))), "%", sep=""), las=2) else axis(2, at=pretty(y), labels=pretty(y)) 
        mtext(paste("(", letters[c(4,7,10,13,16)[i]], ") ", regions[i], sep=""), side=3, line=-1.5, adj=0.05, font=2, cex=0.8)
        
        model1 <- lm(y ~ values(dem.r.i))
        slope1 <- coefficients(model1)[2]
        rSq1 <- summary(model1)$r.squared
        slope1 <- 100*(2^(slope1)-1) # convert slope to % change [not necessary here as residuals are already % change]
        slope1 <- slope1*1000 # convert lapse rate in /m to /km
        abline(model1, col="blue")
        # mtext(paste("slope = lapse rate = ", round(slope,1), "%/100m", sep=""), line=-1.5, adj=0.025, side=3)
        # mtext(paste("R-square =", round(rSq, 2)), line=-2.5, adj=0.025, side=3)
        mtext(paste("lapse = ", round(slope1,0), "%/km", sep=""), line=-1.5, adj=0.95, side=if(e==3) 1 else 3, cex=0.65, col="blue")
        
      }
      
      dev.off()
      
      print(monthcode)
    }
    print(element)
  # }
  print(dataset)
}
