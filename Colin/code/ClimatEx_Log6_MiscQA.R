# Miscellaneous QA scripts to identify issues. 
# Colin Mahony colin.mahony@gov.bc.ca

library(terra)
library(data.table)
library(scales)
library(RColorBrewer)

# --------------------------
# spot check annual cycle over ocean

png(filename=paste("Colin/results/ClimatExEval.T2mOverOcean", "png",sep="."), type="cairo", units="in", width=6.5, height=6, pointsize=10, res=600)
par(mfrow=c(3,1), mgp=c(1.75, 0.25, 0))
for(e in 1:2){
  element=elements[e]
  file <- paste0("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/", element, "_1990_2024.nc")
  wrf <- rast(file)
  
  # pick a locaiton in the pacific ocean
  point <- matrix(c(100,500), 1,2)
  location <- as.vector(unlist(extract(wrf, point)[c(1,2)]))

  wrf <- wrf[[-c(1,2)]]
  wrf_year <- substr(names(wrf), 10, 13) #index of wrf months and years
  wrf_month <- substr(names(wrf), 14, 15) #index of wrf months and years
  wrf_date <- as.Date(paste(wrf_year, wrf_month, "01", sep="-"))
  
  if(e==1){
    par(mar=c(0,3,2,0.1))
  } else {
    par(mar=c(2,3,0.5,0.1))
  }
  ts <- as.vector(unlist(extract(wrf, point)))-273.15
  ts[ts>1000] <- NA
  plot(ts, type="l", xaxt="n", ylab=element.names[e], tck=-0.01, xlab="")
  if(e==1) mtext("ClimatEx-WRF T2m", line=-1.5, adj=0.0125, font=2, cex=1)
  if(e==1) title(main=paste0("Temperature over ocean at (", round(location[1]), "E, ", round(location[2]), "N)"))
  if(e==2) axis(1, at=seq(5, length(wrf_year), 12), label=1990:2024, tck=-0.01)
  # points(ts, bg=c(0,1,0,0,0,0,0,2,0,0,0,0), pch=21, cex=c(0,1,0,0,0,0,0,1,0,0,0,0))
  points(ts, bg=c(0,0,"blue",0,0,1,0,0,0,0,0,2), pch=21, cex=c(0,0,1,0,0,1,0,0,0,0,0,1))
  # points(ts, bg=c(0,1,0,0,0,0,0,2,0,0,0,0), pch=21, cex=c(0,1,0,0,0,0,0,1,0,0,0,0))
  if(e==2) legend("topleft", legend=c("November", "March", "August"), pch=21, pt.bg=c("blue", 1,2), bty="n")
}

# equivalent plot for ERA5
era5 <- rast("C:/Users/CMAHONY/OneDrive - Government of BC/Data/ERA5/ERA5_pr_tas_Monthly_1940_2023.nc")
era5 <- era5[[grep("t2m_expver=1", names(era5))]]
era5_date <- as.Date(substr(time(era5), 1, 10))

era5 <- era5[[era5_date %in% wrf_date]]
era5_date <- as.Date(substr(time(era5), 1, 10))
ts.era5 <- as.vector(unlist(extract(era5, matrix(location, 1,2))))-273.15

par(mar=c(2,3,1.5,0.1))
plot(ts, col="white", xaxt="n", ylab="Tmean", tck=-0.01, xlab="", ylim=range(ts.era5, na.rm=T))
mtext("ERA5 T2m", line=-1.5, adj=0.0125, font=2, cex=1)
axis(1, at=seq(5, length(wrf_year), 12), label=1990:2024, tck=-0.01)
lines(ts.era5)
points(ts.era5, bg=c(0,0,"blue",0,0,1,0,0,0,0,0,2), pch=21, cex=c(0,0,1,0,0,1,0,0,0,0,0,1))
if(e==1) legend("topleft", legend=c("October", "April"), pch=21, pt.bg=c(1,2), bty="n")

dev.off()



# --------------------------
# map the anomalies for suspect months

element="tmin"
file <- paste0("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/", element, "_1990_2024.nc")
wrf <- rast(file)

X.mean <- rast(paste0("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/ClimatExWRF_climatology_", element, "_lambert.tif"))

lim <- 8
breaks <- seq(0-lim,lim,lim/50)
ColScheme <- colorRampPalette(rev(c(brewer.pal(11, "RdBu")[1:5], rep("white", 1), brewer.pal(11, "RdBu")[7:11])))(length(breaks)-1)

par(mfrow=c(3,4), mar=c(0,0,0,0))
for(m in c(7:12, 1:6)){
  year <- c(rep(2018, 6), rep(2017, 6))[m]
  yearmonth <- paste0(year, monthcodes[m]) 
  X <- wrf[[grep(yearmonth, names(wrf))]]
  X.anom <- X-X.mean[[m]]
  X.anom[X.anom>lim] <- lim
  X.anom[X.anom < 0-lim] <- 0-lim
  plot(X.anom, col=ColScheme, breaks=breaks, type="continuous", axes=F, main=paste(month.abb[m], year))
}


# --------------------------
# visualize missing/anomalous data by plotting the spatial mean time series for the full study area

projection="latlon"

library(bcmaps)
bc <- vect(bc_bound())
bc <- project(bc, "EPSG:4326")


e=1
for(e in 1:3){
  element=elements[e]
  
  if(projection == "latlon"){
    file <- paste0("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/", element, "_latlon_1990_2024.nc")
    wrf <- rast(file)
  } else {
    file <- paste0("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/", element, "_1990_2024.nc")
    wrf <- rast(file)
  }
  
  if(projection=="lambert") wrf <- wrf[[-which(names(wrf) %in% c("XLONG", "XLAT"))]]
  
  # plot(wrf[[1]])
  
  library(ncdf4)
  nc <- nc_open(file)
  time_vals <- ncvar_get(nc, if(e==3) "time" else "Times")
  nc_close(nc)
  names(wrf) <- paste(substr(time_vals, 1, 4), substr(time_vals, 5,6), sep="-")
  
  wrf_year <- substr(names(wrf), 1, 4) #index of wrf months and years
  years <- unique(wrf_year)
  
  png(filename=paste("Colin/results/ClimatExEval.GlobalMean", element,  "png",sep="."), type="cairo", units="in", width=6.5, height=4, pointsize=10, res=600)
  
  par(mfrow=c(3,4), mar=c(2,3,2,0), mgp=c(1.75, 0.25, 0))
  
  
  for(m in c(10,11,12,1:9)){
    monthcode <- monthcodes[m]
    
    #MSWX blend anomalies for comparison
    dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/data_climr_blend_monthly_anomalies/clmr_blend_ts_1901_2024/"
    files <- list.files(dir, pattern = paste0(c("tmin", "tmax", "prcp")[e], "_.*.nc"))
    order <- as.integer(sub(".*_mon_([0-9]+)_.*", "\\1", files))
    obs <- rast(paste0(dir, files[which(order==m)]))
    obs <- obs[[which(substr(time(obs), 1, 4) %in% years)]]
    obs <- crop(obs, bc)
    obs <- mask(obs, bc)
    # plot(obs[[1]])
    
    # convert wrf to anomalies
    wrf.m <- wrf[[grep(paste0("-", monthcode), names(wrf))]]
    names(wrf.m) <- names(wrf)[grep(paste0("-", monthcode), names(wrf))]
    baseline <- mean(wrf.m[[which(years%in%1991:2020)]])
    wrf.m <- if(e == 3) wrf.m / baseline else wrf.m - baseline
    
    # aggregate and clip to obs for direct comparison
    wrf.m <- project(wrf.m, obs)
    wrf.m <- mask(wrf.m, obs[[1]])
    # plot(wrf.m[[1]])

    x <- rep(NA, length(years))
    y.wrf <- rep(NA, length(years))
    y.obs <- rep(NA, length(years))
    i=0
    for(year in years){
      i=i+1
      yearmonth <- paste(year, monthcode, sep="-")
      y.wrf[i] <- if(yearmonth%in%names(wrf.m)) as.numeric(unlist(global(wrf.m[[grep(yearmonth, names(wrf.m))]], fun = "mean", na.rm = TRUE))) else NA
      y.obs[i] <- if(yearmonth%in%substr(time(obs), 1, 7)) as.numeric(unlist(global(obs[[grep(yearmonth, time(obs))]], fun = "mean", na.rm = TRUE))) else NA
      # print(year)
    }
    
    # re-calculate mswx anomalies relative to 1991-2020 baseline
    y.obs <- if(e == 3) y.obs / mean(y.obs[which(years%in%1991:2020)], na.rm=T) else y.obs - mean(y.obs[which(years%in%1991:2020)], na.rm=T)
    
    plot(years, y.wrf, 
         main=paste0(month.name[m], " (r = ", round(cor(y.wrf, y.obs, use = "complete.obs"), 2), ")"), 
         ylab=paste("Anomaly:", element.names[e]), xlab="", tck=-0.01, 
         ylim=range(c(y.wrf, y.obs), na.rm=T))
    lines(years, y.obs, col="gray", lwd=2)
    points(years, y.wrf, pch=16)
    print(m)
  }
  dev.off()
  print(element)
}

# are the corrupted values all the same? 
par(mfrow=c(1,1))
yearmonth <- "2014-09"
X <- wrf[[grep(yearmonth, names(wrf))]]
plot(X)
values(X)[which(values(X) < -100)]
h <- hist(X, plot = FALSE)
h

# --------------------------
# map the anomalies for WY2024

element="tmin"
file <- paste0("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/", element, "_monthly_1990_2024.nc")
wrf <- rast(file)

X.mean <- rast(paste0("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/ClimatExWRF_climatology_", element, "_lambert.tif"))

lim <- 8
breaks <- seq(0-lim,lim,lim/50)
ColScheme <- colorRampPalette(rev(c(brewer.pal(11, "RdBu")[1:5], rep("white", 1), brewer.pal(11, "RdBu")[7:11])))(length(breaks)-1)

par(mfrow=c(3,4), mar=c(0,0,0,0))
for(m in c(10:12, 1:9)){
  year <- c(rep(2024, 9), rep(2023, 3))[m]
  yearmonth <- paste0(year, monthcodes[m]) 
  X <- wrf[[grep(yearmonth, names(wrf))]]
  X.anom <- X-X.mean[[m]]
  X.anom[X.anom>lim] <- lim
  X.anom[X.anom < 0-lim] <- 0-lim
  plot(X.anom, col=ColScheme, breaks=breaks, type="continuous", axes=F, main=paste(month.abb[m], year))
}
