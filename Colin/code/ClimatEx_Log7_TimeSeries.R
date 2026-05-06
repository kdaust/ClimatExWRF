# time series comparison figure 
# Colin Mahony colin.mahony@gov.bc.ca

library(terra)
library(data.table)
library(scales)
library(RColorBrewer)
library(bcmaps)
bc <- vect(bc_bound())
bc <- project(bc, "EPSG:4326")


monthcodes <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
month.abb.lowercase <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
monthdays <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 365.25)
elements <- c("tmin", "tmax", "pr")
element.names <- c("Tmin", "Tmax", "Precip.")
element.units <- c("(\u00B0C)", "(\u00B0C)", "")


projection="latlon"

# ------------------------------
# 4 months for three elements

png(filename=paste("Colin/results/ClimatExEval.GlobalMean", "png",sep="."), type="cairo", units="in", width=6.5, height=5, pointsize=10, res=600)

par(mfrow=c(3,4), mar=c(2,3,2,0), mgp=c(1.75, 0.25, 0))

e=1
for(e in 1:3){
  element=elements[e]
  
  if(projection == "latlon"){
    file <- paste0("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/monthly_", element, "_latlon.nc")
    wrf <- rast(file)
  } else {
    file <- paste0("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/monthly_", element, ".nc")
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
  
  
  for(m in c(1,4,7,10)){
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
         main=if(e==1) month.name[m] else "", 
         ylab=if(m==1) paste(element.names[e], "anomaly", element.units[e]) else "", xlab="", tck=-0.01, 
         ylim=range(c(y.wrf, y.obs), na.rm=T))
    lines(years, y.obs, col="gray", lwd=2)
    points(years, y.wrf, pch=16)
    # text(years, y.wrf, years, cex=0.7, pos=2)
    mtext(paste0("r = ", round(cor(y.wrf, y.obs, use = "complete.obs"), 2)), line=-1, side=1, adj=0.95, cex=0.8)
    if(m==1 & e==1) legend("topleft", legend=c("ClimatEx-WRF", "MSWX"), pch=c(16, NA), lwd=c(NA,2), col=c(1, "gray"), bty="n")
    print(m)
  }
  print(element)
}
dev.off()



# ------------------------------
# 12 months for each element

e=1
for(e in 1:3){
  element=elements[e]
  
  png(filename=paste("Colin/results/ClimatExEval.GlobalMean", element, "png",sep="."), type="cairo", units="in", width=6.5, height=5, pointsize=10, res=600)
  
  par(mfrow=c(3,4), mar=c(2,3,2,0), mgp=c(1.75, 0.25, 0))
  
  if(projection == "latlon"){
    file <- paste0("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/monthly_", element, "_latlon.nc")
    wrf <- rast(file)
  } else {
    file <- paste0("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/monthly_", element, ".nc")
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
  
  
  for(m in 1:12){
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
         main=if(e==1) month.name[m] else "", 
         ylab=if(m==1) paste(element.names[e], "anomaly", element.units[e]) else "", xlab="", tck=-0.01, 
         ylim=range(c(y.wrf, y.obs), na.rm=T))
    lines(years, y.obs, col="gray", lwd=2)
    points(years, y.wrf, pch=16)
    # text(years, y.wrf, years, cex=0.7, pos=2)
    mtext(paste0("r = ", round(cor(y.wrf, y.obs, use = "complete.obs"), 2)), line=-1, side=1, adj=0.95, cex=0.8)
    if(m==1 & e==1) legend("topleft", legend=c("ClimatEx-WRF", "MSWX"), pch=c(16, NA), lwd=c(NA,2), col=c(1, "gray"), bty="n")
    print(m)
  }
dev.off()

print(element)
}
