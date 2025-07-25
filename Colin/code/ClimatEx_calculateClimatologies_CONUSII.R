# Process WRF monthly time series into monthly climatologies
# Colin Mahony colin.mahony@gov.bc.ca

library(terra)

# --------------------------
# Initial QC on data gaps

wrf_pr <- rast("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/ppt_latlon_2000_2023.nc")
wrf_tn <- rast("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/tmin_latlon_2000_2023.nc")
wrf_tx <- rast("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/tmax_latlon_2000_2023.nc")

plot(wrf_tx[[202]])
names(wrf_pr)
names(wrf_tn)
names(wrf_tx)
time(wrf_tx)

month_tx <- substr(sub("T2_Times=", "", names(wrf_tx)), 5, 6)
year_tx <- substr(sub("T2_Times=", "", names(wrf_tx)), 1, 4)
paste(year_tx, month_tx, sep="_")

month_tn <- substr(sub("T2_Times=", "", names(wrf_tn)), 5, 6)
year_tn <- substr(sub("T2_Times=", "", names(wrf_tn)), 1, 4)
paste(year_tn, month_tn, sep="_")

# missing: 
# 2005_03 through 2005_07
# 2005_10
# 2010_10 through 2011_04
# 2011_06 through 2011_10
# 201708? (timestamp is anomalous)
# 2022_10 through 2022_12

# --------------------------
# Calculate monthly climatologies

elements <- c("tmin", "tmax", "pr")
monthdays <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
monthcodes <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

projections <- c("lambert", "latlon")

projection = "latlon"
# for(projection in projections){
  
  element = elements[3]
  for(element in elements[1:3]){
    
    if(projection == "latlon"){
      file <- paste0("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_CONUSII/conus2_", element, "_latlon.nc")
      wrf <- rast(file)
    } else {
      file <- paste0("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/", element, "_2000_2023.nc")
      wrf <- rast(file)
    }
    
    # exclude lat and lon fields
    if(projection=="lambert") wrf <- wrf[[-which(names(wrf) %in% c("XLONG", "XLAT"))]]
    
    # plot(wrf[[1]])
    
    library(ncdf4)
    nc <- nc_open(file)
    time_vals <- ncvar_get(nc, "Time")
    time_units <- ncatt_get(nc, "Time", "units")$value
    nc_close(nc)
    origin <- sub(".*since ", "", time_units)
    time_seq <- as.POSIXct(time_vals * 60, origin = origin, tz = "UTC")  # adjust multiplier if units are days

    names(wrf) <- substr(time_seq, 1, 7)
    
    #index of wrf months and years
    wrf_month <- substr(names(wrf), 6, 7)
    wrf_year <- substr(names(wrf), 1, 4)
    table(wrf_month)
    
    X <- NULL # initiate the spatraster
    for(m in 1:12){
      monthcode <- monthcodes[m]
      s <- which(wrf_month==monthcode)
      temp <- wrf[[s]]
      temp <- mean(temp)
      # plot(temp)
      X <- if (is.null(X)) temp else c(X, temp)
      names(X)[m] <- paste(min(wrf_year[s]), max(wrf_year[s]), monthcode, sep="_")
      print(m)
    }
    weights <- monthdays / sum(monthdays)
    X_annual <- if(e==3) sum(X) else weighted.mean(X, w = weights)
    X <- c(X, X_annual)
    names(X)[13] <- paste(min(wrf_year), max(wrf_year), "ann", sep="_")
    
    
    writeRaster(X, paste0("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_CONUSII/conus2_climatology_", element, "_", projection, ".tif"), overwrite=TRUE)
    print(element)
  }
  
  # print(projection)
# }

for(m in 1:13){
  plot(X[[m]], main=names(X)[m])
}

  # --------------------------
  # spot check annual cycle over ocean
  
  
# --------------------------
# spot check annual cycle over ocean

png(filename=paste("results\\ClimatExEval.T2mOverOcean", "png",sep="."), type="cairo", units="in", width=6.5, height=4, pointsize=10, res=600)
par(mfrow=c(2,1), mgp=c(1.75, 0.25, 0))
for(e in 1:2){
  element=elements[e]
  file <- paste0("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_CONUSII/conus2_", element, "_latlon.nc")
  wrf <- rast(file)
  
  if(e==1){
    par(mar=c(0,3,2,0.1))
  } else {
    par(mar=c(2,3,0.5,0.1))
  }
  point <- matrix(c(-137,54), 1,2)
  ts <- as.vector(unlist(extract(wrf, point)))[-c(1,2)]-273.15
  plot(ts, type="l", xaxt="n", ylab=element.names[e], tck=-0.01, xlab="")
  if(e==1) title(main=paste0("Temperature over ocean at (", point[1], "E, ", point[2], "N)"))
  if(e==2) axis(1, at=seq(10, length(ts), 12), label=1995:2015, tck=-0.01)
  points(ts, bg=c(1,0,0,0,0,0,2,0,0,0,0,0), pch=21, cex=c(1,0,0,0,0,0,1,0,0,0,0,0))
  if(e==1) legend("topleft", legend=c("October", "April"), pch=21, pt.bg=c(1,2), bty="n")
}
dev.off()


names(wrf)
# --------------------------
# visualize missing data

e=3
for(e in 1:3){
  element==elements[e]
  png(filename=paste("results\\ClimatExEval.GlobalMean", element,  "png",sep="."), type="cairo", units="in", width=6.5, height=4, pointsize=10, res=600)
  
  par(mfrow=c(3,4), mar=c(2,3,2,0), mgp=c(1.75, 0.25, 0))
  for(m in c(10,11,12,1:9)){
    monthcode <- monthcodes[m]
    x <- rep(NA, 25)
    y <- rep(NA, 25)
    i=0
    for(year in 1999:2023){
      i=i+1
      name <- paste(year, monthcode, sep="-")
      y[i] <- if(name%in%names(wrf)) as.numeric(unlist(global(wrf[[which(names(wrf)==name)]], fun = "mean", na.rm = TRUE))) else 0
      print(year)
    }
    plot(1999:2023, y, ylim=c(0, max(y)), main=month.name[m], ylab=element.names[e], xlab="", pch=16, tck=-0.01)
  }
  dev.off()
  }
