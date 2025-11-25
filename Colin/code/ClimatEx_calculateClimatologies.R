# Process the WRF monthly time series into monthly climatologies
# Colin Mahony colin.mahony@gov.bc.ca

library(terra)

# --------------------------
# Initial QC on data gaps
wrf_pr <- rast("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/ppt_monthly_1990_2024.nc")
wrf_tn <- rast("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/tmin_monthly_1990_2024.nc")
wrf_tx <- rast("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/tmax_monthly_1990_2024.nc")

plot(wrf_pr[[3]])
plot(wrf_tn[[3]])
plot(wrf_tx[[3]])

names(wrf_pr)
names(wrf_tn)
names(wrf_tx)
time(wrf_pr)

month_pr <- substr(sub("T2_Times=", "", names(wrf_pr)), 5, 6)
year_pr <- substr(sub("T2_Times=", "", names(wrf_pr)), 1, 4)
paste(year_pr, month_pr, sep="_")

month_tn <- substr(sub("T2_Times=", "", names(wrf_tn)), 5, 6)
year_tn <- substr(sub("T2_Times=", "", names(wrf_tn)), 1, 4)
paste(year_tn, month_tn, sep="_")

month_tx <- substr(sub("T2_Times=", "", names(wrf_tx)), 5, 6)
year_tx <- substr(sub("T2_Times=", "", names(wrf_tx)), 1, 4)
paste(year_tx, month_tx, sep="_")

library(ncdf4)
nc_pr <- nc_open("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/ppt_monthly_1990_2024.nc")
nc_pr$dim$time$vals  # check if "Time" exists
nc_close(nc_pr)

library(ncdf4)
nc_tx <- nc_open("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/tmax_monthly_1990_2024.nc")
nc_tx$dim$Times$vals  # check if "Time" exists
nc_close(nc_tx)



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
element.names <- c("Tmin (\u00B0C)", "Tmax (\u00B0C)", "precip. (mm)")
monthdays <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
monthcodes <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

projections <- c("lambert", "latlon")

projection = "lambert"
for(projection in projections){
  
  e=3
  # for(e in 1:3){
    element = elements[e]
    
  if(projection == "latlon"){
    file <- paste0("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/", element, "_latlon_1990_2024.nc")
    wrf <- rast(file)
  } else {
    file <- paste0("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/", element, "_monthly_1990_2024.nc")
    wrf <- rast(file)
  }
  
  # exclude lat and lon fields
  if(projection=="lambert") wrf <- wrf[[-which(names(wrf) %in% c("XLONG", "XLAT"))]]
  
  # plot(wrf[[1]])
  
  library(ncdf4)
  nc <- nc_open(file)
  time_vals <- ncvar_get(nc, if(e==3) "time" else "Times")
  nc_close(nc)
  names(wrf) <- paste(substr(time_vals, 1, 4), substr(time_vals, 5,6), sep="-")

  # exclude anomalous months
  if(e==1) wrf <- wrf[[-which(names(wrf)%in%c("2017-10", "2018-04"))]] # there were temperature outliers over ocean in these months; exclude for all
  if(e!=3) wrf <- wrf[[-which(names(wrf)%in%c("2005-09"))]] # corrupted
  if(e==3) wrf <- wrf[[-which(names(wrf)%in%c("1997-06", "1997-07", "1998-07", "1994-08"))]] # corrupted

  # reduce to 1991-2020
  wrf <- wrf[[which(as.numeric(substr(names(wrf), 1, 4))%in% 1991:2020)]] # corrupted
  
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
  
  writeRaster(X, paste0("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/ClimatExWRF_climatology_", element, "_", projection, ".tif"), overwrite=TRUE)
  print(element)
  # }
  
  print(projection)
}

for(m in 1:13){
  plot(X[[m]], main=names(X)[m])
}


# --------------------------
# spot check annual cycle over ocean

png(filename=paste("results\\ClimatExEval.T2mOverOcean", "png",sep="."), type="cairo", units="in", width=6.5, height=4, pointsize=10, res=600)
par(mfrow=c(2,1), mgp=c(1.75, 0.25, 0))
for(e in 1:2){
  element=elements[e]
  file <- paste0("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/", element, "_monthly_1990_2024.nc")
  wrf <- rast(file)

  if(e==1){
    par(mar=c(0,3,2,0.1))
  } else {
    par(mar=c(2,3,0.5,0.1))
  }
  point <- matrix(c(100,500), 1,2)
  ts <- as.vector(unlist(extract(wrf, point)))[-c(1,2)]-273.15
  ts[ts>1000] <- NA
  location <- as.vector(unlist(extract(wrf, point)[c(1,2)]))
  plot(ts, type="l", xaxt="n", ylab=element.names[e], tck=-0.01, xlab="")
  if(e==1) title(main=paste0("Temperature over ocean at (", round(location[1]), "E, ", round(location[2]), "N)"))
  if(e==2) axis(1, at=seq(10, length(wrf_year), 12), label=1990:2024, tck=-0.01)
  points(ts, bg=c(1,0,0,0,0,0,2,0,0,0,0,0), pch=21, cex=c(1,0,0,0,0,0,1,0,0,0,0,0))
  if(e==1) legend("topleft", legend=c("October", "April"), pch=21, pt.bg=c(1,2), bty="n")
}
dev.off()

# --------------------------
# map the anomalies for suspect months

element="tmin"
file <- paste0("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/", element, "_monthly_1990_2024.nc")
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
# isolate corrupted september tmax values

element="tmin"
file <- paste0("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/", element, "_monthly_1990_2024.nc")
wrf <- rast(file)






# --------------------------
# visualize missing/anomalous data

element="tmin"
projection="lambert"

e=1
for(e in 1:3){
  element=elements[e]
  
  file <- paste0("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/", element, "_monthly_1990_2024.nc")
  wrf <- rast(file)

  if(projection=="lambert") wrf <- wrf[[-which(names(wrf) %in% c("XLONG", "XLAT"))]]
  
  # plot(wrf[[1]])
  
  library(ncdf4)
  nc <- nc_open(file)
  time_vals <- ncvar_get(nc, if(e==3) "time" else "Times")
  nc_close(nc)
  names(wrf) <- paste(substr(time_vals, 1, 4), substr(time_vals, 5,6), sep="-")
  
  png(filename=paste("Colin/results/ClimatExEval.GlobalMean", element,  "png",sep="."), type="cairo", units="in", width=6.5, height=4, pointsize=10, res=600)
  
  par(mfrow=c(3,4), mar=c(2,3,2,0), mgp=c(1.75, 0.25, 0))
  
  years <- 1989:2024
  
  for(m in c(10,11,12,1:9)){
    monthcode <- monthcodes[m]
    x <- rep(NA, length(years))
    y <- rep(NA, length(years))
    i=0
    for(year in years){
      i=i+1
      yearmonth <- paste(year, monthcode, sep="-")
      y[i] <- if(yearmonth%in%names(wrf)) as.numeric(unlist(global(wrf[[grep(yearmonth, names(wrf))]], fun = "mean", na.rm = TRUE))) else NA
      print(year)
    }
    plot(years, y, main=month.name[m], ylab=element.names[e], xlab="", pch=16, tck=-0.01)
  }
  dev.off()
}

# are the corrupted values all the same? 
par(mfrow=c(1,1))
yearmonth <- "2005-09"
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
