# Process WRF monthly time series into monthly climatologies
# Colin Mahony colin.mahony@gov.bc.ca

library(terra)

# --------------------------
# Calculate monthly climatologies

elements <- c("tmin", "tmax", "pr")
element.names <- c("Tmin (\u00B0C)", "Tmax (\u00B0C)", "Precip.")
monthdays <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
monthcodes <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

projections <- c("lambert", "latlon")

start <- as.Date("2000-10-01")
end   <- as.Date("2015-09-30")

projection = "latlon"
# for(projection in projections){
  
  element = elements[3]
  for(element in elements[1:3]){
    
      file <- paste0("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_CONUSII/conus2_", element, "_latlon.nc")
      wrf <- rast(file)

    # plot(wrf[[1]])
    
    library(ncdf4)
    nc <- nc_open(file)
    time_vals <- ncvar_get(nc, "Time")
    time_units <- ncatt_get(nc, "Time", "units")$value
    nc_close(nc)
    origin <- sub(".*since ", "", time_units)
    time_seq <- as.POSIXct(time_vals * 60, origin = origin, tz = "UTC")  # adjust multiplier if units are days

    names(wrf) <- substr(time_seq, 1, 7)
    wrf_date <- as.Date(paste0(names(wrf), "-15"))
    
    # reduce to the climatological period
    keep <- which(wrf_date >= start & wrf_date <= end)
    wrf <- wrf[[keep]]
    wrf_date <- wrf_date[keep]
    
    #index of wrf months and years
    wrf_month <- format(wrf_date, "%m")
    wrf_year  <- format(wrf_date, "%Y")
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
    
    if(format(start, "%m")=="10") {
      name <- paste0("WY", as.numeric(format(start, "%Y"))+1, "_", "WY", format(end, "%Y"))
    } else {
      name <- paste0(format(start, "%Y"), "_", format(end, "%Y"))
    }
    writeRaster(X, paste0("C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_CONUSII/conus2_climatology_", name, "_", element, "_", projection, ".tif"), overwrite=TRUE)
    print(element)
  }
  
  # print(projection)
# }

