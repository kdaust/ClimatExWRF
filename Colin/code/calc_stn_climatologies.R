# calculate 2000-2015 climatologies for Aseem's station database, for evaluation of WRF datasets using the USask-WRF period
# USask-WRF period is Oct 1, 2000 to September 29, 2015 (not sure why they missed the last day; sept 30 is there in the other years)

library(data.table)

elements <- c("tmin", "tmax", "pr")
monthdays <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# dir <- "//objectstore2.nrs.bcgov/ffec/bc_stations_database_2025/bc_stns_data_2025/bc_stn_dt_final/"
dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/bc_stn_dt_final/" #local copy for speed

files <- list.files(dir, pattern="*.csv")
files.stn <- sapply(strsplit(files, "_"), "[", 3)

for(e in 1:3){
  stn <- fread(paste0(dir, "bc_obs_stns_metadata_threaded_all.csv"))
  
  stn[, (c(month.abb, "Annual")) := NA_real_]
  
  for(id in stn$stn){
    data <- fread(paste0(dir, files[files.stn==id]))
    data[, year := as.integer(format(date, "%Y"))]
    data[, month := as.integer(format(date, "%m"))]
    for(m in 1:12){
      startyear <- if(m>9) 2000 else 2001
      endyear <- if(m>9) 2014 else 2015
      x <- mean(data[month==m & year >= startyear & year <= endyear, get(elements[e])], na.rm=T)
      if(e==3) x <- round(x*monthdays[m], 1) else round(x, 2)
      
      stn[stn == id, (month.abb[m]) := x]
      # print(month.abb[m])
    }
    if(which(stn$stn==id) %in% seq(0,1000, 10)) print(paste("Processed", which(stn$stn==id), "stations"))
  }

  stn[, Annual :=
        if (e == 3) {
          # precipitation: sum of monthly totals
          rowSums(.SD, na.rm = TRUE)
        } else {
          # temperature: day-weighted mean
          rowSums(.SD * monthdays, na.rm = TRUE) / sum(monthdays)
        },
      .SDcols = month.abb
  ]
  
  # filter out stations without complete data
  stn <- stn[as.integer(format(str_date, "%Y")) <= startyear & as.integer(format(end_date, "%Y")) >= endyear]

  write.csv(stn, paste0(dir, "bc_obs_stns_climatology_WY", startyear+1, "_WY", endyear+1, "_", elements[e], ".csv"), row.names = F)
  print(elements[e])
}

