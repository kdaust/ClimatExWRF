
library(aws.s3)
Sys.setenv(
  "AWS_ACCESS_KEY_ID" = "nr-ffec-prd",
  "AWS_SECRET_ACCESS_KEY" = "GxyY9ljxtj9OVW7AKdEJ6zpk30YgTH2DPbyj48KR",
  "AWS_S3_ENDPOINT" = "nrs.objectstore.gov.bc.ca"
)
save_object("WRF/tmax_latlon_2000_2023.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/tmax_latlon_2000_2023.nc") 
save_object("WRF/tmin_latlon_2000_2023.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/tmin_latlon_2000_2023.nc")
save_object("WRF/ppt_monthly_latlon_jul23.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/pr_latlon_2000_2023.nc")
save_object("WRF/HGT_lonlat_d03.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/HGT_latlon.nc")

save_object("WRF/tmax_monthly_2000_2023.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/tmax_2000_2023.nc") 
save_object("WRF/tmin_monthly_2000_2023.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/tmin_2000_2023.nc")
save_object("WRF/ppt_monthly_jul23.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/pr_2000_2023.nc")
save_object("WRF/HGT_d03.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/HGT.nc")

save_object("Data_Sharing/conusII_HGT_latlon.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_CONUSII/conus2_HGT_latlon.nc")
