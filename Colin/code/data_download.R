
library(aws.s3)
Sys.setenv(
  "AWS_ACCESS_KEY_ID" = "nr-ffec-prd",
  "AWS_SECRET_ACCESS_KEY" = "", #get from Kiri
  "AWS_S3_ENDPOINT" = "nrs.objectstore.gov.bc.ca"
)
save_object("data_WRF_all/climatex_monthly_jan5/monthly_tmax_latlon.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/tmax_latlon_1990_2024.nc") 
save_object("data_WRF_all/climatex_monthly_jan5/monthly_tmin_latlon.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/tmin_latlon_1990_2024.nc")
save_object("data_WRF_all/climatex_monthly_jan5/monthly_precip_latlon.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/pr_latlon_1990_2024.nc")
# save_object("WRF/HGT_lonlat_d03.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/HGT_latlon.nc")

save_object("data_WRF_all/climatex_monthly_jan5/monthly_tmax.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/tmax_1990_2024.nc") 
save_object("data_WRF_all/climatex_monthly_jan5/monthly_tmin.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/tmin_1990_2024.nc")
save_object("data_WRF_all/climatex_monthly_jan5/monthly_precip.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/pr_1990_2024.nc")
# save_object("WRF/HGT_d03.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/HGT.nc")

# save_object("Data_Sharing/conusII_HGT_latlon.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_CONUSII/conus2_HGT_latlon.nc")
