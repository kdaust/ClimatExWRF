
library(aws.s3)
Sys.setenv(
  "AWS_ACCESS_KEY_ID" = "nr-ffec-prd",
  "AWS_SECRET_ACCESS_KEY" = "GxyY9ljxtj9OVW7AKdEJ6zpk30YgTH2DPbyj48KR", #
  "AWS_S3_ENDPOINT" = "nrs.objectstore.gov.bc.ca"
)
save_object("climatex_wrf_monthly_feb20/monthly_tmax_latlon.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/monthly_tmax_latlon.nc") 
save_object("climatex_wrf_monthly_feb20/monthly_tmin_latlon.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/monthly_tmin_latlon.nc")
save_object("climatex_wrf_monthly_feb20/monthly_precip_latlon.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/monthly_pr_latlon.nc")
# save_object("WRF/HGT_lonlat_d03.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/HGT_latlon.nc")

save_object("climatex_wrf_monthly_feb20/monthly_tmax.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/monthly_tmax.nc") 
save_object("climatex_wrf_monthly_feb20/monthly_tmin.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/monthly_tmin.nc")
save_object("climatex_wrf_monthly_feb20/monthly_precip.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/monthly_pr.nc")
# save_object("WRF/HGT_d03.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_ClimatEx/HGT.nc")

# save_object("Data_Sharing/conusII_HGT_latlon.nc", region = "",bucket = 'gmrtde', file = "C:/Users/CMAHONY/OneDrive - Government of BC/Data/WRF_CONUSII/conus2_HGT_latlon.nc")
