from netCDF4 import Dataset
#from wrf import extract_global_attrs
import numpy as np
import sys
import glob


attrs = ["MAP_PROJ","MAP_PROJ_CHAR"]
domain = sys.argv[1]
list_of_wrfuvic = glob.glob("wrfuvic_"+domain+"*00")
list_of_wrfuvic.sort()

if (domain != "d01") and (domain != "d02") and (domain != "d03") and (domain != "d04"):
    print("domain is invalid, must be one of d01, d02, d03 or d04")
    sys.exit(1)

for ind in np.arange(1,len(list_of_wrfuvic)):

    current_file = list_of_wrfuvic[ind]
    prev_file = list_of_wrfuvic[ind-1]

    BUCKET_MM = 100.0
    BUCKET_J = 1.0e9

    cf = Dataset(current_file,"r")
    pf = Dataset(prev_file,"r")

    hourly_precip= (np.array(cf["I_RAINNC"])*BUCKET_MM+np.array(cf["RAINNC"]) + np.array(cf["I_RAINC"])*BUCKET_MM+np.array(cf["RAINC"])
                - (np.array(pf["I_RAINNC"])*BUCKET_MM+np.array(pf["RAINNC"]) + np.array(pf["I_RAINC"])*BUCKET_MM+np.array(pf["RAINC"])))

    hourly_solid = np.array(cf["SNOWNC"]) - np.array(pf["SNOWNC"])

    hourly_graupel = np.array(cf["GRAUPELNC"]) - np.array(pf["GRAUPELNC"])
    hourly_hail = np.array(cf["HAILNC"]) - np.array(pf["HAILNC"])

    hourly_liquid = hourly_precip - hourly_solid
    hourly_snow = hourly_solid - (hourly_graupel + hourly_hail)

    hourlyfile_name = "wrfhourly"+current_file[7:]

    print(hourlyfile_name)

    ncfile = Dataset(hourlyfile_name,mode='w',format='NETCDF4_CLASSIC') 

    south_north = ncfile.createDimension("south_north", hourly_precip.shape[1])
    west_east = ncfile.createDimension("west_east", hourly_precip.shape[2])
    time = ncfile.createDimension("Time", None)

    # Define two variables with the same names as dimensions,
    # a conventional way to define "coordinate variables".
    
    ### PRECIPITATION ###
  
    #cf.MAP_PROJ = extract_global_attrs(cf,"MAP_PROJ")
    #cf.MAP_PROJ_CHAR = extract_global_attrs(cf,"MAP_PROJ_CHAR")

    lat = ncfile.createVariable('XLAT', np.float32, ("Time","south_north","west_east",))
    lat.units = 'degree_north'
    lat.long_name = 'latitude'
    lat.coordinates = 'XLONG XLAT'
    lon = ncfile.createVariable('XLONG', np.float32, ("Time","south_north","west_east",))
    lon.units = 'degree_east'
    lon.long_name = 'longitude'
    lon.coordinates = 'XLONG XLAT' 
   

    hprecipnc = ncfile.createVariable('HPRECIPNC',np.float32,("Time","south_north","west_east",)) # note: unlimited dimension is leftmost
    hprecipnc.units = "mm"
    hprecipnc.coordinates = 'XLONG XLAT XTIME'

    hsolidnc = ncfile.createVariable('HSOLIDNC',np.float32,("Time","south_north","west_east",)) # note: unlimited dimension is leftmost
    hsolidnc.units = "mm"
    hsolidnc.coordinates = 'XLONG XLAT XTIME'

    hgraupelnc = ncfile.createVariable('HGRAUPELNC',np.float32,("Time","south_north","west_east",)) # note: unlimited dimension is leftmost
    hgraupelnc.units = "mm"
    hgraupelnc.coordinates = 'XLONG XLAT XTIME'

    hhailnc = ncfile.createVariable('HHAILNC',np.float32,("Time","south_north","west_east",)) # note: unlimited dimension is leftmost
    hhailnc.units = "mm"
    hhailnc.coordinates = 'XLONG XLAT XTIME'

    hliquidnc = ncfile.createVariable('HLIQUIDNC',np.float32,("Time","south_north","west_east",)) # note: unlimited dimension is leftmost
    hliquidnc.units = "mm"
    hliquidnc.coordinates = 'XLONG XLAT XTIME'

    hsnownc = ncfile.createVariable('HSNOWNC',np.float32,("Time","south_north","west_east",)) # note: unlimited dimension is leftmost
    hsnownc.units = "mm"
    hsnownc.coordinates = 'XLONG XLAT XTIME'

    lat[:] = np.array(cf["XLAT"])
    lon[:] = np.array(cf["XLONG"])

    hprecipnc[:] = hourly_precip[:]
    hsolidnc[:] = hourly_solid[:]
    hgraupelnc[:] = hourly_graupel[:]
    hhailnc[:] = hourly_hail[:]
    hliquidnc[:] = hourly_liquid[:]
    hsnownc[:] = hourly_snow[:]

    ### RADIATION ###
    try:
        #ACSWUPB,ACSWUPBC,ACSWDNB,ACSWDNBC,ACLWUPB,ACLWUPBC,ACLWDNB,ACLWDNBC
        hourly_acswupb  = np.array(cf["I_ACSWUPB"])*BUCKET_J+np.array(cf["ACSWUPB"])   - (np.array(pf["I_ACSWUPB"])*BUCKET_J+np.array(pf["ACSWUPB"]))
        hourly_acswupbc = np.array(cf["I_ACSWUPBC"])*BUCKET_J+np.array(cf["ACSWUPBC"]) - (np.array(pf["I_ACSWUPBC"])*BUCKET_J+np.array(pf["ACSWUPBC"]))
        hourly_acswdnb  = np.array(cf["I_ACSWDNB"])*BUCKET_J+np.array(cf["ACSWDNB"])   - (np.array(pf["I_ACSWDNB"])*BUCKET_J+np.array(pf["ACSWDNB"]))
        hourly_acswdnbc = np.array(cf["I_ACSWDNBC"])*BUCKET_J+np.array(cf["ACSWDNBC"]) - (np.array(pf["I_ACSWDNBC"])*BUCKET_J+np.array(pf["ACSWDNBC"]))
        
        hourly_aclwupb  = np.array(cf["I_ACLWUPB"])*BUCKET_J+np.array(cf["ACLWUPB"])   - (np.array(pf["I_ACLWUPB"])*BUCKET_J+np.array(pf["ACLWUPB"]))
        hourly_aclwupbc = np.array(cf["I_ACLWUPBC"])*BUCKET_J+np.array(cf["ACLWUPBC"]) - (np.array(pf["I_ACLWUPBC"])*BUCKET_J+np.array(pf["ACLWUPBC"]))
        hourly_aclwdnb  = np.array(cf["I_ACLWDNB"])*BUCKET_J+np.array(cf["ACLWDNB"])   - (np.array(pf["I_ACLWDNB"])*BUCKET_J+np.array(pf["ACLWDNB"]))
        hourly_aclwdnbc = np.array(cf["I_ACLWDNBC"])*BUCKET_J+np.array(cf["ACLWDNBC"]) - (np.array(pf["I_ACLWDNBC"])*BUCKET_J+np.array(pf["ACLWDNBC"]))

        hacswupb = ncfile.createVariable('hacswupb',np.float32,("Time","south_north","west_east",))
        hacswupb.units = "J m-2"    
        hacswupb.coordinates = 'XLONG XLAT XTIME'
        hacswupbc = ncfile.createVariable('hacswupbc',np.float32,("Time","south_north","west_east",))
        hacswupbc.units = "J m-2"    
        hacswupbc.coordinates = 'XLONG XLAT XTIME'
        hacswdnb = ncfile.createVariable('hacswdnb',np.float32,("Time","south_north","west_east",))
        hacswdnb.units = "J m-2"    
        hacswdnb.coordinates = 'XLONG XLAT XTIME'
        hacswdnbc = ncfile.createVariable('hacswdnbc',np.float32,("Time","south_north","west_east",))
        hacswdnbc.units = "J m-2"
        hacswdnbc.coordinates = 'XLONG XLAT XTIME'
        
        haclwupb= ncfile.createVariable('haclwupb',np.float32,("Time","south_north","west_east",))
        haclwupb.units = "J m-2"    
        haclwupb.coordinates = 'XLONG XLAT XTIME'
        haclwupbc = ncfile.createVariable('haclwupbc',np.float32,("Time","south_north","west_east",))
        haclwupbc.units = "J m-2"    
        haclwupbc.coordinates = 'XLONG XLAT XTIME'
        haclwdnb = ncfile.createVariable('haclwdnb',np.float32,("Time","south_north","west_east",))
        haclwdnb.units = "J m-2"  
        haclwdnb.coordinates = 'XLONG XLAT XTIME'
        haclwdnbc = ncfile.createVariable('haclwdnbc',np.float32,("Time","south_north","west_east",))
        haclwdnbc.units = "J m-2"
        haclwdnbc.coordinates = 'XLONG XLAT XTIME'

        hacswupb[:] =  hourly_acswupb[:]
        hacswupbc[:] = hourly_acswupbc[:]
        hacswdnb[:] =  hourly_acswdnb[:]
        hacswdnbc[:] = hourly_acswdnbc[:]
        
        haclwupb[:] =  hourly_aclwupb[:]
        haclwupbc[:] = hourly_aclwupbc[:]
        haclwdnb[:] =  hourly_aclwdnb[:]
        haclwdnbc[:] = hourly_aclwdnbc[:]

    except:
        print("Radiation variables may not be available in current file = {} or previous file = {}".format(current_file,prev_file))

    
    cf.close()
    pf.close()
    ncfile.close()



