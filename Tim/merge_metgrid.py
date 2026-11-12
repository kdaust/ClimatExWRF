import numpy as np 
import matplotlib.pyplot as plt
import netCDF4 as nc
import sys

#### IMPORTANT ---> ERA5-LAND MUST **NOT** INCLUDE SNOW_DENSITY ----> ONLY ERA5

allfile_in = sys.argv[1]
landfile_in = sys.argv[2]

allfile = nc.Dataset(allfile_in,"a")
landfile = nc.Dataset(landfile_in,"r")

variables = ["PRES","PSFC","RH","SKINTEMP","SM","SM000007","SM007028","SM028100","SM100289","SNOWH","SNOW",
             "ST","ST000007","ST007028","ST028100","ST100289",
             "TT","UU","VV"]

landsea = np.array(allfile["LANDSEA"])

for each_variable in variables:
    land_variable = np.array(landfile[each_variable])
    land_variable_mag = np.abs(land_variable)

    print(each_variable)

    if each_variable == "RH" or each_variable == "SM" or each_variable == "SM000007" or each_variable == "SM007028" or each_variable == "SM028100" or each_variable == "SM100289" or each_variable == "SNOW" or each_variable == "SNOWH":
        moistflag = True
    else:
        moistflag = False
   
    ndim_land = np.ndim(land_variable)
    if ndim_land == 3:
        if each_variable == "RH":
            land_variable_index = np.where((land_variable_mag < 10**8) & (landsea == 1.0))
        else:
            land_variable_index = np.where((land_variable_mag < 10**8))

        land_variable_subset = land_variable[0,land_variable_index[1],land_variable_index[2]]
        
        allfile_var = allfile[each_variable]
        allfile_numpy = np.array(allfile_var)

        ndim_all = np.ndim(allfile_numpy)

        if ndim_all == 3:
            print(np.max(allfile_numpy))
            print(np.mean(allfile_numpy))
            print(np.min(allfile_numpy))
        else:
            print(np.max(allfile_numpy[0,0,:,:]))
            print(np.mean(allfile_numpy[0,0,:,:]))
            print(np.min(allfile_numpy[0,0,:,:]))


        print("")

        if ndim_all == 3:
            allfile_numpy[0,land_variable_index[1],land_variable_index[2]] = land_variable_subset
        elif ndim_all == 4:
            allfile_numpy[0,0,land_variable_index[1],land_variable_index[2]] = land_variable_subset

        if moistflag:
            #land_variable_subset[land_variable_subset<0] = 0.0
            allfile_numpy[allfile_numpy < 0.0] = 0.0

        allfile_var[:] = allfile_numpy
        allfile_numpy2 = np.array(allfile_var)

        if ndim_all == 3:
            print(np.max(allfile_numpy2))
            print(np.mean(allfile_numpy2))
            print(np.min(allfile_numpy2))
        else:
            print(np.max(allfile_numpy2[0,0,:,:]))
            print(np.mean(allfile_numpy2[0,0,:,:]))
            print(np.min(allfile_numpy2[0,0,:,:]))


        print("")

    elif ndim_land == 4:

        allfile_var = allfile[each_variable]
        allfile_numpy = np.array(allfile_var)

        for k in np.arange(land_variable.shape[1]):

            land_variable_index = np.where(land_variable_mag[0,k,:,:] < 10**8)
            land_variable_subset = land_variable[0,k,land_variable_index[0],land_variable_index[1]]
        
            print(each_variable + " " + str(k))
            print(np.max(allfile_numpy[0,k,:,:]))
            print(np.mean(allfile_numpy[0,k,:,:]))
            print(np.min(allfile_numpy[0,k,:,:]))
            print("")

            allfile_numpy[0,k,land_variable_index[0],land_variable_index[1]] = land_variable_subset

            if moistflag:
                #land_variable_subset[land_variable_subset<0] = 0.0
                allfile_numpy[allfile_numpy<0] = 0.0
            print(np.max(allfile_numpy[0,k,:,:]))
            print(np.mean(allfile_numpy[0,k,:,:]))
            print(np.min(allfile_numpy[0,k,:,:]))
            print("")

        
        allfile_var[:] = allfile_numpy
        allfile_numpy2 = np.array(allfile_var)

    print("")

