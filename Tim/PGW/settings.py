#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
description     Settings namelist for all routines in PGW for ERA5
authors		Before 2022:    original developments by Roman Brogli
                Since 2022:     upgrade to PGW for ERA5 by Christoph Heim 
                2022:           udpates by Jonas Mensch
"""
##############################################################################
##############################################################################

### GENERAL SETTINGS
##############################################################################
# debug output level
i_debug = 2 # [0-2]

# Input and output file naming convention for the HIST and climate delta
# (SCEN-HIST) files from the GCM.
# ({} is placeholder for variable name).
"""
file_name_bases = {
    'SCEN-HIST':    '{}_delta.nc',
    'HIST':         '{}_deg1.nc',
}
"""

file_name_bases = {
    'SCEN-HIST':    '{}_{}_{}_monthly_delta.nc',
    'HIST':         '{}_{}_deg1_monthly.nc',
}

# File naming convention for ERA5 files to be read in and written out.
#era5_file_name_base = 'cas{:%Y%m%d%H}0000.nc'
era5_file_name_base = "ERA5_final_{:%Y%m%d%H}.nc"
#era5_file_name_base = 'caf{:%Y%m%d%H}.nc'

# dimension names in ERA5 file
TIME_ERA        = 'time'
LON_ERA         = 'lon'
LAT_ERA         = 'lat'
#LEV_ERA         = 'level'
LEV_ERA         = 'plev'
HLEV_ERA        = 'level1'
#SOIL_HLEV_ERA   = 'soil1'
SOIL_HLEV_ERA1   = 'depth'
SOIL_HLEV_ERA2   = 'depth_2'
SOIL_HLEV_ERA3   = 'depth_3'
SOIL_HLEV_ERA4   = 'depth_4'

# dimension names in GCM (used for all GCM variables except tos)
TIME_GCM        = 'time'
LON_GCM         = 'lon'
LAT_GCM         = 'lat'
PLEV_GCM        = 'plev'
LEV_GCM         = 'lev'

# dimension names in GCM ocean model (used for tos)
TIME_GCM_OCEAN  = 'time'

LON_GCM_OCEAN   = ['longitude','lon']
LAT_GCM_OCEAN   = ['latitude','lat']

### VARIABLE LIST
##############################################################################
# The names on the left side (dict keys) are CMOR convention names
# The names on the right side (dict values) are the name of the
# respective variables in the ERA5 files (Please adjust to ERA5 format used).
# Not all of these variables are required as climate deltas.
# Only zg,ta,hur,ua,va,tas,tos are required as climate delta (SCEN-HIST)
# while ps is required for the HIST climatology.
var_name_map = {

    ##### climate delta (SCEN-HIST) required
    ####################

    # 3D air temperature
    #'ta'   :'T',
    'ta' : 'var130',
    # 3D lon-wind speed
    #'ua'   :'U',
    'ua' : 'var131',
    # 3D lat-wind speed 
    #'va'   :'V',
    'va' : 'var132',
    # 3D air relative humidity
    #'hur'  :'RELHUM',
    'hur' : 'var157',
    # geopotential
    #'zg'   :'PHI', # used for pressure adjustment only
    'zg' : 'var129',
    # near-surface temperature
    #'tas'  :None, # not modified in ERA5 (auxiliary field for computations)
    'tas' : 'var167',
    
    # near-surface u
    #'tas'  :None, # not modified in ERA5 (auxiliary field for computations)
    #'uas' : 'var167',
    
    # near-surface v
    #'tas'  :None, # not modified in ERA5 (auxiliary field for computations)
    #'tas' : 'var167',


    # near-surface relative humidity
    #'hurs' :None, # not modified in ERA5 (auxiliary field for computations)
    'hurs' : 'hurs',
    'hurs_orig':'hurs_orig',
    # sea-surface temperature (SST)
    #'tos'  :None, # not modified in ERA5 (auxiliary field for computations)
    'tos' : 'var34',

    # dewpoint
    'dpt' : 'var168',

    ##### HIST climatology required
    ####################

    # surface pressure
    #'ps'   :'PS', # auxiliary field for interpolation and pressure adjustm.
    'ps' : 'var134',   
 
    ##### no GCM data required but ERA5 variable used by the code
    ####################

    # air specific humidity
    'hus'  :'QV',
    # surface geopotential
    'zgs'  :'FIS', # used for pressure adjustment
    # surface skin temperature
    #'ts'   :'T_SKIN',
    'ts' : 'var235',
    # soil layer temperature
    #'st'   :'T_SO',
    'st1' : 'var139',
    'st2' : 'var170',
    'st3' : 'var183',
    'st4' : 'var236',
    # land area fraction
    #'sftlf':'FR_LAND',
    'sftlf':'var172',
    # sea-ice area fraction
    #'sic':  'FR_SEA_ICE',
    'sic': 'var31',
}


### 02 PREPROCESS DELTAS 
##############################################################################

### SMOOTHING
####################################

### REGRIDDING
####################################
# depending on whether the xesmf pacakge is installed, it can be used
# for interpolation. Else, an xarray-based method is used.
# the latter should be identical to XESMF
# except for tiny differences that appear to originate from
# numerical precision
i_use_xesmf_regridding = 1

## ## Nan-Ingoring kernel interpolation used for tos climate delta
# maximum kernel radius
# higher values imply that remote lakes (and bays) without GCM SST data will
# receive data from further remote GCM SST grid points instead of falling
# back to the tas (near surface temperature) climate delta

nan_interp_kernel_radius = 1000000 # m
#nan_interp_kernel_radius = 10000000 # m


# sharpness: decrease (increase) for smoother (sharper) interpolation
nan_interp_sharpness = 4
# 4


### SURFACE PRESSURE ADJUSTMENT SETTINGS 
##########################################################################
# reference pressure level
# if set to None, the reference pressure level is chosen locally.
# if the climate deltas have low vertical resolution (e.g. Amon data
# with only 6 vertical levels between 1000-500 hPa), settting
# p_ref_inp = None may help to improve the accuray of the
# pressure adjustment. See publication for more information.
p_ref_inp = 30000 # Pa
#p_ref_inp = None
# surface pressure adjustment factor in the iterative routine
adj_factor = 0.95
# convergence threshold (maximum geopotential error)
# if procedure does not converge, raise this value a little bit.
thresh_phi_ref_max_error = 0.15
# maximum number of iterations before error is raised.
max_n_iter = 20
# re-interpolation turned on/off
i_reinterp = 0
##########################################################################

