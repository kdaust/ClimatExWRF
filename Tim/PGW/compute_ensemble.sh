#!/bin/bash -l

 models=( REGRIDDED/"ACCESS-ESM1-5_2"/tas*delta* REGRIDDED/"BCC-CSM2-MR_2"/tas*delta* REGRIDDED/"CNRM-ESM2-1_2"/tas*delta*
          REGRIDDED/"CanESM5_2"/tas*delta* REGRIDDED/"EC-Earth3_2"/tas*delta* REGRIDDED/"GFDL-ESM4_2"/tas*delta*
          REGRIDDED/"GISS-E2-1-G_2"/tas*delta* REGRIDDED/"INM-CM5-0_2"/tas*delta*
          REGRIDDED/"IPSL-CM6A-LR_2"/tas*delta* REGRIDDED/"MIROC6_2"/tas*delta*
          REGRIDDED/"MPI-ESM1-2-HR_2"/tas*delta* REGRIDDED/"MRI-ESM2-0_2"/tas*delta* 
          REGRIDDED/"UKESM1-0-LL_2"/tas*delta* )


cdo ensmean ${models[@]} REGRIDDED/ens/tas_delta.nc

 


