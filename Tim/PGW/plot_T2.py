import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import netCDF4 as nc
from cartopy import crs
from cartopy.feature import NaturalEarthFeature
import wrf
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords)
import matplotlib.ticker as mticker 
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import sys

# Python script for contour plots on individual files and variables
# Usage: python3 plot_2d_contours.py wrfout_file wrfout_variable
# e.g. python3 plot_2d_contours.py wrfout_d03_2000-01-01_18:00:00 T2

if len(sys.argv) != 4:
    print("")
    print("ERROR - args are: wrfout_file wrfout_variable coord_file")
    print("")
    print(sys.exit(1))

myfile = sys.argv[1]
var = sys.argv[2]
coordfile = sys.argv[3]


print(myfile)
print(var)
directory = myfile.split("/")[0]
if len(directory) <= 1:
    directory = directory = myfile.split("/")[1]

# Open up file
ncfile = nc.Dataset(myfile,"r")
nccoordfile = nc.Dataset(coordfile,"r")

# Get variable
#myvar = getvar(ncfile, var)
try:
    myvar = np.array(ncfile["T2"])[0,:,:]
except:
    myvar = np.array(ncfile["var167"])[0,:,:]


# Get units
#units = myvar.units

# Get the latitude and longitude points
try:
    lats, lons = latlon_coords(myvar)
    cart_proj = get_cartopy(myvar)

except:
    hgt = getvar(nccoordfile,"TSK")
    lats, lons = latlon_coords(hgt)
    #lats = ncfile["lat"]
    #lons = ncfile["lon"]
    cart_proj = get_cartopy(hgt)

# Get the cartopy mapping object
#cart_proj = get_cartopy(myvar)

# Create a figure
fig = plt.figure(figsize=(14,10))
# Set the GeoAxes to the projection used by WRF
ax = plt.axes(projection=cart_proj)

# Download and add the states and coastlines
states = NaturalEarthFeature(category="cultural", scale="50m",
                             facecolor="none",
                             name="admin_1_states_provinces_shp")
ax.add_feature(states, linewidth=1.5, edgecolor="black")
ax.coastlines('10m', linewidth=0.8)

# Filled contours

#levels = np.arange(15,45.25,0.25)
#levels = np.arange(-30,20,0.25)

#myvar = myvar - 273.15
levels = np.arange(15,45.25,0.25)

plt.pcolormesh(to_np(lons), to_np(lats), to_np(myvar)-273.15, vmin=15,vmax=45,
             transform=crs.PlateCarree(),
             cmap=mpl.colormaps["jet"])

# Add a color bar
cbar = plt.colorbar(ax=ax, shrink=.98, extend="both")
cbar.set_label(label=r'$^o$C', size='large', weight='bold')
cbar.ax.tick_params(labelsize=14) 

# Set the map bounds
try:
   ax.set_xlim(cartopy_xlim(myvar))
   ax.set_ylim(cartopy_ylim(myvar))
except:
   ax.set_xlim(cartopy_xlim(hgt))
   ax.set_ylim(cartopy_ylim(hgt))

# Add the gridlines
ax.gridlines(color="black", linestyle="dotted")
ax.plot(-121.5814,  50.2333 , 'om', markeredgecolor="k",markersize=8, transform=crs.PlateCarree())

#ax.scatter(-121.5814,50.2333,s=400,c='red',marker='o',transform=crs.PlateCarree(),alpha=0.8,zorder=100000)

# Get titlestring

myfile = myfile.split("/")[-1]

#datetime = myfile.split("_")[2:]
#dom = myfile.split("_")[1]
#datetime=[0]
#hours = datetime[1].split(":")[0]
#hours="00"
plt.title(r"PGW (+2$^o$C) D03 2-m Temp." +r" ($^o$C)   Valid: 20210629_00 UTC",fontsize=16,fontweight="bold")
#plt.title("D03 (3 km) Terrain Heights",fontsize=16,fontweight="bold")
#plt.savefig("d03_hgt.png",box_inches='tight')
#plt.savefig("./"+directory+"/"+dom+"_"+var+"_"+datetime[0]+"_"+hours+".png",bbox_inches='tight',dpi=400)
#plt.savefig("blah.png",box_inches="tight")
plt.tight_layout()
plt.savefig(myfile+"_T2.png",dpi=400,bbox_inches="tight")


print(max(myvar.flatten()-273.15))
print(np.mean(myvar.flatten()-273.15))
print(min(myvar.flatten()-273.15))


