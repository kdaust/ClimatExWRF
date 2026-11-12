from sys import argv
from netCDF4 import Dataset
from wrf import getvar,interplevel

wrffile = argv[1]
wrfpres = "wrfpres" + wrffile[6:]

wrfin = Dataset(wrffile,"r")

p = getvar(wrfin, "pressure")
ht = getvar(wrfin, "z")

ht_500 = interplevel(ht, p, 500.0)
