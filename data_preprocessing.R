library(sf)
library(raster)  
#install.packages("ncdf4", type = "binary")
#install.packages("CFtime")
library(ncdf4)
library(CFtime)
library(lattice) 
library(RColorBrewer)
#load the ncdf4 package 
library(ncdf4)
    
#open a netCDF file 
ncin <- nc_open("temp_data_SRS.nc")

unzip("geoBoundaries-GBR-ADM1-all.zip", exdir = "GBR_ADM1")
uk <- st_read("GBR_ADM1/geoBoundaries-GBR-ADM1.shp")

unique(uk$shapeName)
scotland <- uk[uk$shapeName == "Scotland", ]

plot(st_geometry(scotland), main = "Scotland")
