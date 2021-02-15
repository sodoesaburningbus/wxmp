#This script tests the merra_dap module in WxMP
#Written by Christopher Phillips

#Import required modules
from wxmp import merra_dap as wm
from datetime import datetime as dt

#Create test dates
start_date = dt.strptime("20160618", "%Y%m%d")
end_date = dt.strptime("20160621", "%Y%m%d")

#Password and username
username = "cphillips574"
password = "V3Betredr"

#Create merra analysis object
merra = wm.MERANL(username, password, start_date, end_date, "atmosphere")

#Print dates
print("Dates loaded: {}".format(merra.dates))

#Print number of datasets
print("Number of dates loaded: {}".format(len(merra.dataset)))

#Print grid shape
print("Shape of grid: {}".format(merra.lats.shape))

#Get the indices of a point
point = (47, -17)
x, y = merra.get_point(point)
print("Extent of domain: {}".format(merra.extent))
print("Point being tested: {}".format(point))
print("Indices of point: {}, {}".format(x, y))
print("Found point: {}, {}".format(merra.lons[y,x], merra.lats[y,x]))

#Set the new extent
region = [-60, -40, -10, 20]
merra.set_anl_region(region)

#Print regional lon/lats
print("New region: {}".format(region))
print("Regional Lons and Lats")
print(merra.rlons)
print(merra.rlats)

#Print variable list
print("Variable list: {}".format(merra.variables))

#Grab a variable
var = merra.get_var(merra.variables[0])
print("All time var shape: {}".format(var.shape))

#Now grab a variable for a specific time
var = merra.get_var(merra.variables[0], filedate=dt.strptime("2016060708", "%Y%m%d%H"))
print("Single time var shape: {}".format(var.shape))

print("TEST SUCCESSFUL")
