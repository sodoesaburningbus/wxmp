#This script tests the HRRR module for WxMP
#Chris Phillips
import matplotlib.pyplot as pp
import pygrib
import wxmp.hrrr as wh

#Location of test file
filep = "/home/christopher/Code/test_data/hrrr"

#Demonstration subset extent (west, east, south, north)
extent = [-100, -85, 25, 45]

#Now open HRRR date set object
hrrr = wh.HRRRANL(filep, hrrrname=None)

obj = pygrib.open(hrrr.files[0])

#Print files and dates in dataset
print("Available files: {}".format(hrrr.files))
print("Available dates: {}".format(hrrr.valid_dates))

#Set analysis region
print("HRRR full extent: {}".format(hrrr.full_extent))
print("Shape of HRRR full grid: {}".format(hrrr.lats.shape))
print("Subsetting HRRR grid to {}".format(extent))
hrrr.set_anl_region(extent)
print("Shape of subsetted grid: {}".format(hrrr.rlats.shape))


#Retrieve a variable and show its shape
temp = hrrr.get_var(name="Temperature", level=850, values=True)
print("Shape of subsetted temperature field: {}".format(temp.shape))

#Make a quick plot
print("Plotting simulated reflectivity...")
fig, ax = hrrr.map_var("Maximum/Composite radar reflectivity", 0, hrrr.start_of_anl)
pp.show()
print("Done plotting, exiting example.")
