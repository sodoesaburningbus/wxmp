#This script tests the HRRR module for WxMP
#Chris Phillips
import matplotlib.pyplot as pp
import pygrib
import wxmp.hrrr as wh

#Location of test file
filep = "/rstor/cphillip/model_exps/grainex/regional/hrrr_analysis/"

#Now open GFS file object
hrrr = wh.HRRRANL(filep)

obj = pygrib.open(hrrr.files[0])

#Set analysis region
print(hrrr.full_extent)
print(hrrr.lats.shape)
hrrr.set_anl_region([-125, -100, 15, 45])
print(hrrr.lats.shape)
print(hrrr.rlats.shape)

#Make a quick plot
fig, ax = hrrr.map_var("Maximum/Composite radar reflectivity", 0, hrrr.start_of_anl)
pp.show()
