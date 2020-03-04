#This script is for testing the GFS module of WxMP
import matplotlib.pyplot as pp
import wxmp.gfs as wg

#Location of test file
filep = "/rstor/cphillip/model_exps/atlantic_dust/input/20170703/gdas1.fnl0p25.2017070300.f00.grib2"

#Creating GFS grid object
grid = wg.GFSgrid("0.25")

#Now open GFS file object
grib = wg.GFSANL(filep)

#Print select messages
messages = grib.get_messages(name="Temperature", level=1000)
#for m in messages:
#    print(m)
 
#for k in messages[0].keys():
#    print(k)

#Grab a subset index array
sub = grib.get_subset((-60, -30, 15, 50))
print(sub.shape)

#Do a quick plot
data = messages[0].values
print(data[sub].shape)
fig, ax = pp.subplots(nrows=1, ncols=1, subplot_kw={"projection":grid.proj})
ax.contourf(grib.lons[sub], grib.lats[sub], data[sub], transform=grid.pcp)
ax.stock_img()
ax.coastlines()
pp.show()