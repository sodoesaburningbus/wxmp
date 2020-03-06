#This script is for testing the GFS module of WxMP
import matplotlib.pyplot as pp
import wxmp.gfs as wg

#Location of test file
filep = "/rstor/cphillip/model_exps/atlantic_dust/input/20170703/gdas1.fnl0p25.2017070300.f00.grib2"

#Creating GFS grid object
grid = wg.GFSgrid("0.25")

#Now open GFS file object
grib = wg.GFSANL(filep)

#Set analysis region
print(grib.lats.shape)
grib.set_anl_region([-135, -45, 15, 50])
print(grib.lats.shape)
print(grib.rlats.shape)

#Print select messages
data = grib.get_messages(name="Temperature", level=1000, values=True)
#for m in messages:
#    print(m)
 
#for k in messages[0].keys():
#    print(k)

#Do a quick plot
print(data.shape)
fig, ax = pp.subplots(nrows=1, ncols=1, subplot_kw={"projection":grid.proj})
ax.contourf(grib.rlons, grib.rlats, data, transform=grid.pcp)
ax.stock_img()
ax.coastlines()
pp.show()

#Now disable sub-setting
grib.disable_subset()
print(grib.lats.shape)
print(grib.rlats.shape)

#Plot analysis region
grib.plot_anl_region()