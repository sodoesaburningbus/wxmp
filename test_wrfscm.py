#Test program for WRFSCM class

#Import supporting modules
import matplotlib.pyplot as pp
import numpy

#Import class
from wrfscm import WRFSCM

#Load the file
scm = WRFSCM("/home/christopher/Code/test_data/wrfscm/wrfout_d01_2006-01-09_00_00_00")

#Grab a variable
temp = scm.get_var(["T2"])["T2"]
print(temp.min())
print(temp.max())

#Set the analysis period
scm.set_anl_period(scm.all_sim_times[24], scm.all_sim_times[48])

#Print the analysis period
print("Start: {}, End: {}".format(scm.start_of_anl.strftime("%Y-%d-%m %H:%M"), scm.end_of_anl.strftime("%Y-%d-%m %H:%M")))

#Plot a meteogram
(fig, ax) = scm.meteogram()
#pp.show()

#Grab a sounding
sounding = scm.get_sounding()
print(sounding.keys())
