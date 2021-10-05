### This script tests the WxMP RAMSSCM module
### Chris Phillips

# Import modules
import numpy
import wxmp.ramsscm as rscm

# Set data directory
fdir = "/rgroup/smdcot/op_model/nwp/rams_20201201_release_6.3.01/run/chris_test/test/"

# Load the data
rams = rscm.RAMS_SCM(fdir)

# Compute geopotential heights
gph = rams.get_gph()
print(gph.shape)
print(numpy.mean(gph, axis=0))

# Compute the model pressure levels
data = rams.get_vars(["PI"])
pres = ((data["PI"]/rams.CP)**(rams.CP/rscm.at.RD))*rams.P00/100.0
print(pres.shape)
print(numpy.mean(pres, axis=0))

