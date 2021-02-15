### This script tests the features of the wxmpi module
### Written by Christopher Phillips
### Univeristy of Alabama in Huntsville
### Atmospheric and Earth Science Department

##### START OPTIONS #####

#Select model dataset
model = "WRF"

#Location of dataset
loc = "/home/christopher/Code/test_data/wrf"

#Timestamp to plot
timestamp = "2017070500"

#####  END OPTIONS  #####

#Importing required modules
import datetime as dt
import matplotlib.pyplot as pp
import wxmp.wxmpi as ww

#Load data set
sim = ww.wxmpi(model, loc)

#var = sim.get_var2D("T", level=850)
#print(var.shape)

#var = sim.get_var3D("T")
#print(var.shape)

#Convert timestamp
date = dt.datetime.strptime(timestamp, "%Y%m%d%H")

#Grab sounding
sounding = sim.get_sounding((0, 0), filedate=date)
pp.plot(sounding["temp"][0,:], sounding["pres"][0,:], color="red")
pp.plot(sounding["dewp"][0,:], sounding["pres"][0,:], color="green")
pp.ylim(1000, 100)
pp.show()

#Plot map
#(fig, ax) = sim.map_var("T", 4, date)
#pp.show()
