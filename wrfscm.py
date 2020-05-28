### This module contains classes for handling Weather Research and Forecasting model 
### grids, inputs, and outputs when configured for single-column operations.
### Unlike other WRF simulations, it is assumed that the SCM output is all contained within a single file.
### This is due to SCM files typically being much smaller than their regional counterparts.
### Written by Christopher Phillips
### University of Alabama in Huntsville
### Atmospheric and Earth Science Department
###
### Requirements:
###   Python 3+
###   Atmos (available at github.com/sodoesaburningbus)
###   Cartopy
###   f90NML
###   Matplotlib
###   NetCDF4
###   Numpy
###   PIL
###   PyMODIS (available at github.com/sodoesaburningbus)

### Importing required modules
import atmos.thermo as at
import datetime as dt
import glob
import f90nml
import matplotlib.patches as patches
import matplotlib.pyplot as pp
import netCDF4 as nc
import numpy
from PIL import Image
import sys
        
############################################################
#---------------------     WRFSCM     ---------------------#
############################################################
### The WRFSCM class provides tools for visualizing WRFSCM
### model output. It provides tools for quick creation
### of common products such as soundings and meteograms.
### Further, it enables the setting of an analysis period
### to quickly subset data to remove spin-up period.
###
### NOTE: This object only supports one data frame per file.
###
### Available Attributes:
###  file - File name
###  lat0 - Center latitude of simulation
###  lon0 - Center longitude of simulation
###  tformat - Format of WRF output file times
###  start_of_sim - Datetime date object with simulation start time
###  end_of_sim - Datetime date object with simulation end time
###  start_of_anl - Datetime date object with start of analysis period
###  end_of_anl - Datetime date object with end of analysis period
###
### Available Methods:
###  get_var(names) - Returns dictionary containing requested variables
###  get_wxc((lon, lat)) - Returns WxChallenge forecast for desired point
###  meteogram(save_location) - Plots a meteogram from simulation output
###  set_period(tstart, tend) - Set beginning and end of analysis period
class WRFSCM:

    ### Method to construct object
    ### Inputs:
    ###  wrfpath, string, path to WRF output files
    ###  wrfpre, string, optional, prefix to WRF files, defaults to "wrfout_"
    def __init__(self, wrfpath, wrfpre=None):
        
        ### Setup WRF file prefix
        ### If statement supports wxmpi compatibility
        if (wrfpre == None):
            wrfpre = "wrfout_"
                
        #Store file name as object attribute
        self.filename = wrfpath
        
        #Go ahead and open the file for access
        self.ncfile = nc.Dataset(self.filename, "r")
    
        #Now construct an array of all times in the file
        self.tformat = "%Y-%m-%d_%H:%M:%S"
        self.all_sim_times = [] #List to hold the times
        times = self.ncfile.variables["Times"][:]
        for i in range(times.shape[0]):
            time = "".join(numpy.array(times[i,:], dtype="str"))
            self.all_sim_times.append(dt.datetime.strptime(time, self.tformat))
        self.all_sim_times = numpy.array(self.all_sim_times)

        #Pull out the start and end times fo the simulation for consistency with other WxMP modules
        self.start_of_sim = self.all_sim_times[0]
        self.end_of_sim = self.all_sim_times[-1]

        #Set period of analysis to simulation period
        self.set_anl_period(self.start_of_sim, self.end_of_sim)

        ### Construct list of available variables
        #Loop across all variables, pulling the 3D fields' level types
        self.type_of_level_list = []
        for k in self.ncfile.variables.keys():
            if ((len(self.ncfile.variables[k].dimensions) == 4) and
                (self.ncfile.variables[k].dimensions[1] not in self.type_of_level_list)):
                self.type_of_level_list.append(self.ncfile.variables[k].dimensions[1])
        
        #Add a level type for 1D and 2D vars
        self.type_of_level_list.append("1D")
        self.type_of_level_list.append("2D")
        
        #Now create dictionaries with lists for the variables
        self.var_list = {}
        self.level_list = {}
        for tl in self.type_of_level_list:
            self.var_list[tl] = []
            self.level_list[tl] = []
                    
        #Now populate variable and level list by level type
        for k in self.ncfile.variables.keys():
            if (len(self.ncfile.variables[k].dimensions) == 4):
                self.var_list[self.ncfile.variables[k].dimensions[1]].append(self.ncfile.variables[k].name)
                self.level_list[self.ncfile.variables[k].dimensions[1]].append(numpy.arange(0,self.ncfile.variables[k].shape[1]))
            elif (len(self.ncfile.variables[k].dimensions) == 3):
                self.var_list["2D"].append(self.ncfile.variables[k].name)
                self.level_list["2D"].append(1)
            else:
                self.var_list["1D"].append(self.ncfile.variables[k].name)
                self.level_list["1D"].append(1)
        
        ### Returning
        return
            
    ### Method to retreive sounding
    ### Inputs:
    ###   point, tuple of floats, (lon, lat) or list of tuples of such points.
    ###   gid, integer, optional, grid id to pull sounding from. Defaults to current working grid.
    ###   period, tuple of date objects, optional, start and end times of desired analysis Defaults to working analysis period.
    ###   filedate, datetime object, optional, date of file to pull. (Only retreives one file's data if set). Defaults to None.
    ###
    ### Outputs:
    ###   sounding, dictionary of lists containing sounding info, or list of such dictionaries with len(points)
    ###     dictionaries are keyed ["temp", "pres", "dewp", "uwind", "vwind"] for temperature (K), pressure (hPa),
    ###     dewpoint (K), zonal wind speed (m/s), and meriodinal wind speed (m/s) respectively.
    ###     Note that winds are rotated to Earth coordinates.
    def get_sounding(self, period=None, anldate=None):
            
        #Determine time period to retreive variables over
        if (period == None):
            tstart = self.start_of_anl
            tend = self.end_of_anl
        else:
            tstart = period[0]
            tend = period[1]
        
        #Pull soundings at given locations
        sounding = []

        #Retreive data
        data = self.get_var(["T", "P", "PB", "QVAPOR", "U", "V", "SINALPHA", "COSALPHA"], period=(tstart, tend), anldate=anldate)

        #Calculate actual temperature and pressure
        ptemp = data["T"]+300.0 #Add WRF base temp to variable
        pres = (data["P"]+data["PB"]) #Pressure in Pa. On return converts to hPa.
        temp = at.poisson(100000.0, pres, ptemp) #Convert potential temperature to temperature

        #Calculate dewpoint
        dewp = at.dewpoint(at.wtoe(pres, data["QVAPOR"]))

        #Rotate winds
        try: #Case for multiple soundings
            uwind = data["U"]*data["COSALPHA"][0]-data["V"]*data["SINALPHA"][0]
            vwind = data["V"]*data["COSALPHA"][0]+data["U"]*data["SINALPHA"][0]
        except: #Case for single sounding
            uwind = data["U"]*data["COSALPHA"][0]-data["V"]*data["SINALPHA"][0]
            vwind = data["V"]*data["COSALPHA"][0]+data["U"]*data["SINALPHA"][0]

        #Append sounding to list
        sounding.append({"temp":numpy.atleast_2d(temp), "pres":numpy.atleast_2d(pres/100.0),
            "dewp":numpy.atleast_2d(dewp), "uwind":numpy.atleast_2d(uwind),
            "vwind":numpy.atleast_2d(vwind)})

        #Return soundings as list
        #or as dictionary if only one.
        if (len(sounding) > 1):
            return sounding
        else:
            return sounding[0]

    ### Method to pull variables from the WRF simulation.
    ### This method retrieves user-requested variables over a time-period
    ### for the requested grid. It eliminates extraneous dimensions in the process.
    ### Inputs:
    ###   var_labels, list of strings, WRF variable names to retrieve.
    ###   period, optional, tuple of datetime objects, (start, end) temporal bounds of plot (inclusive).
    ###     Defaults to current analysis period.
    ###   level, integer, optional, model level to pull (First level is zero). Defaults to all.
    ###   anldate, datetime object, optional, date of analysis to retrieve (a single time stamp)
    ###
    ### Outputs:
    ###   data, dictionary of lists, contains arrays with WRF variables keyed to var_labels (Final array order is (Time, Z).
    ###     or if specific level requested, order is (Time).
    def get_var(self, var_labels, period=None, level=None, anldate=None):
                
        #Determine time period to retreive variables over
        if ((period == None) and (anldate == None)):
            tstart = self.start_of_anl
            tend = self.end_of_anl
        elif ((period != None) and (anldate == None)): #A period is requested
            tstart = period[0]
            tend = period[1]
        else: #A specific date is requested
            tstart = anldate
            tend = tstart
        
        #Find indices of requested dates
        start_ind = numpy.where(self.all_sim_times == tstart)[0][0]
        end_ind = numpy.where(self.all_sim_times == tend)[0][0]

        #Create dictionary to store requested variables
        data = {"date":self.all_sim_times[start_ind:end_ind+1]}
            
        #Grab other variables (handle differently if doing subsetting)
        try:
            for vl in var_labels: #Loop over each variable
                if (len(self.ncfile.variables[vl].shape) == 3): #2D case
                    data[vl] = numpy.squeeze(self.ncfile.variables[vl][start_ind:end_ind+1,0,0])

                elif (len(self.ncfile.variables[vl].shape) == 4): #3D case
                    data[vl] = numpy.squeeze(self.ncfile.variables[vl][start_ind:end_ind+1,:,0,0])        
        except Exception as err:
            print(err)
            raise Exception
    
        #Convert lists to numpy arrays (and pull desired level if so)
        for k in data.keys():
            if (level != None): #Pull a level
                try: #3D case
                    data[k] = numpy.array(data[k][-1][:,level])
                except Exception as err: #Only fails in 2D case
                    pass
            else: #No desired level
                data[k] = numpy.array(data[k])
                
        #Returning
        return data       
    
    ### Method to plot meteogram from WRF simulation
    ### Temperature, dewpoint, wind speed, and solar radiation are plotted for a single
    ### point in a WRF simulation. By default, this is averaged over the whole domain.
    ### If a point is provided by the user, that is used instead.
    ### Inputs:
    ###   point, optional, tuple, (lon, lat) of point to be plotted. Default is average over domain.
    ###   gid, optional, integer, grid to analyze. Default is current grid.
    ###   period, optional, tuple of datetime objects, (start, end) temporal bounds of plot.
    ###     Defaults to current analysis period.
    ###
    ### Outputs:
    ###   (fig, ax) tuple with pyplot figure and axis object.
    def meteogram(self, period=None):
                
        #Determine time period of meteogram
        if (period == None):
            tstart = self.start_of_anl
            tend = self.end_of_anl
        else:
            tstart = period[0]
            tend = period[1]
        
        #Retrieve variables to plot
        var_names = ["T2", "Q2", "U10", "V10", "SWDOWN", "PSFC", "RAINC", "RAINNC"]
        data = self.get_var(var_names, period=period)
                
        #Calculate dewpoint from water vapor mixing ratio and change to Fahrenheit
        dewp = (at.dewpoint(at.wtoe(data["PSFC"], data["Q2"]))-273.15)*9/5+32
        
        #Convert temperature from Kelvin to Fahrenheit
        temp = (data["T2"]-273.15)*9/5+32
        
        #Calculate total wind speed and convert to knots
        wspd = numpy.sqrt(data["U10"]**2+data["V10"]**2)*1.94384
        
        #Calculate accumulated rainfall over analysis period and convert to inches     
        precip = ((data["RAINC"]+data["RAINNC"])-(data["RAINC"][0]+data["RAINNC"][0]))*0.0393701
        
        #Calculate ranges for each variable
        tmin = numpy.floor(numpy.nanmin(dewp)*0.9)
        tmax = numpy.ceil(numpy.nanmax(temp)*1.1)
        wmin = numpy.floor(numpy.nanmin(wspd)*0.9)
        wmax = numpy.ceil(numpy.nanmax(wspd)*1.1)
        pmin = 0
        pmax = max(numpy.ceil(numpy.nanmax(precip)*1.1),2)
        rmin = 0
        rmax = numpy.ceil(numpy.max(data["SWDOWN"])*1.1)
        
        #Compute ticks for each variable
        tticks = numpy.arange(tmin, tmax+4, 4)
        wticks = numpy.arange(wmin, wmax+4, 4)
        pticks = numpy.linspace(pmin, pmax, 11)
        rticks = numpy.linspace(rmin, rmax, 11)        
        
        ### Plot the variables
        #Create figure and axis objects
        fig, ax = pp.subplots(nrows=4, ncols=1, figsize=(6,10))
        
        #Temperature and Dewpoint
        ax[0].plot(data["date"], temp, color="red")
        ax[0].plot(data["date"], dewp, color="green")
        ax[0].set_ylabel("Temperature [K]", fontsize=14)
        ax[0].legend(["Temperature", "Dewpoint"])
        ax[0].tick_params(axis="x", rotation=30)
        ax[0].set_ylim((tmin, tmax))
        ax[0].set_yticks(tticks)
        ax[0].grid()
        
        #Wind Speed
        ax[1].plot(data["date"], wspd, color="black")
        ax[1].set_ylabel("Wind Speed [kt]", fontsize=14)
        ax[1].tick_params(axis="x", rotation=30)
        ax[1].set_ylim((wmin, wmax))
        ax[1].set_yticks(wticks)
        ax[1].grid()
        
        #Precipitation
        ax[2].plot(data["date"], precip, color="dodgerblue")
        ax[2].set_ylabel("Accumulated Precip [inch]", fontsize=14)
        ax[2].tick_params(axis="x", rotation=30)
        ax[2].set_ylim((pmin, pmax))
        ax[2].set_yticks(pticks)
        ax[2].grid()
        
        #Solar radiation
        ax[3].plot(data["date"], data["SWDOWN"], color="black")
        ax[3].set_ylabel("Insolation [Wm$^{-2}$]", fontsize=14)
        ax[3].tick_params(axis="x", rotation=30)
        ax[3].set_ylim((rmin, rmax))
        ax[3].set_yticks(rticks)
        ax[3].grid()
        
    
        #Returning
        return (fig, ax)
        
    ### Method to set working time period
    ### Inputs:
    ###   tstart, datetime object, start of simulation period to be analyzed
    ###   tend, datetime object, end of simulation period to be analyzed
    def set_anl_period(self, tstart, tend):
        
        #Check that times are within bounds of simulation
        if ((tstart < self.start_of_sim) or (tstart > self.end_of_sim) or
            (tend < self.start_of_sim) or (tend > self.end_of_sim)):
            print("ERROR: Requested period outside of simulation period. Exiting...")
            exit()
            
        #Set new analysis bounds
        self.start_of_anl = tstart
        self.end_of_anl = tend
    
        #Returning
        return

    ### Destructor Method
    ### Ensures that reference to file is closed when object is deleted
    def __del__(self):

        #Try closing the file
        try:
            self.ncfile.close()
        except:
            pass #Nothing needed if already closed

        return


