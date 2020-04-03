### WxMPI (Weather Modelling Package Interface)
### This object is the interface bewteen the user and the underlying
### WxMP modules (e.g. wrf.py and gfs.py). While it is possible to import
### and use those modules individually. This object provides a top-level interface
### to simplify the user experience.
### Written by Christopher Phillips
### University of Alabama in Huntsville
### Earth and Atmospheric Sciences Department.
###
### Requirements:
###   Python 3+
###   Atmos (available at github.com/sodoesaburningbus)
###   Cartopy
###   f90NML
###   Matplotlib
###   NetCDF4
###   Numpy
###   PyGRIB
###   PIL
###   PyMODIS (available at github.com/sodoesaburningbus)

### Importing required modules
import atmos.math as am
import numpy
from wxmp import *

### WxMPI object
### Inputs:
###  model, string, model name
###  filedir, string, path to directory containing files for analysis
class wxmpi:
    ### Constructor method
    ### Inputs:
    ###  model, string, model name
    ###  filedir, string, path to directory containing files for analysis
    ###  prefix, string, optional, file prefix for dataset files, default depends on selected model.
    def __init__(self, model, filedir, prefix=None):
    
        ### Select model from wxmp to use for file access
        if (model.upper() == "WRF"): #WRF model
            self.dataset = wrf.WRFANL(filedir, wrfpre=prefix)
            self.dataset_name = model.upper()
            
        elif (model.upper() == "GFS"): #GFS model
            self.dataset = gfs.GFSANL(filedir, gfspre=prefix)
            self.dataset_name = model.upper()
            
        else: #Model not supported
            raise ValueError("{} is unsupported at this time.".format(model.upper()))
                        
        ### Returning
        return

    ### Method to disable subsetting behavior
    def disable_subset(self):
    
        #Call model method
        return disable_subset()

    ### Method to retrieve closest grid point to desired location
    ### Inputs:
    ###   point, tuple of floats, (lon, lat) or list of tuples of such points.
    ###
    ### Outputs:
    ###   (xi, yj), tuple of ints or list of ints, array index of grid point closest to given point
    ###     xi corresponds to lon; yj corresponds to lat.
    def get_point(self, point, gcoord=False):
    
        #Call model method
        return self.dataset.get_point(point, gcoord=gcoord)
    
    ### Method to retrieve soundings from the data set
    ### Inputs:
    ###   point, tuple of floats, (lon, lat) or list of tuples of such points
    ###
    ### Outputs:
    ###   sounding, dictionary of lists containing sounding info, or list of such dictionaries with len(points)
    ###     dictionaries are keyed as ["temp", "pres", "mixr"] for temperature (K), pressure (hPa),
    ###     and mixing ratio (kg/kg) respectively.
    def get_sounding(self, point, **kwords):
        
        #Call model method
        return self.dataset.get_sounding(point, **kwords)

    ### Method to retrieve variables from data set as a 2D field
    ### Inputs:
    ###  var_label, string, name of variable to retreive
    ###  level, integer, pressure level of variable to retrieve (in millibars)
    ### Outputs:
    ###  var, list of 2D arrays of floats, variable at desired level for all time steps
    def get_var2D(self, var_label, level):
    
        ### Retreive variable based on model dataset
        if (self.dataset_name == "WRF"): #WRF
            #Retrieve requested variable
            vars = self.dataset.get_var([var_label, "PB", "P"])
                        
            #Calculate total pressure
            pres = (vars["PB"]+vars["P"])/100.0 #Total pressure Pa -> hPa
            plevels = numpy.mean(pres, axis=(2,3))
                        
            ### Interpolate to desired level
            ### Loop over each time step
            var = []
            for i in range(vars[var_label].shape[0]):
                                                          
                #Locate level just below desired level
                pind = numpy.where(numpy.clip(plevels[i,:]-level, a_min=None, a_max=1) < 0)[0][0]
                
                if (plevels[i,pind] < level): #Closest level is below desired level
                    var.append(am.layer_interp(plevels[i,pind], plevels[i,pind+1], level,
                        vars[var_label][i,pind,:,:], vars[var_label][i,pind+1,:,:]))
                        
                else: #Case in which desired level is present in model
                    var.append(vars[var_label][i,pind,:,:])
                                
        elif (self.dataset_name == "GFS"): #GFS
            #Retreive variables
            try: #Trivial case of level is present in model
                return self.dataset.get_var(name=var_label, level=level, values=True)
                
            except: #Extract variable at all pressure levels
                vars = self.dataset.get_var(name=var_label, typeOfLevel="isobaricInhPa")
            
            #Build list of pressure levels in model
            plevels = []
            for v in vars:
                plevels.append(v.level)
            plevels = numpy.array(plevels, dtype="float")
            
            #Ensure pressure levels are in ascending (height) order
            if (plevels[0] < plevels[1]):
                plevels = plevels[::-1]
                vars = vars[::-1]
            
            #Locate level just below desired level
            pind = numpy.where(numpy.clip(plevels-level, a_min=None, a_max=1) < 0)[-1]
            
            #Interpolate to desired level
            var = am.layer_interp(plevels[pind], plevels[pind+1], level,
                    vars[pind].values, vars[pind+1].values)       
    
        ### Return variable
        return numpy.array(var)
    
    ### Method to retrieve variables from data set as a 3D field
    ### Inputs:
    ###  var_label, string, name of variable to retreive
    ###  leveltype, string, optional, type of model level. Default varies with data set.
    ### Outputs:
    ###  var, 4D arrays of floats, variable at all model levels for all time steps.
    ###   Indexed as (Time, Z, Y, X)
    def get_var3D(self, var_label, leveltype=None):
        
        ### Handle data set differently based on model
        if (self.dataset_name == "WRF"): #WRF
            return self.dataset.get_var([var_label])[var_label]

        elif (self.dataset_name == "GFS"): #GFS
            #Check if using default model level type
            if (leveltype == None):
                leveltype = "isobaricInhPa"
            
            #Loop over each model level
            var = [] #List to hold arrays
            for level in self.dataset.level_list[leveltype]:
                var.append(self.dataset.get_var(name=var_label, level=level, values=True))

            #Convert list to numpy array
            var = numpy.swapaxes(numpy.array(var), 1, 0)

            #Returning
            return var

    ### Method to retreive a WxChallenge forecast
    ### Inputs:
    ###   point, tuple of floats, (lon, lat) of forecast location
    ###   gid, optional, integer, grid to analyze. Default is current grid.
    ###   period, optional, tuple of datetime objects, (start, end) temporal bounds of plot.
    ###     Defaults to current analysis period.
    ### Outputs:
    ###   tmax, float, forecasted maximum temperature [F]
    ###   tmin, float, forecsated minimum temperature [F]
    ###   wmax, float, forecasted maximum wind speed [kt]
    ###   pacc, float, forecasted accumulated precipitation [inch]
    def get_wxc(self, point, period=None):
        
        #Call model method
        #Not all models support this option.
        try:
            return self.dataset.wxchallenge(point, period=period)
        except Exception as err:
            print("Warning: Not all models support this option. Actual error below.")
            print(err)
            raise Exception
    
    ### Method to contour variable at specific model level and time
    ### Inputs:
    ###  var, string, name of variable to map
    ###  level, integer, model level to plot variable on
    ###  date, datetime object, date of file to to plot
    ###
    ### Outputs:
    ### (fig, ax) tuple with pyplot figure and axis object.
    def map_var(self, varname, level, date):
    
        #Call model method
        return self.dataset.map_var(varname, level, date)
    
    ### Method to print a meteogram
    ### Temperature, dewpoint, wind speed, and solar radiation are plotted for a single
    ### point in the dataset. By default, this is averaged over the whole domain.
    ### If a point is provided by the user, that is used instead.
    ### Inputs:
    ###   spath, string, location (with trailing slash) to save meteogram.
    ###   point, optional, tuple, (lon, lat) of point to be plotted. Default is average over domain.
    ###   gid, optional, integer, grid to analyze. Default is current grid.
    ###   period, optional, tuple of datetime objects, (start, end) temporal bounds of plot.
    ###     Defaults to current analysis period.
    def meteogram(self, spath, point=None, period=None):
    
        #Call model method
        #Not all models support this option.
        try:
            return self.dataset.meteogram(spath, point=None, period=None)
        except Exception as err:
            print("Warning: Not all models support this option. Actual error below.")
            print(err)
            raise Exception
    
    ### Method to plot current analysis grid
    ### Inputs:
    ###   gid=, integer, optional, GFS grid to plot. Defaults to 1. There is only 1.
    ###   lc=, boolean, optional, Flag to include land cover imagery. Defaults to False.
    ###   points=, list of tuples, optional, tuples of lon/lat pairs.
    ###
    ### Outputs:
    ###   (fig, ax) tuple with pyplot figure and axis object.
    def plot_grid(self, gid=1, lc=False, points=[]):
        
        #Call model method
        return self.dataset.plot_grid(gid=gid, lc=lc, points=points)
    
    ### Method to set working time period
    ### Inputs:
    ###   tstart, datetime object, start of simulation period to be analyzed
    ###   tend, datetime object, end of simulation period to be analyzed
    def set_anl_period(self, tstart, tend):
    
        #Call model method
        return self.dataset.set_anl_period(tstart, tend)

    ### Method to set analysis region
    ### All data is subset to this region if set
    ### Inputs:
    ###   extent, list of floats, [west lon, east lon, south lat, north lat]
    def set_anl_region(self, extent):
              
        #Call model method
        return self.dataset.set_anl_region(extent)
        
    ### Method to set the current dataset grid being examined
    ### Not supported by all datasets, notably most reanalysis
    ### Inputs:
    ###   gid, integer, grid ID to be set as current working grid
    def set_cg(self, gid):
    
        #Call model method
        #Not all models support this option.
        try:
            return self.dataset.set_cg(gid)
        except Exception as err:
            print("Warning: Not all models support this option. Actual error below.")
            print(err)
            raise Exception
