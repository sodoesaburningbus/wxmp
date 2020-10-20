#This module contains classes for visualizing the various hrrr products.
#Currently, three grids are handled, 0.25, 0.5, and 1.0 degree grids.
#Written by Christopher Phillips, February 2020
#University of Alabama in Huntsville
#Atmospheric and Earth Science Department
#
#Requirements:
#Python 3+
#Cartopy
#Matplotlib
#Numpy
#PyGRIB

#Importing modules
import atmos.thermo as at
import cartopy.crs as ccrs
import cartopy.mpl.gridliner as cmg
import cartopy.mpl.ticker as cticker
import cartopy.feature as cfeature
import glob
import matplotlib.pyplot as pp
import matplotlib.ticker as mticker
import numpy
import pygrib

############################################################
#---------------------     HRRRANL     ---------------------#
############################################################
### This class contains tools for reading hrrr GRIB2 files
### It enables the user to pull variables easily, and plot
### common products such as soundings, and meteograms.
class HRRRANL:

    ### Method to construct hrrr analysis object
    ### Inputs:
    ###  hrrrpath, string, full path to directory containing hrrr files
    ###  hrrrpre, string, optional, prefix to hrrr files, defaults to none
    def __init__(self, hrrrpath, hrrrpre=None):
    
        ### Setup hrrr file prefix
        ### If statement supports wxmpi compatibility
        if (hrrrpre == None):
            hrrrpre = ""
    
        ### Locate all files
        self.files = sorted(glob.glob(hrrrpath+"/{}*".format(hrrrpre)))
                        
        #Extract message table from first file
        grib = pygrib.open(self.files[0])
        string = ""
        for m in grib[0:]:
            string = "{}{}\n".format(string,m)
        self.mtable = string
                
        ### Go ahead and pull the lat/lon grid and convert to 2D grid
        #Exctract first and last points
        self.lons = numpy.reshape(grib[1].longitudes, grib[1].values.shape)-360.0 #West is negative
        self.lats = numpy.reshape(grib[1].latitudes, grib[1].values.shape)

        #Get start time of dataset
        self.start_of_sim = grib[1].validDate
        
        #Now store some basic grid info
        self.nx = grib[1].Nx
        self.ny = grib[1].Ny
        self.dx = grib[1].DxInMetres
        self.dy = grib[1].DyInMetres
        self.full_extent = [self.lons.min(), self.lons.max(), self.lats.min(), self.lats.max()]

        #Set subsetting to off
        self.disable_subset()
        
        #Attach projection information to object
        self.proj_name = "Lambert COnformal"
        self.center_lon = (self.lons.min()+self.lons.max())/2
        self.proj = ccrs.LambertConformal(central_longitude=self.center_lon)
        self.pcp = ccrs.PlateCarree() #This is for transformations within the object
        self.states = cfeature.NaturalEarthFeature(category="cultural",
            name="admin_1_states_provinces_lines",scale="110m", facecolor="none")
        
        ### Build dictionary of variables and levels keyed to type of level
        #Find all types of levels and initialize dictionaries of lists with those keys         
        self.type_of_level_list = list(dict.fromkeys(list(m.typeOfLevel for m in grib[0:])))
        self.var_list = {}
        self.level_list = {}
        for tl in self.type_of_level_list:
            self.var_list[tl] = []
            self.level_list[tl] = []
            
        #Now fill those lists by level type (no duplicates)  
        for m in grib[0:]:
            if m.name not in self.var_list[m.typeOfLevel]:
                self.var_list[m.typeOfLevel].append(m.name)
            if m.level not in self.level_list[m.typeOfLevel]:
                self.level_list[m.typeOfLevel].append(m.level)
                    
        #Close grib file
        grib.close()
        
        #Get end time of dataset
        grib = pygrib.open(self.files[-1])
        self.end_of_sim = grib[1].validDate
        grib.close()
        
        #Set period of analysis to simulation period
        self.set_anl_period(self.start_of_sim, self.end_of_sim)
        
        #Set current working grid
        #Primarily for wxmpi compatibility
        self.cg = 1
        self.ngrids = 1
        
        #Returning
        return
    
    ### Method to disable sub-setting behavior (Off by default)
    ### Doesn't truly disable sub-setting, just sets it to the full domain
    def disable_subset(self):
    
        #Set subsetting to False
        self.subset = False
        
        #Set analysis region to full hrrr domain
        self.extent = self.full_extent
        self.xind1 = 0
        self.xind2 = self.nx
        self.yind1 = 0
        self.yind2 = self.ny
        self.rlats = self.lats[:]
        self.rlons = self.lons[:]
        
        #Returning
        return
    
    ### Method to pull the data for a particular GRIB2 message by number
    ### Inputs:
    ###   mid, integer, number of desired GRIB2 message
    ###   period, tuple of date objects, optional, start and end times of desired analysis
    ### Outputs:
    ###   numpy array containing values indexed as (Time, Y, X)
    def get_message(self, mid, period=None):
        
        #Determine time period to retreive variables over
        if (period == None):
            tstart = self.start_of_anl
            tend = self.end_of_anl
        else:
            tstart = period[0]
            tend = period[1]
        
        #Loop over each file
        vars = []
        for f in self.files:
            #Open grib file
            grib = pygrib.open(f)
            
            #Grab file date
            date = grib[1].validDate
            
            #Skip file if outside analysis period
            if ((date < tstart) or (date > tend)):
                grib.close()
                continue
            
            #Grab variable and close file
            vars.append(grib[mid].values[self.yind1:self.yind2, self.xind1:self.xind2])
            grib.close()
            
        #Return values as numpy array
        return numpy.array(vars)
    

    ### Method to retrieve closest grid point to desired location
    ### Inputs:
    ###   point, tuple of floats, (lon, lat) or list of tuples of such points.
    ###
    ### Outputs:
    ###   (xi, yj), tuple of ints or list of ints, index of grid point closest to given point
    ###     xi corresponds to lon; yj corresponds to lat.
    def get_point(self, point, gcoord=True):
        #Force lon and lats into lists.
        #This enables program consistency and user convenience.
        if not isinstance(point, list):
            point = [point]

        #Lists to hold grid points closest to user points
        xi = []
        yj = []

        #Now for each lon/lat pair find the closest domain point.
        for (ulon, ulat) in point:
            #Transform point to hrrr projection
            (x, y) = self.proj.transform_point(ulon, ulat, self.pcp)
            (x0, y0) = self.proj.transform_point(self.lons[0,0], self.lats[0,0], self.pcp)

            print(x0, x)
            print(y0, y)

            #Check that desired location is within model grid
            if ((ulon < self.extent[0]) or (ulon > self.extent[1]) or (ulat < self.extent[2]) or (ulat > self.extent[3])):
                raise ValueError("Point Lon: {}, Lat: {} is outside the model grid.".format(ulon,ulat))

            #Calculate index of each point
            if gcoord:
                xi.append(int(round(abs(x-x0)/self.dx)))
                yj.append(int(round(abs(y-y0)/self.dy)))
            else:
                xi.append(int(round(abs(x-x0)/self.dx))-1)
                yj.append(int(round(abs(y-y0)/self.dy))-1)

        #Return indices to user (only return lists if user gave list)
        if (len(xi) > 1):
            return (xi, yj)
        else:
            return (xi[0], yj[0])
    
    ### Method to subset data
    ### Generates a mask that can be applied to any variable in the GRIB file
    ### Inputs:
    ###  extent, list of floats, [west lon, east lon, south lat, north lat]
    ###
    ### Outputs,
    ###   subset, tuple of ints (y1, y2, x1, x2), array indices corresponding to subsetted region.
    def get_subset(self, extent):
                    
        #Now verify that coordinates are in increasing order.
        if (extent[0] > extent[1]): #Longitudes
            dummy = extent[0]
            extent[0] = extent[1]
            extent[1] = dummy
        if (extent[2] > extent[3]): #Latitudes
            dummy = extent[2]
            extent[2] = extent[3]
            extent[3] = dummy
        
        #Get indice of boundaries
        print(extent)
        lon1_ind, lat1_ind = self.get_point((extent[0], extent[2]))
        lon2_ind, lat2_ind = self.get_point((extent[1], extent[3]))
                
        #Return indices (Lats firsts because arrays are ordered [lat, lon])
        return (lat1_ind, lat2_ind, lon1_ind, lon2_ind)
    
    ### Method to retreive sounding
    ### Inputs:
    ###   point, tuple of floats, optional, (lon, lat) or list of tuples of such points.
    ###     if absent, the mean sounding for the analysis region is pulled.
    ###   period, tuple of date objects, optional, start and end times of desired analysis Defaults to working analysis period.
    ###   filedate, datetime object, optional, date of file to pull. (Only retreives one file's data if set). Defaults to None.
    ###
    ### Outputs:
    ###   sounding, dictionary of lists containing sounding info, or list of such dictionaries with len(points)
    ###     dictionaries are keyed ["temp", "pres", "dewp", "uwind", "vwind"] for temperature (K), pressure (hPa),
    ###     dewpoint (K), zonal wind speed (m/s), and meriodinal wind speed (m/s) respectively.
    ###     Arrays are (Time, Level) with lowest level first.
    def get_sounding(self, point=None, period=None, filedate=None):
        
        #Variable name list
        var_names = ["Temperature", "U component of wind", "V component of wind", "Relative humidity"]
        dict_keys = ["temp", "uwind", "vwind", "dewp", "pres"]

        #Create list to hold soundings
        sounding = []

        #Loop over each point
        if (point != None):
            #Force lon and lats into lists.
            #This enables program consistency and user convenience.
            if not isinstance(point, list):
                point = [point]            
                
            for p in point:

                #Retrieve grid indices
                [xind, yind] = self.get_point(p, gcoord=False)

                #Create dictionary to hold sounding
                data = {}

                #Retrieve messages from files
                for [vn, dk] in zip(var_names, dict_keys):
                    messages = self.get_var(name=vn, typeOfLevel="isobaricInhPa")
                    
                    #Loop over time and layers
                    data[dk] = []
                    dummy = numpy.zeros(messages.shape)
                    for i in range(messages.shape[0]):
                        for j in range(messages.shape[1]):
                            dummy[i,j] = messages[i,j].values[yind, xind]
                    data[dk].append(dummy)
                                    
                #Now force everything into arrays
                for k in data.keys():
                    data[k] = numpy.atleast_2d(numpy.squeeze(numpy.array(data[k])))

                #Grab pressure levels
                data["pres"] = numpy.atleast_2d(list(m.level for m in messages[0,:]))

                #Calculate dewpoint
                data["dewp"] = at.dewpoint(at.sat_vaporpres(data["temp"])*(data["dewp"]/100))

                #Reverse levels if pressure not start at surface
                if (data["pres"][0,0] < data["pres"][0,-1]):
                    for k in data.keys():
                        data[k] = numpy.flip(data[k], axis=1)

                #Append sounding to list
                sounding.append(data)

        #Compute composite sounding over analysis region if no point given.
        else:

            #Create dictionary to hold sounding
            data = {}

            #Retrieve messages from files
            for [vn, dk] in zip(var_names, dict_keys):
                messages = self.get_var(name=vn, typeOfLevel="isobaricInhPa")
                
                #Loop over time and layers
                data[dk] = []
                dummy = numpy.zeros(messages.shape)
                for i in range(messages.shape[0]):
                    for j in range(messages.shape[1]):
                        dummy[i,j] = numpy.mean(messages[i,j].values[self.yind1:self.yind2, self.xind1:self.xind2])
                data[dk].append(dummy)

            #Now force everything into arrays
            for k in data.keys():
                data[k] = numpy.atleast_2d(numpy.squeeze(numpy.array(data[k])))

            #Grab pressure levels
            data["pres"] = numpy.atleast_2d(list(m.level for m in messages[0,:]))

            #Calculate dewpoint
            data["dewp"] = at.dewpoint(at.sat_vaporpres(data["temp"])*(data["dewp"]/100))

            #Reverse levels if pressure not start at surface
            if (data["pres"][0,0] < data["pres"][0,-1]):
                for k in data.keys():
                    data[k] = numpy.flip(data[k], axis=1)

            #Append sounding to list
            sounding.append(data)          

        #Return soundings as list
        #or as dictionary if only one.
        if (len(sounding) > 1):
            return sounding
        else:
            return sounding[0]


    ### Method to retrieve messages by name and level
    ### Inputs: (either is optional, but at least one is required)
    ###   name=, string, optional, name of variable to retrieve
    ###   level=, integer, optional, level of variable to retrieve
    ###   values=False, boolean, optional, if only one message whether to return the values
    ###     instead of the full message. Defaults to False.
    ###   period, tuple of date objects, optional, start and end times of desired analysis
    ###   filedate, datetime object, optional, date of file to pull. (Only retreives one file's data if set).
    ### Outputs:
    ###   vars, array of messages that match given criteria for each file in dataset. First index is Time.
    def get_var(self, name=None, level=None, values=False, period=None, filedate=None, **kwords):
        
        #Determine time period to retreive variables over
        if (period == None):
            tstart = self.start_of_anl
            tend = self.end_of_anl
        else:
            tstart = period[0]
            tend = period[1]
        
        #Loop over each file in hrrr dataset
        vars = []
        for f in self.files:
            #Open file
            grib = pygrib.open(f)
            
            #Grab file date
            date = grib[1].validDate
            
            #Skip file if outside analysis period
            if (((date < tstart) or (date > tend)) and (filedate == None)):
                grib.close()
                continue
            elif ((filedate != None) and (date != filedate)):
                grib.close()
                continue
        
            #First test that appropriate inputs are given
            if ((name == None) and (level == None)): #No name or level
                print("ERROR: Must pass at one of name or level.")
                raise ValueError("Must pass at least name or level.")
            elif ((name != None) and (level == None)): #Only name
                vars.append(grib.select(name=name, **kwords))
            elif ((name == None) and (level != None)): #Only level
                vars.append(grib.select(level=level, **kwords))
            elif values: #Name, level, and requesting values only
                vars.append(grib.select(name=name, level=level, **kwords)[0].values[self.yind1:self.yind2, self.xind1:self.xind2])
            else: #Name and level but returning the messages themselves
                vars.append(grib.select(name=name, level=level, **kwords))
                
            #Close file
            grib.close()
            
        #Return data
        return numpy.squeeze(numpy.array(vars))
    
    ### Method to contour variable at specific model level and time
    ### Inputs:
    ###  var, string, name of variable to map
    ###  level, integer, hrrr level to plot variable on
    ###  date, datetime object, date of file to to plot
    ###
    ### Outputs:
    ### (fig, ax) tuple with pyplot figure and axis object.
    def map_var(self, varname, level, date):
    
        #Grab selected variable
        data = self.get_var(name=varname, level=level, filedate=date, values=True)
        
        ### Make the map
        #Create figure and axis
        fig, ax = pp.subplots(nrows=1, ncols=1, subplot_kw={"projection":self.pcp})
        
        #Contour
        cont = ax.contourf(self.rlons, self.rlats, data, transform=self.pcp)
        cb = pp.colorbar(cont)
        cb.set_label(varname, fontsize=14)
        
        #Add map features
        ax.coastlines()
        ax.add_feature(cfeature.BORDERS, edgecolor="black")
        ax.add_feature(self.states, edgecolor="black")
        gl = ax.gridlines(crs=self.pcp, draw_labels=True, linewidth=1, color="black", alpha=0.6, linestyle="--")
        gl.xlabels_top = False
        gl.ylabels_right = False
        
        #Label plot
        ax.set_title("Date: {}    Level: {}".format(date, level), fontsize=14, horizontalalignment="left", loc="left")
        
        #Return figure and axis objects
        return (fig, ax)
    
    ### Method to plot analysis region
    ### Inputs:
    ###   gid=, integer, optional, hrrr grid to plot. Defaults to 1. There is only 1.
    ###   lc=, boolean, optional, Flag to include land cover imagery. Defaults to False.
    ###   points=, list of tuples, optional, tuples of lon/lat pairs.
    ###
    ### Outputs:
    ###   (fig, ax) tuple with pyplot figure and axis object.
    def plot_grid(self, gid=1, lc=False, points=[]):
        
        ### Plot subset region
        #Create figure and axis objects
        fig, ax = pp.subplots(nrows=1, ncols=1, subplot_kw={"projection":self.pcp})
        
        #Add continents and states
        ax.coastlines()
        ax.add_feature(self.states)
        ax.stock_img()
                
        #Set extent and make gridlines
        gl = ax.gridlines(crs=self.pcp, draw_labels=False, linewidth=1, linestyle=":", color="grey")
        if self.subset: #If subsetting is on
            ax.set_extent(self.extent, crs=self.pcp)
            gl.xlabels_bottom = True
            gl.ylabels_left = True
            gl.xformatter = cmg.LONGITUDE_FORMATTER
            gl.yformatter = cmg.LATITUDE_FORMATTER
            
        else: #If no analysis region is set (i.e. working with full hrrr domain)
            #Format meridions
            gl.xlocator = mticker.FixedLocator(numpy.linspace(-180, 180, 13))
            ax.set_xticks(numpy.linspace(-150, 180, 12), crs=self.pcp)
            ax.set_xticklabels(numpy.linspace(-150, 180, 12))
            ax.xaxis.set_major_formatter(cticker.LongitudeFormatter())

            #Format parallels
            gl.ylocator = mticker.FixedLocator(numpy.linspace(-90, 90, 7))
            ax.set_yticks(numpy.linspace(-90, 90, 7), crs=self.pcp)
            ax.set_yticklabels(numpy.linspace(-90, 90, 7))
            ax.yaxis.set_major_formatter(cticker.LatitudeFormatter())
        
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
            print("Requested period outside of simulation period. Returning.")
            return -1
            
        #Set new analysis bounds
        self.start_of_anl = tstart
        self.end_of_anl = tend
    
        #Returning
        return
    
    ### Method to set analysis region
    ### All data is subset to this region if set
    ### Inputs:
    ###   extent, list of floats, [west lon, east lon, south lat, north lat]
    def set_anl_region(self, extent):
        
        #First set subset flag to True. This turns on subsetting behavior for data functions.
        self.subset = True
        
        #Now store the analysis extent
        self.extent = extent
        
        #Now pull the necessary indices for subsetting
        [self.yind1, self.yind2, self.xind1, self.xind2] = self.get_subset(extent)
        
        #Create lats and lats for analysis region
        self.rlons = self.lons[self.yind1:self.yind2, self.xind1:self.xind2]
        self.rlats = self.lats[self.yind1:self.yind2, self.xind1:self.xind2]
        
        #Returning
        return