#This module contains classes for visualizing the various GFS products.
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
#---------------------     GFSANL     ---------------------#
############################################################
### This class contains tools for reading GFS GRIB2 files
### It enables the user to pull variables easily, and plot
### common products such as soundings, and meteograms.
class GFSANL:

    ### Method to construct GFS analysis object
    ### Inputs:
    ###  gfspath, string, full path to directory containing GFS files
    ###  gfspre, string, optional, prefix to GFS files, defaults to none
    def __init__(self, gfspath, gfspre=None):
    
        ### Setup GFS file prefix
        ### If statement supports wxmpi compatibility
        if (gfspre == None):
            gfspre = ""
    
        ### Locate all files
        self.files = sorted(glob.glob(gfspath+"/{}*".format(gfspre)))
                        
        #Extract message table from first file
        grib = pygrib.open(self.files[0])
        string = ""
        for m in grib[0:]:
            string = "{}{}\n".format(string,m)
        self.mtable = string
        
        ### Go ahead and pull the lat/lon grid and convert to 2D grid
        #Exctract first and last points
        lat1 = grib[1].latitudeOfFirstGridPointInDegrees
        lon1 = grib[1].longitudeOfFirstGridPointInDegrees
        lat2 = grib[1].latitudeOfLastGridPointInDegrees
        lon2 = grib[1].longitudeOfLastGridPointInDegrees
        
        #Get start time of dataset
        self.start_of_sim = grib[1].validDate
        
        #Now store some basic grid info
        self.nx = grib[1].values.shape[1]
        self.ny = grib[1].values.shape[0]
        self.res = (lon2-lon1)/(self.nx-1)
        
        #Build up vectors for grid
        lats = numpy.linspace(lat1, lat2, self.ny)
        lons = numpy.linspace(lon1, lon2, self.nx)
                
        #Finally build and store 2D lat/lon arrays as attributes
        self.lons, self.lats = numpy.meshgrid(lons, lats)
        
        #Set subsetting to off
        self.disable_subset()
        
        #Attach projection information to object
        self.proj_name = "Global Latitude/Longitude"
        self.center_lon = (self.lons.min()+self.lons.max())/2
        self.proj = ccrs.PlateCarree(central_longitude=self.center_lon)
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
        
        #Set analysis region to full GFS domain
        self.extent = [self.lons.min(), self.lons.max(), self.lats.min(), self.lats.max()]
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
            #Transform point to GFS projection
            (x, y) = self.proj.transform_point(ulon, ulat, self.pcp)

            #Check that desired location is within model grid
            if ((x < self.extent[0]) or (x > self.extent[1]) or (y < self.extent[2]) or (y > self.extent[3])):
                print("Point Lon: {}, Lat: {}) is outside the model grid. Returning...".format(ulon,ulat))
                return None

            #Calculate index of each point
            if gcoord:
                xi.append(int(round((x-self.extent[0])/self.res)))
                yj.append(int(round((self.extent[3]-y)/self.res)))
            else:
                xi.append(int(round((x-self.extent[0])/self.res))-1)
                yj.append(int(round((self.extent[3]-y)/self.res))-1)

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
    ###   subset, boolean array, subset indices corresponding to subsetted region.
    def get_subset(self, extent):
        
        #First add 360 to any longitudes that are negative (b/c GFS goes from [0, 360])
        if (extent[0] < 0):
            extent[0] = extent[0] + 360
        if (extent[1] < 0):
            extent[1] = extent[1] + 360
            
        #Now verify that coordinates are in increasing order.
        if (extent[0] > extent[1]): #Longitudes
            dummy = extent[0]
            extent[0] = extent[1]
            extent[1] = dummy
        if (extent[2] > extent[3]): #Latitudes
            dummy = extent[2]
            extent[2] = extent[3]
            extent[3] = dummy
                
        #Calculate coordinates of each bounding side
        lon1_ind = int((extent[0]-self.lons[0,0])/self.res)
        lon2_ind = int((extent[1]-self.lons[0,0])/self.res)
        lat1_ind = int((self.lats[0,0]-extent[3])/self.res) #Do north first because GFS stores north first.
        lat2_ind = int((self.lats[0,0]-extent[2])/self.res)
        
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
                print(point)
                print(p)
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
        
        #Loop over each file in GFS dataset
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
    ###  level, integer, GFS level to plot variable on
    ###  date, datetime object, date of file to to plot
    ###
    ### Outputs:
    ### (fig, ax) tuple with pyplot figure and axis object.
    def map_var(self, varname, level, date):
    
        #Grab selected variable
        data = self.get_var(name=varname, level=level, filedate=date, values=True)
        
        ### Make the map
        #Create figure and axis
        fig, ax = pp.subplots(nrows=1, ncols=1, subplot_kw={"projection":self.proj})
        
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
    ###   gid=, integer, optional, GFS grid to plot. Defaults to 1. There is only 1.
    ###   lc=, boolean, optional, Flag to include land cover imagery. Defaults to False.
    ###   points=, list of tuples, optional, tuples of lon/lat pairs.
    ###
    ### Outputs:
    ###   (fig, ax) tuple with pyplot figure and axis object.
    def plot_grid(self, gid=1, lc=False, points=[]):
        
        ### Plot subset region
        #Create figure and axis objects
        fig, ax = pp.subplots(nrows=1, ncols=1, subplot_kw={"projection":self.proj})
        
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
            
        else: #If no analysis region is set (i.e. working with full GFS domain)
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
        
    
#############################################################
#---------------------     GFSgrid     ---------------------#
#############################################################
### This class contains tools for using GFS grids.
### It allows for easy visualization of a domain, and mapping user-defined
### points to the grid.
### Inputs:
###   gid, string, the GFS grid you want to visualize ("1.0", "0.5", "0.25")
###
### Attributes:
###  self.lons - GFS longitudes
###  self.slons - longitude grid clipped to [-180,180]
###  self.lats - GFS latitudes
###  self.proj_name - Name of GFS projection
###  self.proj - Cartopy projection object for GFS grid
###  self.pcp - Cartopy PlateCarree projection object (for internal use)
###  self.center_lon - central longitude of GFS projection
###  self.extent - Extent of GFS domain in GFS projection coordinates
###  self.nx - Number of x grid points
###  self.ny - Number of y grid points
###  self.dx - x grid spacing
###  self.dy - y grid spacing
###  self.lat1 - First latitude on grid
###  self.lat2 - Last latitude on grid
###  self.lon1 - First longitude on grid
###  self.lon2 - Last longitude on grid
###
### Methods:
###  __init__(self, gid) - Object constructor
###  g1p0(self) - Provides data for GFS 1.0 deg grid
###  g0p5(self) - Provides data for GFS 0.5 deg grid
###  g1p25(self) - Provides data for GFS 0.25 deg grid
###  find_point(self, lon, lat) - Finds grid index corresponding to given lon/lat
###  plot_grid(self, [points], [extent]) - Plots model domain
class GFSgrid:
    ### Funtion to construct object
    def __init__(self, gid):
        #Call different function based on grid id.
        #These attach basic grid info to object.
        setup = {"1.0":self.g1p0, "0.5":self.g0p5, "0.25":self.g0p25}
        try:
            setup[gid]()
        except:
            print('That was not a valid grid. Options are: "1.0", "0.5", "0.25"')
            print("Deconstructing object.")
            del(self)

        #Calculate lon/lats of the grid.
        #These are 1D arrays forming the rows and columns of the grid.
        self.lons = numpy.linspace(self.lon1, self.lon2, self.nx)
        self.lats = numpy.linspace(self.lat1, self.lat2, self.ny)

        #Shift longitude grid for when [-180, 180] is necessary
        dummy_lons = numpy.linspace(self.lon1, self.lon2, self.nx)
        dummy_lons[dummy_lons>180] -= 360
        self.slons = dummy_lons

        #Attach projection information to object
        self.proj_name = "Global Latitude/Longitude"
        self.center_lon = (self.lon1+self.lon2)/2
        self.proj = ccrs.PlateCarree(central_longitude=self.center_lon)
        self.pcp = ccrs.PlateCarree() #This is for transformations within the object

        #Calculate full domain extent (left, right, bottom, top)
        (x0, y0) = self.proj.transform_point(self.lon1, self.lat2, self.pcp)
        (xf, yf) = self.proj.transform_point(self.lon2, self.lat1, self.pcp)
        self.extent = (x0,xf,y0,yf) #Note that this is in GFS projection coords.

        #Returning
        return

    ### Function to setup 1.0' grid.
    def g1p0(self):
        #Attach basic grid info to object
        self.nx = 360     #Number of x grid points
        self.ny = 181     #Number of y grid points
        self.dx = 1.0     #x grid spacing
        self.dy = 1.0     #y grid spacing
        self.lat1 = 90.0  #First latitude on grid
        self.lat2 = -90.0 #Last latitude on grid
        self.lon1 = 0.0   #First longitude on grid
        self.lon2 = 359.0 #Last longitude on grid

        #Returning
        return

    ### Function to setup 0.5' grid.
    def g0p5(self):
        #Attach basic grid info to object
        self.nx = 720     #Number of x grid points
        self.ny = 361     #Number of y grid points
        self.dx = 0.5     #x grid spacing
        self.dy = 0.5     #y grid spacing
        self.lat1 = 90.0  #First latitude on grid
        self.lat2 = -90.0 #Last latitude on grid
        self.lon1 = 0.0   #First longitude on grid
        self.lon2 = 359.5 #Last longitude on grid

        #Returning
        return

    ### Function to setup 0.25' grid
    def g0p25(self):
        #Attach basic grid info to object
        self.nx = 1440     #Number of x grid points
        self.ny = 721     #Number of y grid points
        self.dx = 0.25     #x grid spacing
        self.dy = 0.25     #y grid spacing
        self.lat1 = 90.0  #First latitude on grid
        self.lat2 = -90.0 #Last latitude on grid
        self.lon1 = 0.0   #First longitude on grid
        self.lon2 = 359.75 #Last longitude on grid

        #Returning
        return

    ### Funtion to retrieve closest grid point to desired location
    ### Inputs:
    ###   lon, float or list of floats, longitude of desired point
    ###   lat, float or list of floats, latitude of desired point
    ###
    ### Outputs:
    ###   (xi, yj), tuple of ints or list of ints, index of grid point closest to given point
    ###     xi corresponds to lon; yj corresponds to lat.
    def get_point(self, lon, lat):
        #Force lon and lats into lists.
        #This enables program consistency and user convenience.
        if not isinstance(lon, list):
            lon = [lon]
            lat = [lat]

        #Lists to grid points closest to user points
        xi = []
        yj = []

        #Now for each lon/lat pair find the closest domain point.
        for (ulon, ulat) in zip(lon, lat):
            #Transform point to GFS projection
            (x, y) = self.proj.transform_point(ulon, ulat, self.pcp)

            #Check that desired location is within model grid
            if ((x < self.extent[0]) or (x > self.extent[1]) or (y < self.extent[2]) or (y > self.extent[3])):
                print("Point Lon: {}, Lat: {}) is outside the model grid. Returning...".format(ulon,ulat))
                return None

            #Calculate index of each point
            xi.append(int(round((x-self.extent[0])/self.dx)))
            yj.append(int(round((self.extent[3]-y)/self.dy)))

        #Return indices to user (only return lists if user gave list)
        if (len(xi) > 1):
            return (xi, yj)
        else:
            return (xi[0], yj[0])

    ### Function to plot model domain
    ### Also allows user to determine a sub-extent and plot points
    ### Inputs:
    ###  extent=, optional, tuple of floats, (west lon, east lon, south lat, north lat)
    ###     Note that extent cannot cross the map boundary (i.e. long=0).
    ###  points=, optional, list of tuples corresponding to lon/lat pairs. (Note the order)
    def plot_grid(self, extent=None, points=[]):

        #Create figure and axis objects. Also use stock image from cartopy
        fig, ax = pp.subplots(nrows=1, ncols=1, subplot_kw={"projection":self.proj})
        ax.stock_img()

        #Add gridlines
        gl = ax.gridlines(crs=self.pcp, draw_labels=False, linewidth=1, linestyle=":", color="grey")

        if (extent == None): #Only do fancy gridline format for default plot.
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
        else:
            gl.xlabels_bottom = True
            gl.ylabels_left = True
            gl.xformatter = cmg.LONGITUDE_FORMATTER
            gl.yformatter = cmg.LATITUDE_FORMATTER

        #Add user defined points
        for pt in points:
            #Now plot the point
            ax.scatter(pt[0], pt[1], color="red", transform=self.pcp)

        #Change map extent if requested
        if (extent != None):
            #First must transform extent corners
            llc = self.proj.transform_point(extent[0], extent[2], self.pcp)
            urc = self.proj.transform_point(extent[1], extent[3], self.pcp)

            #Now can set extent
            ax.set_extent((llc[0],urc[0],llc[1],urc[1]), crs=self.proj)

        #Display image
        pp.show()

        #Returning
        return
