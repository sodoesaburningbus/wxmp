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
import cartopy.crs as ccrs
import cartopy.mpl.gridliner as cmg
import cartopy.mpl.ticker as cticker
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
    ###   filepath, string, full path to GFS file
    def __init__(self, filepath):
    
        #First store some basic attributes about the file
        self.name = filepath           #Full path of file
        self.date = filepath[-20:-10]  #Timestamp on file (YYYYMMDDHH)
                
        #Now open the file
        self.grib = pygrib.open(filepath) #Store GRIB2 file object as attribute
        
        #Extract message table
        string = ""
        for m in self.grib[0:]:
            string = "{}{}\n".format(string,m)
        self.mtable = string
        
        ### Go ahead and pull the lat/lon grid and convert to 2D grid
        #Exctract first and last points
        lat1 = self.grib[1].latitudeOfFirstGridPointInDegrees
        lon1 = self.grib[1].longitudeOfFirstGridPointInDegrees
        lat2 = self.grib[1].latitudeOfLastGridPointInDegrees
        lon2 = self.grib[1].longitudeOfLastGridPointInDegrees
        
        #Now store some basic grid info
        self.nx = self.grib[1].values.shape[1]
        self.ny = self.grib[1].values.shape[0]
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
        
        #Returning
        return
    
    ### Method to pull the data for a particular GRIB2 message by number
    ### Inputs:
    ###   mid, integer, number of desired GRIB2 message
    ### Outputs:
    ###   numpy array containing message values
    def get_value(self, mid):
        
        #Return selected values
        return self.grib[mid].values[self.yind1:self.yind2, self.xind1:self.xind2]
        
    ### Method to retrieve messages by name and level
    ### Inputs: (either is optional, but at least one is required)
    ###   name=, string, optional, name of variable to retrieve
    ###   level=, integer, optional, level of variable to retrieve
    ###   values=False, boolean, optional, if only one message whether to return the values
    ###     instead of the full message. Defaults to False.
    ### Outputs:
    ###   messages, list of messages that match given criteria.
    def get_messages(self, name=None, level=None, values=False):
        
        #First test that appropriate inputs are given
        if ((name == None) and (level == None)):
            print("ERROR: Must pass at one of name or level.")
            raise ValueError("Must pass at least name or level.")
        elif ((name != None) and (level == None)):
            return self.grib.select(name=name)
        elif ((name == None) and (level != None)):
            return self.grib.select(level=level)
        elif values:
            return self.grib.select(name=name, level=level)[0].values[self.yind1:self.yind2, self.xind1:self.xind2]
        else:
            return self.grib.select(name=name, level=level)
            
    
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
    
    ### Method to plot analysis region
    ### Inputs:
    ###   spath=, string, optional, place to save plot. Displays if not given.
    def plot_anl_region(self, spath=None):
        
        ### Plot subset region
        #Create figure and axis objects
        fig, ax = pp.subplots(nrows=1, ncols=1, subplot_kw={"projection":self.proj})
        
        #Add continents
        ax.coastlines()
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
    
        #Showing or saving plot
        if (spath == None):
            pp.show()
        else:
            pp.savefig(spath+"/gfs__analysis_domain.png")
    
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
        
    
    #------ Methods below this line are primarily for internal use ------#
    
    ### Method to close file on object destruction
    def __del__(self):
        try:
            self.grib.close()
        except:
            pass

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
    def find_point(self, lon, lat):
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
