#This module contains a class for visualizing OPENDAP MERRA files
#made available by UCAR.
#Written by Christopher Phillips, February 2020
#University of Alabama in Huntsville
#Atmospheric and Earth Science Department
#
#Requirements:
#Python 3+
#Cartopy
#NetCDF4
#Numpy
#Pydap

#Importing required modules
from pydap.client import open_url
from pydap.cas.urs import setup_session
import numpy
import netCDF4 as nc
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime as dt
from datetime import timedelta as td

############################################################
#---------------------     MERANL     ---------------------#
############################################################
### This class contains tools for reading OpenDAP MERRA 2
### files. It features variable retrieval for 2D, 3D, and 1D
### fields with automatic spatial and temporal subsetting.
class MERANL:

    ### The constructor method
    ### Inputs:
    ###  username, string, username for UCAR MERRA2 OpenDAP server
    ###  password, string, password for UCAR MERRA2 OpenDAP server
    ###  start_date, datetime object, start time of analysis
    ###  end_date, datetime object, end time of analysis
    ###  dataset_name, string, type of mERRA files to be retreived.
    ###    options are: "Atmosphere", "Aerosol", "Ocean"
    def __init__(self, username, password, start_date, end_date, dataset_name):

        #Store username and password
        self.username = username
        self.password = password

        #Set analysis period
        self.set_anl_period(start_date, end_date)

        #Attach base url of dataset
        self.base_url = "https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/"

        #Create dictionary of available datasets
        self.available_datasets = {"atmosphere":("M2T1NXSLV.5.12.4", "slv"),
                                   "aerosol":("M2T1NXAER.5.12.4", "aer"),
                                   "ocean":("M2T1NXOCN.5.12.4", "ocn")}

        #Load files
        self.load_dataset(dataset_name)
        
        #Set spatial subsetting to off
        self.disable_subset()
        
        #Set current working grid
        #Primarily for wxmpi compatibility
        self.cg = 1
        self.ngrids = 1

        ### Returning
        return

    ### Method to disable sub-setting behavior (Off by default)
    ### Doesn't truly disable sub-setting, just sets it to the full domain
    def disable_subset(self):
    
        #Set subsetting to False
        self.subset = False
        
        #Set analysis region to full GFS domain
        self.extent = self.full_extent
        self.xind1 = 0
        self.xind2 = self.nx
        self.yind1 = 0
        self.yind2 = self.ny
        self.rlats = self.lats[:]
        self.rlons = self.lons[:]
        
        #Returning
        return

    ### Method to load the dataset
    ### This method pulls grid information, sets available times, etc.
    ### Inputs:
    ###  dataset_name, string, type of mERRA files to be retreived.
    ###    options are: "Atmosphere", "Aerosol", "Land", "Ocean"
    def load_dataset(self, dataset_name):

        #For dataset into lowercase and attach to object
        self.dataset_name = dataset_name.lower()

        #Start session
        session = setup_session(self.username, self.password,
            check_url=self.base_url+self.available_datasets[self.dataset_name][0])

        #Now load the appropriate dataset
        self.dataset = []
        #Loop over dates
        for d in self.dates:

            #Loop over file versions (there were changes over the years; start with most recent)
            for i in [4, 3, 2, 1]:
                #Construct url
                url = ("{0}{1}/{2}/{3:02d}/MERRA2_{4}00.tavg1_2d_{5}_Nx.{2}{3:02d}{6:02d}.nc4"
                    .format(self.base_url, self.available_datasets[self.dataset_name][0],
                    d.year, d.month, i, self.available_datasets[self.dataset_name][1], d.day))
                
                #Load data
                try:                
                    self.dataset.append(open_url(url, session=session))
                    #data = open_url(url, session=session)
                except:
                    continue

                #Break to next date on successful load
                break

        #Store variables
        self.var_list = sorted(list(self.dataset[0].keys()))

        #Now pull grid info
        lons = numpy.array(self.dataset[0]["lon"][:])
        lats = numpy.array(self.dataset[0]["lat"][:])
        self.lons, self.lats = numpy.meshgrid(lons, lats)
        self.full_extent = [numpy.min(self.lons), numpy.max(self.lons),
            numpy.min(self.lats), numpy.max(self.lats)]
        self.nx = len(self.lons)
        self.ny = len(self.lats)
        self.dx = float(self.dataset[0].attributes["HDF5_GLOBAL"]["LongitudeResolution"])
        self.dy = float(self.dataset[0].attributes["HDF5_GLOBAL"]["LatitudeResolution"])

        #Now create projection information using first file in the dataset
        self.proj_name = "Global Latitude/Longitude"
        self.center_lon = (self.full_extent[0]+self.full_extent[1])/2
        self.proj = ccrs.PlateCarree(central_longitude=self.center_lon)
        self.pcp = ccrs.PlateCarree() #This is for transformations within the object
        self.states = cfeature.NaturalEarthFeature(category="cultural",
            name="admin_1_states_provinces_lines",scale="110m", facecolor="none")

    ### Method to retrieve closest grid point to desired location
    ### Inputs:
    ###   point, tuple of floats, (lon, lat) or list of tuples of such points.
    ###   gcoord=, optional, boolean, dpes nothing. Here for backwards compatibility.
    ###
    ### Outputs:
    ###   (xi, yj), tuple of ints or list of ints, index of grid point closest to given point
    ###     xi corresponds to lon; yj corresponds to lat.
    def get_point(self, point, gcoord=None):
        #Force lon and lats into lists.
        #This enables program consistency and user convenience.
        if not isinstance(point, list):
            point = [point]

        #Lists to hold grid points closest to user points
        xi = []
        yj = []

        #Now for each lon/lat pair find the closest domain point.
        for (x, y) in point:

            #Check that desired location is within model grid
            if ((x < self.extent[0]) or (x > self.extent[1]) or (y < self.extent[2]) or (y > self.extent[3])):
                print("Point Lon: {}, Lat: {}) is outside the model grid. Returning...".format(ulon,ulat))
                return None

            #Calculate index of each point
            xi.append(int(round((x-self.extent[0])/self.dx)))
            yj.append(int(round((y-self.extent[2])/self.dy)))

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
                
        #Calculate coordinates of each bounding side
        lon1_ind = int((extent[0]-self.lons[0,0])/self.dx)
        lon2_ind = int((extent[1]-self.lons[0,0])/self.dx)
        lat1_ind = int((extent[2]-self.lats[0,0])/self.dy) #Do south first because MERRA stores south first.
        lat2_ind = int((extent[3]-self.lats[0,0])/self.dy)
        
        #Return indices (Lats firsts because arrays are ordered [lat, lon])
        return (lat1_ind, lat2_ind, lon1_ind, lon2_ind)

    ### Method to retrieve variables by name
    ### Inputs:
    ###   name, string, name of variable to retrieve
    ###   period, tuple of date objects, optional, start and end times of desired analysis
    ###   filedate, datetime object, optional, date of file to pull. (Only retreives one file's data if set).
    ### Outputs:
    ###   var, array containing variable. Indices are: (Time, Y, X)
    def get_var(self, name, period=None, filedate=None, **kwords):
        
        #Determine time period to retreive variables over
        if (period == None):
            tstart = self.start_of_anl
            tend = self.end_of_anl
        else:
            tstart = period[0]
            tend = period[1]

        #Test if a specific time has been requested
        #Set requested hour if so
        if (filedate != None):
            thour = int(filedate.hour)
     
        #Loop over each file in dataset
        var = []
        for [f, date] in zip(self.dataset, self.dates):
            
            #Skip file if outside analysis period
            if (((date < tstart) or (date > tend)) and (filedate == None)):
                continue
            elif ((filedate != None) and (date != dt.strptime(filedate.strftime("%Y%m%d"), "%Y%m%d"))):
                continue
        
            #First test that appropriate inputs are given
            if (filedate != None):
                var.append(numpy.array(f[name][thour, self.yind1:self.yind2, self.xind1:self.xind2]))
            else:
                var.append(numpy.array(f[name][:, self.yind1:self.yind2, self.xind1:self.xind2]))

        #Convert list to array
        var = numpy.array(var)

        #Re-arrange array if extra time dimension
        #This arises from each file containing a whole day.
        if (len(var.shape) == 4):
            dummy = list(var[i,:,:,:] for i in range(var.shape[0]))
            var = numpy.concatenate(dummy, axis=0)

        #Return data
        return numpy.squeeze(var)

    ### Method to set analysis period
    ### Method to set working time period
    ### Inputs:
    ###   tstart, datetime object, start of simulation period to be analyzed
    ###   tend, datetime object, end of simulation period to be analyzed
    def set_anl_period(self, tstart, tend):
                    
        #Set new analysis bounds
        self.start_of_anl = tstart
        self.end_of_anl = tend

        #Create datetime objects for each day within bounds
        duration = (self.end_of_anl-self.start_of_anl)
        duration = int(numpy.ceil((duration.days*24+duration.seconds*3600)/24))
        self.dates = list(self.start_of_anl+td(days=i) for i in range(duration+1))
    
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
