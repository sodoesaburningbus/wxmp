### This module contains classes for handling Weather Research and Forecasting model
### grids, inputs, and outputs.
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
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import datetime as dt
import glob
import f90nml
from matplotlib.patches import Rectangle
import matplotlib.pyplot as pp
import netCDF4 as nc
import numpy
from PIL import Image
from pymodis.modislc import MCD12C1
import sys

############################################################
#----------------------    WPSNL     ----------------------#
############################################################
### The WPSNL class provides tools for creating the WPS
### namelist. It first reads in a template and allows the
### user to make changes. It can also generate a matching WRF
### namelist provided a template.
###
### Available Attributes:
###  nml - F90NML namelist object for WPS namelist
###
### Available Methods:
###  create_wrfnml(template_location, save_location) - Create WRF namelist based on WPS namelist
###  create_wpsnml(save_location) - Write out WPS namelist (call set_nml first)
###  set_nml(share={}, geogrid={}, ungrib={}, metgrid={}) - Change WPS namelist 
class WPSNL:

    ### Method to construct object
    ### Inputs:
    ###   namelist, string, path to WRF namelist template
    def __init__(self, namelist):
    
        #Create namelist object
        try:
            self.nml = f90nml.read(namelist)
        except:
            print("{} not found.".format(namelist))
            sys.exit(2)
        
        #Pull in namelist information
        self.read_nml()
        
        #Returning
        return
        
    ### Method to create WRF namelist from WPS namelist
    ### This method builds a namelist using the WPS settings
    ### Users can override the defaults using optional dictionary inputs.
    ### By default grids with dx > 4km have CU parameterization turned on.
    ### Inputs:
    ###   template, string, path to WRF namelist template
    ###   output, string, path to write new WRF namelist
    ###   time_control, dictionary, optional, user-dseired changes to time_control section
    ###     of namelist. Format are {"option_name":value}. Applies to inputs below too.
    ###   domains, dictionary, optional, user-dseired changes to domains section
    ###   physics, dictionary, optional, user-dseired changes to physics section
    ###   fdda, dictionary, optional, user-dseired changes to fdda section
    ###   dynamics, dictionary, optional, user-desired changes to dynamics section
    ###   bdy_control, dictionary, optional, user-dseired changes to bdy_control section
    ###   grib2, dictionary, optional, user-dseired changes to grib2 section
    ###   namelist_quilt, dictionary, optional, user-dseired changes to namelist_quilt section
    def create_wrfnml(self, template, output, time_control={}, domains={}, physics={},
        fdda={}, dynamics={}, bdy_control={}, grib2={}, namelist_quilt={}):
    
        #Create WRF namelist object
        try:
            wrfnml = f90nml.read(template)
        except:
            print("{} not found.".format(template))
            sys.exit(2)
            
        #Update WPS namelist information
        self.read_nml()

        ### Setup new time controls
        #Start times
        wrfnml["time_control"]["start_year"] = list(int(date.strftime("%Y")) for date in self.sdate)
        wrfnml["time_control"]["start_month"] = list(int(date.strftime("%m")) for date in self.sdate)
        wrfnml["time_control"]["start_day"] = list(int(date.strftime("%d")) for date in self.sdate)
        wrfnml["time_control"]["start_hour"] = list(int(date.strftime("%H")) for date in self.sdate)
        wrfnml["time_control"]["start_minute"] = list(int(date.strftime("%M")) for date in self.sdate)
        wrfnml["time_control"]["start_second"] = list(int(date.strftime("%S")) for date in self.sdate)

        #End times
        wrfnml["time_control"]["end_year"] = list(int(date.strftime("%Y")) for date in self.edate)
        wrfnml["time_control"]["end_month"] = list(int(date.strftime("%m")) for date in self.edate)
        wrfnml["time_control"]["end_day"] = list(int(date.strftime("%d")) for date in self.edate)
        wrfnml["time_control"]["end_hour"] = list(int(date.strftime("%H")) for date in self.edate)
        wrfnml["time_control"]["end_minute"] = list(int(date.strftime("%M")) for date in self.edate)
        wrfnml["time_control"]["end_second"] = list(int(date.strftime("%S")) for date in self.edate)
        
        #Boundary file interval
        wrfnml["time_control"]["interval_seconds"] = self.interval
        
        ### Setup new domain specifications
        wrfnml["domains"]["time_step"] = int(numpy.ceil(self.dx[0]/1000*6))
        wrfnml["domains"]["max_dom"] = self.ngrids
        wrfnml["domains"]["e_we"] = list(self.nx)
        wrfnml["domains"]["e_sn"] = list(self.ny)
        wrfnml["domains"]["dx"] = list(self.dx)
        wrfnml["domains"]["dy"] = list(self.dy)
        wrfnml["domains"]["grid_id"][1:] = list(numpy.array(self.nml["geogrid"]["parent_id"], dtype="int")+1)[1:]
        wrfnml["domains"]["parent_id"] = self.nml["geogrid"]["parent_id"]
        wrfnml["domains"]["i_parent_start"] = self.nml["geogrid"]["i_parent_start"]
        wrfnml["domains"]["j_parent_start"] = self.nml["geogrid"]["j_parent_start"]
        wrfnml["domains"]["parent_grid_ratio"] = list(self.grid_ratio)
        wrfnml["domains"]["parent_time_step_ratio"] = list(self.grid_ratio) #Keep space/time constant

        ### Physics
        #Cumulus (Option 3 over US to match NCAR)
        cu = [] #list to hold CU settings
        for dx in self.dx:
            if (dx > 4000): #Limit CU to spacing > 4km
                cu.append(3)
            else:
                cu.append(0)
        wrfnml["physics"]["cu_physics"] = cu #Set new CU parameterization
        
        ### Write user-defined options
        #Setup lists with namelist options
        namelist_vars = [time_control, domains, physics, fdda, dynamics,
            bdy_control, grib2, namelist_quilt]
        namelist_names = ["time_control", "domains", "physics", "fdda", "dynamics",
            "bdy_control","grib2", "namelist_quilt"]
        
        #Set new variables
        for i in range(len(namelist_vars)):
            for k in namelist_vars[i].keys():
                wrfnml[namelist_names[i]][k] = namelist_vars[i][k]
                    
        ### Write new WRF namelist
        wrfnml.write(output, force=True)

        #Deleting WRF namelist object after writing
        del(wrfnml)
            
        #Returning
        return
    
    ### Method to write WPS namelist
    ### Inputs:
    ###  output, string, location to save new namelist
    def create_wpsnml(self, output):
    
        #Write namelist
        self.nml.write(output, force=True)
        
        #Returning
        return
    
    ### Method to change WPS namelist
    ### Inputs:
    ###   share, dictionary, contains WPS namelist variables within the share group
    ###     must be keyed with variable name corresponding to namelist, e.g. "wrf_core"
    ###   geogrid, dictionary, contains WPS namelist variables within the geogrid group
    ###   ungrib, dictionary, contains WPS namelist variables within the ungrib group
    ###   metgrid, dictionary, contains WPS namelist variables within the metgrid group
    def set_nml(self, share={}, geogrid={}, ungrib={}, metgrid={}):
        
        #Set new share variables
        for k in share.keys():
            self.nml["share"][k] = share[k]
            
        #Set new geogrid variables
        for k in geogrid.keys():
            self.nml["geogrid"][k] = geogrid[k]
            
        #Set new ungrib variables
        for k in ungrib.keys():
            self.nml["ungrib"][k] = ungrib[k]
            
        #Set new metgrid variables
        for k in metgrid.keys():
            self.nml["metgrid"][k] = metgrid[k]
            
        #Returning
        return    

        
    #------ Methods below this line are primarily for internal use ------#
    
    ### Method to read in WPS namelist attributes
    ### and derive further grid information
    def read_nml(self):

        ### Pull time information
        self.sdate = list(dt.datetime.strptime(date, "%Y-%m-%d_%H:%M:%S") for date in self.nml["share"]["start_date"])
        self.edate = list(dt.datetime.strptime(date, "%Y-%m-%d_%H:%M:%S") for date in self.nml["share"]["end_date"])
        self.interval = self.nml["share"]["interval_seconds"] #Boundary file interval

        ### Pull map information
        self.lat0 = self.nml["geogrid"]["ref_lat"]       #Domain center
        self.lon0 = self.nml["geogrid"]["ref_lon"]       #Domain center
        self.slat1 = self.nml["geogrid"]["truelat1"]     #Standard parallel 1
        self.slat2 = self.nml["geogrid"]["truelat2"]     #Standard parallel 2
        self.slon1 = self.nml["geogrid"]["stand_lon"]    #Standard longitude
        self.proj_name = self.nml["geogrid"]["map_proj"] #Map projection name
                
        ### Pull WRF Grid information                
        
        #Number of grids in domain
        self.ngrids = self.nml["share"]["max_dom"]
        
        #Number of points in x and y directions
        self.nx = numpy.array(self.nml["geogrid"]["e_we"])[:self.ngrids]
        self.ny = numpy.array(self.nml["geogrid"]["e_sn"])[:self.ngrids]
        self.dx = numpy.ones(self.ngrids)*self.nml["geogrid"]["dx"]
        self.dy = numpy.ones(self.ngrids)*self.nml["geogrid"]["dy"]

        #Ratio of each grid to its parent and corresponding grid spacings.
        self.grid_ratio = numpy.array(self.nml["geogrid"]["parent_grid_ratio"])[:self.ngrids]
        
        #Compute grid spacing
        ratio = self.grid_ratio[0]
        for i in range(1, len(self.grid_ratio)):
            ratio *= self.grid_ratio[i]
            self.dx[i] /= ratio
            self.dy[i] /= ratio
            
        #Returning
        return

        
############################################################
#---------------------     WRFANL     ---------------------#
############################################################
### The WRFANL class provides tools for visualizing WRF
### model output. It rquires only the directory where the
### output is stored. It provides tools for quick creation
### of common products such as radar loops and provides file-level
### access to all variables.
### Further, it enables the setting of a working grid and analysis period
### to quickly subset data to desired period and domain.
###
### NOTE: This object only supports one data frame per file.
###
### Available Attributes:
###  ngrids - Number of grids in simulation
###  files - List of list of files, separated by grid
###  proj_name - Name of grid map projection
###  lat0 - Center latitude of simulation
###  lon0 - Center longitude of simulation
###  slat1 - Standard latitude 1 for map projection
###  slat2 - Standard latitude 2 for map projection
###  proj - Cartopy projection object for grid map projection
###  pcp - Cartopy PlateCarree projection object
###  tformat - Format of WRF output file times
###  start_of_sim - Datetime date object with simulation start time
###  end_of_sim - Datetime date object with simulation end time
###  start_of_anl - Datetime date object with start of analysis period
###  end_of_anl - Datetime date object with end of analysis period
###  cg - Current working grid
###
### Available Methods:
###  get_point((lon, lat)) - Returns grid coordinates of (lon, lat) pair
###  get_var(names) - Returns dictionary containing requested variables
###  get_wxc((lon, lat)) - Returns WxChallenge forecast for desired point
###  meteogram(save_location) - Plots a meteogram from simulation output
###  radar_loop(save_location) - Creates a radar loop (gif format)
###  set_cg(grid_number) - Set current working grid
###  set_period(tstart, tend) - Set beginning and end of analysis period
class WRFANL:

    ### Method to construct object
    ### Inputs:
    ###  wrfpath, string, path to WRF output files
    ###  wrfpre, string, optional, prefix to WRF files, defaults to "wrfout_"
    def __init__(self, wrfpath, wrfpre=None):
        ### Setup WRF file prefix
        ### If statement supports wxmpi compatibility
        if (wrfpre == None):
            wrfpre = "wrfout_"
        
        ### Locate all of the output files for each domain
        #Find all files and sort them
        all_files = sorted(glob.glob(wrfpath+"/{}*".format(wrfpre)))

        #Determine number of domains
        dfile = nc.Dataset(all_files[-1], "r")
        self.ngrids = dfile.__dict__["GRID_ID"]
        
        #Seperate list of files into a list for each domain
        self.files = [] #List to hold file lists for each domain
        for gid in range(1, self.ngrids+1): #Loop over grids
            dummy = [] #Dummy list to hold files for each domain
            for f in all_files: #Loop over files
                if (int(f[f.rfind("_d")+2:f.rfind("_d")+4]) == gid):
                    dummy.append(f) #Store files in dummy list if the match domain
            self.files.append(dummy) #Store domain list
            
        ### Determine model map information
        #Create dictionary with possible projections
        proj_list = {"Mercator":self.mercator, "Lambert Conformal":self.lambert,
            "Polar Stereographic":self.polar, "Lat-Lon":self.latlon}
        
        #Pull basic map info
        self.proj_name = dfile.__dict__["MAP_PROJ_CHAR"]
        self.lat0 = dfile.__dict__["MOAD_CEN_LAT"]
        self.lon0 = dfile.__dict__["STAND_LON"]
        self.slat1 = dfile.__dict__["TRUELAT1"]
        self.slat2 = dfile.__dict__["TRUELAT2"]   
        self.proj = proj_list[self.proj_name]()
        self.pcp = ccrs.PlateCarree()
        self.states = cfeature.NaturalEarthFeature(category="cultural",
            name="admin_1_states_provinces_lines",scale="110m", facecolor="none")

        #Pull simulation start and end times
        self.tformat = "%Y-%m-%d_%H:%M:%S"
        time = dfile.__dict__["SIMULATION_START_DATE"]
        self.start_of_sim = dt.datetime.strptime(time, self.tformat)        
        time = "".join(numpy.array(dfile.variables["Times"][0,:], dtype="str"))
        self.end_of_sim = dt.datetime.strptime(time, self.tformat)
        
        ### Construct list of available variables
        #Loop across all variables, pulling the 3D fields' level types
        self.type_of_level_list = []
        for k in dfile.variables.keys():
            if ((len(dfile.variables[k].dimensions) == 4) and
                (dfile.variables[k].dimensions[1] not in self.type_of_level_list)):
                self.type_of_level_list.append(dfile.variables[k].dimensions[1])
        
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
        for k in dfile.variables.keys():
            if (len(dfile.variables[k].dimensions) == 4):
                self.var_list[dfile.variables[k].dimensions[1]].append(dfile.variables[k].name)
                self.level_list[dfile.variables[k].dimensions[1]].append(numpy.arange(0,dfile.variables[k].shape[1]))
            elif (len(dfile.variables[k].dimensions) == 3):
                self.var_list["2D"].append(dfile.variables[k].name)
                self.level_list["2D"].append(1)
            else:
                self.var_list["1D"].append(dfile.variables[k].name)
                self.level_list["1D"].append(1)
        
        #Pull grid information for each domain
        self.get_grids()
        
        ### Set the current working grid to grid 1 (the outermost grid)
        self.cg = 1
        
        #Set period of analysis to simulation period
        self.set_anl_period(self.start_of_sim, self.end_of_sim)
        
        ### Grab latitude and longitudes
        self.subset = False #Necessary here to disable subsetting check in get_var
        lonlats = self.get_var(["XLONG", "XLAT"])
        self.lons = lonlats["XLONG"][0]
        self.lats = lonlats["XLAT"][0]
        
        ### Destroy now unnecessary dummy file object
        dfile.close()
                                
        ### Disable subsetting by default
        self.disable_subset()
        
        ### Returning
        return
    
    ### Method to disable sub-setting behavior (Off by default)
    def disable_subset(self):
    
        #Set subsetting to False
        self.subset = False
        
        #Set extent to none
        self.anl_extent = None
        
        #Reset regional lat/lons
        self.rlons = self.lons
        self.rlats = self.lats
                
        #Returning
        return
    
    ### Method to locate nearest point on the WRF grid given (lon, lat).
    ### Note that the mass (i.e. unstaggered) grid coordinates are used.
    ### Inputs:
    ###   point, tuple of floats, (lon, lat)
    ###   gid, optional, integer, grid to find point on, defautls to current grid
    ###   gcoord, optional, boolean, flag to use grid coordinates (default) or array indices.
    ###     array indices are grid coordinates minus one.
    ### Outputs:
    ###   (xi, yj), tuple of ints or list of ints, index of grid point closest to given point
    ###     xi corresponds to lon; yj corresponds to lat.
    def get_point(self, point, gid=None, gcoord=True):
                
        #Set default grid if necessary
        if (gid == None):
            gid = self.cg
        
        #Force lon and lats into lists.
        #This enables program consistency and user convenience.
        if not isinstance(point, list):
            point = [point]

        #Lists to grid points closest to user points
        xi = []
        yj = []
        
        #Now for each lon/lat pair find the closest domain point.
        for (ulon, ulat) in point:
        
            #Convert points to grid coordinates (including grid center)
            x, y = self.proj.transform_point(ulon, ulat, self.pcp)
            cx, cy = self.proj.transform_point(self.clon[gid-1], self.clat[gid-1], self.pcp)
            
            #Calculate change in distance in grid points
            dxind = int(round((x-cx)/self.dx[gid-1]))
            dyind = int(round((y-cy)/self.dy[gid-1]))
            
            #Calculate desired point index using center indices
            xind = int(self.nx[gid-1]/2+dxind)
            yind = int(self.ny[gid-1]/2+dyind)
            
            #Check that points are still valid.
            if ((xind < 1) or (yind < 1) or (xind > self.nx[gid-1]) or (yind > self.ny[gid-1])):
                print("ERROR: Point outside grid extent. Exiting...")
                exit()
            
            #Store indices in lists
            if gcoord:
                xi.append(xind)
                yj.append(yind)
            else:
                xi.append(xind-1)
                yj.append(yind-1)
                
            
        #Return point
        if (len(xi) > 1):
            return (xi, yj)
        else:
            return (xi[0], yj[0])      
    
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
    def get_sounding(self, point, gid=None, period=None, filedate=None):

        #Force lons and lats into lists.
        #This enables program consistency and user convenience.
        if not isinstance(point, list):
            point = [point]

        #Set default grid if necessary
        if (gid == None):
            gid = self.cg
            
        #Determine time period to retreive variables over
        if (period == None):
            tstart = self.start_of_anl
            tend = self.end_of_anl
        else:
            tstart = period[0]
            tend = period[1]
        
        #Pull soundings at given locations
        sounding = []
        for p in point:
            #Retreive data
            data = self.get_var(["T", "P", "PB", "QVAPOR", "U", "V", "SINALPHA", "COSALPHA"], point=p, gid=gid, period=(tstart, tend), filedate=filedate)

            #Calculate actual temperature and pressure
            ptemp = data["T"]+300.0 #Add WRF base temp to variable
            pres = (data["P"]+data["PB"]) #Pressure in Pa. On return converts to hPa.
            temp = at.poisson(100000.0, pres, ptemp) #Convert potential temperature to temperature

            #Calculate dewpoint
            dewp = at.dewpoint(at.wtoe(pres, data["QVAPOR"]))

            #Rotate winds
            uwind = data["U"]*data["COSALPHA"]-data["V"]*data["SINALPHA"]
            vwind = data["V"]*data["COSALPHA"]+data["U"]*data["SINALPHA"]

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
    ###   point, optional, tuple, (lon, lat) of point to retrieve variables for.
    ###     Defaults to entire grid.
    ###   gid, optional, integer, grid to retrieve variables for. Defaults to current grid.
    ###   period, optional, tuple of datetime objects, (start, end) temporal bounds of plot.
    ###     Defaults to current analysis period.
    ###   level, integer, optional, model level to pull (First level is zero). Defaults to all.
    ###   filedate, datetime object, optional, date of file to read in (for single files only)
    ###
    ### Outputs:
    ###   vars, dictionary of lists, contains arrays with WRF variables keyed to var_labels (Final array order is (Time, Z, Y, X).
    ###     of if specific level requested, order is (Time, Y, X).
    def get_var(self, var_labels, point=None, gid=None, period=None, level=None, filedate=None):
    
        #Set default grid if necessary
        if (gid == None):
            gid = self.cg
            
        #Determine time period to retreive variables over
        if (period == None):
            tstart = self.start_of_anl
            tend = self.end_of_anl
        else:
            tstart = period[0]
            tend = period[1]
        
        #Calulcate user point if given
        if (point != None):
            (xind, yind) = self.get_point(point, gid, gcoord=False)
        
        #Create dictionary to store requested variables
        vars = {"date":[]}
        for vl in var_labels:
            vars[vl] = []
        
        #Loop over WRF output files
        for f in self.files[gid-1]:
        
            #Create netCDF object
            fn = nc.Dataset(f, "r")
            
            #Pull date
            date = "".join(numpy.array(fn.variables["Times"][0,:], dtype="str"))
            date = dt.datetime.strptime(date, self.tformat)
            
            #Skip file if outside analysis period
            if (((date < tstart) or (date > tend)) and (filedate == None)):
                fn.close()
                continue
            elif ((filedate != None) and (date != filedate)):
                fn.close()
                continue
                
            #Store date
            vars["date"].append(date)
                        
            #Grab other variables (handle differently if doing subsetting)
            try:
                if self.subset: #Subset the data
                    for vl in var_labels: #Loop over each variable
                            if (point == None): #Grab variable over region
                                if (fn.variables[vl].ndim == 3): #2D case (because time axis)
                                    vars[vl].append(numpy.squeeze(fn.variables[vl][:,self.yind1:self.yind2,self.xind1:self.xind2]))
                                else: #3D case
                                    vars[vl].append(numpy.squeeze(fn.variables[vl][:,:,self.yind1:self.yind2,self.xind1:self.xind2]))
                            else: #Grab var at a point
                                if (fn.variables[vl].ndim == 3): #2D case (because time axis)
                                    vars[vl].append(numpy.squeeze(fn.variables[vl][:,yind,xind]))
                                else: #3D case
                                    vars[vl].append(numpy.squeeze(fn.variables[vl][:,:,yind,xind]))

                else: #No subsetting of data
                    for vl in var_labels: #Loop over each variable
                            if (point == None): #Grab variable over region
                                vars[vl].append(numpy.squeeze(fn.variables[vl]))
                            else: #Grab var at a point
                                if (fn.variables[vl].ndim == 3): #2D case (because time axis)
                                    vars[vl].append(numpy.squeeze(fn.variables[vl][:,yind,xind]))
                                else: #3D case #3D case
                                    try:
                                        vars[vl].append(numpy.squeeze(fn.variables[vl][:,:,yind,xind]))
                                    except: #1D case like XTIME
                                        vars[vl].append(numpy.squeeze(fn.variables[vl][:]))
            
            except Exception as err:
                fn.close()
                raise Exception(err)
    
            #Close file
            fn.close()
    
        #Convert lists to numpy arrays (and pull desired level if so)
        for k in vars.keys():
            if (level != None): #Pull a level
                try: #3D case
                    vars[k] = numpy.squeeze(numpy.array(vars[k][-1][level,:,:]))
                except Exception as err: #Only fails in 2D case
                    pass
            else: #No desired level
                vars[k] = numpy.squeeze(numpy.array(vars[k]))
                
        #Returning
        return vars
    
    ### Method to retrieve WxChallenge forecast
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
    def get_wxc(self, point, gid=None, period=None):
    
        #Set default grid if necessary
        if (gid == None):
            gid = self.cg
            
        #Determine time period of forecast
        if (period == None):
            tstart = self.start_of_anl
            tend = self.end_of_anl
        else:
            tstart = period[0]
            tend = period[1]
    
        #Retrieve forecast variables
        var_names = ["T2", "U10", "V10", "RAINC", "RAINNC"]
        data = self.get_var(var_names, point=point, period=period, gid=gid)
    
        #Convert temperature from Kelvin to Fahrenheit
        temp = (data["T2"]-273.15)*9/5+32
        
        #Calculate total wind speed and convert to knots
        wspd = numpy.sqrt(data["U10"]**2+data["V10"]**2)*1.94384
        
        #Calculate accumulated rainfall over analysis period and convert to inches
        precip = ((data["RAINC"]+data["RAINNC"])-(data["RAINC"][0]+data["RAINNC"][0]))*0.0393701
        
        #Calculate forecast
        tmax = round(temp.max(), 2)
        tmin = round(temp.min(), 2)
        wmax = round(wspd.max(), 2)
        pacc = round(precip[-1], 2)
        
        #Returning
        return [tmax, tmin, wmax, pacc]
    
    ### Method to contour variable at specific model level and time
    ### Inputs:
    ###  var, string, name of variable to map
    ###  level, integer, WRF level to plot variable on
    ###  date, datetime object, date of file to to plot
    ###
    ### Outputs:
    ### (fig, ax) tuple with pyplot figure and axis object.
    def map_var(self, varname, level, date):
    
        #Grab selected variable
        data = self.get_var([varname], level=level, filedate=date)
        
        ### Make the map
        #Create figure and axis
        fig, ax = pp.subplots(nrows=1, ncols=1, subplot_kw={"projection":self.proj})
        
        #Contour
        print(data[varname].shape)
        cont = ax.contourf(self.rlons, self.rlats, data[varname], transform=self.pcp)
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
    def meteogram(self, point=None, gid=None, period=None):
    
        #Set default grid if necessary
        if (gid == None):
            gid = self.cg
            
        #Determine time period of meteogram
        if (period == None):
            tstart = self.start_of_anl
            tend = self.end_of_anl
        else:
            tstart = period[0]
            tend = period[1]
        
        #Retrieve variables to plot
        var_names = ["T2", "Q2", "U10", "V10", "SWDOWN", "PSFC", "RAINC", "RAINNC"]
        data = self.get_var(var_names, point=point, period=period, gid=gid)
        
        #Average variables over domain
        if (point == None):
            for k in data.keys():
                if (k == "date"): #Skip the date
                    continue
                
                data[k] = numpy.mean(data[k], axis=(1,2))
        
        #Calculate dewpoint from water vapor mixing ratio and change to Fahrenheit
        dewp = (at.dewpoint(at.wtoe(data["PSFC"], data["Q2"]))-273.15)*9/5+32
        
        #Convert temperature from Kelvin to Fahrenheit
        temp = (data["T2"]-273.15)*9/5+32
        
        #Calculate total wind speed and convert to knots
        wspd = numpy.sqrt(data["U10"]**2+data["V10"]**2)*1.94384
        
        #Calculate accumulated rainfall over analysis period and convert to inches
        precip = ((data["RAINC"]+data["RAINNC"])-(data["RAINC"][0]+data["RAINNC"][0]))*0.0393701
        
        #Calculate ranges for each variable
        tmin = numpy.floor(numpy.min(dewp)*0.9)
        tmax = numpy.ceil(numpy.max(temp)*1.1)
        wmin = numpy.floor(numpy.min(wspd)*0.9)
        wmax = numpy.ceil(numpy.max(wspd)*1.1)
        pmin = 0
        pmax = max(numpy.ceil(numpy.max(precip)*1.1),2)
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
        ax.set_extent(self.extent, crs=self.pcp)
        gl.xlabels_bottom = True
        gl.ylabels_left = True
        gl.xformatter = cmg.LONGITUDE_FORMATTER
        gl.yformatter = cmg.LATITUDE_FORMATTER
                    
        #Returning
        return (fig, ax)
    
    ### Method to create radar animation
    ### Animates composite reflectivity from a single WRF domain.
    ### Inputs:
    ###   spath, string, Directory in which to save radar loop (both gif and pngs)
    ###   gid, optional, integer, grid to create radar loop for.
    ###     Defaults to current grid
    ###   period, optional, tuple of datetime objects, (start, end) temporal bounds of animation.
    ###     Defaults to current analysis period.
    ###   point, optional, tuple of lon/lat (lon/lat) to highlight on map
    def radar_loop(self, spath, gid=None, period=None, point=None):
                
        #Set default grid if necessary
        if (gid == None):
            gid = self.cg
            
        #Determine time period of radar loop
        if (period == None):
            tstart = self.start_of_anl
            tend = self.end_of_anl
        else:
            tstart = period[0]
            tend = period[1]
                
        #Plot each frame of the gif
        for f in self.files[gid-1]:
        
            #Create netCDF object
            fn = nc.Dataset(f, "r")
            
            #Pull date
            date = "".join(numpy.array(fn.variables["Times"][0,:], dtype="str"))
            date = dt.datetime.strptime(date, self.tformat)
            
            #Skip file if outside analysis period
            if ((date < tstart) or
                (date > tend)):
                fn.close()
                continue
            
            #Pull data
            lats = numpy.squeeze(fn.variables["XLAT"])
            lons = numpy.squeeze(fn.variables["XLONG"])
            try:
                refl = numpy.max(numpy.squeeze(fn.variables["REFL_10CM"]), axis=0)
            except:
                print("Reflectivity not found. Check WRF namelist for do_radar_ref=1\nExiting...")
                exit()
                            
            ### Plot the data
            #Create fig and axis objects
            fig, ax = pp.subplots(nrows=1, ncols=1, subplot_kw={"projection":self.proj})
            
            #Contour the reflectivity
            cont =  ax.contourf(lons, lats, refl, cmap="gist_ncar", extend="min",
                levels=numpy.linspace(5, 60, 12), transform=self.pcp)
            cont.cmap.set_under("white")
            cb = pp.colorbar(cont, orientation="horizontal")
            cb.set_clim(5, 60)
            cb.set_label("Composite Reflectivity [dBZ]", fontsize=14)
            cb.set_ticks(numpy.linspace(5, 60, 12))
            
            #Add map features
            ax.coastlines()
            ax.add_feature(cfeature.BORDERS, edgecolor="grey")
            ax.add_feature(self.states, edgecolor="gray")
            if (point != None):
                ax.scatter(point[0], point[1], transform=self.pcp, color="black")
            
            #Title the plot
            ax.set_title(date.strftime("%Y-%m-%d   %H:%M:%S"), fontsize=16)
            
            #Save figure as png
            pp.savefig(spath+"radar_{}.png".format(date.strftime("%Y-%m-%d_%H-%M-%S")))
            pp.close(fig)
            
            #Close file
            fn.close()
            
        ### Once all pngs are saved, make gif.
        #Locate pngs just saved
        pngs = sorted(glob.glob(spath+"radar_*.png"))
        
        #Store each frame
        frames = []
        for img in pngs:
            frames.append(Image.open(img))
            
        #Make gif
        frames[0].save(spath+"radar_loop.gif", format="GIF", append_images=frames[1:],
            save_all=True, duration=500, loop=0)
            
        #Returning
        return
    
    ### Method to set working time period
    ### Inputs:
    ###   tstart, datetime object, start of simulation period to be analyzed
    ###   tend, datetime object, end of simulation period to be analyzed
    def set_anl_period(self, tstart, tend):
        
        #Check that times are within bounds of simulation
        if ((tstart < self.start_of_sim) or (tstart > self.end_of_sim) or
            (tend < self.start_of_sim) or (tend > self.end_of_sim)):
            print("Requested period outside of simulation period. Exiting.")
            exit()
            
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
        
        #Set subset flag to True. This turns on subsetting behavior for data functions.
        self.subset = True
          
        #Store the analysis extent
        self.anl_extent = extent
        
        #Calculate indices of analysis region (uses nearest point to boundary)
        self.xind1, self.yind1 = self.get_point((extent[0], extent[2]), gcoord=False)
        self.xind2, self.yind2 = self.get_point((extent[1], extent[3]), gcoord=False)
                        
        #Create lats and lats for analysis region
        lonlats = self.get_var(["XLONG", "XLAT"])
        self.rlons = lonlats["XLONG"][0]
        self.rlats = lonlats["XLAT"][0]
                        
        #Returning
        return
    
    ### Method to set the current WRF grid being examined
    ### Inputs:
    ###   gid, integer, WRF grid ID to be set as current working grid
    def set_cg(self, gid):
    
        #Check that new grid is not larger than number of grids in model
        if (gid > self.ngrids): #Print warning and return without modifying curent grid
            print("Grid {} non-existent. Max grid is: {}".format(gid, self.ngrids))
            return
                
        #Set new current grid
        self.cg = gid
        
        #Reset extent if current grid changes
        if self.subset:
            try:
                self.set_anl_region(self.anl_extent)
            except:
                self.disable_subset()
                print("Warning: Current grid does not include full extent of analysis region. Subsetting disabled.")
        
        #Returning
        return
        
    #------ Methods below this line are primarily for internal use ------#
    
    ### Method for retreiving grid information for each domain
    def get_grids(self):
    
        #Initializing arrays to store grid info
        self.nx = numpy.zeros(self.ngrids, dtype="int")  #x-dim
        self.nxs = numpy.zeros(self.ngrids, dtype="int") #staggered x-dim
        self.ny = numpy.zeros(self.ngrids, dtype="int")  #y-dim
        self.nys = numpy.zeros(self.ngrids, dtype="int") #staggered y-dim
        self.nz = numpy.zeros(self.ngrids, dtype="int")  #z-dim
        self.nzs = numpy.zeros(self.ngrids, dtype="int") #staggered z-dim
        self.ns = numpy.zeros(self.ngrids, dtype="int")  #staggered soil layers
        self.clat = numpy.zeros(self.ngrids, dtype="float")  #grid center lat
        self.clon = numpy.zeros(self.ngrids, dtype="float")  #grid center lon
        self.dx = numpy.zeros(self.ngrids, dtype="float")  #grid center lon
        self.dy = numpy.zeros(self.ngrids, dtype="float")  #grid center lon
        
        #Use the last file as our template since that timestamp will have all grids
        ftemplate = self.files[-1][-1]
                
        #Loop over each domain
        for i in range(self.ngrids):
            
            #Open appropriate file to retrieve grid
            fname = ftemplate[:ftemplate.rfind("_d")+2]+"{:02d}".format(i+1)+ftemplate[ftemplate.rfind("_d")+4:]
            dfile = nc.Dataset(fname, "r")
            
            #Retrieve grid variables
            self.nx[i] = dfile.dimensions["west_east"].size
            self.nxs[i] = dfile.dimensions["west_east_stag"].size
            self.ny[i] = dfile.dimensions["south_north"].size
            self.nys[i] = dfile.dimensions["south_north_stag"].size
            self.nz[i] = dfile.dimensions["bottom_top"].size
            self.nzs[i] = dfile.dimensions["bottom_top_stag"].size
            self.ns[i] = dfile.dimensions["soil_layers_stag"].size
            self.clat[i] = dfile.__dict__["CEN_LAT"]
            self.clon[i] = dfile.__dict__["CEN_LON"]
            self.dx[i] = dfile.__dict__["DX"]
            self.dy[i] = dfile.__dict__["DY"]
            
        #Returning
        return
     
    ### Methods for setting the WRF grid map projection
    #Lambert Conformal projection
    def lambert(self):
    
        return ccrs.LambertConformal(central_longitude=self.lon0, central_latitude=self.lat0,
            standard_parallels=(self.slat1, self.slat2))
        
    #Lat/lon porjection
    def latlon(self):
    
        print("WARNING: lat-lon option not finished. Use at own risk.")
        return ccrs.PlateCarree()
    
    #Mercator projection
    def mercator(self):
    
        return ccrs.Mercator(latitude_true_scale=self.slat1)
    
    #Polar Stereographic projection
    def polar(self):
    
        return ccrs.NorthPolarStereo(central_longitude=self.lon0, true_scale_latitude=self.slat1)
        
        
############################################################
#--------------------     WRFgrid     ---------------------#
############################################################
### The WRFgrid class provides tools for visualizing the WRF
### grid before the model is even run. It requires only a WPS
### namelist to provide the grid information.
###
### Available Attributes:
###  nml - f90nml namelist object
###  proj_name - Name of model grid projection
###  proj - Cartopy projection object for model grid
###  pcp - Cartopy PlateCarree projection object
###  lat0 - Domain center latitude
###  lon0 - Domain center longitude
###  slat1 - Standard parallel 1
###  slat2 - Standard parallel 2
###  slon1 - Standard longitude
###  clat - Center latitude of each grid
###  clon - Center longitude of each grid
###  extent - Bounds of each grid (in grid projection coordinates)
###  gwidth - Width of each grid (x span)
###  gheight - Height of each grid (y span)
###  cg - Current grid being handled
###  ngrids - Number of grids in domain
###  nx - Number of x points on each grid
###  ny - Number of y points on each grid
###  grid_ratio - Ratio of grid spacing between each grid and its parent
###  dx - x grid spacing on each grid
###  dy - y grid spacing on each grid
###  lcfile - File path of landcover dataset used in plotting
###
### Available Methods:
###  find_point(lon, lat) - Compute nearest grid location to a given lon/lat
###  plot_grid() - Plot grid and all interior grids
###  set_cg(gid) - Set current working grid
###  set_lcfile(filepath) - Set location of landcover data file used in plotting

class WRFgrid:
    ### Method to construct object
    ### Inputs:
    ###   namelist, string, path to WPS namelist file to be read in
    def __init__(self, namelist):
        
        #First thing is to open WRF namelist file
        try:
            self.nml = f90nml.read(namelist)
        except:
            print("{} not found.".format(namelist))
            sys.exit(2)
                            
        ### Identify WRF grid map projection
        #Create dictionary with possible projections
        proj_list = {"mercator":self.mercator, "lambert":self.lambert, "polar":self.polar,
            "lat-lon":self.latlon}
            
        #Pull map information
        self.lat0 = self.nml["geogrid"]["ref_lat"]       #Domain center
        self.lon0 = self.nml["geogrid"]["ref_lon"]       #Domain center
        self.slat1 = self.nml["geogrid"]["truelat1"]     #Standard parallel 1
        self.slat2 = self.nml["geogrid"]["truelat2"]     #Standard parallel 2
        self.slon1 = self.nml["geogrid"]["stand_lon"]    #Standard longitude
        self.proj_name = self.nml["geogrid"]["map_proj"] #Map projection name
        
        #Setup map projection
        self.proj = proj_list[self.proj_name]() #Model projection
        self.pcp = ccrs.PlateCarree()           #Lat/lon projection for interal use
        
        ### Compute WRG Grid information                
        #Set current grid
        self.cg = 1 #Start with the parent grid
        
        #Grab domain information
        #Number of grids in domain
        self.ngrids = self.nml["share"]["max_dom"]
        
        #Number of points in x and y directions
        self.nx = numpy.atleast_1d(self.nml["geogrid"]["e_we"])[:self.ngrids]
        self.ny = numpy.atleast_1d(self.nml["geogrid"]["e_sn"])[:self.ngrids]

        #Ratio of each grid to its parent and corresponding grid spacings.
        self.grid_ratio = numpy.atleast_1d(self.nml["geogrid"]["parent_grid_ratio"])[:self.ngrids]
        self.dx = numpy.ones(self.ngrids)*self.nml["geogrid"]["dx"] #Only starting value
        self.dy = numpy.ones(self.ngrids)*self.nml["geogrid"]["dy"] #Only starting value
        
        #Compute grid spacing
        ratio = self.grid_ratio[0]
        for i in range(1, len(self.grid_ratio)):
            ratio *= self.grid_ratio[i]
            self.dx[i] /= ratio
            self.dy[i] /= ratio
        
        #Calculate grid widths and heights
        self.gwidth = (self.nx-1)*self.dx
        self.gheight = (self.ny-1)*self.dy
        
        #Calculate grid centers
        (self.clat, self.clon) = self.get_centers()
        
        #Calculate grid extents
        self.extent = list(self.get_extent(gid=i) for i in range(1, self.ngrids+1))

        #Attach land cover dataset as attribute.
        self.lcfile = "MCD12C1.A2018001.006.2019200161458.hdf"

        #Returning
        return
    
    ### Method to retrieve closest grid point to desired location
    ### Inputs:
    ###   lon, float or list of floats, longitude of desired point
    ###   lat, float or list of floats, latitude of desired point
    ###   gid=, integer, optional, grid to find nearest point for. Defaults to current grid.
    ###
    ### Outputs:
    ###   (xi, yj), tuple of ints or list of ints, index of grid point closest to given point
    ###     xi corresponds to lon; yj corresponds to lat.
    def find_point(self, lon, lat, gid=None):
        #Set gid to default if required
        if (gid == None):
            gid = self.cg
        
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
            if ((x < self.extent[gid-1][0]) or (x > self.extent[gid-1][1]) or
                (y < self.extent[gid-1][2]) or (y > self.extent[gid-1][3])):
                
                print("Point Lon: {}, Lat: {}) is outside model grid {}. Exiting...".format(ulon,ulat,gid))
                sys.exit(2)

            #Calculate index of each point
            xi.append(int(round((x-self.extent[gid-1][0])/self.dx[gid-1])))
            yj.append(int(round((self.extent[gid-1][3]-y)/self.dy[gid-1])))

        #Return indices to user (only return lists if user gave list)
        if (len(xi) > 1):
            return (xi, yj)
        else:
            return (xi[0], yj[0])
        
    ### Method to plot WRF grids
    ### Given a grid ID, it plots that grid and all within.
    ### Inputs:
    ###   gid=, integer, optional, Outer most grid ID to plot. If not provided, defaults to current grid.
    ###   lc=, boolean, optional, Flag to include land cover imagery. Defaults to False.
    ###   points=, list of tuples, optional, tuples of lon/lat pairs.
    ###   savep=, string, optional, location to save domain plot. Displays if omitted.
    def plot_grid(self, gid=None, lc=False, points=[], spath=None):
               
        #Set gid to default if required
        if (gid == None):
            gid = self.cg
               
        #Check that requested grid is not larger than number of grids in model
        if (gid > self.ngrids): #Print warning
            print("Grid {} non-existent. Plotting outermost grid".format(gid))
            gid = 1
        
        #Setup figure and axis objects
        fig, ax = pp.subplots(nrows=1, ncols=1, subplot_kw={"projection":self.proj})
        
        #Set plot extent
        ax.set_extent(self.extent[gid-1], self.proj)
        
        #Add background image
        if lc:
            self.plot_landcover(ax, gid)
        else:
            ax.stock_img()
            
        #Add map features
        ax.coastlines()
        ax.add_feature(cfeature.STATES)
        ax.add_feature(cfeature.BORDERS)
        gl = ax.gridlines(color="black", linestyle="--", draw_labels=False)
        
        #Now add interior domains
        for i in range(gid-1, self.ngrids):
            box = Rectangle((self.extent[i][0], self.extent[i][2]), self.gwidth[i],
                self.gheight[i], color="red", fill=False, transform=self.proj)
            ax.add_patch(box)
        
        #Add any given points
        for p in points:
            ax.scatter(p[0], p[1], color="black", transform=self.pcp)
        
        #Display the image
        pp.tight_layout()
        if (spath == None):
            pp.show()
        else:
            pp.savefig(spath)
        
        #Returning
        return
        
    ### Method to set the current WRF grid being examined
    ### Inputs:
    ###   gid, integer, WRF grid ID to be set as current working grid
    def set_cg(self, gid):
    
        #Check that new grid is not larger than number of grids in model
        if (gid > self.ngrids): #Print warning and return without modifying curent grid
            print("Grid {} non-existent. Max grid is: {}".format(gid, self.ngrids))
            return
        
        #Set new current grid
        self.cg = gid
        
    ### Method to specify location of land cover data file
    ### By default, object searches within the local directory.
    ### Inputs:
    ###   filepath, string, full path to MODIS land cover data file.
    def set_lcfile(self, filepath):
        
        #Set new file location
        self.lcfile = filepath
        
        #Returning
        return
    
    #------ Methods below this line are primarily for internal use ------#
    
    ### Method for calculating grid center lat/lons
    ### Outputs:
    ###   (clat, clon), tuple of arrays of floats, Grid center latitudes, grid center longitudes
    def get_centers(self):
       
        #Initialize local arrays for calculations
        clat = numpy.ones(self.ngrids)*self.lat0 #Multiply by lat0 sets the first grid center to the proper value
        clon = numpy.ones(self.ngrids)*self.lon0 #Multiply by lon0 sets the first grid center to the proper value
                
        #Calculate coordinates of lower-left corner of outermost domain
        (xc, yc) = self.proj.transform_point(self.lon0, self.lat0, self.pcp) #Get center
        x0 = xc-self.gwidth[0]/2  #Identify left boundary
        y0 = yc-self.gheight[0]/2 #Identify bottom boundary
               
        #Grab starting points of each domain
        i0 = self.nml["geogrid"]["i_parent_start"]
        j0 = self.nml["geogrid"]["j_parent_start"]
                        
        #Now calculate grid centers
        for i in range(1, self.ngrids):
            #Compute lower left point of current domain
            x0 += (i0[i]-1)*self.dx[i-1]
            y0 += (j0[i]-1)*self.dy[i-1]
            
            #Calculating grid center
            xc = x0+self.gwidth[i]/2
            yc = y0+self.gheight[i]/2
            
            #Convert to lat/lon
            (clon[i], clat[i]) = self.pcp.transform_point(xc, yc, self.proj)
            
        #Return grid centers
        return (clat, clon)
            
    ### Method for calculating grid extent
    ### Inputs:
    ###   gid=, integer, optional, grid ID number (outermost grid is grid number 1). Defaults to current grid.
    ###
    ### Outputs:
    ###   (w, e, s, n), tuple of floats, grid extent in projection coordinates. Order is west, east, south, north.
    def get_extent(self, gid=None):
        #Set gid to default if required
        if (gid == None):
            gid = self.cg
        
        #Transform grid center into projection coordinates
        (xc, yc) = self.proj.transform_point(self.clon[gid-1], self.clat[gid-1], self.pcp)
        
        #Return grid extent
        return (xc-self.gwidth[gid-1]/2, xc+self.gwidth[gid-1]/2,
            yc-self.gheight[gid-1]/2, yc+self.gheight[gid-1]/2)
    
    ### Method to plot landcover if requested
    ### Inputs:
    ###   ax, pyplot axis object, axis upon which to plot land cover
    ###   gid, integer, optional, Outer most grid ID on plot.
    def plot_landcover(self, ax, gid):
        #Read in land cover data
        lco = MCD12C1(self.lcfile)
        landcover = lco.get("Majority_Land_Cover_Type_1")

        #Colorbar to use for plotting
        cmap = lco.lc1_cmap

        #Data levels
        levs = range(len(lco.lc1_legend)+1)

        #Populate 2D location grids
        (lon_grid, lat_grid) = numpy.meshgrid(lco.lons, lco.lats)
        
        #Convert lat/lons to projection coordinates
        tgrid = self.proj.transform_points(self.pcp, lon_grid, lat_grid)
        
        #Split 3D tgrid into its constituents
        x = tgrid[:,:,0]
        y = tgrid[:,:,1]
        
        #Limit extent to model domain for faster plotting
        mask1 = x < self.extent[gid-1][0]
        mask2 = x > self.extent[gid-1][1]
        mask3 = y < self.extent[gid-1][2]
        mask4 = y > self.extent[gid-1][3]
        mask = mask1 | mask2 | mask3 | mask4
        x = numpy.ma.masked_array(x, mask=mask)
        y = numpy.ma.masked_array(y, mask=mask)
        landcover = numpy.ma.masked_array(landcover, mask=mask)

        #Creating plot
        cont = ax.contourf(x, y, landcover, transform=self.proj, cmap=cmap, levels=levs)

        #Colorbar
        ratio = x.shape[1]/x.shape[0]
        cb = pp.colorbar(cont, ax=ax, ticks=range(len(lco.lc1_legend)), fraction=0.045*ratio, pad=0.04)
        cb.ax.locator_params(nbins=len(lco.lc1_legend))
        cb.ax.set_yticklabels(lco.lc1_legend, rotation=-30)
        
        #Clean up data object
        del(lco)
        
        #Returning
        return
    
    ### Methods for setting the WRF grid map projection
    #Lambert Conformal projection
    def lambert(self):
    
        return ccrs.LambertConformal(central_longitude=self.lon0, central_latitude=self.lat0,
            standard_parallels=(self.slat1, self.slat2))
        
    #Lat/lon porjection
    def latlon(self):
    
        print("WARNING: lat-lon option not finished. Use at own risk.")
        return ccrs.PlateCarree()
    
    #Mercator projection
    def mercator(self):
    
        return ccrs.Mercator(latitude_true_scale=self.slat1)
    
    #Polar Stereographic projection
    def polar(self):
    
        return ccrs.NorthPolarStereo(central_longitude=self.slon1, true_scale_latitude=self.slat1)
    
    ### Method to cleanly destroy object
    def __del__(self):
        
        #Clean up namelist object to free memory
        try:
            del(self.nml)
        except:
            pass
            
        #Returning
        return
