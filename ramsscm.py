### This object is used to handle SCM runs from the RAMS modeling system.
### Note: The RAMS SCM is a 5x5 grid of points with cyclic boundary conditions.
###   This object is not compatible with the optional SCM output for larger simulations.
###
### All RAMS output files should be located in a single directory with one time stamp
### per file.
###
### Written by Christopher Phillips

############################################################
#---------------------     RAMS_SCM     ---------------------#
############################################################
### The RAMS_SCM class provides tools for visualizing RAMS SCM
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

### Importing required modules
import atmos.thermo as at
import atmos.math as am
from datetime import datetime
from glob import glob
import h5py
import matplotlib.pyplot as pp
import numpy

class RAMS_SCM:

    ### Method to construct object
    ### Inputs:
    ###   fdir, string, directory containing RAMS output
    ###   pref, optional, string, file prefix as set in the namelist,
    ###     defaults to "a"
    def __init__(self, fdir, pref="a"):
    
        # Set some basic constants
        self.P00 = 100000.0 # Reference pressure in Pa (based on the RAMS code)
        self.CP = 1004.0    # Specific heat at constant pressure as set in RAMS (J/kg/K)
    
        # Locate the files and construct list of valid times
        self.fdir = fdir+"/" #Add final slash just in case
        self.files = numpy.array(sorted(glob("{}/{}*.h5".format(self.fdir, pref))))
        self.valid_dates = numpy.array([datetime.strptime(f[-23:-6], "%Y-%m-%d-%H%M%S") for f in self.files])
        self.start_of_sim = self.valid_dates[0]
        self.end_of_sim = self.valid_dates[-1]
        
        # Extract grid data
        self.file = h5py.File(self.files[0], "r")
        self.lat = self.file["GLAT"][2,2]
        self.lon = self.file["GLON"][2,2]
        self.nx = 1
        self.ny = 1
        self.nz = self.file["UC"].shape[0]
        self.topo = self.file["TOPT"][2,2]
        
        # Set the analysis period (start with full simulation)
        self.set_anl_period(self.start_of_sim, self.end_of_sim)
        
        # Return
        return
    
    
    ### Method to compute model geopotential height levels at each analysis time.
    ### Heights are above mean sea-level
    ### Inputs:
    ###  period, optional, tuple of datetime objects, the analysis period to retrieve
    ###    defaults to the current analysis period
    ###  base_state, optional, boolean, whether to compute the reference state geopotential height,
    ###    defaults to False (i.e. the actual geopotential height). Base state may not work properly yet.
    ### Outputs:
    ###  gph_out, array of floats, shape is (Time, Z), the geopotential height above mean-sea level.
    def get_gph(self, period=None, base_state=False):
    
        ### Retrieve necessary variables
        if base_state:
            exner = "PI0"
        else:
            exner = "PI"
            
        data = self.get_vars([exner, "THETA", "RV"], period=period)
        
        ### Compute virtual temperature and pressure
        pres = ((data[exner]/self.CP)**(self.CP/at.RD))*self.P00
        temp = data[exner]/self.CP*data["THETA"]
        vtemp = at.virt_temp(temp, data["RV"])
        
        ### Now integrate the hypsometric equation for each time step
        gph_out = []
        for t in range(pres.shape[0]):
            gph = [self.topo]
            for k in range(vtemp.shape[1]-1):
                tvbar = am.layer_average(pres[t,k:k+2], vtemp[t,k:k+2])
                gph.append(gph[-1]+at.RD*tvbar/at.G*numpy.log(pres[t,k]/pres[t,k+1]))
            gph_out.append(gph)
        
        ### Convert to numpy array and return
        return numpy.array(gph_out)
    
    ### Method to retrieve a list of variables
    ### They are returned as a dictionary of arrays keyed to the variable names
    ### Inputs:
    ###  names, list of strings, the variable names to retrieve
    ###  period, optional, tuple of datetime objects, the analysis period to retrieve
    ###    defaults to the current analysis period
    def get_vars(self, names, period=None):
    
        ### Set the temporary analysis period if necessary
        if (period != None):
            
            # Check that times are within bounds of simulation
            if ((tstart < self.start_of_sim) or (tstart > self.end_of_sim) or
                (tend < self.start_of_sim) or (tend > self.end_of_sim)):
                print("WARNING: Requested period outside of simulation period. Returning without changing analysis period.")
                return -1
            
            # Grab file indices that match those bounds
            ind1 = numpy.arange(0, self.files.size, dtype="int")[self.valid_dates >= period[0]][0]
            ind2 = numpy.arange(0, self.files.size, dtype="int")[self.valid_dates <= period[1]][-1]
            
        else:
        
            ind1 = self.start_of_anl_ind
            ind2 = self.end_of_anl_ind
        
        ### Create dictionary
        data = {"time":[]}
        for n in names:
            data[n] = []
        
        ### Loop across the files for the analysis period
        ### and pull the data
        for f, d in zip(self.files[ind1:ind2+1], self.valid_dates[ind1:ind2+1]):
        
            # Try closing any previously open files
            try:
                self.file.close()
            except:
                pass
                
            # Open file
            self.file = h5py.File(f, "r")
        
            # Store time
            data["time"].append(d)
            
            # Get the other variables
            for n in names:
            
                #Check that name is in file
                if (n not in self.file.keys()):
                    raise ValueError("{} not in RAMS output file!".format(n))
            
                try:
                    data[n].append(self.file[n][:,:,2,2])
                except:
                    try:
                        data[n].append(self.file[n][:,2,2])
                    except:
                        data[n].append(self.file[n][2,2])
                
        # Convert lists to arrays
        for k in data.keys():
            data[k] = numpy.array(data[k])
                
        # Return the dictionary
        return data
            
        
    ### Method to set working time period
    ### Inputs:
    ###   tstart, datetime object, start of simulation period to be analyzed
    ###   tend, datetime object, end of simulation period to be analyzed
    def set_anl_period(self, tstart, tend):

        # Check that times are within bounds of simulation
        if ((tstart < self.start_of_sim) or (tstart > self.end_of_sim) or
            (tend < self.start_of_sim) or (tend > self.end_of_sim)):
            print("WARNING: Requested period outside of simulation period. Returning without changing analysis period.")
            return -1

        # Set new analysis bounds
        self.start_of_anl = tstart
        self.end_of_anl = tend
        
        # Grab file indices that match those bounds
        self.start_of_anl_ind = numpy.arange(0, self.files.size, dtype="int")[self.valid_dates >= tstart][0]
        self.end_of_anl_ind = numpy.arange(0, self.files.size, dtype="int")[self.valid_dates <= tend][-1]

        #Returning
        return