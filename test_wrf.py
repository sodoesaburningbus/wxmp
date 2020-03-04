#This script tests features of the WRF grid object
import matplotlib
matplotlib.use("Agg")

from wxmp.wrf import WRFgrid

#Namelist to read in
namelist = "/rstor/cphillip/model_exps/seasia_bb/smoke/wps/namelist.wps"

#Creating WRF grid object
grid = WRFgrid(namelist)

#Print some basic grid properties
print("Namelist: {}".format(namelist))
print("Number of grids: {}".format(grid.ngrids))
print("Projection: {}".format(grid.proj_name))
print("Domain center ({:.02f}, {:.02f})".format(grid.lat0, grid.lon0))
print("Outermost domain width: {:.0f}m".format(grid.gwidth[0]))
print("Outermost domain height: {:.0f}m".format(grid.gheight[0]))
print("Grid spacing in x-direction for each grid: {}".format(grid.dx))
print("Grid spacing in y-direction for each grid: {}".format(grid.dy))

#Make list of grid centers
clist = list((grid.clon[i], grid.clat[i]) for i in range(grid.ngrids))

#Find grid point corresponding to each center
for i in range(1, grid.ngrids+1):
    grid.set_cg(i)
    (x, y) = grid.find_point(clist[i-1][0], clist[i-1][1])
    print("Center of grid {}: xi = {} yj = {}".format(i,x,y))

#Display grids with centers and landcover
grid.set_lcfile("/rhome/cphillip/datasets/modis/landcover/MCD12C1.A2018001.006.2019200161458.hdf")
grid.plot_grid(points=clist, lc=True, gid=1)