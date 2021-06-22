import netCDF4 as nc
import numpy as np
import glob
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("TkAgg")

files = glob.glob("FORCING_*.nc")
nens = len(files)
nx,ny = 20,20

with nc.Dataset(files[0],'r') as f:
    npoints = f.dimensions["Number_of_points"].size
    ntimes = f.dimensions["time"].size



keys = ["Rainf", "DIR_SWdown", "LWdown"]
var = {}

for i in range(len(keys)):
    var[keys[i]] = np.zeros((nens,ntimes,npoints))

for i in range(nens):
    with nc.Dataset(files[i]) as f:
        for j in range(len(keys)):
            var[keys[j]][i, :, :] = f.variables[keys[j]][:]

for key in keys:
    plt.figure()
    plt.plot(var[key][:,:,0].transpose(),)
    plt.title(key)