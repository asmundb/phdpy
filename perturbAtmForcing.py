"""
#################################################
Python program for perturbation of forcing files
Created 11.8.17 JB
11.8.17: AR(1) process included 
25.8.17: AR(1) with perturbed correlations, ref Stefan S, tair = alpha*tair + eps_tair - 0.8*SW + eps_SW?
##################################################
"""
import netCDF4 as nc
import numpy as np
#import matplotlib.pylab as plt
import numpy.ma as ma
import shutil
import glob

def readSurfaceVariables():
    """
    Function reading the surface variables we wish 
    to perturb. 
    1. Precipitation
    2. Long wave downwards radiation
    3. Short wave downwards radiation
    4. Air temperature
    """

    # Read forcing file

    forcing = nc.Dataset('FORCING.nc', mode='r')

    # List of variables
    rainf = forcing.variables['Rainf'][:]
    lwdown = forcing.variables['LWdown'][:]
    swdown = forcing.variables['DIR_SWdown'][:]
#    tair = forcing.variables['Tair'][:]    
    forcing.close()

    return rainf, lwdown, swdown #, tair


def readNoise():
    noisefiles = glob.glob("noise_*.nc")
    nens = len(noisefiles)
    with nc.Dataset(noisefiles[0],'r') as f:
        N = f.dimensions["Number_of_points"].size
        nt = f.dimensions["time"].size

    q_rainf = np.zeros((nens,nt,N))
    q_sw = np.zeros((nens,nt,N))
    q_lw = np.zeros((nens,nt,N))
    #q_tair = np.zeros((nens,nt,N))

    for ens in range(nens):
        with nc.Dataset(noisefiles[ens],mode='r') as noise:
            rainftmp = noise.variables['Rainf'][:]
            swtmp = noise.variables['DIR_SWdown'][:]
            lwtmp = noise.variables['LWdown'][:]
    #        tairtmp = noise.variables['Tair'][:]

        q_rainf[ens,:,:] = rainftmp[:,:]
        q_sw[ens,:,:] = swtmp[:,:]
        q_lw[ens,:,:] = lwtmp[:,:]
#        q_tair[ens,:,:,:] = tairtmp[:,:,:]
    return q_rainf, q_sw, q_lw #, q_tair


def pertForcing(rainf, lwdown, swdown, q_rainf, q_sw, q_lw):
    """
    Function doing univariate perturbation of atmospheric forcing
    Using eq. 31 in Evensen 2003
    """

    nens = q_rainf.shape[0]
    nt = q_rainf.shape[1]
    N = q_rainf.shape[2]

#    lw_pert_mean = np.zeros((nens,nt,N))
#    sw_pert_mean = np.zeros((nens,nt,N))
#    rainf_pert_mean = np.zeros((nens,nt,N))
#    tair_pert_mean = np.zeros((nens,nt,ni,nj))

    lwdown_pert = np.zeros((nens,nt,N))
#    tair_pert = np.zeros((nens,nt,ni,nj))
    swdown_pert = np.zeros((nens,nt,N))
    rainf_pert = np.zeros((nens,nt,N))

#    origrainf = rainf
#    origsw = swdown
#    origlw = lwdown

#    q_sw[:,:,:,:] = 0.3*q_sw[:,:,:,:]
#    q_rainf[:,:,:,:] = 0.5*q_rainf[:,:,:,:]
    
    for ens in range(nens):
        # Additive error
        lwdown_pert[ens,:,:] = lwdown[:,:] + 30.0*q_lw[ens,:,:]
#        tair_pert[ens,:,:,:] = tair[:,:,:] + q_tair[ens,:,:,:]

        # Multiplicative error
        swdown_pert[ens,:,:] = swdown[:,:]*np.exp(0.3*q_sw[ens,:,:] - 0.5*(0.3*0.3))
        rainf_pert[ens,:,:] = rainf[:,:]*np.exp(0.5*q_rainf[ens,:,:] - 0.5*(0.5*0.5))

#    tair_pert[-1,:,:,:] = tair[:,:,:]

    for ens in range(nens):
        writeSurface(rainf_pert[ens,:,:], lwdown_pert[ens,:,:], swdown_pert[ens,:,:], ens)

def writeSurface(rainf_pert, lwdown_pert, swdown_pert, ens):
    """
    Function writing perturbed variables to netCDF file
    """

    print( ens)
    print( "!!!!!!!! Writing to file !!!!!!!!!")
    outfile = 'FORCING_%02d.nc'%(ens+1)
    shutil.copyfile("FORCING.nc",outfile)
    forcing_out = nc.Dataset(outfile, mode='r+')

    forcing_out.variables['Rainf'][:] = rainf_pert[:]
    forcing_out.variables['LWdown'][:] = lwdown_pert[:]
    forcing_out.variables['DIR_SWdown'][:] = swdown_pert[:]
#    forcing_out.variables['Tair'][:] = tair_pert[:]

    forcing_out.close()

def main():

    rainf, lwdown, swdown = readSurfaceVariables()
    q_rainf, q_sw, q_lw, = readNoise()

    pertForcing(rainf, lwdown, swdown, q_rainf,q_sw,q_lw)

            
if __name__ == '__main__':     
    main() 
