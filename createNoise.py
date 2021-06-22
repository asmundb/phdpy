#!/usr/bin/env python3

"""
#################################################
White noise atmospheric forcing
Created 6.11.17 JB
##################################################
"""

from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import time
from datetime import timedelta, date, datetime
import os,glob,sys

#from multiprocessing import Process

class Dimensions:

    def __init__(self,ncForcing):
        self.npoints = ncForcing.dimensions["Number_of_points"].size
        self.ntimes = ncForcing.dimensions["time"].size
        self._params = {"tau": 24.0, "dt": 1.0, "nvar": 3, "nens": 1}
        tmfmt = "hours since %Y-%m-%d %H:%M:%S 0:00"
        self.timeUnit = ncForcing.variables["time"].units
        self.startTime = datetime.strptime(self.timeUnit, tmfmt)
        self.endTime = self.startTime + timedelta(hours=float(ncForcing.variables["time"][-1]))

    @property
    def nens(self):
        return self._params["nens"]

    @property
    def nvar(self):
        return self._params["nvar"]

    @property
    def tau(self):
        return self._params["tau"]

    @property
    def dt(self):
        return self._params["dt"]

    @nens.setter
    def nens(self,nens):
        self._params["nens"] = nens



def datespan(startDate, endDate, delta=timedelta(days=1)):
    currentDate = startDate
    while currentDate < endDate:
        yield currentDate
        currentDate += delta


def createNoise(dims):

    # NLDAS
    #outpath = '/home/sm_josbl/hm_home/programs/python/pertForcing/'

    nens = dims.nens           # number of ensemble members - 1
    nvar = dims.nvar            # number of variables, rainf, sw, lw
    T = dims.ntimes              # number of timesteps in forcing file
    dt = dims.dt            # forcing timestep
    tau = dims.tau          # decorrelation time

 # AA specific variables

    N = dims.npoints
#
    alpha = 1 - (dt/tau)
    beta = np.sqrt(1 - alpha*alpha)

    corr = [[1, -0.8,0.5],[-0.8, 1, -0.5],[0.5,-0.5,1]]
    stdm = [1.0, 1.0, 1.0]
    stdm = np.array(stdm)
    corr = np.array(corr)
    stdm = np.diag(stdm)
    covm = np.dot(stdm, np.dot(corr,stdm))
    L = np.linalg.cholesky(covm)

    #ntot = 0
    #for timestamp in datespan(datetime(2019, 12, 15, 0, 00),
    #                        datetime(2019, 12, 15, 3, 00),
    #                        delta=timedelta(hours=4)):
    #    ntot += 1

    # Loop through dates

    rainf_pert_prev = np.zeros((nens,N))
    sw_pert_prev = np.zeros((nens,N))
    lw_pert_prev = np.zeros((nens,N))

    for timestamp in datespan(datetime(2019, 12, 15, 0, 0),
                            datetime(2019, 12, 15, 3, 0),
                            delta=timedelta(hours=4)):

       #t0 = time.time()
       #date = str(timestamp)
       #year = date[2:4]
       #month = date[5:7]
       #day = date[8:10]
       #hour = date[11:13]

       rainf_pert = np.zeros((nens,T,N))
       sw_pert = np.zeros((nens,T,N))
       lw_pert = np.zeros((nens,T,N))

       q_rainf = np.zeros((nens,T,N))
       q_sw = np.zeros((nens,T,N))
       q_lw = np.zeros((nens,T,N))

       #print(date)

       U = np.zeros((nens, nvar, T, N))
       E = np.zeros((nens, nvar, T, N))

 #      for i in xrange(N):

       q_rainf[:,0,:] = rainf_pert_prev[:,:]
       q_sw[:,0,:] = sw_pert_prev[:,:]
       q_lw[:,0,:] = lw_pert_prev[:,:]


       for ens in range(nens):
           for idx in range(1,T):
               rainf_rand = beta*np.random.normal(0,1,N)
               sw_rand = beta*np.random.normal(0,1,N)
               lw_rand = beta*np.random.normal(0,1,N)
               # AR(1) auto-correlation:
               q_rainf[ens,idx,:] = alpha*q_rainf[ens,idx-1,:] + rainf_rand[:]
               q_sw[ens,idx,:] = alpha*q_sw[ens,idx-1,:] + sw_rand[:]
               q_lw[ens,idx,:] = alpha*q_lw[ens,idx-1,:] + lw_rand[:]

       for i in range(N):
           for ens in range(nens):
               U[ens,:,:,i] =  [ q_rainf[ens,:,i], q_sw[ens,:,i], q_lw[ens,:,i] ]
               E[ens,:,:,i] = L.dot(U[ens,:,:,i])

           rainf_pert_prev[:,i] = q_rainf[:,-1,i]
           sw_pert_prev[:,i] = q_sw[:,-1,i]
           lw_pert_prev[:,i] = q_lw[:,-1,i]

           rainf_pert[:,:,i] = E[:,0,:,i]
           sw_pert[:,:,i] = E[:,1,:,i]
           lw_pert[:,:,i] = E[:,2,:,i]

       for ens in range(nens):

 #          forcingpath = 'noise_%s%s%s%s_%d.nc'%(year,month,day,hour, ens+1)
           forcingpath = 'noise_%02d.nc'%(ens+1)
           saveForcing2File(rainf_pert[ens,:,:], sw_pert[ens,:,:], lw_pert[ens,:,:], forcingpath,dims)

def saveForcing2File(rainf_pert, sw_pert, lw_pert, forcingpath,dims):

    timespan = dims.ntimes
    N = dims.npoints
    nens = dims.nens

    met_path = "."

#    date_val = forcingpath.split('/')[-1].split('_')[-1].split('.')[0]

    surfex = Dataset(forcingpath, 'w', format='NETCDF3_CLASSIC')
    
    # Dim the dimensions of NetCDF
    surfex.createDimension('time', timespan)
    surfex.createDimension('Number_of_points', N)


    FRC_STP = surfex.createVariable('FRC_TIME_STP', 'f4')
    FRC_STP.units = "s"

        
#    time_string = "hours since " + date_val[0:4]+"-"+date_val[4:6]+"-"+date_val[6:] + " " + "00:00:00"

    TIME = surfex.createVariable('time','f4',('time',))    
    TIME.units = dims.timeUnit #time_string
    time = range(0,timespan,1)
    TIME[:] = time[:]
    
   
#    TAIR = surfex.createVariable('Tair','f4',('time', 'LAT','LON'),fill_value=1.00000002004e+20)
#    TAIR.missing_value = 1.00000002004e+20
#    TAIR.long_name = "TAIR"
#    TAIR.history = met_path
#    TAIR.coordinates = "LAT LON"
#    TAIR[:,:,:] = tair_pert[:,:,:]
        
    
    LW_DOWN = surfex.createVariable('LWdown', 'f4',('time','Number_of_points'),fill_value=1.00000002004e+20)
    LW_DOWN.missing_value = 1.00000002004e+20
    LW_DOWN.long_name = "Longwave down"
    LW_DOWN.history = met_path
    LW_DOWN.coordinates = "Number_of_points"
    LW_DOWN[:,:] = lw_pert[:,:]
    
    RAIN_SURF = surfex.createVariable('Rainf', 'f4',('time', 'Number_of_points'),fill_value=1.00000002004e+20)
    RAIN_SURF.missing_value = 1.00000002004e+20
    RAIN_SURF.long_name = "Rain surface"
    RAIN_SURF.history = met_path  
    RAIN_SURF.coordinates = "Number_of_points"
    RAIN_SURF[:,:] = rainf_pert[:,:]
    
    
    SW_DOWN = surfex.createVariable('DIR_SWdown', 'f4',('time','Number_of_points'),fill_value=1.00000002004e+20)
    SW_DOWN.missing_value = 1.00000002004e+20
    SW_DOWN.long_name = "Shortwave down"
    SW_DOWN.history = met_path 
    SW_DOWN.coordinates = "Number_of_points"
    SW_DOWN[:,:] = sw_pert[:,:]
            
    FRC_STP[:] = 3600                                                           # Forcing timestep    
    surfex.close()
 
      
def main():
    ncfile = "/home/asmundb/Projects/H2O/wp3/pert_forcing/FORCING.nc"
    nens = 10

    with Dataset(ncfile,'r') as f:
        dims = Dimensions(f)
    dims.nens = nens
    createNoise(dims)


if __name__ == '__main__':
    #inFile = sys.argv[1]
    #nens = sys.ergv[2]
    dims = main()

