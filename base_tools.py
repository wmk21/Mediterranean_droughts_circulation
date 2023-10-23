"""
   Some additional funcions. 
   
"""

import numpy as np
from netCDF4 import Dataset

def save_netcdf_time(pc, titles):  #time, neof, types
    
    nt=pc.shape[0]
    ntype=pc.shape[1]
    time_m=np.arange(0,nt)
    
    dataset =Dataset(titles,'w',format='NETCDF3_CLASSIC') 
    time=dataset.createDimension('time',None)
    level=dataset.createDimension('level',ntype)
    
    pca=dataset.createVariable('pc',np.double,('time','level'),fill_value=-9999.)
    
    times=dataset.createVariable('time',np.float32,('time',),fill_value=-999.)
    levels=dataset.createVariable('level',np.float32,('level',),fill_value=-999.)
    
    pca.units='unitless'
    times.units='year since 0850-01-01 00:00:00'
    levels.long_name='number of PC'
    
    levels[:]=np.arange(1,ntype+1)
    pca[:]=pc
    times[0:nt]=time_m
    dataset.close()
	


	
def save_netcdf_var(eof_m, expvar, lat_m, lon_m, titles):
    
    nt=eof_m.shape[0]
    ntype=eof_m.shape[1]
    time_m=np.arange(0,nt)
    nlat=len(lat_m)
    nlon=len(lon_m)
    
    dataset =Dataset(titles,'w',format='NETCDF3_CLASSIC') 
    time=dataset.createDimension('time',None)
    level=dataset.createDimension('level',ntype)
    lat=dataset.createDimension('lat',nlat)
    lon=dataset.createDimension('lon',nlon)
    
    eof=dataset.createVariable('eof',np.double,('time','level','lat','lon'),fill_value=-9999.)
    pc_variance=dataset.createVariable('pc_variance',np.double,('time','level'),fill_value=-9999.)
    
    times=dataset.createVariable('time',np.float32,('time',),fill_value=-999.)
    levels=dataset.createVariable('level',np.float32,('level',),fill_value=-999.)
    lons=dataset.createVariable('lon',np.float32,('lon',),fill_value=-999.)
    lats=dataset.createVariable('lat',np.float32,('lat',),fill_value=-999.)
    
    eof.units='unitless'
    eof.long_name='Empirical orthogonal function'
    pc_variance.units='%'
    pc_variance.long_name='Percentage of variance explained by PCs'
    times.units='year since 0850-01-01 00:00:00'
    
    lats.units='degree_east'
    lons.units='degree_north'
    levels.long_name='type of EOF. In order: all droughts, initiation, transition, termination, termination+1'
    
    lats[:]=lat_m
    lons[:]=lon_m 
    levels[:]=np.arange(1,ntype+1)
    eof[:]=eof_m
    pc_variance[:]=expvar
    times[0:nt]=time_m
    dataset.close()



def save_netcdf_time(pc, titles):  #time, neof, types
    
    nt=pc.shape[0]
    ntype=pc.shape[1]
    time_m=np.arange(0,nt)
    
    dataset =Dataset(titles,'w',format='NETCDF3_CLASSIC') 
    time=dataset.createDimension('time',None)
    level=dataset.createDimension('level',ntype)
    
    pca=dataset.createVariable('pc',np.double,('time','level'),fill_value=-9999.)
    
    times=dataset.createVariable('time',np.float32,('time',),fill_value=-999.)
    levels=dataset.createVariable('level',np.float32,('level',),fill_value=-999.)
    
    pca.units='unitless'
    times.units='year since 0850-01-01 00:00:00'
    levels.long_name='number of PC'
    
    levels[:]=np.arange(1,ntype+1)
    pca[:]=pc
    times[0:nt]=time_m
    dataset.close()



def save_netcdf_simple(zval, sval, lat_m, lon_m, titles):
    
    nt=zval.shape[0]
    time_m=np.arange(0,nt)
    nlat=len(lat_m)
    nlon=len(lon_m)
    dataset =Dataset(titles,'w',format='NETCDF3_CLASSIC')
    time=dataset.createDimension('time',None)
    lat=dataset.createDimension('lat',nlat)
    lon=dataset.createDimension('lon',nlon)
    
    S=dataset.createVariable('soil',np.double,('time','lat','lon'),fill_value=-9999.)
    Z=dataset.createVariable('zg',np.double,('time','lat','lon'),fill_value=-9999.)
    times=dataset.createVariable('time',np.float32,('time',),fill_value=-999.)
    
    lons=dataset.createVariable('lon',np.float32,('lon',),fill_value=-999.)
    lats=dataset.createVariable('lat',np.float32,('lat',),fill_value=-999.)
    times.units='year since 0850-01-01 00:00:00'
    times.long_name='drought year'
    S.units='mm'
    Z.units='gpm'
    lats.units='degree_east'
    lons.units='degree_north'
    lats.axis='Y'
    lons.axis='X'
    times.axis='Z'

    lats[:]=lat_m
    lons[:]=lon_m
    S[:]=tval
    Z[:]=zval
    times[0:nt]=time_m
    dataset.close()




def weighting(zvar):
    nt=zvar.shape[0]
    nlon=zvar.shape[-1]
    
    weight_lat=np.sqrt(np.cos(np.pi*lats_nh/180.))
    z_weighted=np.zeros(zvar.shape)
    for t in range(nt):
        for j in range(nlon):
            z_weighted[t,:,j]=zvar[t,:,j]*weight_lat
        
    return(z_weighted)


def dec(x):
    return(np.round(x,2))
    
