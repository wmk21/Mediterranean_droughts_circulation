"""
    Step 2 in the analysis:
    PC/EOF analysis for each models
	
	=== Input file ===
	
	- Z500 anomalies
	- Drought indices (depending on the region)  
	
"""

import numpy as np
from netCDF4 import Dataset
import copy

from sklearn.decomposition import PCA
from sklearn import preprocessing

import base_tools as bt

"""
   Input information
"""

varname = 'zg'
model = sys.argv[1] # model name
region = sys.argv[2] # east or west
std = 'yes'
neof = 15

maindir = './'
savedir=filedir


#==== Reading files

medi_index=np.load(filedir+ model + '_' + region + '_ymean_0850-1849_sorted_index_lm.npy',allow_pickle=True)[0]

idx_all = [x for y in medi_index for x in y]

data_z = Dataset(maindir+'ymean.anomaly.'+model+'.z500.0850-1849.nc','r')
tmp = data_z.variables[varname][:]
z = np.squeeze(tmp)

lats=data_z.variables['lat'][:]
lons=data_z.variables['lon'][:]

#==== Selecting the area to apply PC in NH 
    
lons_new=lons.copy()
    
for n in range(len(lons_new)):
    if lons[n]>180.:
        lons_new[n]=lons[n]-360.

nh_lat=np.where((lats>=20.) & (lats<=85))[0]
nh_lon=np.where((lons_new>=-70.) & (lons_new<=70.))[0]

lats_nh=lats[nh_lat]
lons_nh=lons[nh_lon]
nlat=len(lats_nh)
nlon=len(lons_nh)

z_dr = z_sel[i_all]


#=== Defining EOF function. 

class pc_number():
    
    def __init__(self, data,neof,std='no'):
	
        lenlat=data.shape[1]
        lenlon=data.shape[2]
        lentim=data.shape[0]
        self.data=data
        if std=='no':
            self.data_reshape=np.reshape(data, (lentim, lenlat* lenlon), order='F')
        elif std=='yes':
            data_shape = np.reshape(data, (lentim, lenlat* lenlon), order='F')
            scaler  = preprocessing.StandardScaler()
            scaler_data = scaler.fit(data_shape)
            self.data_reshape = scaler_data.transform(data_shape)
            
        self.neof=neof
        
        pca = PCA()
        pca.fit(self.data_reshape)
        pca_variance = pca.explained_variance_ratio_
        PCs = pca.transform(self.data_reshape)
        PCs = PCs[:,:self.neof]
        
        EOFs = pca.components_
        EOFs = EOFs[:self.neof,:]
        EOF_pattern = np.reshape(EOFs, (self.neof, lenlat, lenlon), order='F')
        
        self.pc_variance_all = pca_variance
        self.pc_variance = pca_variance[0:self.neof]*100
        self.pc = PCs
        self.eof = EOF_pattern
    


z_all = bt.weighting(z_dr)

eof_all_data = pc_number(z_all, neof ,std)
pc_all = eof_all_data.pc_variance
pc_time_all = eof_all_data.pc
eof_all = eof_all_data.eof

title_pc = savedir + model+'.PC_all_zg.'+region+'_annual_lm.nc'
bt.save_netcdf_time(pc_time_all, title_pc)

title_eof = savedir + model+'.EOF_all_zg.'+region+'_annual_lm.nc'
bt.save_netcdf_var(eof_all, pc_time_all, lats_nh, lons_nh, title_eof)

