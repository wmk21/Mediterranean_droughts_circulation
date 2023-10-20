"""
   Step 1 in the analysis: get drought events, and save the indices to *.npy
   
   Input files: soil moisture anomalies, and land masking files over western and eastern Mediterranean. 
   Requires a pre-process soil moisture data, interpolated at 70cm, anomalies calculated, and separated by west and east Mediterranean. 
   
"""

import numpy as np
from netCDF4 import Dataset

from drought_estimate import drought_calculation, drought_sorting

maindir = './'
savedir = './'

"""
   Input information
"""

model='GISS-E2-R'
region='wmedi'
mean_type='ymean'  ## read smean or ymean
period = '0850-1949'
varname='mrlsl'


"""
   Reading soil moisture and land mask for west and east regions. 
"""

#==== Read the soil file

# West
reg_all_west = 'west.'+mean_type+'.anomaly.'+model+'.SOILIQ.'+period+'.nc'

print(reg_all_west)

reg_all_east = 'east.'+mean_type+'.anomaly.'+model+'.SOILIQ.'+period+'.nc'
print(reg_all_east)

var_tmp_west = Dataset(maindir + reg_all_west,'r')
series_all_west = var_tmp_west.variables[varname][:]
series_all_west = np.squeeze(series_all_west)

reg_mean_west = 'regmean.' + reg_all     # this is the regional averaged time series
var_tmp2_west = Dataset(maindir + reg_mean_west, 'r')
series_reg_west = var_tmp2_west.variables[varname][:]
series_reg_west = np.squeeze(series_reg_west)

# East
var_tmp_east = Dataset(maindir + reg_all_east,'r')
series_all_east = var_tmp_east.variables[varname][:]
series_all_east = np.squeeze(series_all_east)

reg_mean_east = 'regmean.' + reg_all     # this is the regional averaged time series
var_tmp2_east = Dataset(maindir + reg_mean_east, 'r')
series_reg_east = var_tmp2_east.variables[varname][:]
series_reg_east = np.squeeze(series_reg_east)


#==== Reading the land mask file

landfile_west = maindir + model + 'west.landfile.nc'
var_land_west = Dataset(landfile_west,'r')
landmass_west = var_land_west.variables['land'][:]
landmass_west = np.squeeze(landmass_west)

landfile_east = maindir + model + 'east.landfile.nc'
var_land_east = Dataset(landfile_east, 'r')
landmass_east = var_land_east.variables['land'][:]
landmass_east = np.squeeze(landmass_east)


#==== Drought estimation

Drought_west = drought_calculation(series_reg_west, series_all_west, landmass_west, reg_per=60.)
Drought_east = drought_calculation(series_reg_east, series_all_east, landmass_east, reg_per=60.)

sorted_indices = drought_sorting(Drought_east[0], Drought_east[-2], Drought_west[0], Drought_west[-2])

emedi_drought = sorted_indices[0]
wmedi_drought = sorted_indices[1]


np.save(savedir+ model +'_east_ymean_0850-1849_sorted_index_lm.npy',emedi_drought)

np.save(savedir+ model +'_west_ymean_0850-1849_sorted_index_lm.npy',wmedi_drought)



