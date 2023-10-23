"""
  Step 3 in the analysis:
  
  Clustering with k over N EOF/PC fields for each of the regions: 

  Testing which N and k to choose
  Each region and period needs to be run separately

  Gather Z and S based on cluster labels and generate means of these values for each cluster. 
  
  ====== Input files ===
  
  (They need to be produced ahead)
  1) PC of Z500 during droughts (needs to have from1 to n>N PCs)  (.nc)
  2) Z500 anomaly  (.nc)
  3) Soil anomaly (.nc)
  4) Drought indices (.npy)
  5) strongs for: model, region, number of eof

"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import seaborn as sns

import matplotlib.cm as cm

from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score

import base_tools as bt

"""
   Input information from the terminal
   model name, region, number of pc/eof
"""

model= sys.argv[1]   # Insert model name & other information from the termninal
region = sys.argv[2] 
neof = int(sys.argv[3])
period = 'all' #sys.argv[4]


filedir = './'
savedir = './'

figtype='.png'


if neof==5:
    ner=2
elif neof==6:
    ner=3
elif neof==7:
    ner=4
elif neof==8:
    ner=5


#==== Read file : PC save from the EOF (calculation_step2)

#=== PC 
m_data=Dataset(filedir + model+'.PC_all_zg.'+region+'_annual_lm.nc')
pc_all_mmedi=m_data.variables['pc'][:]

# ZG
data_z = Dataset(maindir+'ymean.anomaly.'+model+'.z500.0850-1849.nc','r')
z_all_mmedi=mmedi_data.variables['zg'][:]

# SOIL
data_s = Dataset(maindir+'ymean.anomaly.'+model+'.soil.0850-1849.nc','r')
s_all_mmedi=mmedi_data.variables['soil'][:]


# Drought indices
medi_index=np.load(filedir+ model + '_' + region + '_ymean_0850-1849_sorted_index_lm.npy',allow_pickle=True)[0]
idx_all = [x for y in medi_index for x in y]

lat_all = data_z.variables['lat'][:]
lon_all = data_z.variables['lon'][:]

nlat = len(lat_all)
nlon = len(lon_all)

lat_s = data_s.variables['lat'][:]
lon_s = data_s.variables['lon'][:]

nlat_s = len(lat_t)
nlon_s = len(lon_t)


#----------------------------------------
#  Cluster calculation
#----------------------------------------

class cluster_number():
    
    def __init__(self, data, ncluster, std='no'):
        
        self.data=data
        
        cluster = KMeans(n_clusters=ncluster, random_state=10, max_iter=1200)
        z_cluster = cluster.fit(self.data)
        z_labels = z_cluster.labels_
        z_center = cluster.cluster_centers_
        
        cluster_frac=[]
        for i in range(ncluster):
            cluster_frac.append((z_labels==i).sum()/(len(z_labels)*1.0))
        
        self.ncluster = ncluster
        self.clabel = z_labels
        self.centroid = z_center
        self.fraction = cluster_frac
        self.distortion = cluster.inertia_



#-------------------------------------------
#  Silhouette analysis from https://scikit-learn.org/stable/auto_examples/cluster/plot_kmeans_silhouette_analysis.html
#-------------------------------------------

#-- calculating silhouette_score for a reange of PCs (neof) and cluster number k (ncl). 

ncluster_range = range(3,10)
ncl = len(ncluster_range)

neof_range = range(3,11)
ne = len(neof_range)

sil_avg = np.empty((ncl,ne+1))*np.nan

for i in range(ncl):
    cvalue=ncluster_range[i]
    for j in range(ne):
        evalue=neof_range[j]
        z_cluster = cluster_number(pc_all_mmedi[:,0:evalue],cvalue)
        avg = silhouette_score(pc_all_mmedi[:,0:evalue], z_cluster.clabel)
        sil_avg[i,j] = avg


print(sil_avg)



#--- Maximum in EOF=5 (depending on the model/region)

sil_eof = sil_avg_mm[:,ner]
ic = np.argmax(sil_eof)

#ic=1
cluster_max = ncluster_range[ic]

print(sil_eof)


#--------------------------------------
#  Clustering 
#-------------------------------------- 

ncluster=cluster_max

z_all_cluster_mmedi = cluster_number(pc_all_mmedi[:,0:neof],ncluster)
z_all_clabel_mmedi = z_all_cluster_mmedi.clabel
z_all_centroid_mmedi = z_all_cluster_mmedi.centroid
z_all_fraction_mmedi = z_all_cluster_mmedi.fraction
z_all_distance = z_all_cluster_mmedi.distortion

zmean_mmedi = np.zeros((ncluster,nlat,nlon))
smean_mmedi = np.zeros((ncluster,nlat_t,nlon_t))

for i in range(ncluster):
    zmean_mmedi[i] = np.nanmean(z_all_mmedi[z_all_clabel_mmedi==i],axis=0)
    tmean_mmedi[i] = np.nanmean(s_all_mmedi[z_all_clabel_mmedi==i],axis=0)

#-- save clusters in netcdf

titles=savedir + model+'.'+region+'ZS.gathered.cluster_ordered.nc'

bt.save_netcdf_simple(zmean_mmedi, tmean_mmedi, lat_all, lon_all, titles)

np.save(savedir +model+'_Clabel_'+region+'_labels.npy',z_all_clabel_mmedi)


