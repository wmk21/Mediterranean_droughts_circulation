"""
    Step 4 and 5 in the analysis:
    Clustering, then grouping based on the correlation coefficients
	
	=== Input file ===
	
	- Z500 of clusters, gathered all together (needs to have saved the data-labels (with model, region and period).)
	
"""

import xarray as xr
import numpy as np
from scipy import stats 

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score


#=== Reading all files:

filedir = './'

data_zg = xr.open_dataset(filedir + model+'.' +region+ 'ZS.gathered.cluster_ordered.nc')
ntot = data_only_zg.cluster.size #71
zg = data_zg.zg

#------------------
#   Clustering Funcion
#------------------

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


#----------------------
#   Shilhouette coefficient
#----------------------

ncluster_range = range(3,50)
ncl = len(ncluster_range)

sil_avg = np.empty((ncl))
sil_avg[:]=np.NaN


for i in range(ncl):
    cvalue=ncluster_range[i]
    z_cluster = cluster_number(zg, cvalue)
    avg = silhouette_score(data_only_zg, z_cluster.clabel)
    sil_avg[i]=avg


#----------------------
#   Clustering
#----------------------

xmax = np.argmax(sil_avg)
k = ncluster_range[xmax]

clustered = cluster_number(data_only_zg, k)
clabel = clustered.clabel

len_val=[]
for i in range(k):
    len_val.append(len(clabel[clabel==i]))



#------------------
#   Correlation coefficient process
#------------------


index=np.arange(0,ntot)

def matrix_detection(k):
    z_val = zg[clabel==k]
    index_c = index[clabel==k]
    nz = z_val.shape[0]
    
    r_mat = np.empty((nz,nz))*np.nan
    p_mat = np.empty((nz,nz))*np.nan
    i_mat = np.empty((nz,nz))*np.nan
	
    for i in range(nz):
        for j in range((nz)):
            if i<=j:
                r = stats.pearsonr(z_val[i], z_val[j])
                r_mat[i,j] = r[0]
                p_mat[i,j] = r[1]
                i_mat[i,j] = j
    
    tmp_r = r_mat.flatten()
    tmp_p = p_mat.flatten()
    r_sort = tmp_r[~np.isnan(tmp_r)]
    p_sort = tmp_p[~np.isnan(tmp_p)]
    tmp_i = i_mat.flatten()
    i_sort = tmp_i[~np.isnan(tmp_i)]
    
    #--- Excluded variables
    selval = i_sort[(r_sort<0.) | (p_sort>0.05)]  # Excluding negative and non sig. 
    excluded0 = np.unique(selval).astype(int).tolist()
    
    #--- included 
    selval_i = i_sort[(r_sort>0.) & (p_sort<=0.05)]
    included_tmp = np.unique(selval_i).astype(int).tolist()
    included0 = list(set(included_tmp) ^ set(excluded0))
    
    #--- valid R values. 
    selval_valid_r = r_sort[(p_sort<=0.05) & (r_sort>0.)]
    r0 = np.mean(selval_valid_r[selval_valid_r<1.])
    print('Cluster '+str(int(k))+' Excluded: '+str(excluded0))
    print('Cluster '+str(int(k))+' Included: '+str(included0))
    print('Cluster '+str(int(k))+' mean R all: '+str(np.mean(r_sort)))
    
    return [index_c[excluded0], index_c[included0],r0]


#=== r_threshold: set the mean correlation from the previous step. 


def side_detection(ind, r_thred=0.1)):
    z_val = zg[ind]
    index_c = ind
    nz = z_val.shape[0]
    
    r_mat = np.empty((nz,nz))*np.nan
    p_mat = np.empty((nz,nz))*np.nan
    i_mat = np.empty((nz,nz))*np.nan
	
    for i in range(nz):
        for j in range((nz)):
            if i<=j:
                r = stats.pearsonr(z_val[i], z_val[j])
                r_mat[i,j] = r[0]
                p_mat[i,j] = r[1]
                i_mat[i,j] = j
    
    
    tmp_r = r_mat.flatten()
    tmp_p = p_mat.flatten()
    r_sort = tmp_r[~np.isnan(tmp_r)]
    p_sort = tmp_p[~np.isnan(tmp_p)]
    tmp_i = i_mat.flatten()
    i_sort = tmp_i[~np.isnan(tmp_i)]
    
    #--- Excluded variables
    selval = i_sort[(r_sort<r_thred) | (p_sort>0.05)]  # Excluding negative and non sig. 
    excluded0 = np.unique(selval).astype(int).tolist()
    
    #--- included 
    selval_i = i_sort[(r_sort>=r_thred) & (p_sort<=0.05)]
    included_tmp = np.unique(selval_i).astype(int).tolist()
    included0 = list(set(included_tmp) ^ set(excluded0))
    
    #--- valid R values. 
    selval_valid_r = r_sort[(p_sort<=0.05) & (r_sort>=r_thred)]
    r0 = np.mean(selval_valid_r[selval_valid_r<1.])
    print('Cluster '+str(int(k))+' Excluded: '+str(excluded0))
    print('Cluster '+str(int(k))+' Included: '+str(included0))
    print('Cluster '+str(int(k))+' mean R all: '+str(np.mean(r0)))
    
    if len(excluded0)==0:
        e1=[]
    else:
        e1=np.array(index_c)[excluded0]
    if len(included0)==0:
        e2=[]
    else:
        e2=np.array(index_c)[included0]
    return [e1, e2, r0]



def corr_r(z_val):  
    # simply calculation without selecting anything.
    nz=z_val.shape[0]
    r_mat = np.empty((nz,nz))*np.nan
    p_mat = np.empty((nz,nz))*np.nan
	
    for i in range(nz):
        for j in range((nz)):
            if i<j:
                r = stats.pearsonr(z_val[i], z_val[j])
                p_mat[i,j]=r[1]
                r_mat[i,j] = r[0]
    
    tmp_r = r_mat.flatten()
    r_sort = tmp_r[~np.isnan(tmp_r)]
    tmp_p = p_mat.flatten()
    p_sort = tmp_p[~np.isnan(tmp_r)]
    
    #--- valid R values. 
    selval_valid_r = r_sort[(p_sort<=0.05) & (r_sort>=0)]
    r0 = np.mean(selval_valid_r[selval_valid_r<1.])
    
    return r0
   


#===== This process needs to be repeated until no more gathering is possible. Use "side_detection" in the next step. 

nc = np.max(clabel)

valid_clusters=[]
excluded_clusters=[]
rvalid_clusters=[]

for i in range(nc+1):
    if len_val[i]>1:
        tmp = matrix_detection(i)
        valid_clusters.append(tmp[1].tolist())
        excluded_clusters.append(tmp[0].tolist())
        rvalid_clusters.append(tmp[2])
    else:
        excluded_clusters.append(index[clabel==i].tolist())


