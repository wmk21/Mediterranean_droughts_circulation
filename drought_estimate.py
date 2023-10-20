"""

    Estimate drought events defined over Mediterranean region
    
    Two functions: drought_calculation, drought_sorting

"""


"""
    1) drought_calculation:
        First drought estimation setting conditions for droughts of two negative soil moisture anomalies, and spatial threshold of reg_per. 
       Input: 4 inputs
        - Spatially averaged time series of Soil moisture anomalies [dim=time], 
        - Spatial soil moisture anomalies [dim= time x lat x lon]
        - landmask files (0=ocean, 1=land)
        - a spatial threshold for negative soil moisture anomalies (50-60%)
"""


import numpy as np

def drought_calculation(series, series_reg, landmass, reg_per=0.):
    """
    #---- Drought starts with two consecutive negative years, then it stop with two wet period It should be continuously negative to be considered as droughts. 
    """
    nt=len(series)
    series_neg=np.where(series > 0., np.zeros(nt), 1)
    
    #--- Take continuous negative anomalies. 
    x=0.
    series_long=np.zeros((nt))
    series_ind=[]
    series_int=[]
    
    for i in range(nt):
        if series_neg[i]==1.:
            x=x+1
            series_long[i]=x
            series_int.append(i)
        else:
            x=0.
            if len(series_int)!=0:        
                series_ind.append(series_int)
            series_int=[]

    ns=len(series_ind)
    
    #--- Calculate percentages of land grid under negative soil. 
    def reg_percentage(regional_value):
        reg_val_tmp=regional_value.flatten()
        land_tmp=landmass.flatten()
        reg_val_land=reg_val_tmp[land_tmp==1]
        nland=float(len(land_tmp[land_tmp==1]))
        reg_val=reg_val_land[~np.isnan(reg_val_land)]
        reg_per=float(len(reg_val[reg_val<0]))/nland
        return(reg_per*100.)

    #--- Sort out only valid drought clusters. --> initial two values need to be negative. 
    #--- Calculate Duration, all values, regional percentages of each drought clusters.
    
    duration=[]   # list of the clusters with duration of events
    drought_val=[]   # list of the clusters with indices
    drought_per=[]   # list of the percentages of land coverages
    drought_index=[]  # list of the indices of clusters

    for k in range(ns):       #cluster numbers
        kval=series_ind[k]    #indices
        index_value=series[kval]  
        if len(index_value)>2:   # at least two continuous negative values=Drought. 
            drought_val.append(index_value)
            duration.append(max(series_long[kval])) 
            reg_tmp=series_reg[kval]
            reg_v=[]
            for i in range(len(kval)):
                rval=reg_percentage(reg_tmp[i])
                reg_v.append(rval)
                del(rval)
            drought_per.append(reg_v)
            drought_index.append(kval)
            del(reg_tmp,reg_v)  
        del(kval,index_value)
    
    #--- Setting a threshold 
    if reg_per!=0:
        duration_final=[]   
        drought_val_final=[]  
        drought_per_final=[]  
        drought_index_final=[] 
        n=len(drought_per)  

        for i in range(n):
            Dp=drought_per[i]
            per_iden=np.array(Dp)[np.array(Dp)<reg_per]
            if len(per_iden)==0.:
                Dv=drought_val[i]
                Du=duration[i]
                Di=drought_index[i]
                
                drought_val_final.append(Dv)
                drought_per_final.append(Dp)  
                drought_index_final.append(Di)
                duration_final.append(Du)
    
    elif reg_per==0.:
        duration_final=duration 
        drought_val_final=drought_val  
        drought_per_final=drought_per  
        drought_index_final=drought_index

    # Return the final drought indices, the mean duration of each clusters, the values of duration, the percentage of duration, cumulative duration. 
    return([drought_index_final, duration_final, drought_val_final, drought_per_final, series_long])
    

"""

    2) drought_sorting: 
    
    Set the next threshold for drought: droughts are terminated when they encounter 2 consecutive years of positive anomalies. And it identifies west and east regional droughts. 
    
    Uses the Output from drought_estimate as Input (drought_index_final, and drought_per_final)
    
     It requires 4 inputs, 2 for each region (Western and Eastern Mediterranean regions): 
    - Original drought indices for eastern Mediterranean
    - Percentage of drought coverages for eastern Mediterranean
    - Original drought indices for western Mediterranean
    - Percentage of drought coverages for western Mediterranean
    
"""

def drought_sorting(emedi_all_index, emedi_all_perc, wmedi_all_index, wmedi_all_perc):

    emedi_sorted_index=[]
    emedi_sorted_perc=[]

    for i,e in enumerate(ensemble):
        drought_e=emedi_all_index[i]
        drought_p=emedi_all_perc[i]
        nd=len(drought_e)
        d_i=[]
        p_i=[]
        for n in range(nd):
            drought_i=drought_e[n]
            drought_ip=drought_p[n]
            n_1=drought_i[-1]+1
            n_2=drought_i[-1]+2
            if n_1>=ny or n_2>=ny:
                d_i.append(drought_i)
                p_i.append(drought_ip)
            else:
                if soil_re[i,n_1]>0 and soil_re[i,n_2]>0:
                    d_i.append(drought_i)
                    p_i.append(drought_ip)
            del(drought_i,drought_ip)
        emedi_sorted_index.append(d_i)
        emedi_sorted_perc.append(p_i)
        
    
    # WMEDI
    wmedi_sorted_index=[]
    wmedi_sorted_perc=[]

    for i,e in enumerate(ensemble):
        drought_e=wmedi_all_index[i]
        drought_p=wmedi_all_perc[i]
        nd=len(drought_e)
        d_i=[]
        p_i=[]
        for n in range(nd):
            drought_i=drought_e[n]
            drought_ip=drought_p[n]
            n_1=drought_i[-1]+1
            n_2=drought_i[-1]+2
            if n_1>=ny or n_2>=ny:
                d_i.append(drought_i)
                p_i.append(drought_ip)
            else:
                if soil_rw[i,int(n_1)]>0 and soil_rw[i,int(n_2)]>0:
                    d_i.append(drought_i)
                    p_i.append(drought_ip)
            del(drought_i,drought_ip)
        wmedi_sorted_index.append(d_i)
        wmedi_sorted_perc.append(p_i)
    
    """
     Finding indices that are coherent between wmedi and emedi. 
    
      Criteria: Match of at least 3 years (approximately >=60-80% of drought periods. This value needs to be change depending on the mean % of overlapped periods of one's interests.). 
      A negative soil moisture anomaly condition can starts earliear or can last one year longer, but no longer than that. The final index can be obtained by merging both indices (of each region). 
    
     wonly and eonly: when they show no match between them. 
    (These are the regional drought values to use for each region.)
    """

    wmedi_only=[]
    emedi_only=[]

    ew_coh=[]

    wperc_only=[]
    eperc_only=[]

    for i in range(ne):
        # wmedi
        wd_i=wmedi_sorted_index[i]
        wp_i=wmedi_sorted_perc[i]
        nwm=len(wd_i)
    
        # emedi
        ed_i=emedi_sorted_index[i]
        ep_i=emedi_sorted_perc[i]
        nem=len(ed_i)
    
        dwco_all=[]
        deco_all=[]
        dco_3=[]
    
        # for emedi first
        for j in range(nwm):
            wcluster=wd_i[j]
            for k in range(nem):
                ecluster=ed_i[k]
                coherence=list(set(wcluster).intersection(ecluster))
                # one way to define +-1 year
                if len(coherence)>=3:   # set N years = 3 (around )
                    if wcluster[0]==ecluster[0] or wcluster[0]==ecluster[0]+1 or wcluster[0]==ecluster[0]-1:
                        if wcluster[-1]==ecluster[-1] or wcluster[-1]==ecluster[-1]+1 or wcluster[-1]==ecluster[-1]-1:
                            all_d=list(dict.fromkeys(wcluster+ecluster))
                            dco_3.append(all_d)
                            del(all_d)
                #-----------------------
                if len(coherence)>=1: # 2 or 1 depending on my choices
                    deco_all.append(ecluster)
                del(ecluster)
            del(wcluster)
    
        #---- emedi
        for j in range(nem):
            ecluster=ed_i[j]
            for k in range(nwm):
                wcluster=wd_i[k]
                coherence=list(set(wcluster).intersection(ecluster))
                if len(coherence)>=1:
                    dwco_all.append(wcluster)
                del(wcluster)
            del(ecluster)
    
        #----Eliminate repeated indices
        dwonly=[]
        pwonly=[]
        for l in range(nwm):
            x=wd_i[l]
            y=wp_i[l]
            if x not in dwco_all:
                dwonly.append(x)
                pwonly.append(y)
            del(x,y)
    
        deonly=[]
        peonly=[]
        for l in range(nem):
            x=ed_i[l]
            y=ep_i[l]
            if x not in deco_all:
                deonly.append(x)
                peonly.append(y)
            del(x,y)
    
        wmedi_only.append(dwonly)
        emedi_only.append(deonly)
        wperc_only.append(pwonly)
        eperc_only.append(peonly)
        ew_coh.append(dco_3)
        
    #--- Return sorted final indices for eastern and western regions. 
    #--- ew_coh is the indices with pan-mediterranean drought
    return [[emedi_only,eperc_only], [wmedi_only,wperc_only], ew_coh)
    
    

