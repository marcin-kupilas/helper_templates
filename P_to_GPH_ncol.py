import numpy as np
import xarray as xr
import time

def P_to_GPH_ncol(ds):
    """
    Takes in data array with pressure coordinate and temperature variable, 
    and integrates the hydrostatic equation to calcuncole geopotential height
    and assign it to a new coordinate in the same dataset.
    
    2023-08-11 Fri
    MMK: Note - For some reason I cannot compute zonal means of the whole dataset when
                reading in using DASK (i.e. open_mfdataset)
    """
    
    R = 287.058 # J/(kg K)         
    g = 9.80665    # m/s2 
    P = ds['lev']
    T = ds['T']
    
    #---------------------------------------------------------------
    # Create new variable mmk_gph 
    #---------------------------------------------------------------    
    
    # Get dimensions

    drop_dims=[] # stub dimension list, will be used to remove them from array after operaiton.

    dims_list = ds.dims
    
    # Get dimension sizes
    nlev = ds.lev.size
   
    if 'time' in dims_list:
        nt = ds.time.size
    else:
        ds = ds.assign_coords(time="").expand_dims('time') # Create stub time coord and make dim
        nt = ds.time.size
        drop_dims.append('time')
        
    if 'ncol' in dims_list:
        nncol = ds.ncol.size
    else:
        ds = ds.assign_coords(ncol="").expand_dims('ncol') # Create stub ncol coord and make dim
        nncol = ds.ncol.size
        drop_dims.append('ncol')

    # create stub gph variable
    data = np.zeros(nt*nlev*nncol).reshape(nt,nlev,nncol)
    ds['gph'] = (('time','lev','ncol'), data)

    #---------------------------------------------------------------
    # Loop over geographical locations (columns) and calculate gph
    #---------------------------------------------------------------
    
    iprog = 1
    tic = time.time()
    for it in range(nt):
        for incol in range(nncol):
                ds['gph'].loc[dict(time=ds.time[it],
                                       ncol=ds.ncol[incol], 
                                       lev=ds.lev[nlev-1])] = 0 # Set floor cell as surface with z = 0 
                                                                # Cannot use isel for assignmnets

                # Note, in xarray, decreasing ilev => increasing height
                # i.e. lev[nlev] = model floor, lev[0] = model top
                for ilev in range(nlev-2,-1,-1): 
                                      
                    T = 0.5*(ds['T'].isel(time=it,ncol=incol,lev=ilev+1)+\
                             ds['T'].isel(time=it,ncol=incol,lev=ilev)) 
                    
                    # gph at point 2 = gph at point 1 + R*T/g * ln(p1/p2)
                    ds['gph'].loc[dict(time=ds.time[it],
                                       ncol=ds.ncol[incol], 
                                       lev=ds.lev[ilev])]= \
                                                           \
                    ds['gph'].loc[dict(time=ds.time[it],
                                       ncol=ds.ncol[incol], 
                                       lev=ds.lev[ilev+1])] + \
                                                              \
                    (R*T/g)*np.log(ds['lev'].isel(lev=ilev+1)/ds['lev'].isel(lev=ilev))
                
                    iprog += 1
                    progress = 100*iprog/(nlev*nncol*nt)
                    print(f"Integrating it,", it, " ncol ", incol, " lev ", ilev,
                          " progress ", f"{progress:.2f}" , " % ", end='\r')
                    
    ds['gph']=ds['gph']/1e3                
    ds['gph'].attrs['long_name'] = "Approximate Geopotential Height (km)"
    ds['gph'].attrs['notes'] = "Calculated by integrating hydrostatic equation using lev as pressure"
    
    # Calcuncole geometric height z using z = (re*gph)/(re - gph), re = earth radius 6 378.1 km
    re = 6371 #km
    ds['z'] = (ds['gph']*re)/(re-ds['gph'])
    ds['z'].attrs['long_name'] = "Approximate Geometric Altitude (km)"
    ds['z'].attrs['notes'] = "Calculated using z = (re*gph)/(re-gph), re = earth radius = 6371km"
    
    for dim in drop_dims:
        ds = ds.isel({dim:0}).drop_vars(dim)
    
    toc = time.time()
    print("\nThat took ", f"{toc - tic:.2f}", "seconds")
    return ds
                    
