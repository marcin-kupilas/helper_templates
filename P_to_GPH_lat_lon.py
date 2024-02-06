# Generated Z3 via integrating https://en.wikipedia.org/wiki/Hypsometric_equation
import numpy as np
import xarray as xr
import time

def P_to_GPH_lat_lon(ds):
    """
    Takes in data array with pressure coordinate and temperature variable, 
    and integrates the hydrostatic equation to calculate geopotential height
    and assign it to a new coordinate in the same dataset.
    
    See: https://en.wikipedia.org/wiki/Hypsometric_equation
    
    Author: Marcin Kupilas
    2023-08-11 Fri
    Marcin Kupilas: Note - For some reason I cannot compute zonal means of the whole dataset when
                reading in using DASK (i.e. open_mfdataset)
    """
    print("Using P to GPH for structured lat lon grid")
    print("WARNING \n\
    - 99% accurate - off by about -500m in the MLT from WACCM Z3 \n\
    - Takes a long time for large datasets \n\
    - Any masked dimensions will be removed \n\
      e.g. ds.sel(lat=90...) removes lat as dim and makes it coord. \n\
      This script drops all coordinates that aren't dimensions.")

    R = 287.058 # J/(kg K)         
    g = 9.80665    # m/s2 
    P = ds['lev']
    
    #---------------------------------------------------------------
    # Create new variable mmk_gph 
    #---------------------------------------------------------------   
    
    # Get dimensions
    dims_list = ds.dims
    # Get dimension sizes
    nlev = ds.lev.size
    
    drop_dims=[] # stub dimension list, will be used to remove them from array after operaiton.
    if 'time' in dims_list:
        nt = ds.time.size
    else:
        ds = ds.assign_coords(time="").expand_dims('time') # Create stub time coord and make dim
        nt = ds.time.size
        drop_dims.append('time')
        
    if 'lat' in dims_list:
        nlat = ds.lat.size
    else:
        ds = ds.assign_coords(lat="").expand_dims('lat') # Create stub lat coord and make dim
        nlat = ds.lat.size
        drop_dims.append('lat')

    if 'lon' in dims_list:
        nlon = ds.lon.size
    else:
        ds = ds.assign_coords(lon="").expand_dims('lon') # Create stub lon coord and make dim
        nlon = ds.lon.size
        drop_dims.append('lon')

    # create stub gph variable
    data = np.zeros(nt*nlev*nlat*nlon).reshape(nt,nlev,nlat,nlon)
    ds['gph'] = (('time','lev','lat','lon'), data)

    #---------------------------------------------------------------
    # Loop over geographical locations (lat lon) and calculate gph
    #---------------------------------------------------------------
    
    iprog = 1
    tic = time.time()
    for it in range(nt):
        for ilat in range(nlat):
            for ilon in range(nlon):
                ds['gph'].loc[dict(time=ds.time[it],
                                       lat=ds.lat[ilat], 
                                       lon=ds.lon[ilon],                   
                                       lev=ds.lev[nlev-1])] = 0 # Set floor cell as surface with z = 0 
                                                                # Cannot use isel for assignmnets

                # Note, in xarray, decreasing ilev => increasing height
                # i.e. lev[nlev] = model floor, lev[0] = model top
                for ilev in range(nlev-2,-1,-1): 
                                      
                    T = 0.5*(ds['T'].isel(time=it,lat=ilat,lon=ilon,lev=ilev+1)+\
                             ds['T'].isel(time=it,lat=ilat,lon=ilon,lev=ilev)) 
                    
                    # gph at point 2 = gph at point 1 + R*T/g * ln(p1/p2)
                    ds['gph'].loc[dict(time=ds.time[it],
                                       lat=ds.lat[ilat], 
                                       lon=ds.lon[ilon],                   
                                       lev=ds.lev[ilev])]= \
                                                           \
                    ds['gph'].loc[dict(time=ds.time[it],
                                       lat=ds.lat[ilat], 
                                       lon=ds.lon[ilon],                   
                                       lev=ds.lev[ilev+1])] + \
                                                              \
                    (R*T/g)*np.log(ds['lev'].isel(lev=ilev+1)/ds['lev'].isel(lev=ilev))
                
                    iprog += 1
                    progress = 100*iprog/(nlev*nlon*nlat*nt)
                    print(f"Integrating it,", it, " lat ", ilat, " lon ", 
                          ilon, " lev ", ilev,
                          " progress ", f"{progress:.2f}" , " % ", end='\r')
                    
    ds['gph']=ds['gph']/1e3                
    ds['gph'].attrs['long_name'] = "Approximate Geopotential Height (km)"
    ds['gph'].attrs['notes'] = "Calculated by integrating hydrostatic equation using lev as pressure"
    
    # Calculate geometric height z using z = (re*gph)/(re - gph), re = earth radius 6 378.1 km
    re = 6371 #km
    ds['z'] = (ds['gph']*re)/(re-ds['gph'])
    ds['z'].attrs['long_name'] = "Approximate Geometric Altitude (km)"
    ds['z'].attrs['notes'] = "Calculated using z = (re*gph)/(re-gph), re = earth radius = 6371km"
    
    for dim in drop_dims:
        ds = ds.isel({dim:0}).drop_vars(dim)
    
    toc = time.time()
    print("\nThat took ", f"{toc - tic:.2f}", "seconds")
    return ds
                    

    
"""
Notes:

Can create list of all dim names using:
a = [item for item in ds1.coords.keys()]

Can drop dims using

    
#     for dim in drop_dims:
#         ds = ds.isel({dim:0}).drop_vars(dim)

for some reason drop_dims(dim) method gets rid of data, so haven't used that.

"""