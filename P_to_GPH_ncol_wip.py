import numpy as np
import xarray as xr

def P_to_GPH(ds):
    """
    Takes in data array with pressure coordinate and temperature variable, 
    and integrates the hydrostatic equation to calculate geopotential height
    and assign it to a new coordinate in the same dataset.
    
    2023-08-11 Fri
    MMK: Note - For some reason I cannot compute zonal means of the whole dataset when
                reading in using DASK (i.e. open_mfdataset)
    """
    
    R = 287.058 # J/(kg K)         
    g = 9.80665    # m/s2 
    P = ds['lev']
    T = ds['T']
    
    ## create new variable mmk_gph 
    
    # Get dimensions
    dims_list = ds.dims
    # Get dimension sizes
    nlev = ds.lev.size
    
    if 'time' in dims_list:
        nt = ds.time.size
    else:
        ds = ds.assign_coords(time="").expand_dims('time') # Create stub time coord and make dim
        nt = ds.time.size
        
    if 'lat' in dims_list:
        nlat = ds.lat.size
    else:
        ds = ds.assign_coords(lat="").expand_dims('lat') # Create stub lat coord and make dim
        nlat = ds.lat.size

    if 'lon' in dims_list:
        nlon = ds.lon.size
    else:
        ds = ds.assign_coords(lon="").expand_dims('lon') # Create stub lon coord and make dim
        nlon = ds.lon.size
    
    # create stub mmk_gph variable
    data = np.zeros(nt*nlev*nlat*nlon).reshape(nt,nlev,nlat,nlon)
    ds['mmk_gph'] = (('time','lev','lat','lon'), data)

#### TODO - MMK -  below is untrustworthy!

T = lambda it, incol, inlev : ds['T'].isel(time=it).isel(ncol=incol).isel(lev=inlev) 


# time stamp loop
for it in range(nt):
    # column loop
    for incol in range(ncol): 

        ds['mmk_z3'].loc[dict(time=ds.time[it], ncol=ds.ncol[incol], lev=ds.lev[nlev-1])] = 0 # Set floor cell as surface with z = 0 
        
        # level loop
        for inlev in range(nlev-2,-1,-1): 
            
            # Virtual temperature
            T_virt = 0.5*(T(it,incol,inlev) + T(it,incol,inlev+1))
            
            # Calculate z3
            ds['mmk_z3'].loc[dict(time=ds.time[it], ncol=ds.ncol[incol], lev=ds.lev[inlev])] = \
            ds['mmk_z3'].loc[dict(time=ds.time[it], ncol=ds.ncol[incol], lev=ds.lev[inlev+1])] + \
            ((R*T_virt/g))* \
            np.log(ds['lev'].isel(lev=inlev+1)/ds['lev'].isel(lev=inlev))

ds['mmk_z3'].to_netcdf(writepath + writefile)
