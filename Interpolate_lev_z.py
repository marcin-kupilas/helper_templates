import numpy as np
import xarray as xr
import time

def interpolate_lev_to_z_ncol(ds):
    toc = time.time()
    nlev = ds.lev.size
    ntime = ds.time.size
    
    createnew = False
    # create ncol dim if not in ds
    if "ncol" not in ds.dims: 
        print("ncol not in dims")
        print("nlev, ntime", nlev, ntime)
        createnew = True

    else:
        nncol = ds.ncol.size
        print("nlev, ntime, nncol", nlev, ntime, nncol)

    
    z_interp = np.arange(120,50,-0.25)
    ds_interp_icol_step = []
    ds_interp_itime_step = []
    counter = 0
    if createnew == False:
        for itime in range(ntime):
            toc1 = time.time()
            for icol in range(nncol):
                counter += 1
                # print("Interpolating, ", "{:.3f}".format(100*counter/(ntime*nncol)) + " %", end="\r")
                
                # Set up temp dataset replacing lev with z as indexed dimension
                ds_temp = ds.isel(time=itime,ncol=icol).swap_dims({'lev':'z'}).reset_coords("lev")
                # Interpolate lev to z 
                ds_temp = ds_temp.interp(z=z_interp)
                # Append 
                ds_interp_icol_step.append(ds_temp)
                
            # Concat
            ds_interp_icol_step = xr.concat(ds_interp_icol_step, dim="ncol")
            # Append
            ds_interp_itime_step.append(ds_interp_icol_step)
            # Reset
            ds_interp_icol_step = []
        
            print("Interpolating time", itime, "took", "{:.4f}".format(time.time() - toc1), "seconds")
            print("Completed", "{:.3f}".format(100*counter/(ntime*nncol)) + " %")
    else:
        for itime in range(ntime):
            print("Interpolating time", itime)
            toc1 = time.time()
            counter += 1
            # print("Interpolating, ", "{:.3f}".format(100*counter/(ntime*nncol)) + " %", end="\r")
            
            # Set up temp dataset replacing lev with z as indexed dimension
            ds_temp = ds.isel(time=itime).swap_dims({'lev':'z'}).reset_coords("lev")
            # Interpolate lev to z 
            ds_temp = ds_temp.interp(z=z_interp)
            # Append 
            ds_interp_itime_step.append(ds_temp)
    
        print("That took", "{:.4f}".format(time.time() - toc1), "seconds")
        print("Completed", "{:.3f}".format(100*counter/(ntime)) + " %")

        
        
    # Concat itime step
    ds_interpolated = xr.concat(ds_interp_itime_step, dim="time")
    print("Interpolating took", "{:.4f}".format(time.time() - toc), "seconds")
          
    return ds_interpolated

def interpolate_lev_to_z_latlon(ds): # In Progress
    
    toc = time.time()
    nlev = ds.lev.size
    ntime = ds.time.size
    
    if "lat" not in ds.dims: 
        print("lat not in dims, createing lat")
        #     print("lat not in dims")
        #     print("nlev, ntime", nlev, ntime)
        ds['lat'] = 1
        ds = ds.assign_coords({"lat":ds.lat}).expand_dims(dim={"lat":1})
        create_lat = True

    if "lon" not in ds.dims:
        print("lon not in dims, creating lon")
        ds['lon'] = 1
        ds = ds.assign_coords({"lon":ds.lon}).expand_dims(dim={"lon":1})
        create_lon = True
        
    # TODO create time?
    # TODO create lev?


    nlat = ds.lat.size
    nlon = ds.lon.size
        
    print("nlev, ntime, nlat, nlon " + str(nlev) + " " + str(ntime) + " " + str(nlat) + " " + str(nlon))

    print()

    
    z_interp = np.arange(120,50,-1)
    ds_interp_itime_step = []
    ds_interp_ilat_step = []
    ds_interp_ilon_step = []
    
    counter = 0
    for itime in range(ntime):
        toc1 = time.time()
        for ilat in range(nlat):
            for ilon in range(nlon):
                counter += 1
                
                # Set up temp dataset replacing lev with z as indexed dimension
                ds_temp = ds.isel(time=itime,lat=ilat,lon=ilon).swap_dims({'lev':'z'}).reset_coords("lev")
                # Interpolate lev to z 
                ds_temp = ds_temp.interp(z=z_interp)
                # Append 
                ds_interp_ilon_step.append(ds_temp)
            
            # Concat ilon step
            ds_interp_ilon_step = xr.concat(ds_interp_ilon_step, dim="lon")
            # Append to ilat step
            ds_interp_ilat_step.append(ds_interp_ilon_step)
            # Reset ilon step
            ds_interp_ilon_step = []
        
        # Concat ilat step
        ds_interp_ilat_step = xr.concat(ds_interp_ilat_step, dim="lat")
        # Append to itime step
        ds_interp_itime_step.append(ds_interp_ilat_step)
        # Reset ilat step
        ds_interp_ilat_step = []
    
        print("Interpolating time", itime, "took", "{:.4f}".format(time.time() - toc1), "seconds")
        print("Completed", "{:.3f}".format  (100*counter/(ntime*nlat*nlon)) + " %")
    
    # Concat itime step
    ds_interpolated = xr.concat(ds_interp_itime_step,dim="time")
         
    print("Interpolating took", "{:.4f}".format(time.time() - toc), "seconds")
    
    if create_lat == True:
        ds_interpolated = ds_interpolated.isel(lat=0).drop("lat")
    if create_lon == True:
        ds_interpolated = ds_interpolated.isel(lon=0).drop("lon")
          
    return ds_interpolated


def interpolate_z_to_lev_latlon(ds): 
    
    toc = time.time()
    nz = ds.z.size
    ntime = ds.time.size
    
    if "lat" not in ds.dims: 
        print("lat not in dims, createing lat")
        ds['lat'] = 1
        ds = ds.assign_coords({"lat":ds.lat}).expand_dims(dim={"lat":1})
        create_lat = True

    if "lon" not in ds.dims:
        print("lon not in dims, creating lon")
        ds['lon'] = 1
        ds = ds.assign_coords({"lon":ds.lon}).expand_dims(dim={"lon":1})
        create_lon = True
        
    # TODO create time?

    nlat = ds.lat.size
    nlon = ds.lon.size
        
    print("nz, ntime, nlat, nlon " + str(nz) + " " + str(ntime) + " " + str(nlat) + " " + str(nlon))
    
    lev_interp = np.logspace(0,-5,51).tolist()
    ds_interp_itime_step = []
    ds_interp_ilat_step = []
    ds_interp_ilon_step = []
    
    counter = 0
    for itime in range(ntime):
        toc1 = time.time()
        for ilat in range(nlat):
            for ilon in range(nlon):
                counter += 1
                
                # Set up temp dataset replacing lev with z as indexed dimension
                ds_temp = ds.isel(time=itime,lat=ilat,lon=ilon).swap_dims({'z':'lev'}).reset_coords("z")
                # Interpolate lev to z 
                ds_temp = ds_temp.interp(lev=lev_interp)
                # Append 
                ds_interp_ilon_step.append(ds_temp)
            
            # Concat ilon step
            ds_interp_ilon_step = xr.concat(ds_interp_ilon_step, dim="lon")
            # Append to ilat step
            ds_interp_ilat_step.append(ds_interp_ilon_step)
            # Reset ilon step
            ds_interp_ilon_step = []
        
        # Concat ilat step
        ds_interp_ilat_step = xr.concat(ds_interp_ilat_step, dim="lat")
        # Append to itime step
        ds_interp_itime_step.append(ds_interp_ilat_step)
        # Reset ilat step
        ds_interp_ilat_step = []
    
        print("Interpolating time", itime, "took", "{:.4f}".format(time.time() - toc1), "seconds")
        print("Completed", "{:.3f}".format  (100*counter/(ntime*nlat*nlon)) + " %")
    
    # Concat itime step
    ds_interpolated = xr.concat(ds_interp_itime_step,dim="time")
         
    print("Interpolating took", "{:.4f}".format(time.time() - toc), "seconds")
    
    if create_lat == True:
        ds_interpolated = ds_interpolated.isel(lat=0).drop("lat")
    if create_lon == True:
        ds_interpolated = ds_interpolated.isel(lon=0).drop("lon")
          
    return ds_interpolated