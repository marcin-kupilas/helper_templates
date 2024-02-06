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
                
            # Merge
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

        
        
    # Grand append    
    ds_interpolated = xr.concat(ds_interp_itime_step, dim="time")
    print("Interpolating took", "{:.4f}".format(time.time() - toc), "seconds")
          
    return ds_interpolated

def interpolate_z_to_lev_latlon(ds): # TODO
    
    return ds_interpolated