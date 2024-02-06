import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import time
import nc_time_axis

# Prepare data for vertical velocity perturbation calculations
def prepare_data_omega(ds_raw, ds_T):
   
    print("Converting omega (Pa/s) to vertical velocity (m/s)")  
    # ds_raw : xr dataset containing 'lev' and 'OMEGA'
    # ds_T : xr dataset containing "Temperature"

    # Conversion constants
    R = 287.058 # J/ kg^-1 K^-1 => m^2 s^-2 K^-1
    to_Pa = 100 # convert from hPa to Pa. Pa => kg m^-1 s^-2
    g = 9.80665 # m s^-2
    
    # Calculate mass density: kg m^-3    
    rho = ds_raw['lev']*to_Pa/(R*ds_T) # Units: kg m^-1 s^-2 * m^-2 s^2 K * K^-1 => kg m^-3
   
    # Convert velocity
    ds = -ds_raw/(rho*g) # Units: kg m^-1 s^-2 * s^-1 * kg^-1 m^3 * m^-1 s^2 => m s^-1
    ds.attrs['long_name'] = "Vertical velocity"
    ds.name = "w"
    
    # mmk debug
    print("Debug Notes")
    print("ds_T.isel(time=0,ncol=0).values\n", ds_T.isel(time=0,ncol=0).values)
    print("ds_raw.isel(time=0,ncol=0).values\n", ds_raw.isel(time=0,ncol=0).values)
    print("ds.isel(time=0,ncol=0).values\n", ds.isel(time=0,ncol=0).values)
    
    return ds
 
def prepare_data_constituent(ds_raw, ds_T, var):
    print("Converting mixing ratio to concentration (n cm^-3)")
    
    # ds_raw : xr dataset containing 'lev' and 'constituent mixing ratio'
    # ds_T : xr dataset containing "Temperature"

    ### Calculate number density: particles cm^-3 (take care with units!)
    # Conversion constants
    k = 1.380649e-23 # Boltzmann constant J K^-1 =  kg m^2 s^-2 K^-1
    to_cm3 = 1e6 # convertion from m^3 to cm^3
    to_Pa = 100 # convert from hPa to Pa. Pa = kg m^-1 s^-2
    n = ds_raw['lev']*to_Pa/(to_cm3*k*ds_T) # Units: kg m^-1 s^-2 * kg^-1 m^-2 s^2 K * K^-1 = m^-3 | m^-3 * 1e-6 = cm^-3

    ds = ds_raw*n
    ds.attrs['long_name'] = var + " number density (particles cm$^{-3}$)"
    ds.name = var + "_n"
    
    return ds

def perform_linear_regression(x,y):
    A = np.vstack([x, np.ones(len(x))]).T
    m, c = np.linalg.lstsq(A, y, rcond=None)[0]    
    return m, c
    
def calculate_perturbations_12hlinear(ds, file, var):
    """
    Using 12 h linear fit, calculates perturbations according to Appendix A: Gardner and Liu (2007)
    
    """
    debug=False # mmk debug 
    #-------------------------------------------------------------
    # Step 1 - calculate perturbations using 12 h linear fit 
    #-------------------------------------------------------------
    print("#################")
    print("Starting perturbation calculation step 1")
    tic = time.time()
    
    step1_ds = ds
    step1_result = ds*0
    
    ncol = step1_result['ncol'].size
    nlev = step1_result['lev'].size
    ntime = step1_result['time'].size
    
    # Loop over columns
    for icol in range(ncol):
        print("computing icol", icol, "\r")
        # loop over levels
        for ilev in range(nlev):
        
            # Extract timeseries data for a single column and level
            data = step1_ds.isel(ncol=icol,lev=ilev).values
            x = np.arange(step1_ds.time.dt.hour.values[0],step1_ds.time.dt.hour.values[-1]+1)

            # perform linear regression to calculate slope and intercept
            slope, intercept = perform_linear_regression(x, data)
            
            # Calculate fitted data
            fitted_data = slope * x + intercept
            
            # Calculate residuals
            residuals = data - fitted_data
            
            # Add residuals to resulting dataset
            lev_temp = step1_result.lev[ilev]
            ncol_temp = step1_result.ncol[icol]
            step1_result.loc[dict(lev=lev_temp, ncol=ncol_temp)] = residuals

            # Turn on/off at top of function
            if icol == ncol-1 and ilev == 0:
                if (debug):
                    var_debug = ds.long_name
                    fig, (ax1, ax2) = plt.subplots(2,figsize=(6,8))
                    
                    ax1.plot(x,data, marker="x", linewidth=1, label="step1_ds")
                    ax1.plot(x,fitted_data, label="step1_mean")
        
                    residuals = step1_result.isel(lev=ilev,ncol=icol)
                    
                    ax2.plot(x,residuals, marker="x",label="step1 perturbations")
                   
                    ax1.set_title("Step 1 calculating " + var_debug + " perturbations \n" \
                                 + " at icol ilev " + str(icol) + " " + str(ilev))
                    ax1.set_ylabel(var_debug)
                    ax1.set_xlabel("UT hour")
                    ax1.legend()
                    
                    ax2.set_title(" ")
                    ax2.set_ylabel(ds.name + "'")
                    ax2.legend()
                    plt.show()
                    
                
    toc = time.time()
    print("Step 1 took", "{:.4f}".format(toc - tic), "seconds\n")
    
    #-------------------------------------------------------------
    # Step 2 - subtract data exceeding 3*std away from average 
    #-------------------------------------------------------------
#    print("#################")
#    print("Starting perturbation calculation step 2")
#    toc = time.time()
#    
#    # Loop over columns and levels
#    # Calculate standard deviation over time window
#    std = step1_result.std(axis=0)
#    
#    # Loop over columns and levels
#    # remove data from calculation window exceeding 3*std of average 
#
#    # Loop over columns and levels
#    for icol in range(ncol):
#        tic = time.time()
#        for ilev in range(nlev):
#            toc_1 = time.time()
#            icol_ilev_mean = np.mean(step1_result.isel(ncol=icol, lev=ilev))
#
#            # Identify outliers - returns list of boolean, True = outlier, False = not outlier
#            delta = np.abs(step1_result.isel(ncol=icol, lev=ilev) - icol_ilev_mean)
#            outliers = delta > 3 * std.isel(lev=ilev, ncol=icol)
#            # Convert outlier array to positional index array (e.g. [True,False,True] returns [0,2]
#            condition_indices = np.where(outliers==True)[0]
#            #print(condition_indices)
#
#            # Get array data at location and add nan where they should be
#            array_data = step1_result.isel(lev=ilev,ncol=icol).values
#            array_data[condition_indices] = np.nan
#            # Set level data with updated values
#            step1_result.loc[dict(lev=step1_result.lev[ilev],\
#                                  ncol=step1_result.ncol[icol])] = array_data
#           #
#           #if (debug):
#           #    print("col", icol+1, "/", ncol,\
#           #          "lev", ilev+1, "/", nlev,\
#           #          "took", "{:.4f}".format(time.time() - toc_1), "seconds")
#            print("Computing ncol", ncol, "ilev", ilev, "took", "{:.4f}".format(time.time() - toc_1), "seconds\n")
#
#
    step2_result = step1_result
#    
#    # mmk: You now have the dataset step2_result with nan values where
#    # data in time window at a single icol and ilev 
#    # exceeds 3*std of step1_result average
#    
#    tic = time.time() 
#    print("removed", np.count_nonzero(np.isnan(step2_result.values)), "values to nan")
#    print("Step 2 took", "{:.4f}".format(tic - toc), "seconds\n")
#
    #-------------------------------------------------------------
    # Step 3 - calculate perturbations from Step 2
    #-------------------------------------------------------------
    print("#################")
    print("Starting perturbation calculation step 3")
    tic = time.time()
    
    step3_ds = step2_result
    step3_result = step2_result*0

    # Loop over columns
    for icol in range(ncol):
        print("Computing icol", icol)
        # loop over levels
        for ilev in range(nlev):
        
            # Extract timeseries data for a single column and level
            data = step3_ds.isel(ncol=icol,lev=ilev).values
            x = np.arange(step3_ds.time.dt.hour.values[0],step3_ds.time.dt.hour.values[-1]+1)

            # perform linear regression to calculate slope and intercept
            slope, intercept = perform_linear_regression(x, data)
            
            # Calculate fitted data
            fitted_data = slope * x + intercept
            
            # Calculate residuals
            residuals = data - fitted_data
            
            # Add residuals to resulting dataset
            lev_temp = step3_result.lev[ilev]
            ncol_temp = step3_result.ncol[icol]
            step3_result.loc[dict(lev=lev_temp, ncol=ncol_temp)] = residuals

            # Turn on/off at top of function
            if icol == ncol-1 and ilev == 0:
                if (debug):
                    var_debug = ds.long_name
                    fig, (ax1, ax2) = plt.subplots(2,figsize=(6,8))
                    
                    ax1.plot(x,data, marker="x", linewidth=1, label="step3_ds")
                    ax1.plot(x,fitted_data, label="step3_mean")
        
                    residuals = step3_result.isel(lev=ilev,ncol=icol)
                    
                    ax2.plot(x,residuals, marker="x",label="step3 perturbations")
                   
                    ax1.set_title("Step 3 calculating " + var_debug + " perturbations \n" \
                                 + " at icol ilev " + str(icol) + " " + str(ilev))
                    ax1.set_ylabel(var_debug)
                    ax1.set_xlabel("UT hour")
                    ax1.legend()
                    
                    ax2.set_title(" ")
                    ax2.set_ylabel(ds.name + "'")
                    ax2.legend()
                    plt.show()
                    
                
    toc = time.time()
    print("Step 3 took", "{:.4f}".format(toc - tic), "seconds\n")
 
    #-------------------------------------------------------------
    # Step 4 - subtract data exceeding 3*std from step 3
#    #-------------------------------------------------------------
#    print("#################")
#    print("Starting perturbation calculation step 4")
#    toc = time.time()
#
#    # Loop over columns and levels
#    # Calculate standard deviation over time window
#    std = step3_result.std(axis=0)
#    
#    # Loop over columns and levels
#    # remove data from calculation window exceeding 3*std of average 
#
#    # Loop over columns and levels
#    for icol in range(ncol):
#        print("computing icol", icol)
#        for ilev in range(nlev):
#            toc_1 = time.time()
#            icol_ilev_mean = np.mean(step3_result.isel(ncol=icol, lev=ilev))
#            
#            # Identify outliers - returns list of boolean, True = outlier, False = not outlier
#            delta = np.abs(step3_result.isel(ncol=icol, lev=ilev) - icol_ilev_mean)
#            outliers = delta > 3 * std.isel(lev=ilev, ncol=icol)
#            #print(outliers)
#            # Convert outlier array to positional index array (e.g. [True,False,True] returns [0,2]
#            condition_indices = np.where(outliers==True)[0]
#
#            # Get array data at location and add nan where they should be
#            array_data = step3_result.isel(lev=ilev,ncol=icol).values
#            array_data[condition_indices] = np.nan
#
#            # Set level data with updated values
#            step3_result.loc[dict(lev=step3_result.lev[ilev],\
#                                  ncol=step3_result.ncol[icol])] = array_data
#           #
#           #if (debug):
#           #    print("col", icol+1, "/", ncol,\
#           #          "lev", ilev+1, "/", nlev,\
#           #          "took", "{:.4f}".format(time.time() - toc_1), "seconds")
#
    step4_result = step3_result
#    
#    # mmk: You now have the dataset step2_result with nan values where
#    # data in time window at a single icol and ilev 
#    # exceeds 3*std of step3_result average
#    
#    tic = time.time() 
#    print("removed", np.count_nonzero(np.isnan(step4_result.values)), "values to nan")
#
#    print("Step 4 took", "{:.4f}".format(tic - toc), "seconds\n")
#
    #----------------------------------------------------------
    # Step 5 - Subtract vertical mean from perturbation profile.
    #----------------------------------------------------------
    toc = time.time()
    print("#################")
    print("Starting perturbation calculation step 5 - calculating vertical mean")
    
    val_mean_tot = 0 # total of weighted values in pressure slice
    dlev_tot = 0 # total pressure thickness

    for ilev in range(nlev-1): 
        # mmk: nlev-1 because final loop iteration needs access to nlev-1 and nlev.
        # If loop was in range(nlev), the final calculation would try to access
        # nlev and nlev+1

        # Isolate pressure and variable values
        lev2 = step4_result.isel(lev=ilev)['lev'] # Upper level, smaller value than lev1
        lev1 = step4_result.isel(lev=ilev+1)['lev'] # Lower level, larger value than lev2
        val2 = step4_result.isel(lev=ilev)  # Value at upper level
        val1 = step4_result.isel(lev=ilev+1) # Value at lower level

        print("Calculating average for layer between lev2, lev1", lev2.values, lev1.values)

        # Calculate layer thickness between grid midpoints
        dlev = lev1 - lev2 
        # Calculate average value at center of layer with thickness dlev
        val_mean_dlev = dlev*0.5*(val1 + val2)  

        # Add to totals
        val_mean_tot += val_mean_dlev
        dlev_tot += dlev

    val_mean = val_mean_tot/dlev_tot
    varname = var + "_prime"
    step5_result = xr.Dataset()
    step5_result[varname] = step4_result - val_mean
 
    return step5_result
    
    
def calculate_perturbations_12hlinear_1_only(ds, file, var):
    """
    Using 12 h linear fit, calculates perturbations. Only 1 cleaning step.
    
    """
    debug=False # mmk debug 
    #-------------------------------------------------------------
    # Step 1 - calculate perturbations using 12 h linear fit 
    #-------------------------------------------------------------
    print("#################")
    print("Starting perturbation calculation step 1")
    tic = time.time()
    
    step1_ds = ds
    step1_result = ds*0
    
    ncol = step1_result['ncol'].size
    nlev = step1_result['lev'].size
    ntime = step1_result['time'].size
    
    # Loop over columns
    for icol in range(ncol):
        #print("computing icol", icol, "\r")
        # loop over levels
        for ilev in range(nlev):
        
            # Extract timeseries data for a single column and level
            data = step1_ds.isel(ncol=icol,lev=ilev).values
            x = np.arange(step1_ds.time.dt.hour.values[0],step1_ds.time.dt.hour.values[-1]+1)

            # perform linear regression to calculate slope and intercept
            slope, intercept = perform_linear_regression(x, data)
            
            # Calculate fitted data
            fitted_data = slope * x + intercept
            
            # Calculate residuals
            residuals = data - fitted_data
            
            # Add residuals to resulting dataset
            lev_temp = step1_result.lev[ilev]
            ncol_temp = step1_result.ncol[icol]
            step1_result.loc[dict(lev=lev_temp, ncol=ncol_temp)] = residuals

            # Turn this part of code on/off at top of function
            if icol == ncol-1 and ilev == 0:
                if (debug):
                    var_debug = ds.long_name
                    fig, (ax1, ax2) = plt.subplots(2,figsize=(6,8))
                    
                    ax1.plot(x,data, marker="x", linewidth=1, label="step1_ds")
                    ax1.plot(x,fitted_data, label="step1_mean")
        
                    residuals = step1_result.isel(lev=ilev,ncol=icol)
                    
                    ax2.plot(x,residuals, marker="x",label="step1 perturbations")
                   
                    ax1.set_title("Step 1 calculating " + var_debug + " perturbations \n" \
                                 + " at icol ilev " + str(icol) + " " + str(ilev))
                    ax1.set_ylabel(var_debug)
                    ax1.set_xlabel("UT hour")
                    ax1.legend()
                    
                    ax2.set_title(" ")
                    ax2.set_ylabel(ds.name + "'")
                    ax2.legend()
                    plt.show()
                    
                
    toc = time.time()
    print("Step 1 took", "{:.4f}".format(toc - tic), "seconds\n")

    varname = var + "_prime"
    step5_result = xr.Dataset()
    step5_result[varname] = step1_result
 
    return step5_result

   
    
def calculate_perturbations_12hlinear_1_5_only(ds, file, var):
    """
    Using 12 h linear fit, calculates perturbations. Only 1 cleaning step.
    
    """
    debug=False # mmk debug 
    #-------------------------------------------------------------
    # Step 1 - calculate perturbations using 12 h linear fit 
    #-------------------------------------------------------------
    print("#################")
    print("Starting perturbation calculation step 1")
    tic = time.time()
    
    step1_ds = ds
    step1_result = ds*0
    
    ncol = step1_result['ncol'].size
    nlev = step1_result['lev'].size
    ntime = step1_result['time'].size
    
    # Loop over columns
    for icol in range(ncol):
        #print("computing icol", icol, "\r")
        # loop over levels
        for ilev in range(nlev):
        
            # Extract timeseries data for a single column and level
            data = step1_ds.isel(ncol=icol,lev=ilev).values
            x = np.arange(step1_ds.time.dt.hour.values[0],step1_ds.time.dt.hour.values[-1]+1)

            # perform linear regression to calculate slope and intercept
            slope, intercept = perform_linear_regression(x, data)
            
            # Calculate fitted data
            fitted_data = slope * x + intercept
            
            # Calculate residuals
            residuals = data - fitted_data
            
            # Add residuals to resulting dataset
            lev_temp = step1_result.lev[ilev]
            ncol_temp = step1_result.ncol[icol]
            step1_result.loc[dict(lev=lev_temp, ncol=ncol_temp)] = residuals

            # Turn on/off at top of function
            if icol == ncol-1 and ilev == 0:
                if (debug):
                    var_debug = ds.long_name
                    fig, (ax1, ax2) = plt.subplots(2,figsize=(6,8))
                    
                    ax1.plot(x,data, marker="x", linewidth=1, label="step1_ds")
                    ax1.plot(x,fitted_data, label="step1_mean")
        
                    residuals = step1_result.isel(lev=ilev,ncol=icol)
                    
                    ax2.plot(x,residuals, marker="x",label="step1 perturbations")
                   
                    ax1.set_title("Step 1 calculating " + var_debug + " perturbations \n" \
                                 + " at icol ilev " + str(icol) + " " + str(ilev))
                    ax1.set_ylabel(var_debug)
                    ax1.set_xlabel("UT hour")
                    ax1.legend()
                    
                    ax2.set_title(" ")
                    ax2.set_ylabel(ds.name + "'")
                    ax2.legend()
                    plt.show()
                    
                
    toc = time.time()
    print("Step 1 took", "{:.4f}".format(toc - tic), "seconds\n")

    #----------------------------------------------------------
    # Step 5 - Subtract vertical mean from perturbation profile.
    #----------------------------------------------------------
    toc = time.time()
    print("#################")
    print("Starting perturbation calculation step 5 - calculating vertical mean")
    
    val_mean_tot = 0 # total of weighted values in pressure slice
    dlev_tot = 0 # total pressure thickness

    for ilev in range(nlev-1): 
        # mmk: nlev-1 because final loop iteration needs access to nlev-1 and nlev.
        # If loop was in range(nlev), the final calculation would try to access
        # nlev and nlev+1

        # Isolate pressure and variable values
        lev2 = step1_result.isel(lev=ilev)['lev'] # Upper level, smaller value than lev1
        lev1 = step1_result.isel(lev=ilev+1)['lev'] # Lower level, larger value than lev2
        val2 = step1_result.isel(lev=ilev)  # Value at upper level
        val1 = step1_result.isel(lev=ilev+1) # Value at lower level

        print("Calculating average for layer between lev2, lev1", lev2.values, lev1.values)

        # Calculate layer thickness between grid midpoints
        dlev = lev1 - lev2 
        # Calculate average value at center of layer with thickness dlev
        val_mean_dlev = dlev*0.5*(val1 + val2)  

        # Add to totals
        val_mean_tot += val_mean_dlev
        dlev_tot += dlev

    val_mean = val_mean_tot/dlev_tot
    varname = var + "_prime"
    step5_result = xr.Dataset()
    step5_result[varname] = step1_result - val_mean
 
    return step5_result
 


def calculate_perturbations_rolling(ds, file, var):
    """
    IN:
    ds - dataset for perturbation calculations
        
    Following 5 steps from Appendix A: Gardner and Liu (2007) 
    Seasonal variations of the vertical fluxes of heat and horizontal momentum 
    in the mesopause region at Starfire Optical Range, New Mexico
    
    file - filename for debug purposes only for png filename
    
    var - variable name - for debug purposes only for png filename
    
    OUT:
    
    ds_prime - perturbation dataset.
    
    NOTES:
    
    Designed to work with 15 hourly readings, reducing it a dataset with 7 readings.
    
    This script has two main stages - Calculating perturbations, and debug.

    """
    print("Calculating perturbations")
    debug=False # mmk debug
    
    # mmk dev - returned by function for inspection by notebook
    step1_result = xr.Dataset()
    step2_result = xr.Dataset()
    step3_result = xr.Dataset()
    step4_result = xr.Dataset()
    ds_prime = xr.Dataset()

    
    #----------------------------------------------------------
    # Step 1 - Calculate perturbations using 5 day rolling mean
    #----------------------------------------------------------
    
    toc = time.time()
    
    step1_ds = ds
    step1_rolling = step1_ds.rolling(time=5, center=True).mean() # First and last two values become nans
    step1_result = step1_ds - step1_rolling 
    
    step1_rolling_3 = step1_ds.rolling(time=3, center=True).mean() # First and last two values become nans
    step1_rolling_7 = step1_ds.rolling(time=7, center=True).mean() # First and last two values become nans

    step1_result_3 = step1_ds - step1_rolling_3
    step1_result_7 = step1_ds - step1_rolling_7

    tic = time.time()
    print("#################")
    print("Step 1 took", "{:.4f}".format(tic - toc), "seconds\n")
        
    #----------------------------------------------------------
    # Step 2 - Remove data from calculation window that exceeds 3*std of average in window
    #----------------------------------------------------------
   
    toc = time.time()

    # n columns, n levels, n time stamps
    ncol = step1_result['ncol'].size
    nlev = step1_result['lev'].size
    ntime = step1_result['time'].size
    
    # Loop over columns and levels
    # Calculate standard deviation over time window
    std = step1_result.std(axis=0)
    
    # Loop over columns and levels
    # remove data from calculation window exceeding 3*std of average 

    # Loop over columns and levels
    for icol in range(ncol):
        for ilev in range(nlev):
            toc_1 = time.time()
            icol_ilev_mean = np.mean(step1_result.isel(ncol=icol, lev=ilev))

            # Identify outliers - returns list of boolean, True = outlier, False = not outlier
            delta = np.abs(step1_result.isel(ncol=icol, lev=ilev) - icol_ilev_mean)
            outliers = delta > 3 * std.isel(lev=ilev, ncol=icol)

            # Convert outlier array to positional index array (e.g. [True,False,True] returns [0,2]
            condition_indices = np.where(outliers==True)[0]

            # Get array data at location and add nan where they should be
            array_data = step1_result.isel(lev=ilev,ncol=icol).values
            array_data[condition_indices] = np.nan

            # Set level data with updated values
            step1_result.loc[dict(lev=step1_result.lev[ilev],\
                                  ncol=step1_result.ncol[icol])] = array_data
            
            if (debug):
                print("col", icol+1, "/", ncol,\
                      "lev", ilev+1, "/", nlev,\
                      "took", "{:.4f}".format(time.time() - toc_1), "seconds")

    step2_result = step1_result
    
    # mmk: You now have the dataset step2_result with nan values where
    # data in time window at a single icol and ilev 
    # exceeds 3*std of step1_result average
    
    tic = time.time() 
    print("#################")
    print("Step 2 took", "{:.4f}".format(tic - toc), "seconds\n")
    
    #----------------------------------------------------------
    # Step 3 - Calculate perturbations using 5 day rolling mean of data from Step 2
    #----------------------------------------------------------
    
    toc = time.time()
    
    step3_ds = step1_result.copy()
    step3_rolling = step3_ds.rolling(time=5, center=True).mean() # First and last two values become nans
    step3_result = step3_ds - step3_rolling 
    
    tic = time.time()
    print("#################")
    print("Step 2 took", "{:.4f}".format(tic - toc), "seconds\n")
    
    #----------------------------------------------------------
    # Step 4 - Remove data from calculation windows that exceeds 3*std of average in window
    #----------------------------------------------------------
  
    toc = time.time()

    # n columns, n levels, n time stamps
    ncol = step3_result['ncol'].size
    nlev = step3_result['lev'].size
    ntime = step3_result['time'].size
    
    # Loop over columns and levels
    # Calculate standard deviation over time window
    std = step3_result.std(axis=0)
    
    # Loop over columns and levels
    # remove data from calculation window exceeding 3*std of average 

    # Loop over columns and levels
    for icol in range(ncol):
        for ilev in range(nlev):
            toc_1 = time.time()
            icol_ilev_mean = np.mean(step3_result.isel(ncol=icol, lev=ilev))

            # mmk debug
            if icol == 0 and ilev == 10:
                print("WARNING! DEBUG TEST VALUES!")
                step3_result.loc[dict(lev=step3_result.lev[10],\
                                      ncol=step3_result.ncol[0],\
                                      time=step3_result.time[8])] = 999
            if icol == 0 and ilev == 109:      
                print("WARNING! DEBUG TEST VALUES!")
                step3_result.loc[dict(lev=step3_result.lev[109],\
                                      ncol=step3_result.ncol[0],\
                                      time=step3_result.time[5])] = 999

            # Identify outliers - returns list of boolean, True = outlier, False = not outlier
            delta = np.abs(step3_result.isel(ncol=icol, lev=ilev) - icol_ilev_mean)
            outliers = delta > 3 * std.isel(lev=ilev, ncol=icol)

            # Convert outlier array to positional index array (e.g. [True,False,True] returns [0,2]
            condition_indices = np.where(outliers==True)[0]

            # Get array data at location and add nan where they should be
            array_data = step3_result.isel(lev=ilev,ncol=icol).values
            array_data[condition_indices] = np.nan

            # Set level data with updated values
            step3_result.loc[dict(lev=step3_result.lev[ilev],\
                                  ncol=step3_result.ncol[icol])] = array_data
            
            if (debug):
                print("col", icol+1, "/", ncol,\
                      "lev", ilev+1, "/", nlev,\
                      "took", "{:.4f}".format(time.time() - toc_1), "seconds")

    step4_result = step3_result
    # mmk: You now have the dataset step4_result with nan values where
    # data in time window at a single icol and ilev 
    # exceeds 3*std of step3_result average
    
    tic = time.time() 
    print("#################")
    print("Step 4 took", "{:.4f}".format(tic - toc), "seconds\n")

    #------------------------------------------------
    # mmk debug
    #------------------------------------------------
    if (debug):
        print("Debug true for calculate_perturbations")
        print("# step1_ds", step1_ds)
        print("# step1_rolling", step1_rolling)
        print("# step1_result", step1_result)
        
        fig, (ax1,ax2,ax3) = plt.subplots(3,figsize=(10,15))
        icol = 0
        ilev = 5
        
        step1_ds.isel(ncol=icol,lev=ilev).plot(ax=ax1, marker="x", label="step1_ds")
        step1_rolling_3.isel(ncol=icol,lev=ilev).plot(ax=ax1, marker="x", label="3 hour rolling mean")
        step1_rolling.isel(ncol=icol,lev=ilev).plot(ax=ax1, marker="x", label="5 hour rolling mean")
        step1_rolling_7.isel(ncol=icol,lev=ilev).plot(ax=ax1, marker="x", label="7 hour rolling mean")

        step1_result_3.isel(ncol=icol,lev=ilev).plot(ax=ax2, marker="x", label="perturbations from 3 hour mean")
        step3_ds.isel(ncol=icol,lev=ilev).plot(ax=ax2, marker="x", label="perturbations from 5 hour mean")
        step1_result_7.isel(ncol=icol,lev=ilev).plot(ax=ax2, marker="x", label="perturbations from 7 hour mean")
       # step3_rolling.isel(ncol=icol,lev=ilev).plot(ax=ax2, marker="x", label="step3_rolling")
        
        step4_result.isel(ncol=icol,lev=ilev).plot(ax=ax3, marker="x", label="step4_result")
        #step5_result.isel(ncol=icol,lev=ilev).plot(ax=ax3, marker="x", label="step5_result")
        
        ax1.legend()
        ax2.legend()
        ax3.legend()
        
        filename = var + "_" + file + ".png"
        plt.show()
        #plt.savefig("./debug/" + filename)
    else:
        print("Debug false for calculate_perturbations")
        
    return step5_result