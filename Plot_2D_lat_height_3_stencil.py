# %load "../Helper_Templates/Plot_2D_lat_height_3_stencil.py"
months= ['Jan']
# variables = ['T','CO2','CO','H2O','NO','O','O3']
variables = ['O']

for month in months:
    for var in variables:
        datapath = "./Processed Data/"
        filename_rr =
        filename_nonrr = 

        ## Read in data
        ds_nonrr_var = xr.open_dataset("./Processed Data/"+ month + "_" + var + "_ZM_lon_mask_210_310_Non-RR.nc")[var]
        ds_rr_var = xr.open_dataset("./Processed Data/"+ month + "_" + var + "_ZM_lon_mask_210_310_INTERP_RR.nc")[var]

        ## Set-up arg_dict
        arg_dict = { 
            # dataset
            'ds':[ds_nonrr_var, ds_rr_var],
            'var': var,
            'diff_type': "abs",
            # 'diff_type': "%",
            'y_r_p': [1e-2,1e-5],  
            'label_size':10,
            'GPH':False,
            'exp':10**0, # Scale Plot
            'log':False,
            'ncont':10,
            'qr':None, # [start, stop, nsteps]
        }

        ## Plot
        Plot_2D_lat_height(**arg_dict)
