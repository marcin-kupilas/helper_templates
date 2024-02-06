#!/usr/bin/env python
# coding: utf-8
# %%

### Module import ###
import numpy as np
import math
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.colors as colors
import matplotlib.cm as cm
  

def Plot_1D_profiles_monthly(ds=None, var=None, gph=False, z=False, log=False, concentration=False, diff_type=None, 
                       lat_r=[-90,90], y_r=None, nxticks=None, exp=None, qr=None,diffqr=None,
                       cont=True, ncont=10, fill=True, nfill=100, diffcont=10, fontsize=10,cont_label_size=15, savefig=False, filename=None,suptitle=None,std=False,xlabel=None,zero_tick=False,legend=None):
    
    """
    
    zero_tick - controls if a zero tick is plotted, most useful for diff plots
    legend - list containing a row and column index for what subplot should contain legend
    """
   
    

        
    #========================================================================
    #===== Error check and pass input values to class-accessible values=====
    #========================================================================
    print("Start error check")
    ## Currently only allowing 3 plots

    ## ds check 
    
    ## Dimensions check
    
    ## xarray, lat, height check
    
    ## Assert no diff if size(ds)=1
    
    
    #=== End Error check and pass input values to class-accessible values===
    #========================================================================


    #=======================================================================
    #============================ Initial Setup ============================
    #=======================================================================          
 
    ## Initialise datasets for plotting format
    ds_var = []
    
    # Convert to concentration
    # Note - if want to compare two datasets, one which already is
    # in concencentration units, and one that isn't, I've had to 
    # implement a workaround to ignore files that don't need converting
    # via the "skip_convert" attribute
    # TODO REMOVE THIS - shouldn't have a function that does more than one thing!
    for i in range(len(ds)):
        if concentration==True:
            if "skip_convert" not in ds[i].attrs:
                print("Converting ", var, " to concentration for ", ds[i].attrs['ds_tag'])
                
                # Conversion constants
                k=1.380649e-23 # Boltzmann constant J K^-1= kg m^2 s^-2 K^-1
                to_cm3=1e6 # convertion from m^3 to cm^3
                to_Pa=100 # convert from hPa to Pa. Pa=kg m^-1 s^-2
                
                # Get pressure data (same for both datasets)
                P=ds[i]['lev']
                # Get temperature data
                T=ds[i]['T']
                
                # Calculate number density, ds0, ds1
                # Units: kg m^-1 s^-2 * kg^-1 m^-2 s^2 K * K^-1=m^-3 | m^-3 * 1e-6=cm^-3
                n=P*to_Pa/(to_cm3*k*T) 
                # Convert molar mixing ratio to number density
                ds_var.append(ds[i][var]*n)
                
            else:
                print("skip_convert in ds",i)
                ds_var.append(ds[i][var])
        else:
            ds_var.append(ds[i][var])   
    
#     ## Compute diffs
#     if len(ds)==2:
#         if (diff_type=='abs'): 
#             ds_diff=ds_var[1] - ds_var[0]
#         if (diff_type=='%'): 
#             ds_diff=100*(ds_var[1] - ds_var[0])/ds_var[0]
#         if (diff_type=='frac'):
#             ds_diff=ds_var[1]/ds_var[0]
        
#     # Diff operation removes gph if manually calculated.
#     # Need to add again.  
#     if diff_type!=None and gph==True and len(ds)==2:
#         print("\nWARNING: currently using arithmetic mean \n\
#         of gph between ds0['gph'] and d1['gph'] as ds_diff['gph'] \n ")
#         gph_diff = 0.5*(ds_var[0]['gph']+ds_var[1]['gph'])
#         ds_diff['gph'] = gph_diff.assign_coords(gph=gph_diff)
    
    ## Scale vars and other variable controls
             
    for i in range(len(ds)):
        ds_var[i]=ds_var[i]*exp
        if log==True:
            ds_var[i]=np.log10(ds_var[i])
            
    ## Compute var range for i months
    
    # default
    if y_r==None:
        if gph==True:
            y_r=[0,150]
        else:
            y_r=[1e3,1e-6]  
            

    if qr==None:  

        v_min = 1e50
        v_max = -1e50
        for i in range(len(ds)):
            if gph==True:
                minval = np.min(ds_var[i].where((ds_var[i]['gph']>=y_r[0]) & (ds_var[i]['gph']<=y_r[1]),drop=True)).values.item()
                v_min = min(v_min,minval)
                
                maxval = np.max(ds_var[i].where((ds_var[i]['gph']>=y_r[0]) & (ds_var[i]['gph']<=y_r[1]),drop=True)).values.item()
                v_max = max(v_max,maxval)

                print("ds_var",i, "v_min, maxval",minval, maxval)

            elif z==True:    
                minval = np.min(ds_var[i].where((ds_var[i]['z']>=y_r[0]) & (ds_var[i]['z']<=y_r[1]),drop=True)).values.item()
                v_min = min(v_min,minval)
                
                maxval = np.max(ds_var[i].where((ds_var[i]['z']>=y_r[0]) & (ds_var[i]['z']<=y_r[1]),drop=True)).values.item()
                v_max = max(v_max,maxval)

                print("ds_var",i, "minval, maxval",minval, maxval)
                if std==True:
                    ds_min = ds_var[i] - 2*ds[i][var+"_std"]
                    minval = np.min(ds_min.where((ds_var[i]['z']>=y_r[0]) & (ds_var[i]['z']<=y_r[1]),drop=True)).values.item()
                    v_min = min(v_min,minval)
                    
                    ds_max = ds_var[i] + 2*ds[i][var+"_std"]
                    maxval = np.max(ds_max.where((ds_var[i]['z']>=y_r[0]) & (ds_var[i]['z']<=y_r[1]),drop=True)).values.item()
                    v_max = max(v_max,maxval)
    
            else:
                minval = np.min(ds_var[i].where((ds_var[i]['lev']<=y_r[0]) & (ds_var[i]['lev']>=y_r[1]),drop=True)).values.item()
                v_min = min(v_min,minval)
                
                maxval = np.max(ds_var[i].where((ds_var[i]['lev']<=y_r[0]) & (ds_var[i]['lev']>=y_r[1]),drop=True)).values.item()
                v_max = max(v_max,maxval)

                print("ds_var",i, "minval, maxval",minval, maxval)
                
        print("ds_var",i, "v_min, v_max",v_min, v_max)

                   
        qr=[]
        qr.append(v_min)
        qr.append(v_max)
        qr.append(4) # number of x ticks given a limit between v_min and v_max

    else:    
        v_min=qr[0]
        v_max=qr[1]
        
    
    print("qr0 qr1 qr2", qr[0], qr[1], qr[2])

    ## Compute diff range
#     if diff_type!=None:
#         if diffqr==None:
#             if gph==True:
#                 diff_min=np.min(ds_diff.where((ds_diff['gph']>=y_r[0]) & (ds_diff['gph']<=y_r[1])))
#                 diff_max=np.max(ds_diff.where((ds_diff['gph']>=y_r[0]) & (ds_diff['gph']<=y_r[1])))
#             else:
#                 diff_min=np.min(ds_diff.where((ds_diff['lev']<=y_r[0]) & (ds_diff['lev']>=y_r[1])))
#                 diff_max=np.max(ds_diff.where((ds_diff['lev']<=y_r[0]) & (ds_diff['lev']>=y_r[1])))

#             diff_min=diff_min.values
#             diff_max=diff_max.values

#             diffqr=[]
#             diffqr.append(diff_min)
#             diffqr.append(diff_max)
#             diffqr.append(ncont)

#         else:

#             diff_min=diffqr[0]
#             diff_max=diffqr[1]
           
        
        

    
    print("Finished initial set-up\n")
    #========================== End Initial Setup ==========================
    #=======================================================================
    
    #=======================================================================
    #=============================== Plot===================================
    #=======================================================================      
    
    print("Start plotting script\n")

    # Set-up figure
    num_cols = math.ceil(len(ds_var[0].time.values) / 2)
    num_rows = 2 if len(ds_var[0].time.values) >= 2 else 1

    fig_width = num_cols * 2.5
    fig_height = num_rows * 4
        
    fig, axes = plt.subplots(nrows=num_rows,ncols=num_cols, figsize=(fig_width*1.2,fig_height*1.2))
#     fig.set_dpi(300)
    
#     if suptitle!=None:
#         fig.suptitle(suptitle)
    
    # Set up coordinates
    if gph==True:
        axis_types = {"y":"gph"}
    elif z==True:
        axis_types = {"y":"z"}
    else:
        axis_types = {"y":"lev","yscale":"log"}

    # Set up sub-plot line labels
    for i in range(len(ds)):
        ds_var[i].attrs['ds_label'] = ds[i].attrs['ds_label']

    # colors
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
]

    # debug
    print(axes.flatten())
    print("ds_var[0].time", ds_var[0].time.values)
    print("ds_var[0].time", ds_var[0].time.values)
    #

    ### Create figures and plot data
    for i in range(num_rows*num_cols):
        if i+1 > len(ds_var[0].time.values):
            break
        row = i // num_cols 
        col = i % num_cols
        print("i", i+1)
        print("len dsvar0.time.values", len(ds_var[0].time.values))
        
        ## Loop over different data ii per month i
        for ii, ds_plot in enumerate(ds_var):
            print("plotting month i var",i,ds_plot.attrs['ds_label'])
           
            ds_plot_var = ds_plot.isel(time=i)

            ds_plot_var.plot(ax=axes[row,col],ylim=[y_r[0],y_r[1]],xlim=[qr[0],qr[1]], color=colors[ii],label=ds_plot_var.attrs['ds_label'],**axis_types)

            if std==True:
                if "MIPAS" in ds[ii].attrs['ds_label']:
                    y_std = ds[ii][axis_types['y']]
                    var_std = ds[ii][var + "_std"].isel(time=i)
                    var_std_pve = ds[ii][var].isel(time=i) + 2*var_std
                    var_std_nve = ds[ii][var].isel(time=i) - 2*var_std
                    axes[row,col].fill_betweenx(y_std,var_std_pve,var_std_nve,alpha=0.3)                   
                else:
                    y_std = ds[ii][axis_types['y']].isel(time=i)
                    var_std = ds[ii][var + "_std"].isel(time=i)
                    var_std_pve = ds[ii][var].isel(time=i) + 2*var_std
                    var_std_nve = ds[ii][var].isel(time=i) - 2*var_std
                    axes[row,col].fill_betweenx(y_std,var_std_pve,var_std_nve,alpha=0.3)

            ## Subplot titles
            title_label=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
            axes[row,col].set_title(title_label[i])
            ii +=1

        ##########################################
        ## Final subplot controls
        ##########################################   
         
        if gph==True:
            ylabel="GPH (km)"
        elif z==True:
            ylabel="Geometric altitude (km)"
        else:
            ylabel="Pressure (hPa)"

        # Axis legend
        if legend==None:
            axes[row,col].legend()
        else:
            axes[legend[0],legend[1]].legend()

        ## Axis labels
        axes[row,col].set_ylabel(ylabel,fontsize=fontsize)

        # x label
        if var!='T':            
            if concentration==True:
                # Calculate scale str based on exponent
                for i in range(2,18):
                    str_exp = "{:.0e}".format(exp)
                    iexp = 10**-i
                    iexp = "{:.0e}".format(iexp)
                    if iexp==str_exp:
                        scalestr=fr" (x10$^{i})$"
                        break
                    else:
                        scalestr=""

                if log==True:
                    axes[row,col].set_xlabel("log10(" +var + r") cm$^{{-3}}$ " + scalestr,fontsize=fontsize)
                else:
                    axes[row,col].set_xlabel(var +  r" cm$^{{-3}}$" + scalestr,fontsize=fontsize)

                # if diff_type!=None: axes[row,col].set_xlabel(diff_type + " difference",fontsize=fontsize)

            else: # if concentration
                
                # Calculate scale str based on exponent
                for i in range(2,18):
                    str_exp = "{:.0e}".format(exp)
                    iexp = 10**i
                    iexp = "{:.0e}".format(iexp)
                    if iexp==str_exp:
                        scalestr=f" (x10^-{i})"
                        break
                    else:
                        scalestr=""
                        
                # if ppm or ppb
                if exp==1e6: scalestr="ppm"
                if exp==1e9: scalestr="ppb"

                if log==True:
                    axes[row,col].set_xlabel("log10(" +var + ") mmr " + scalestr,fontsize=fontsize)
                else:
                    axes[row,col].set_xlabel(var + " mmr " + scalestr,fontsize=fontsize)
                    
                if xlabel!=None:
                    axes[row,col].set_xlabel(xlabel)
                    
                # if diff_type!=None: axes[row,col].set_xlabel(diff_type + " difference",fontsize=fontsize)

        else:
            axes[row,col].set_xlabel(var ,fontsize=fontsize)
            # if diff_type!=None:  axes[row,col].set_title(diff_type + " difference",fontsize=fontsize)
            if xlabel!=None:
                axes[row,col].set_xlabel(xlabel)

        ## End Axis labels                

        # Tick controls
        axes[row,col].set_xticks(np.linspace(qr[0],qr[1],qr[2]))
        x_tick_labels = ['{:.1f}'.format(label) for label in np.linspace(qr[0], qr[1], qr[2])]
        axes[row,col].set_xticklabels(x_tick_labels)
        
        if zero_tick==True:
            axes[row,col].axvline(x=0,linestyle="--",c="k",linewidth=0.75)
            axes[row,col].text(0,y_r[0],"0")
            
        axes[row,col].tick_params(which="major",labelsize=fontsize,width=1,length=6,direction="out",right=False)
        axes[row,col].tick_params(which="minor",labelsize=fontsize,width=0,length=0,direction="out",right=False)

#     if diff_type!=None: ax[2].tick_params(which="major",labelsize=fontsize,width=1,length=6,direction="out",right=False)
    
        # Grid controls
        axes[row,col].grid(True)
 
        ## Axis suptitle todo

        ##########################################
        ## End Final subplot controls
        ##########################################   

    # Remove empty subplots        
    if len(ds_var[0].time.values) < num_rows * num_cols:
        for i in range(len(ds_var[0].time.values), num_rows * num_cols):
            fig.delaxes(axes.flatten()[i])

        
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
       
#     # savefig=False
#     if savefig:
#         plt.savefig("/glade/scratch/mmkupilas/Analysis/Plots/"+ filename + ".png")
        
#     return ds_diff
    #============================= End Plot================================
    #=======================================================================

