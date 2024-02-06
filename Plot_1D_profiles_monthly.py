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
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.ticker as ticker
  

def Plot_1D_profiles_monthly(ds=None, variables=None, gph=False, z=False, log=False, concentration=False, diff_type=None, 
                       lat_r=[-90,90], y_r=None, nxticks=None, exp=None, qr=None,diffqr=None,
                       cont=True, ncont=10, fill=True, nfill=100, diffcont=10, fontsize=10,cont_label_size=15, savefig=False, filename=None,suptitle=None,std=False,xlabel=None,ylabel=None,zero_tick=False,legend=None, xsize=3,ysize=4):
    
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
    ds_var = ds.copy()
    
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
            
    ## Calculate var range 
    if qr==None:  

        v_min = 1e50
        v_max = -1e50

        ## Loop over datasets
        for i in range(len(ds_var)):
            ## Loop over variables in dataset i
            ## Note - algorithm considers all time stamps for var in dataset i
            for var in variables:
                if var in ds_var[i].data_vars:
                    if gph==True:
                        minval = np.min(ds_var[i][var].where((ds_var[i]['gph']>=y_r[0]) & (ds_var[i]['gph']<=y_r[1]),drop=True)).values.item()
                        v_min = min(v_min,minval)
                        
                        maxval = np.max(ds_var[i][var].where((ds_var[i]['gph']>=y_r[0]) & (ds_var[i]['gph']<=y_r[1]),drop=True)).values.item()
                        v_max = max(v_max,maxval)

                        print("ds_var",i, "v_min, maxval",minval, maxval)

                    elif z==True:    
                        minval = np.min(ds_var[i][var].where((ds_var[i]['z']>=y_r[0]) & (ds_var[i]['z']<=y_r[1]),drop=True)).values.item()
                        v_min = min(v_min,minval)
                        
                        maxval = np.max(ds_var[i][var].where((ds_var[i]['z']>=y_r[0]) & (ds_var[i]['z']<=y_r[1]),drop=True)).values.item()
                        v_max = max(v_max,maxval)

                        print("ds_var",i, "minval, maxval",minval, maxval)
                        if std==True:
                            if var + "_std" in ds_var[i].data_vars:
                                ds_min = ds_var[i][var] - 2*ds[i][var+"_std"]
                                minval = np.min(ds_min.where((ds_var[i]['z']>=y_r[0]) & (ds_var[i]['z']<=y_r[1]),drop=True)).values.item()
                                v_min = min(v_min,minval)
                                
                                ds_max = ds_var[i][var] + 2*ds[i][var+"_std"]
                                maxval = np.max(ds_max.where((ds_var[i]['z']>=y_r[0]) & (ds_var[i]['z']<=y_r[1]),drop=True)).values.item()
                                v_max = max(v_max,maxval)
                
                    else:
                        minval = np.min(ds_var[i][var].where((ds_var[i]['lev']<=y_r[0]) & (ds_var[i]['lev']>=y_r[1]),drop=True)).values.item()
                        v_min = min(v_min,minval)
                        
                        maxval = np.max(ds_var[i][var].where((ds_var[i]['lev']<=y_r[0]) & (ds_var[i]['lev']>=y_r[1]),drop=True)).values.item()
                        v_max = max(v_max,maxval)

                        print("ds_var",i, "minval, maxval",minval, maxval)
                    
            print("ds_var",i, "v_min, v_max",v_min, v_max)

                    
            qr=[]
            qr.append(v_min)
            qr.append(v_max)
            qr.append(4) ## number of x ticks given a limit between v_min and 
            qr.append(1) ## Number of minor ticks
            qr.append('{:.0f}') ## Number of decimal points
            v_max

    else:    
        v_min=qr[0]
        v_max=qr[1]
        qr[4] = "{:." + str(qr[4]) + "f}"
        
    
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
    if "time" not in ds_var[0].dims:
        time_size = 1
        num_cols = 1
        num_rows = 1
    else:    
        time_size = len(ds_var[0].time.values)
        # Set-up figure
        num_cols = math.ceil(time_size / 2)
        num_rows = 2 if time_size >= 2 else 1

    fig_width = num_cols * 3
    fig_height = num_rows * 4

    if "time" not in ds_var[0].dims:
        fig, axes = plt.subplots(nrows=num_rows,ncols=num_cols, figsize=(xsize,ysize),dpi=150)
    else:
        fig, axes = plt.subplots(nrows=num_rows,ncols=num_cols, figsize=(fig_width,fig_height),dpi=300)
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


    # debug
    # print(axes.flatten())
    # print("ds_var[0].time", ds_var[0].time.values)
    # print("ds_var[0].time", ds_var[0].time.values)
    #

    ### Create figures and plot data

    for i in range(num_rows*num_cols):
        if i+1 > time_size:
            break
        if "time" not in ds_var[0].dims:
            row = 0
            col = 0
        else:  
            row = i // num_cols 
            col = i % num_cols
            print("i", i+1)
            print("len dsvar0.time.values", time_size)
        
        ## There are len(ds) datasets 
        ## Each dataset has len(variables) variables
        ## Each dataset has len(ds[ii].time.values) time steps
        ## ii is just index of dataset from ds
        ## Loop over different datasets ii per month i
        if "time" not in ds_var[0].dims:
            for ii, dataset in enumerate(ds_var):
                print("ds",ii)
                jj = 0
                for var in variables:
                    if var in dataset.data_vars:
                        # print(var, "in", dataset.attrs['ds_label'])
    
                        ## Line style
                        line_style = ['-', '--', ':', '-.',(0, (5, 10)),(0, (1, 10)),(0, (3, 10, 1, 10))]
    
                        ## colors
                        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
                        '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf' ]
    
                        ## Plot all variables contained in each dataset ii
                        print("plotting ", dataset.attrs['ds_label'], "var", var)
                        
                        ## Plot obs as a bolder and more visible line
                        print(dataset.attrs['ds_props'])
                        if "OBS" in dataset.attrs['ds_props']:
                            dataset[var].plot(ax=axes,ylim=[y_r[0],y_r[1]],xlim=[qr[0],qr[1]], \
                                              label=dataset.attrs['ds_label'], linewidth="1.5", marker="o", \
                                              markersize=4, color="k",linestyle=line_style[jj], **axis_types)
                        else:
                            dataset[var].plot(ax=axes,ylim=[y_r[0],y_r[1]],xlim=[qr[0],qr[1]], \
                                              color=colors[ii],label=dataset.attrs['ds_label'],linestyle=line_style[jj],**axis_types)
    
                        if std==True:
                            if var + "_std" in dataset.data_vars:
                                if "MIPAS" in dataset.attrs['ds_label']:
                                    y_std = dataset[axis_types['y']]
                                    var_std = dataset[var + "_std"]
                                    var_std_pve = dataset[var] + 2*var_std
                                    var_std_nve = dataset[var] - 2*var_std
                                    axes.fill_betweenx(y_std,var_std_pve,var_std_nve,color="k",alpha=0.075)   
                                    print("#### ok for dataset", dataset.attrs['ds_label'], "var", var + "_std")
                        
                                else:
                                    y_std = dataset[axis_types['y']]
                                    var_std = dataset[var + "_std"]
                                    var_std_pve = dataset[var] + 2*var_std
                                    var_std_nve = dataset[var] - 2*var_std
                                    axes.fill_betweenx(y_std,var_std_pve,var_std_nve,color=colors[ii],alpha=0.15)
                                    print("#### ok for dataset", dataset.attrs['ds_label'], "var", var + "_std")
    
                        jj += 1
    
    
                ## Subplot titles
                title_label=[suptitle]
                axes.set_title(title_label[i])
        else:
            for ii, dataset in enumerate(ds_var):
                print("ds",ii)
                jj = 0
                for var in variables:
                    if var in dataset.data_vars:
                        # print(var, "in", dataset.attrs['ds_label'])
    
                        ## Line style
                        line_style = ['-', '--', ':', '-.',(0, (5, 10)),(0, (1, 10)),(0, (3, 10, 1, 10))]
    
                        ## colors
                        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
                        '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf' ]
    
                        ## Plot all variables contained in each dataset ii
                        print("plotting month i", i+1, "for dataset", dataset.attrs['ds_label'], "var", var)
                        
                        ## Plot obs as a bolder and more visible line
                        print(dataset.attrs['ds_props'])
                        if "OBS" in dataset.attrs['ds_props']:
                            dataset[var].isel(time=i).plot(ax=axes[row,col],ylim=[y_r[0],y_r[1]],xlim=[qr[0],qr[1]],label=dataset.attrs['ds_label'],linewidth="1.5", marker="o",markersize=4, color="k",linestyle=line_style[jj],**axis_types)
                        else:
                            dataset[var].isel(time=i).plot(ax=axes[row,col],ylim=[y_r[0],y_r[1]],xlim=[qr[0],qr[1]], color=colors[ii],label=dataset.attrs['ds_label'],linestyle=line_style[jj],**axis_types)
    
                        if std==True:
                            if var + "_std" in dataset.data_vars:
                                if "MIPAS" in dataset.attrs['ds_label']:
                                    y_std = dataset[axis_types['y']]
                                    var_std = dataset[var + "_std"].isel(time=i)
                                    var_std_pve = dataset[var].isel(time=i) + 2*var_std
                                    var_std_nve = dataset[var].isel(time=i) - 2*var_std
                                    axes[row,col].fill_betweenx(y_std,var_std_pve,var_std_nve,color="k",alpha=0.075)   
                                    print("#### ok for dataset", dataset.attrs['ds_label'], "var", var + "_std")
                        
                                else:
                                    y_std = dataset[axis_types['y']]
                                    var_std = dataset[var + "_std"].isel(time=i)
                                    var_std_pve = dataset[var].isel(time=i) + 2*var_std
                                    var_std_nve = dataset[var].isel(time=i) - 2*var_std
                                    axes[row,col].fill_betweenx(y_std,var_std_pve,var_std_nve,color=colors[ii],alpha=0.15)
                                    print("#### ok for dataset", dataset.attrs['ds_label'], "var", var + "_std")
    
                        jj += 1
    
    
                ## Subplot titles
                title_label=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
                axes[row,col].set_title(title_label[i])

        ##########################################
        ## Final subplot controls
        ##########################################   
         
        # plt.tight_layout()
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
                    if "time" not in ds_var[0].dims:
                        axes.set_xlabel("log10(" +var + r") cm$^{{-3}}$ " + scalestr,fontsize=fontsize)
                    else:
                        axes[row,col].set_xlabel("log10(" +var + r") cm$^{{-3}}$ " + scalestr,fontsize=fontsize)
                else:
                    if "time" not in ds_var[0].dims:
                        axes.set_xlabel(var +  r" cm$^{{-3}}$" + scalestr,fontsize=fontsize)
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
                    if "time" not in ds_var[0].dims:
                        axes.set_xlabel("log10(" +var + ") mmr " + scalestr,fontsize=fontsize)
                    else:
                        axes[row,col].set_xlabel("log10(" +var + ") mmr " + scalestr,fontsize=fontsize)
                else:
                    if "time" not in ds_var[0].dims:
                        axes.set_xlabel(var +  " mmr " + scalestr,fontsize=fontsize)
                    else:
                        axes[row,col].set_xlabel(var +  " mmr " + scalestr,fontsize=fontsize)

                if 'flux' in list(ds_var[0].keys()):
                    axes.set_xlabel("flux")

                    
                if xlabel!=None:
                    if "time" not in ds_var[0].dims:
                        axes.set_xlabel(xlabel)
                    else:
                        axes[row,col].set_xlabel(xlabel)
                    
                # if diff_type!=None: axes[row,col].set_xlabel(diff_type + " difference",fontsize=fontsize)

        else:
            if "time" not in ds_var[0].dims:
                axes.set_xlabel(var ,fontsize=fontsize)
            else:
                axes[row,col].set_xlabel(var ,fontsize=fontsize)
            # if diff_type!=None:  axes[row,col].set_title(diff_type + " difference",fontsize=fontsize)
            if xlabel!=None:
                if "time" not in ds_var[0].dims:
                    axes.set_xlabel(xlabel)
                else:
                    axes[row,col].set_xlabel(xlabel)
        ###############
        # y label
        ###############

        if ylabel==None:
            if gph==True:
                ylabel="GPH (km)"
            elif z==True:
                ylabel="Geometric altitude (km)"
            else:
                ylabel="Pressure (hPa)"
        else:
            ylabel=ylabel

        # Axis legend
        if legend==None:
            if "time" not in ds_var[0].dims:
                axes.legend(loc='lower right', bbox_to_anchor=(1.25,0),framealpha=1)
            else:
                axes[row,col].legend(loc='lower right', bbox_to_anchor=(1.25,0),framealpha=1)
        else:
            if "time" not in ds_var[0].dims:
                axes.legend(framealpha=legend[2])
            else:
                axes[legend[0],legend[1]].legend(framealpha=legend[2])

        ## Axis labels
        if "time" not in ds_var[0].dims:
            axes.set_ylabel(ylabel,fontsize=fontsize)

        else:
            axes[row,col].set_ylabel(ylabel,fontsize=fontsize)

     
        ###############
        # End y label
        ###############
        
        ################################
        # End Axis labels               
        ################################

        ################################
        # Tick and grid controls
        ################################
        if "time" not in ds_var[0].dims:
            axes.grid(True)
            axes.grid(which='minor', alpha=0.3,linestyle='--')
            axes.set_xticks(np.linspace(qr[0],qr[1],qr[2]))
            x_tick_labels = [qr[4].format(label) for label in np.linspace(qr[0], qr[1], qr[2])]
            axes.set_xticklabels(x_tick_labels)
            axes.xaxis.set_minor_locator(MultipleLocator(float(qr[3])))
    
            # y ticks
            axes.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=15))
    
            if gph==False or z==False:
                axes.tick_params(which="minor",length=10)
                axes.yaxis.set_major_locator(MultipleLocator(y_r[2]))
            else:
                axes.yaxis.set_minor_locator(MultipleLocator(y_r[2]))
    
    
    
            if zero_tick==True:
                axes.axvline(x=0,linestyle="--",c="k",linewidth=0.75)
                axes.text(0,y_r[0],"0")
                
            axes.tick_params(which="major",labelsize=fontsize,width=1,length=6,direction="out",right=False)
            axes.tick_params(which="minor",labelsize=fontsize,width=1,length=3,direction="out",right=False)
            
            # Grid controls
            axes.grid(True)
            axes.grid(which='minor', alpha=0.2,linestyle='--')
    
    #     if diff_type!=None: ax[2].tick_params(which="major",labelsize=fontsize,width=1,length=6,direction="out",right=False)
        else:
            axes[row,col].grid(True)
            axes[row,col].grid(which='minor', alpha=0.3,linestyle='--')
            axes[row,col].set_xticks(np.linspace(qr[0],qr[1],qr[2]))
            x_tick_labels = [qr[4].format(label) for label in np.linspace(qr[0], qr[1], qr[2])]
            axes[row,col].set_xticklabels(x_tick_labels)
            axes[row,col].xaxis.set_minor_locator(MultipleLocator(float(qr[3])))
    
            # y ticks
            axes[row,col].yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=15))
    
            if gph==False or z==False:
                axes[row,col].tick_params(which="minor",length=10)
            else:
                axes[row,col].yaxis.set_minor_locator(MultipleLocator(y_r[2]))
    
    
    
            if zero_tick==True:
                axes[row,col].axvline(x=0,linestyle="--",c="k",linewidth=0.75)
                axes[row,col].text(0,y_r[0],"0")
                
            axes[row,col].tick_params(which="major",labelsize=fontsize,width=1,length=6,direction="out",right=False)
            axes[row,col].tick_params(which="minor",labelsize=fontsize,width=1,length=3,direction="out",right=False)
    
    #     if diff_type!=None: ax[2].tick_params(which="major",labelsize=fontsize,width=1,length=6,direction="out",right=False)
        
            # Grid controls
            axes[row,col].grid(True)
            axes[row,col].grid(which='minor', alpha=0.2,linestyle='--')

        #############################
        # End tick and grid controls
        ############################
        ## Axis suptitle todo

        ##########################################
        ## End Final subplot controls
        ##########################################   

    

    # Remove empty subplots        
    if "time" in ds_var[0].dims:
        if time_size < num_rows * num_cols:
            for i in range(time_size, num_rows * num_cols):
                fig.delaxes(axes.flatten()[i])

        
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
       
#     # savefig=False
#     if savefig:
#         plt.savefig("/glade/scratch/mmkupilas/Analysis/Plots/"+ filename + ".png")
        
#     return ds_diff
    #============================= End Plot================================
    #=======================================================================

