#!/usr/bin/env python
# coding: utf-8
# %%
'''
Plot_2D_lat_height.py
this code is designed for plotting CESM output 
that has been reduced to dimensions of lat and height.
TODO: Allow for Pressure and gph

'''


### Module import ###
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import calendar
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


# set the colormap and centre the colorbar

import matplotlib.colors as colors

class MidpointNormalize(colors.Normalize):
 
    """
    Normalise the colorbar so that colorbar is centered on a prescribed midpoint value
    e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    """

    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint=midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)
        
    def __call__(self, value, clip=None):

        # I'm ignoring masked values and all kinds of edge cases to make a simple example
        x, y=[self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

# ========================================================================
# ===== Function Plot_2D_lat_height
# ========================================================================
    
def Plot_2D_lat_height(ds=None, var=None, gph=False, log=False, concentration=False, diff_type=None, 
                       lat_r=[-90,90], y_r=None, nxticks=6, exp=None, qr=None,diffqr=None,
                       cont=True, ncont=10, fill=True, nfill=100, diffcont=10, fontsize=10,cont_label_size=15, savefig=False, filename=None,suptitle=None,z=None,ds1_title=None,ds2_title=None,diff_title=None,qrstr=None,decimal_places=1,decimal_places_diff=1,dx_xminor=0,dy_yminor=0):
   
    
    '''
    NAME:
           Plot_2D_lat_height

    PURPOSE:
           2D plot of lat vs height 

    INPUTS:
           ds: List of datasets (max 2) to be plotted 
           (must contain T in Kelvin if concentration=True, Z in km if gph=true)
           
           var: Variable to plot
           
           concentration=False (default - only applies to species variables)
           
           gph: False (default)

           log: True (default) - log plot for var
           
           diff_type: None (default), if 2 ds given difference type between ds[1] and ds[0]
           
           lat_r: [-90,90] (default)           
           
           y_r=[1e3,1e-6] (default if gph=False)
           
           y_r_z=[0,140] (default if gph=True)
           
           cont=True (default) - switch cont for all plots
           
           ncont=10 (default) - number of cont on all plots
           
           labelsize=15 (default)
           
           savefig=False (default)
           
           fill=True (default) - plots filled contours on all plots, False plots raw data
           
           nfill=number of filled contours on all plots
           
           z = geometric coordinate (bool)
           
           qrstr = bool to control if axis ranges is printed on plot
           
           decimal_places = 1 (default) - integer that controls number of decimal places to display on main plots
            
           decimal_places_diff = 1 (default) - integer that controls number of decimal places to display on diff plots
    '''
        
        
        

        
    # ========================================================================
    # ===== Error check and pass input values to class-accessible values =====
    # ========================================================================
    print("Start error check")
    ## Currently only allowing 3 plots
    assert len(ds)==2, "ATM must have 2 datasets for ds0 and ds1"
    assert diff_type!=None, "ATM must have diff plot"

    ## ds check 
    
    ## Dimensions check
    
    ## xarray, lat, height check
    
    ## Assert no diff if size(ds)=1
    
    # 
    # === End Error check and pass input values to class-accessible values ===
    # ========================================================================


    # =======================================================================
    # ============================ Initial Setup ============================
    # =======================================================================      
    print("Start initial set-up")
    
    ## Set-up datasets
    if ds!=None:
        if len(ds)==2:
            ds0 = ds[0]
            ds1 = ds[1]
            ds0_var=ds[0][var]
            ds1_var=ds[1][var]
            # set_up_plot_data(ds0,ds1) todo
        else:
            ds0=ds[0]
            ds0_var=ds[0][var]
            # set_up_plot_data(ds0) 
    
    ## Convert to concentration
    if concentration==True:
        if var!='T':
            print("Converting ", var, " to concentration")

            # Conversion constants
            k=1.380649e-23 # Boltzmann constant J K^-1= kg m^2 s^-2 K^-1
            to_cm3=1e6 # convertion from m^3 to cm^3
            to_Pa=100 # convert from hPa to Pa. Pa=kg m^-1 s^-2

            # Get pressure data (same for both datasets)
            P=ds0['lev']

            # Get temperature data
            T0=ds0['T']
            T1=ds1['T']

            # Calculate number density, ds0, ds1
            # Units: kg m^-1 s^-2 * kg^-1 m^-2 s^2 K * K^-1=m^-3 | m^-3 * 1e-6=cm^-3
            n0=P*to_Pa/(to_cm3*k*T0) 
            n1=P*to_Pa/(to_cm3*k*T1)

            # Convert molar mixing ratio to number density
            ds0_var=ds0_var*n0
            ds1_var=ds1_var*n1    
    
    ## Compute diffs
    if (diff_type=='abs'): 
        ds_diff=ds1_var - ds0_var
        ds_diff = ds_diff*exp
        
    if (diff_type=='%'): 
        ds_diff=100*(ds1_var - ds0_var)/ds0_var
        
    if (diff_type=='frac'):
        ds_diff=ds1_var/ds0_var
    
    ## Scale vars and other variable controls            
    ds0_var=ds0_var*exp
    ds1_var=ds1_var*exp

           
    ## Compute log of var if not T
    if (var!='T') and (log==True):
        ds0_var=np.log10(ds0_var)
        ds1_var=np.log10(ds1_var)
            
    
    ## Compute var range 
    
    # default
    if y_r==None:
        
        if gph==True:
            y_r=[1e3,1e-6]
        else:
            y_r=[0,140]
    
    if qr==None:

        if gph==True: 
            v0_min=np.min(ds0_var.where((ds0['gph']>=y_r[0]) & (ds0['gph']<=y_r[1])).values)
            v0_max=np.max(ds0_var.where((ds0['gph']>=y_r[0]) & (ds0['gph']<=y_r[1])).values)     
            v1_min=np.min(ds1_var.where((ds1['gph']>=y_r[0]) & (ds1['gph']<=y_r[1])).values)
            v1_max=np.max(ds1_var.where((ds1['gph']>=y_r[0]) & (ds1['gph']<=y_r[1])).values)
        elif z==True:
            v0_min=np.nanmin(ds0_var.where((ds0['z']>=y_r[0]) & (ds0['z']<=y_r[1])).values)
            v0_max=np.nanmax(ds0_var.where((ds0['z']>=y_r[0]) & (ds0['z']<=y_r[1])).values)     
            v1_min=np.nanmin(ds1_var.where((ds1['z']>=y_r[0]) & (ds1['z']<=y_r[1])).values)
            v1_max=np.nanmax(ds1_var.where((ds1['z']>=y_r[0]) & (ds1['z']<=y_r[1])).values)
        else:
            v0_min=np.min(ds0_var.where((ds0['lev']<=y_r[0]) & (ds0['lev']>=y_r[1])).values)
            v0_max=np.max(ds0_var.where((ds0['lev']<=y_r[0]) & (ds0['lev']>=y_r[1])).values) 
            v1_min=np.min(ds1_var.where((ds1['lev']<=y_r[0]) & (ds1['lev']>=y_r[1])).values)
            v1_max=np.max(ds1_var.where((ds1['lev']<=y_r[0]) & (ds1['lev']>=y_r[1])).values)

        v0_min=v0_min
        v0_max=v0_max
        v1_min=v1_min
        v1_max=v1_max 
        
        print(" v0_min ",  v0_min)
        print(" v0_max ",  v0_max)
        print(" v1_min ",  v1_min)
        print(" v1_max ",  v1_max)
        
        # for diagnostic printing on figure
        v0_min_print = v0_min 
        v0_max_print = v0_max
        v1_min_print = v1_min
        v1_max_print = v1_max
         
        ## Choose max and min between datasets - set as limits for all plots
        v0_max=max(v0_max,v1_max)
        v0_min=min(v0_min,v1_min)
        v1_max=v0_max
        v1_min=v0_min
        
        qr=[]
        qr.append(v1_min)
        qr.append(v1_max)
        qr.append(ncont)
        qr.append(5)

        
    else:    
        v0_min=qr[0]
        v0_max=qr[1]     
        v1_min=qr[0]
        v1_max=qr[1]
        
        # for diagnostic printing on figure
        v0_min_print = v0_min 
        v0_max_print = v0_max
        v1_min_print = v1_min
        v1_max_print = v1_max
    
    ## Compute diff range
    if diffqr==None:
        if gph==True:
            diff_min=np.min(ds_diff.where((ds_diff['gph']>=y_r[0]) & (ds_diff['gph']<=y_r[1])))
            diff_max=np.max(ds_diff.where((ds_diff['gph']>=y_r[0]) & (ds_diff['gph']<=y_r[1])))
        elif z==True:
            diff_min=np.min(ds_diff.where((ds_diff['z']>=y_r[0]) & (ds_diff['z']<=y_r[1])))
            diff_max=np.max(ds_diff.where((ds_diff['z']>=y_r[0]) & (ds_diff['z']<=y_r[1])))
        else:
            diff_min=np.min(ds_diff.where((ds_diff['lev']<=y_r[0]) & (ds_diff['lev']>=y_r[1])))
            diff_max=np.max(ds_diff.where((ds_diff['lev']<=y_r[0]) & (ds_diff['lev']>=y_r[1])))
   
        diff_min=diff_min.values
        diff_max=diff_max.values
        
        diffqr=[]
        diffqr.append(diff_min)
        diffqr.append(diff_max)
        diffqr.append(ncont)
        diffqr.append(1)
        
    
    else:
        
        diff_min=diffqr[0]
        diff_max=diffqr[1]
           
        
        

    
    ### Debugging Initial Setup
    print("qr0 qr1 qr2", qr[0], qr[1], qr[2], "\n")
    print("diffqr vals", diffqr[0], diffqr[1],"\n")
    # ========================== End Initial Setup ==========================
    # =======================================================================
    
    # =======================================================================
    # =============================== Plot ==================================
    # =======================================================================      
    
    print("Start plotting script")
    
    ## Set-up axes - todo - currently not working as script requires diff - see set-up
    if diff_type!=None:
        icols=3
    else:
        icols=np.size(ds)      
    fig, ax=plt.subplots(nrows=1,ncols=icols, figsize=(7.5*icols,6))
    fig.set_dpi(150)
    
    if gph==True:
        axis_types={"y":"gph"}
    elif z==True:
        axis_types={"y":"z"}
    else:
        axis_types={"y":"lev","yscale":"log"}
        
    print("Finish setting up axes")
    
    
    ##########################################
    ## ds0 plot
    ##########################################
    
    print("Start ds0 plot")
    print(v0_min,v0_max)
    
    if fill==True:
        plot0=ds0_var.plot.contourf(ax=ax[0], **axis_types,
                        ylim=[y_r[0],y_r[1]],
                        cmap="turbo",
                        vmin=v1_min,vmax=v1_max,levels=nfill,
                        add_colorbar=False)
    else:
        plot0=ds0_var.plot(ax=ax[0], **axis_types,
                ylim=[y_r[0],y_r[1]],
                cmap="turbo",
                vmin=v1_min,vmax=v1_max,
                add_colorbar=False)
     
    # Create a colorbar
    print("Start ds0 colorbar")
    cbar=plt.colorbar(plot0, ax=ax[0])

    # Set cbar ticks
    custom_ticks=np.linspace(qr[0],qr[1],qr[2])
    cbar.set_ticks(custom_ticks)
    cbar.ax.yaxis.set_minor_locator(MultipleLocator(float(qr[3])))
    cbar.set_ticklabels(['{:.{}f}'.format(tick, decimal_places) for tick in custom_ticks])
    cbar.ax.tick_params(labelsize=fontsize)
    
    # Plot cont  
    if cont:
        if gph==True:  
            # lined cont
            cs_ds0=ds0_var.plot.contour(**axis_types,
                                                vmin=v0_min,vmax=v0_max,
                                                ylim=[y_r[0],y_r[1]],
                                                levels=ncont,
                                                colors="k",
                                                linewidths=1,
                                                ax=ax[0],
                                                # norm=log_norm, # uncomment for log scale - not working properly for cont!
                                                )
        else:
            # lined cont
            cs_ds0=ds0_var.plot.contour(**axis_types,
                                                vmin=v0_min,vmax=v0_max,
                                                ylim=[y_r[0],y_r[1]],
                                                levels=ncont,
                                                colors="k",
                                                linewidths=1,
                                                ax=ax[0],
                                                # norm=log_norm, # uncomment for log scale - not working properly for cont!
                                                )
            

        # color label controls
        fmt={}
        for l in cs_ds0.levels:
            if var!='T':
                # l=l*exp
                fmt[l] ='{:.{}f}'.format(l, decimal_places)  # uncomment for nolog
            else:
                fmt[l] ='{:.{}f}'.format(l, decimal_places) 

            # fmt[l] =f"{l:.3f}" # uncomment for log
        clabels=ax[0].clabel(cs_ds0, cs_ds0.levels, inline=True, fmt=fmt, fontsize=cont_label_size)
        [txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0)) for txt in clabels]

    print("Finish _var0 plot")
    
    ##########################################
    ## End ds0 plot
    ##########################################    
    
    
    ##########################################
    ## ds1 plot
    ########################################## 
    if fill==True:
        plot1=ds1_var.plot.contourf(ax=ax[1],**axis_types,
                                  ylim=[y_r[0],y_r[1]],
                                  cmap="turbo",
                                  vmin=v1_min,vmax=v1_max,levels=nfill,
                                  add_colorbar=False)
    else:
        plot1=ds1_var.plot(ax=ax[1],**axis_types,
                                  ylim=[y_r[0],y_r[1]],
                                  cmap="turbo",
                                  vmin=v1_min,vmax=v1_max,
                                  add_colorbar=False)
    # Create a colorbar
    cbar=plt.colorbar(plot1, ax=ax[1])

    # Set cbar ticks
    custom_ticks=np.linspace(qr[0],qr[1],qr[2])
    cbar.set_ticks(custom_ticks)
    cbar.ax.yaxis.set_minor_locator(MultipleLocator(float(qr[3])))
    cbar.set_ticklabels(['{:.{}f}'.format(tick, decimal_places) for tick in custom_ticks])
    cbar.ax.tick_params(labelsize=fontsize)    
    
    # Plot cont
    if cont:
        # lined cont
        cs_ds1=ds1_var.plot.contour(**axis_types,
                                       vmin=v1_min,vmax=v1_max,
                                       ylim=[y_r[0],y_r[1]],
                                       levels=ncont,
                                       colors="k",
                                       linewidths=1,
                                       ax=ax[1],
                                        # norm=log_norm, # uncomment for log scale - not working properly for cont!
                                       )

        # color label controls
        fmt={}
        for l in cs_ds1.levels:
            if var!='T':
                # l=l*exp
                fmt[l] ='{:.{}f}'.format(l, decimal_places) # uncomment for nolog
            else:
                fmt[l] ='{:.{}f}'.format(l, decimal_places)
        clabels=ax[1].clabel(cs_ds1, cs_ds1.levels, inline=True, fmt=fmt, fontsize=cont_label_size)
        [txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0)) for txt in clabels]

    ##########################################
    ## End ds1 plot
    ##########################################
    
    
    ##########################################
    ## diff plot
    ##########################################
    
    if diff_type=="%" or diff_type=="abs": 
        norm=MidpointNormalize(midpoint=0)
    if diff_type=="frac": 
        norm= MidpointNormalize(midpoint=1) 
        
    if fill==True:
        plot2=ds_diff.plot.contourf(ax=ax[2],
                              **axis_types,
                              ylim=[y_r[0],y_r[1]],
                              cmap="seismic",
                              vmin=diff_min,vmax=diff_max,
                              levels=nfill,
                              add_colorbar=False,
                              norm=norm,
                              )
    else:
        plot2=ds_diff.plot(ax=ax[2],
                              **axis_types,
                              ylim=[y_r[0],y_r[1]],
                              cmap="seismic",
                              vmin=diff_min,vmax=diff_max,
                              add_colorbar=False,
                              norm=norm,
                              )

    # Create a colorbar
    cbar=plt.colorbar(plot2, ax=ax[2])

    # Set cbar ticks
    custom_ticks=np.linspace(diffqr[0],diffqr[1],diffqr[2])
    cbar.set_ticks(custom_ticks)
    cbar.ax.yaxis.set_minor_locator(MultipleLocator(float(diffqr[3])))

    cbar.set_ticklabels(['{:.{}f}'.format(tick, decimal_places_diff) for tick in custom_ticks])   
    cbar.ax.tick_params(labelsize=fontsize)
    
    # Plot cont
    if cont:
        # lined cont
        cs_diff=ds_diff.plot.contour(**axis_types,
                                        vmin=diff_min,vmax=diff_max,
                                        ylim=[y_r[0],y_r[1]],
                                        levels=diffcont,
                                        colors="k",
                                        linewidths=1,
                                        ax=ax[2],)

        # color label controls
        fmt={}
        for l in cs_diff.levels:
            # l=l # uncomment for log
            fmt[l] ='{:.{}f}'.format(l, decimal_places_diff) # uncomment for nolog
            # fmt[l] =f"{l:.3f}" # uncomment for log
        clabels=ax[2].clabel(cs_diff, cs_diff.levels, inline=True, fmt=fmt, fontsize=cont_label_size)
        [txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0)) for txt in clabels]


        # Plot contour at 0 value    
        cs_diff_0=ds_diff.plot.contour(**axis_types,
                                        vmin=0,vmax=0,
                                        ylim=[y_r[0],y_r[1]],
                                        levels=1,
                                        colors="white",
                                        linewidths=2,
                                        ax=ax[2],
                                        )


        # color label controls
        # fmt={}
        # for l in cs_diff_0.levels:
        #     # l=l # uncomment for log
        #     fmt[l] ='{:.{}f}'.format(l, decimal_places_diff)# uncomment for nolog
        #     # fmt[l] =f"{l:.3f}" # uncomment for log
        # clabels_0=ax[2].clabel(cs_diff_0, cs_diff_0.levels, inline=True, fmt=fmt, fontsize=cont_label_size)
        # [txt.set_bbox(dict(facecolor='k', edgecolor='none', pad=0)) for txt in clabels_0]
    
    ##########################################
    ## End diff plot
    ##########################################


    
    ##########################################
    ## Final plot controls
    ##########################################    
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.25)  # Adjust the value as needed
    
    ## Axis labels
    if qrstr==True:
        qrstr = (
        "v0_min  " + str(v0_min_print) + "\n"
        "v0_max  " + str(v0_max_print) + "\n"
        "v1_min  " + str(v1_min_print) + "\n"
        "v1_max  " + str(v1_max_print) + "\n"
        "qr0 qr1  " + str(qr[0]) + " " + str(qr[1]) + "\n"
        "diffqr vals  " + str(diffqr[0]) + " " + str(diffqr[1])
        )
    
    if qrstr==True:        
        ax[0].set_xlabel("latitude (degrees north)" + "\n" + qrstr,fontsize=fontsize)
    else:
        ax[0].set_xlabel("latitude (degrees north)",fontsize=fontsize)
    ax[1].set_xlabel("latitude (degrees north)",fontsize=fontsize)
    ax[2].set_xlabel("latitude (degrees north)",fontsize=fontsize)
    
    if gph==True:       
        ax[0].set_ylabel("GPH (km)",fontsize=fontsize)
        ax[1].set_ylabel("GPH (km)",fontsize=fontsize)
        ax[2].set_ylabel("GPH (km)",fontsize=fontsize)
        
    elif z==True:
        ax[0].set_ylabel("geometric height (km)",fontsize=fontsize)
        ax[1].set_ylabel("geometric height (km)",fontsize=fontsize)
        ax[2].set_ylabel("geometric height (km)",fontsize=fontsize)
    
    else:
        ax[0].set_ylabel("Pressure (hPa)",fontsize=fontsize)
        ax[1].set_ylabel("Pressure (hPa)",fontsize=fontsize)
        ax[2].set_ylabel("Pressure (hPa)",fontsize=fontsize)

    
    # Axis ticks
    ax[0].set_xticks(np.linspace(lat_r[0],lat_r[1],nxticks))
    ax[1].set_xticks(np.linspace(lat_r[0],lat_r[1],nxticks))
    ax[2].set_xticks(np.linspace(lat_r[0],lat_r[1],nxticks))
    
    ax[0].tick_params(which="major",labelsize=fontsize,width=1,length=7,direction="out",right=False)
    ax[1].tick_params(which="major",labelsize=fontsize,width=1,length=7,direction="out",right=False)
    ax[2].tick_params(which="major",labelsize=fontsize,width=1,length=7,direction="out",right=False)
    
    ax[0].xaxis.set_minor_locator(MultipleLocator(dx_xminor))
    ax[1].xaxis.set_minor_locator(MultipleLocator(dx_xminor))
    ax[2].xaxis.set_minor_locator(MultipleLocator(dx_xminor))

    ax[0].yaxis.set_minor_locator(MultipleLocator(dy_yminor))
    ax[1].yaxis.set_minor_locator(MultipleLocator(dy_yminor))
    ax[2].yaxis.set_minor_locator(MultipleLocator(dy_yminor))    
    
    ax[0].tick_params(which="minor",length=4)
    ax[1].tick_params(which="minor",length=4)
    ax[2].tick_params(which="minor",length=4)


    
    # Axis titles
    if var!='T':            
        if concentration==True:
            
            # Calculate scale str based on exponent
            for i in range(2,18):
                str_exp = "{:.0e}".format(exp)
                iexp = 10**-i
                iexp = "{:.0e}".format(iexp)
                if iexp==str_exp:
                    scalestr = fr" (x10$^{{{i}}}$)"
                    break
                else:
                    scalestr=""
                    
            if log==True:
                ax[0].set_title("log10(" +var + ") concentration ds0",fontsize=fontsize)
                ax[1].set_title("log10(" +var + ") concentration ds1",fontsize=fontsize)
            else:
                ax[0].set_title(var +  " concentration ds0" + scalestr,fontsize=fontsize)
                ax[1].set_title(var + " concentration ds1" + scalestr,fontsize=fontsize)
        else: #if concentration          
                
            # Calculate scale str based on exponent
            for i in range(2,18):
                str_exp = "{:.0e}".format(exp)
                iexp = 10**i
                iexp = "{:.0e}".format(iexp)
                if iexp==str_exp:
                    scalestr = fr" (x10$^{{-{i}}}$)"
                    break
                else:
                    scalestr=""
            print("scalestr", scalestr)
            # if ppm or ppb
            if exp==1e6: scalestr="ppm"
            if exp==1e9: scalestr="ppb"

            if log==True:
                ax[0].set_title("log10(" +var + ") mixing ratio ds0 " + scalestr,fontsize=fontsize)
                ax[1].set_title("log10(" +var + ") mixing ratio ds1 " + scalestr,fontsize=fontsize)
            else:
                ax[0].set_title(var + " mixing ratio ds0 " + scalestr,fontsize=fontsize)
                ax[1].set_title(var + " mixing ratio ds1 " + scalestr,fontsize=fontsize)
                
            ax[2].set_title(diff_type + " difference",fontsize=fontsize)

    else:
        ax[0].set_title(var + " ds0",fontsize=fontsize)
        ax[1].set_title(var + " ds1",fontsize=fontsize)
        ax[2].set_title(diff_type + " difference",fontsize=fontsize)
        
    if ds1_title!=None:
        ax[0].set_title(ds1_title,fontsize=fontsize)
    if ds2_title!=None:
        ax[1].set_title(ds2_title,fontsize=fontsize)
    if diff_title!=None:
        ax[2].set_title(diff_title,fontsize=fontsize)



    # text box with additional information
    #fig.text(0, -0.2, qrstr,horizontalalignment='left', fontsize=fontsize,wrap=True ) 
    
    # Define the position and style of the text box
    fig.suptitle(filename, x=0.5, y=1.05, fontsize=fontsize, ha='center')  
    print("filename", filename)
    
    # savefig=False
    if savefig:
        filanem = filename + ".png"
        plt.savefig("/glade/scratch/mmkupilas/Analysis/Plots/"+ filename,dpi=150, bbox_inches="tight")
    # ============================= End Plot  ===============================
    # =======================================================================
    
# ========================================================================
# ===== Function Plot_2D_lat_height
# ========================================================================
    

