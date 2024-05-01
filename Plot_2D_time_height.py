#!/usr/bin/env python
# coding: utf-8
# %%
'''
Plot_2D_time_height.py
this code is designed for plotting CESM output 
that has been reduced to dimensions of time and height.
TODO: Allow for Pressure and gph

'''


### Module import ###
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import calendar
import matplotlib.colors as colors
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as ticker

# set the colormap and centre the colorbar


class MidpointNormalize(colors.Normalize):
 
    """
    Normalise the colorbar so that colorbar is centered on a prescribed midpoint value
    e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    """

    def __init__(self, vmin=None, vmax=None, midpoint=None, clip =False):
        self.midpoint=midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)
        
    def __call__(self, value, clip =None):

        # I'm ignoring masked values and all kinds of edge cases to make a simple example
        x, y=[self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))
    

def Plot_2D_time_height(ds=None, var=None, gph=False, log=False, concentration=False, diff_type=None, 
                       lat_r=[-90,90], y_r=None, nxticks=6, exp=None, qr=None,diffqr=None,
                       cont=True, ncont=10, fill=True, nfill=100, diffcont=10, fontsize=10,cont_label_size=15, savefig=False, filename=None,suptitle=None,x_size=6.5,y_size=8):
   
    
    '''
    NAME:
           Plot_2D_time_height

    PURPOSE:
           2D plot of time vs height 

    INPUTS:
           ds: List of datasets (max 2) to be plotted 
           (must contain T in Kelvin if concentration=True, Z in km if gph=true)
           
           var: Variable to plot
           
           concentration=False (default - only applies to species variables)
           
           gph: False (default)

           log: True (default) - log plot for var
           
           diff_type: None (default), if 2 ds given difference type between ds[1] and ds[0]
           
           lat_r: [-90,90] (default)           
           
           y_r=None, (y axis range - default calculated below)
           
           cont=True (default) - switch cont for all plots
           
           ncont=10 (default) - number of cont on all plots
           
           labelsize=15 (default)
           
           savefig=False (default)
           
           fill=True (default) - plots filled contours on all plots, False plots raw data
           
           nfill=number of filled contours on all plots
           
    NOTES:
    2023-08-14 Mon: MMK - currently titles set to ds0 and ds1.
    
    '''
    print("WARNING - see notes")    
        
        

        
    #========================================================================
    #===== Error check and pass input values to class-accessible values=====
    #========================================================================
    print("Start error check")
    ## Currently only allowing 3 plots
    assert len(ds)==2, "ATM must have 2 datasets for ds0 and ds1"

    ## ds check 
    
    ## Dimensions check
    
    ## xarray, lat, height check
    
    ## Assert no diff if size(ds)=1
    
    
    #=== End Error check and pass input values to class-accessible values===
    #========================================================================


    #=======================================================================
    #============================ Initial Setup ============================
    #=======================================================================      
    print("Start initial set-up")
    ## Set-up datasets
    if ds!=None:
        if len(ds)==2:
            ds0=ds[0][var]
            ds1=ds[1][var]
            # set_up_plot_data(ds0,ds1) todo
        else:
            ds0=ds[0][var]
            # set_up_plot_data(ds0) 
    
    ## Convert to concentration
    # Note - if want to compare two datasets, one which already is
    # in concencentration units, and one that isn't, I've had to 
    # implement a workaround to ignore files that don't need converting
    if concentration==True:
        for data in ds:
            if "skip_convert" not in data.attrs:
                print("Converting ", var, " to concentration for ", data.attrs['ds_tag'])

                # Conversion constants
                k=1.380649e-23 # Boltzmann constant J K^-1= kg m^2 s^-2 K^-1
                to_cm3=1e6 # convertion from m^3 to cm^3
                to_Pa=100 # convert from hPa to Pa. Pa=kg m^-1 s^-2

                # Get pressure data (same for both datasets)
                P=data['lev']

                # Get temperature data
                T=data['T']

                # Calculate number density, ds0, ds1
                # Units: kg m^-1 s^-2 * kg^-1 m^-2 s^2 K * K^-1=m^-3 | m^-3 * 1e-6=cm^-3
                n=P*to_Pa/(to_cm3*k*T) 

                # Convert molar mixing ratio to number density
                if data.attrs['ds_tag']=="ds0":
                    ds0=data[var]*n
                else:
                    ds1=data[var]*n
            else:
                if data.attrs['ds_tag']=="ds0":
                    ds0=data[var]
                else:
                    ds1=data[var]                

            
    else:
        ds0=ds0[var]
        ds1=ds1[var]
            
    ## Scale vars and other variable controls
             
    ds0=ds0*exp
    ds1=ds1*exp
    
    ## Compute diffs
    if diff_type!=None:
        if (diff_type=='abs'): 
            ds_diff=ds1 - ds0
            print("Calculating abs ds_diff")

        if (diff_type=='%'): 
            ds_diff=100*(ds1 - ds0)/ds0

        if (diff_type=='frac'):
            ds_diff=ds1/ds0
    else:
        ds_diff = None
        
    # Diff operation removes gph if manually calculated.
    # Need to add again.  
    if diff_type!=None and gph==True:
        print("\nWARNING: currently using arithmetic mean \n\
        of gph between ds0['gph'] and d1['gph'] as ds_diff['gph'] \n ")
        gph_diff = 0.5*(ds0['gph']+ds1['gph'])
        ds_diff['gph'] = gph_diff.assign_coords(gph=gph_diff)
    
    
    ## Compute log of var if not T
    if log==True:
        ds0=np.log10(ds0)
        ds1=np.log10(ds1)
            
    
    ## Compute var range
    
    # default
    if y_r==None:
        if gph==True:
            y_r=[0,150]
        else:
            y_r=[1e3,1e-6]    
            
    if qr==None:    
        if gph==True: 
            v0_min=np.min(ds0.where((ds0['gph']>=y_r[0]) & (ds0['gph']<=y_r[1])))
            v0_max=np.max(ds0.where((ds0['gph']>=y_r[0]) & (ds0['gph']<=y_r[1])))     
            v1_min=np.min(ds1.where((ds1['gph']>=y_r[0]) & (ds1['gph']<=y_r[1])))
            v1_max=np.max(ds1.where((ds1['gph']>=y_r[0]) & (ds1['gph']<=y_r[1])))

        else:
            v0_min=np.min(ds0.where((ds0['lev']<=y_r[0]) & (ds0['lev']>=y_r[1])))
            v0_max=np.max(ds0.where((ds0['lev']<=y_r[0]) & (ds0['lev']>=y_r[1]))) 
            v1_min=np.min(ds1.where((ds1['lev']<=y_r[0]) & (ds1['lev']>=y_r[1])))
            v1_max=np.max(ds1.where((ds1['lev']<=y_r[0]) & (ds1['lev']>=y_r[1])))


        v0_min=v0_min.values
        v0_max=v0_max.values
        v1_min=v1_min.values
        v1_max=v1_max.values  

        
        ## Choose max and min between datasets - set as limits for all plots
        v0_max=max(v0_max,v1_max)
        v0_min=min(v0_min,v1_min)
        v1_max=v0_max
        v1_min=v0_min
        
        qr=[]
        qr.append(v1_min)
        qr.append(v1_max)
        qr.append(ncont)
        
    else:    
        v0_min=qr[0]
        v0_max=qr[1]     
        v1_min=qr[0]
        v1_max=qr[1]
    print("v0_min, v0_max", v0_min,v0_max)

    ## Compute diff range
    if diff_type!=None:
        if diffqr==None:
            if gph==True:
                diff_min=np.min(ds_diff.where((ds_diff['gph']>=y_r[0]) & (ds_diff['gph']<=y_r[1])))
                diff_max=np.max(ds_diff.where((ds_diff['gph']>=y_r[0]) & (ds_diff['gph']<=y_r[1])))
            else:
                diff_min=np.min(ds_diff.where((ds_diff['lev']<=y_r[0]) & (ds_diff['lev']>=y_r[1])))
                diff_max=np.max(ds_diff.where((ds_diff['lev']<=y_r[0]) & (ds_diff['lev']>=y_r[1])))

            diff_min=diff_min.values
            diff_max=diff_max.values

            diffqr=[]
            diffqr.append(diff_min)
            diffqr.append(diff_max)
            diffqr.append(ncont)

        else:

            diff_min=diffqr[0]
            diff_max=diffqr[1]
           
        
        

    
    ### Debugging Initial Setup
    print("qr0 qr1 qr2", qr[0], qr[1], qr[2], "\n")
    if diff_type!=None:
        print("diff_min diff_max ", diff_min, diff_max,"\n")
    #========================== End Initial Setup ==========================
    #=======================================================================
    
    #=======================================================================
    #=============================== Plot===================================
    #=======================================================================      
    
    print("Start plotting script\n")
    
    ## Set-up axes - todo - currently not working as script requires diff - see set-up
    if diff_type!=None:
        icols=3
    else:
        icols=len(ds)  
        
    fig, ax=plt.subplots(nrows=1,ncols=icols, figsize=(x_size*icols,y_size)) 
    fig.set_dpi(300)
    
    if suptitle!=None:
        fig.suptitle(suptitle)
    
    # height coordinate
    if gph==True:
        axis_types = {"y":"gph", "x":"time"}
    else:
        axis_types = {"y":"lev","yscale":"log", "x":"time"}

    print("Finished setting up axes")
    
    
    ##########################################
    ## ds0 plot
    ##########################################
    
    print("Start ds0 plot")
    print(v0_min,v0_max)

    if fill==True:
        plot0=ds0.plot.contourf(ax=ax[0],
                                ylim=[y_r[0],y_r[1]],
                                cmap ="turbo",
                                vmin=v1_min,vmax=v1_max,levels=nfill,
                                add_colorbar=False,**axis_types)
    else:
        plot0=ds0.plot(ax=ax[0],
                       ylim=[y_r[0],y_r[1]],
                       cmap ="turbo",
                       vmin=v1_min,vmax=v1_max,
                       add_colorbar=False,**axis_types)
        
    
    # Create a colorbar
    print("Start ds0 colorbar")
    cbar0=plt.colorbar(plot0, ax=ax[0], aspect=40, pad=0.015)

    # Set cbar ticks
    custom_ticks=np.linspace(qr[0],qr[1],qr[2])
    cbar0.set_ticks(custom_ticks)
    cbar0.set_ticklabels([f'{tick:.1f}' for tick in custom_ticks])
    cbar0.ax.tick_params(labelsize=fontsize)
    cbar0.ax.tick_params(which="minor",length=0)
    
    # Plot cont    
    if cont:
        # lined cont
        cs_nonrr=ds0.plot.contour(vmin=v0_min,vmax=v0_max,
                                  ylim=[y_r[0],y_r[1]],
                                  levels=ncont,
                                  colors="k",
                                  linewidths=1,
                                  ax=ax[0],
                                  **axis_types
                                 )

        # color label controls
        fmt={}
        for l in cs_nonrr.levels:
            fmt[l]=f"{l:.2f}"


        clabels=ax[0].clabel(cs_nonrr, cs_nonrr.levels, inline=True, fmt=fmt, fontsize=cont_label_size)
        [txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0)) for txt in clabels]

    print("Finish ds0 plot")
    
    ##########################################
    ## End ds0 plot
    ##########################################    
    
    
    ##########################################
    ## ds1 plot
    ########################################## 
    if fill==True:
        plot1=ds1.plot.contourf(ax=ax[1],
                                ylim=[y_r[0],y_r[1]],
                                cmap ="turbo",
                                vmin=v1_min,vmax=v1_max,levels=nfill,
                                add_colorbar=False,**axis_types)
    else:
        plot1=ds1.plot(ax=ax[1],
                       ylim=[y_r[0],y_r[1]],
                       cmap ="turbo",
                       vmin=v1_min,vmax=v1_max,
                       add_colorbar=False,**axis_types)
    # Create a colorbar
    cbar1=plt.colorbar(plot1, ax=ax[1], aspect=40,pad = 0.015)

    # Set cbar ticks
    custom_ticks=np.linspace(qr[0],qr[1],qr[2])
    cbar1.set_ticks(custom_ticks)
    cbar1.set_ticklabels([f'{tick:.1f}' for tick in custom_ticks])
    cbar1.ax.tick_params(labelsize=fontsize)    
    cbar1.ax.tick_params(which="minor",length=0)
    
    # Plot cont
    if cont:
        # lined cont
        cs_rr=ds1.plot.contour(vmin=v1_min,vmax=v1_max,
                               ylim=[y_r[0],y_r[1]],
                               levels=ncont,
                               colors="k",
                               linewidths=1,
                               ax=ax[1],
                               **axis_types)

        # color label controls
        fmt={}
        for l in cs_rr.levels:
            fmt[l]=f"{l:.2f}"
  
  
        clabels=ax[1].clabel(cs_rr, cs_rr.levels, inline=True, fmt=fmt, fontsize=cont_label_size)
        [txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0)) for txt in clabels]

    ##########################################
    ## End ds1 plot
    ##########################################
    
    
    ##########################################
    ## diff plot
    ##########################################
    if diff_type!=None:
        if diff_type=="%" or diff_type=="abs": 
            norm=MidpointNormalize(midpoint=0)
        if diff_type=="frac": 
            norm= MidpointNormalize(midpoint=1) 

        if fill==True:
            plot2=ds_diff.plot.contourf(ax=ax[2],
                                        ylim=[y_r[0],y_r[1]],
                                        cmap ="seismic",
                                        vmin=diff_min,vmax=diff_max,
                                        levels=nfill,
                                        add_colorbar=False,
                                        norm=norm,
                                        **axis_types)
        else:
            plot2=ds_diff.plot(ax=ax[2],
                               ylim=[y_r[0],y_r[1]],
                               cmap ="seismic",
                               vmin=diff_min,vmax=diff_max,
                               add_colorbar=False,
                               norm=norm,
                               **axis_types)

        # Create a colorbar
        cbar2=plt.colorbar(plot2, ax=ax[2],aspect=40, pad = 0.015)

        # Set cbar ticks
        custom_ticks=np.linspace(diffqr[0],diffqr[1],diffqr[2])
        cbar2.set_ticks(custom_ticks)
        cbar2.set_ticklabels([f'{tick:.1f}' for tick in custom_ticks])   
        cbar2.ax.tick_params(labelsize=fontsize)
        cbar2.ax.tick_params(which="minor",length=0)

        # Plot cont
        if cont:
            # lined cont
            cs_diff=ds_diff.plot.contour(vmin=diff_min,vmax=diff_max,
                                         ylim=[y_r[0],y_r[1]],
                                         levels=diffcont,
                                         colors="k",
                                         linewidths=1,
                                         ax=ax[2],**axis_types,norm=norm)

            # color label controls
            fmt={}
            for l in cs_diff.levels:
                fmt[l]=f"{l:.2f}"

            clabels=ax[2].clabel(cs_diff, cs_diff.levels, inline=True, fmt=fmt, fontsize=cont_label_size)
            [txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0)) for txt in clabels]


            # Plot contour at 0 value    
            cs_diff_0=ds_diff.plot.contour(vmin=0,vmax=0,
                                           ylim=[y_r[0],y_r[1]],
                                           levels=1,
                                           colors="white",
                                           linewidths=2,
                                           ax=ax[2],
                                           **axis_types)


            # color label controls
            fmt={}
            for l in cs_diff_0.levels:
                fmt[l]=f"{l:.2f}" 

            clabels_0=ax[2].clabel(cs_diff_0, cs_diff_0.levels, inline=True, fmt=fmt, fontsize=cont_label_size)
            [txt.set_bbox(dict(facecolor='k', edgecolor='none', pad=0)) for txt in clabels_0]

    ##########################################
    ## End diff plot
    ##########################################


    
    ##########################################
    ## Final plot controls
    ##########################################    
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.25)  # Adjust the value as needed
    
    # Axis labels    
    ax[0].set_xlabel("Month",fontsize=fontsize)
    ax[1].set_xlabel("Month",fontsize=fontsize)
    if diff_type!=None: ax[2].set_xlabel("Month",fontsize=fontsize)
    
    if gph==True:
        ylabel="GPH (km)"
    else:
        ylabel="Pressure (hPa)"
        
    ax[0].set_ylabel(ylabel,fontsize=fontsize)
    ax[1].set_ylabel(ylabel,fontsize=fontsize)
    if diff_type!=None: ax[2].set_ylabel(ylabel,fontsize=fontsize)

    
    # Axis ticks

    # x ticks
    time_values0=ds0.time.values
    time_values1=ds1.time.values
    
    ax[0].set_xticks(time_values0)
    ax[1].set_xticks(time_values1)
    
    if diff_type!=None: ax[2].set_xticks(time_values1)
    
    tick_labels=['J','F','M','A','M','J','J','A','S','O','N','D']
    
    ax[0].set_xticklabels(tick_labels,rotation=0)
    ax[1].set_xticklabels(tick_labels,rotation=0)    
    if diff_type!=None: ax[2].set_xticklabels(tick_labels,rotation=0)

    # y ticks
    if gph==False:
        ax[0].tick_params(which="minor",length=3)
        ax[1].tick_params(which="minor",length=3)
    else:
        ax[0].tick_params(which="minor",length=0)
        ax[1].tick_params(which="minor",length=0)
        
    # diff ticks
    if diff_type!=None: 
        ax[2].tick_params(which="major",labelsize=fontsize,width=1,length=6,direction="out",right=False)
        if gph==False:
            ax[2].tick_params(which="minor",length=3)
        else:
            ax[2].tick_params(which="minor",length=0)

    # Global tick params
    ax[0].tick_params(which="major",labelsize=fontsize,width=1,length=6,direction="out",right=False)
    ax[1].tick_params(which="major",labelsize=fontsize,width=1,length=6,direction="out",right=False)
    if diff_type!=None:
        ax[2].tick_params(which="major",labelsize=fontsize,width=1,length=6,direction="out",right=False)


    # Axis titles
    if var!='T':            
        if concentration==True:
            
            # Calculate scale str based on exponent
            for i in range(2,18):
                str_exp = "{:.0e}".format(exp)
                iexp = 10**-i
                iexp = "{:.0e}".format(iexp)
                if iexp==str_exp:
                    scalestr=f" (x10^{i})"
                    break
                else:
                    scalestr=""

            if log==True:
                ax[0].set_title("log10(" +var + ") concentration ds0",fontsize=fontsize)
                ax[1].set_title("log10(" +var + ") concentration ds1",fontsize=fontsize)
            else:
                ax[0].set_title(var +  " concentration ds0" + scalestr,fontsize=fontsize)
                ax[1].set_title(var + " concentration ds1" + scalestr,fontsize=fontsize)
                
            if diff_type!=None: ax[2].set_title(diff_type + " difference",fontsize=fontsize)

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
                ax[0].set_title("log10(" +var + ") mixing ratio ds0 " + scalestr,fontsize=fontsize)
                ax[1].set_title("log10(" +var + ") mixing ratio ds1 " + scalestr,fontsize=fontsize)
            else:
                ax[0].set_title(var + " mixing ratio ds0 " + scalestr,fontsize=fontsize)
                ax[1].set_title(var + " mixing ratio ds1 " + scalestr,fontsize=fontsize)
                
            if diff_type!=None: ax[2].set_title(diff_type + " difference",fontsize=fontsize)

    else:
        ax[0].set_title(var + " " + ds[0][var].attrs['ds_tag'],fontsize=fontsize)
        ax[1].set_title(var + " " + ds[1][var].attrs['ds_tag'],fontsize=fontsize)
        if diff_type!=None: ax[2].set_title(diff_type + " difference",fontsize=fontsize)
        
    # savefig=False
    if savefig:
        plt.savefig("/glade/scratch/mmkupilas/Analysis/Plots/"+ filename + ".png")
        
    return ds_diff
    #============================= End Plot================================
    #=======================================================================

