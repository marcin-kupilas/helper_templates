 print("Start ds0 plot")
    print(v0_min,v0_max)


    if fill==True:
        plot0=ds0.plot.contourf(ax=ax[0],
                                ylim=[y_r[0],y_r[1]],
                                cmap ="seismic",
                                vmin=v1_min,vmax=v1_max,levels=nfill, norm=norm,
                                add_colorbar=False,**axis_types)
    else:
        plot0=ds0.plot(ax=ax[0],
                       ylim=[y_r[0],y_r[1]],
                       cmap ="seismic",
                       vmin=v1_min,vmax=v1_max, norm=norm,
                       add_colorbar=False,**axis_types)
        
    
    # Create a colorbar
    print("Start ds0 colorbar")
    cbar0=plt.colorbar(plot0, ax=ax[0], aspect=40, pad=0.015)

    # Set cbar ticks
    custom_ticks=np.linspace(qr[0],qr[1],qr[2])
    cbar0.set_ticks(custom_ticks)
    cbar0.set_ticklabels([f'{tick:.2f}' for tick in custom_ticks])
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
    