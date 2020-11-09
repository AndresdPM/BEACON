"""
This is the plotting library for BEACON
"""

import sys, os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats  
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import libraries as lib
from matplotlib import colorbar

def plot_data_maps(galaxy, stars, plotting_options, File_name = './Data.pdf'):

    fig = plt.figure(0, figsize=(5.5,5.5))
    ax_Data = fig.add_subplot(111)
    ax_Data.grid(True)
    
    if plotting_options.wcs:
      Plotting_Coo = stars.Sky_Coo
    else:
      Plotting_Coo = stars.Gref_car_Coo 

    ax_Data.plot(Plotting_Coo[:,0], Plotting_Coo[:,1], '.', ms = 3., color = [1.0, 0.2, 0.2], zorder=0)

    ax_Data.plot(galaxy.Sky_Coo[0], galaxy.Sky_Coo[1], "x", markersize=plotting_options.Markers_size, color=plotting_options.Markers_color, zorder=20)

    try:
       ra_core, dec_core = lib.Create_Ellipse(galaxy.Sky_Coo, galaxy.RH/60., galaxy.E, galaxy.PA, Spherical = True)[0:2]
       ax_Data.plot(ra_core ,dec_core, '--', color=plotting_options.Markers_color, linewidth = plotting_options.Markers_linewidth, zorder=20)
    except:
       pass
    try:
       ra_tidal, dec_tidal = lib.Create_Ellipse(galaxy.Sky_Coo, galaxy.rt/60., galaxy.E, galaxy.PA, Spherical = True)[0:2]
       ax_Data.plot(ra_tidal ,dec_tidal, '--', color=plotting_options.Markers_color, linewidth = plotting_options.Markers_linewidth, zorder=20)
    except:
       pass

    ax_Data.set_xlabel("RA [deg]")
    ax_Data.set_ylabel("Dec [deg]")
    plt.subplots_adjust(left = 0.18, bottom = 0.12, right=0.99, top = 0.88)

    plt.savefig(File_name, bbox_inches="tight")        


def plot_streams(galaxy, stars, plotting_options, KL, RAW_histogram = None, Indices_list = None, File_name = './FIGURAS_PAPER/Circular_Streams.pdf'):

    if Indices_list is None:
        print 'PLOTTING INFO: All stars are being used for quantities calculation'
        Indices_list = [list(np.arange(len(stars.Sky_Coo)))]           
                   
    color1 = [0.25,0.25,0.25]
    color2 = [0.75,0.75,0.75]
    color3 = [.9,0.2,0.2]
    
    norma_kl = 0.5*np.amax(stars.Gref_car_Coo)/np.amax(KL)

    radec_KL = lib.xy2wcs(KL[:,0:2]*norma_kl, galaxy.Sky_Coo, galaxy.Distance)
 
    lin_FeH = np.linspace(plotting_options.Auxiliar_axes_limits[0], plotting_options.Auxiliar_axes_limits[1], 64)
   
    if plotting_options.wcs:
      Plotting_Coo = stars.Sky_Coo
      Plotting_Limits = np.amax(Sky_Coo)
    else:
      Plotting_Coo = stars.Gref_car_Coo 
      Plotting_Limits = [-plotting_options.Maps_plotting_radius, plotting_options.Maps_plotting_radius]
   
    if plotting_options.Draw_markers:
       try:
          ra_core, dec_core = lib.Create_Ellipse(galaxy.Sky_Coo, galaxy.RH/60., galaxy.E, galaxy.PA, Spherical = True)[0:2]
       except:
          pass
       try:
          ra_tidal, dec_tidal = lib.Create_Ellipse(galaxy.Sky_Coo, galaxy.rt/60., galaxy.E, galaxy.PA, Spherical = True)[0:2]
       except:
          pass
    """
    NUMBER OF PLOTS CALCULATION FOR FIGURE SIZE
    """
    N_plots = len(Indices_list)
    figure_Streams_size = [11., 3.59*N_plots]
    params = {'figure.figsize': figure_Streams_size}

    plt.rcParams.update(params)

    fig_Streams = plt.figure()

    for ii, (Indices, kl)  in enumerate(zip(Indices_list, radec_KL)):
        Indices = list(Indices)
        ax_Velocity = plt.subplot(N_plots, 2, 2*ii+1)      
        ax_Auxiliar1 = plt.subplot(N_plots,2 ,2*ii+2)
        ax_Velocity.grid(True)
        if plotting_options.wcs:
           ax_Velocity.set_ylabel('Dec')
        else:
           ax_Velocity.set_ylabel('Y')
        ax_Auxiliar1.xaxis.grid(True)
        ax_Auxiliar1.yaxis.set_ticklabels([])

        if ii < len(Indices_list)-1:
            ax_Auxiliar1.xaxis.set_ticklabels([])
            ax_Velocity.xaxis.set_ticklabels([])

        if ii == 0:
            host_axes = ax_Velocity

        Velocity = ax_Velocity.scatter(Plotting_Coo[Indices,0], Plotting_Coo[Indices,1], s = 10, c = stars.Gref_car_Vel[Indices,2], vmin = plotting_options.Maps_vlos_limits[0], vmax = plotting_options.Maps_vlos_limits[1], cmap = plt.cm.spectral)
        
        ax_Velocity.set_ylim(Plotting_Limits)
        ax_Velocity.set_xlim(Plotting_Limits)
       
        try:
            ax_Velocity.quiver(Plotting_Coo[Indices,0], Plotting_Coo[Indices,1], stars.Gref_car_Vel[Indices,0], stars.Gref_car_Vel[Indices,1])
        except:
            pass

        if plotting_options.Draw_markers:
           ax_Velocity.plot(galaxy.Sky_Coo[0], galaxy.Sky_Coo[1], 'x', color ='r', markeredgewidth=2, ms=5, zorder=30)
           try:
              ax_Velocity.plot(ra_core ,dec_core, '-', color=plotting_options.Markers_color, linewidth = plotting_options.Markers_linewidth, zorder=20)
           except:
              pass
           try:
              ax_Velocity.plot(ra_tidal ,dec_tidal, '--', color=plotting_options.Markers_color, linewidth = plotting_options.Markers_linewidth, zorder=20)
           except:
              pass

        ax_Velocity.annotate("", xytext=(galaxy.Sky_Coo[0], galaxy.Sky_Coo[1]), xy=(kl[0], kl[1]), arrowprops=dict(arrowstyle="->", color='0.1', linewidth = 2.), zorder = 20)

        if RAW_histogram:
           ax_Auxiliar1.step(RAW_histogram[1][0:-1], RAW_histogram[0], where='post', color = [0.1,0.3,0.8], zorder = 1)

        try:
           Histo_metal = stars.FeH[Indices]
        except:
           try:
              Histo_metal = stars.Extra_Coo[Indices]
           except:
              pass
        
        Histo_metal = np.array([item for sublist in Histo_metal for item in sublist])

        count_ax_Auxiliar1, bins_ax_Auxiliar1, patches_ax_Auxiliar1 = ax_Auxiliar1.hist(Histo_metal, bins=50, range = [plotting_options.Auxiliar_axes_limits[0], plotting_options.Auxiliar_axes_limits[1]], color = color3, normed=1, zorder = 0)

        kernel_PDF = stats.gaussian_kde(Histo_metal)

        try:
           kernel_PDF = stats.gaussian_kde(Histo_metal)
           ax_Auxiliar1.plot(lin_FeH, kernel_PDF(lin_FeH)/1.5, '-', linewidth = 2., color = color3)
        except:
              print 'Impossible to derive metallicity kernel.'
        if ii == len(Indices_list)-1:
              if plotting_options.Auxiliar_axes_label:
                 ax_Auxiliar1.set_xlabel(plotting_options.Auxiliar_axes_label)
              else:
                 ax_Auxiliar1.set_xlabel(r'[Fe/H]')

        ax_Auxiliar1.set_xlim(plotting_options.Auxiliar_axes_limits[0], plotting_options.Auxiliar_axes_limits[1])
        
        ax_Auxiliar1.set_ylim([0, 2.])

    plt.subplots_adjust(hspace = .05, wspace = .05)#, top = .95)
    
    cax_Velocity = inset_axes(host_axes,
                   width="100%",
                   height="5%",
                   bbox_transform=host_axes.transAxes,
                   bbox_to_anchor=(0.05, 0.00001, 1.0, 1.12),
                   loc= 1)
      
    norm_Velocity = plt.Normalize(vmin=plotting_options.Maps_vlos_limits[0], vmax=plotting_options.Maps_vlos_limits[1])      
    Velocity_cb = colorbar.ColorbarBase(cax_Velocity, 
                  cmap=plt.cm.spectral, norm=norm_Velocity, 
                  orientation='horizontal', ticklocation='top')
  
    Velocity_cb.set_label(r'$v_{los}$ [km s$^{-1}$]')
    
    if plotting_options.wcs:
        ax_Velocity.set_ylabel('RA')
    else:
        ax_Velocity.set_xlabel('X')
    
    plt.savefig(File_name, bbox_inches="tight")

    return [count_ax_Auxiliar1, bins_ax_Auxiliar1, kernel_PDF]

def plot_angular_momenta(galaxy, KL, MM, plotting_options, File_name = './FIGURAS_PAPER/Angular_momentum.eps'):
    COLORMAP = plt.cm.RdYlBu
    
    if plotting_options.Plot_extra_coo:
       COLORLABEL='Extra coo'
    else:
       COLORLABEL=r'[Fe/H]'
      
    METAL_colorscale = 1./(plotting_options.Auxiliar_axes_limits[1]-plotting_options.Auxiliar_axes_limits[0])*(MM[:,0] - plotting_options.Auxiliar_axes_limits[0])
    COLORSCALE = COLORMAP(METAL_colorscale)    
    COLORMAPABLE = plt.cm.ScalarMappable(cmap=COLORMAP, norm=plt.Normalize(vmin=plotting_options.Auxiliar_axes_limits[0], vmax=plotting_options.Auxiliar_axes_limits[1]))
    COLORMAPABLE._A = []

 
    fig_LK = plt.figure(figsize=(5.5, 5.5))
    ax = fig_LK.add_subplot(111, polar=True)
        
    try:
        """
        MINOR AN MAJOR AXIS PLOT
        """  
	PA_rad = galaxy.PA/180.*np.pi
	ax.plot([PA_rad, PA_rad], [0, 1e10], linestyle = '-.', color = [0.,.5,0.], linewidth = 2, alpha = 0.9)
	ax.plot([PA_rad+np.pi/2., PA_rad+np.pi/2], [0, 1e10], linestyle = '--',  color = [.5,0.,0.], linewidth = 2, alpha = 0.9)
	ax.plot([PA_rad+np.pi, PA_rad+np.pi], [0, 1e10], linestyle = '-.', color = [0.4, 1., 0.4], linewidth = 2, alpha = 0.9)
	ax.plot([PA_rad+3*np.pi/2., PA_rad+3*np.pi/2], [0, 1e10], linestyle = '--',  color = [1.,0.4,0.4], linewidth = 2, alpha = 0.9)
    except:
        pass
    
    Angle = np.arctan2(KL[:,1], KL[:,0])
    Modulus = np.sqrt(np.sum(KL[:,0:2]**2, axis = 1))
   
    for ang, mod, color in zip(Angle, Modulus, COLORSCALE):
        ax.annotate("", xytext=(0.,0.), xy=(np.pi-ang, mod),
            arrowprops=dict(facecolor=color, alpha=0.9, edgecolor = '0.1', linewidth = 0.75))
        
    ax.set_ylim(0,np.amax(Modulus))
   
    ax.set_xticklabels(['W', 'NW', 'N', 'NE', 'E', 'SE', 'S', 'SW'])
    ax.set_yticklabels([])    

    cb = plt.colorbar(COLORMAPABLE, ax=ax, shrink=.75, pad=.1, aspect=30, orientation='horizontal')
    cb.set_label(COLORLABEL)

    plt.tight_layout()
    plt.savefig(File_name)
 
"""
Andres del Pino Molina
"""
