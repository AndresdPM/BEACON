"""
BEACON
"""

import sys, os
import numpy as np
import libraries as lib

def parameters_generator(Filename='Options'):
    import imp

    options = imp.load_source(os.path.splitext(os.path.basename(Filename))[0], Filename)

    class galaxy:
        def __init__(self):
            self.Name = options.Galaxy_Name
            self.Sky_Coo = np.array(options.Galaxy_Coordinates)
            self.Distance =  options.Galaxy_Distance
            self.RH = options.Galaxy_rh
            self.RT = options.Galaxy_rt
            self.PA = options.Galaxy_pa
            self.E =  options.Galaxy_e
            self.Sky_Vel =  options.Galaxy_Velocity
            
    class read_data_options:
        def __init__(self):
            self.Data_file = options.Galaxy_Data_file
            self.Clean_sample = options.Clean_sample
            self.Admisible_range_Vel = options.Admisible_range_Velocity
            self.Admisible_max_error_Vel = options.Admisible_max_error_Velocity
            self.Admisible_max_error_Extra_Coo = options.Admisible_max_error_Extra_Coo
   
    class beacon_options:
        def __init__(self):
            self.Smoothing = options.Smoothing
            self.Min_cluster_size = options.Min_cluster_size
            self.Maxima_ratio = options.Maxima_ratio
            self.Epsilon = options.Epsilon
            self.Standarization_Method = options.Standarization_Method
            self.Standarization_Weights = np.array(options.Standarization_Weights)
            self.Coordinate_system = options.Coordinate_system
            self.Velocity_system = options.Velocity_system
            self.Uniqueness_Method = options.Uniqueness_Method
            self.NStars_normed = options.NStars_normed
     
    class sref_options:
        def __init__(self):
            self.WCS = options.Galaxy_WCS
            self.Extra_Coo_is_FeH = options.Extra_Coo_is_FeH
            self.Use_elliptical_coordinates = options.Use_elliptical_coordinates
            self.Find_centre_sky_radius = options.Find_centre_sky_radius
            self.Find_centre_sky_resolution = options.Find_centre_sky_resolution
            self.Find_centre_vhel_radius = options.Find_centre_vhel_radius
            self.Find_centre_vhel_resolution = options.Find_centre_vhel_resolution

    class plotting_options:
        def __init__(self):
            self.Plot_maps = options.Plot_maps
            self.wcs = options.Galaxy_WCS
            self.Maps_plotting_radius = options.Maps_plotting_radius
            self.Maps_vlos_limits = options.Maps_vlos_limits
            self.Auxiliar_axes_limits = options.Auxiliar_axes_limits
            self.Draw_markers = options.Draw_markers
            self.Markers_color = options.Markers_color
            self.Markers_size = options.Markers_size
            self.Markers_linewidth = options.Markers_linewidth
            self.Plot_angular_momenta = options.Plot_angular_momenta
            self.Plot_extra_coo = options.Plot_extra_coo
            self.Plot_elliptical_radius = options.Use_elliptical_coordinates         
            if options.Auxiliar_axes_label:            	
               self.Auxiliar_axes_label = options.Auxiliar_axes_label

    return galaxy(), read_data_options(), beacon_options(), sref_options(), plotting_options()
  
def stars_generator(Rest_frame_coo):
    class stars_data(object):
        def __init__(self):
            self.Sky_Coo = Rest_frame_coo[0]
            self.Gref_car_Coo = Rest_frame_coo[1]
            self.Gref_car_Vel = Rest_frame_coo[2]
            self.Gref_car_eVel = Rest_frame_coo[3]
            self.Gref_pol_Coo = Rest_frame_coo[4]
            self.Gref_pol_Vel = Rest_frame_coo[5]    
            self.Gref_pol_eVel = Rest_frame_coo[6]
            try:
               self.Z = Rest_frame_coo[7]
               self.eZ = Rest_frame_coo[8]
               self.FeH = Rest_frame_coo[9]
               self.eFeH = Rest_frame_coo[10]
            except:
               pass
            try:
               self.Extra_Coo = Rest_frame_coo[7]
               self.eExtra_Coo = Rest_frame_coo[8]
            except:
               pass

    return stars_data

def read_data(read_data_options):
  
    Stars_data = np.loadtxt(read_data_options.Data_file, skiprows=1)

    try:
       Sky_Coo = Stars_data[:,0:2]
       Sky_Vel = Stars_data[:,2:7:2]
       Sky_eVel = Stars_data[:,3:8:2]
       Extra_Coo = Stars_data[:,-2:7:-2]
       Extra_eCoo = Stars_data[:,-1:7:-2]
    except:
       print 'Bad format in input table.\n'
       print 'Format shoud be:\n'
       print 'RA Dec VRA eVRA VDec eVDec Vlos eVlos Extra_Coo eExtra_Coo Extra_Coo_2 eExtra_Coo_2 ...'
       print 'Quitting now...'
       sys.exit(0)

    Sky_eVel = np.where(Sky_eVel < 1e-5, np.average(Sky_eVel), Sky_eVel)
    Extra_eCoo = np.where(Extra_eCoo < 1e-5, np.average(Extra_eCoo), Extra_eCoo)

    Std_stars_data = [Sky_Coo, Sky_Vel, Sky_eVel, Extra_Coo, Extra_eCoo]

    return Std_stars_data

def ref_frame(Std_data, galaxy, sref_options):

    """
    Sky, WCS and physical coordinates
    """
    if sref_options.WCS:
       Sky_Coo = Std_data[0]
       Gref_car_Coo = lib.wcs2xy(Sky_Coo, galaxy.Sky_Coo, galaxy.Distance)

    else:
       Gref_car_Coo = Std_data[0]
       Sky_Coo = lib.xy2wcs(Gref_car_Coo, galaxy.Sky_Coo, galaxy.Distance)

    """
    Galactic velocity rest frame
    """
    Gref_car_Vel = Std_data[1] - galaxy.Sky_Vel
    Gref_car_eVel = Std_data[2]

    if sref_options.Use_elliptical_coordinates:
       Gref_car_Coo = lib.deproject_ellip(Gref_car_Coo, Ellipticity = galaxy.E, PA = galaxy.PA)
       try:
          Gref_car_Vel[:,0:2] = lib.deproject_ellip(Gref_car_Vel[:,0:2], Ellipticity = galaxy.E, PA = galaxy.PA)
          Gref_car_eVel[:,0:2] = lib.deproject_ellip(Gref_car_eVel[:,0:2], Ellipticity = galaxy.E, PA = galaxy.PA)
       except:
          pass
    """
    Galactic elliptical coordinates
    """
    Gref_pol_Coo = lib.polar_coo(Gref_car_Coo)
    Gref_pol_Vel = lib.polar_vel(Gref_car_Coo, Gref_car_Vel)
    Gref_pol_eVel = np.abs(lib.polar_vel(Gref_car_Coo, Gref_car_eVel))
    
    print 'Right now beacon accepts either FeH and eFeH, or whatever linear coordinates user may provide, but NOT both (FeH and linear coordinates) at the same time!'
    try:
       Extra_Coo, Extra_eCoo = Std_data[3], Std_data[4]
       if sref_options.Extra_Coo_is_FeH:
          Z, eZ = lib.log2lin_metal(Extra_Coo, Extra_eCoo)
          FeH, eFeH = Extra_Coo, Extra_eCoo
          Rest_frame_coo = [Sky_Coo, Gref_car_Coo, Gref_car_Vel, Gref_car_eVel, Gref_pol_Coo, Gref_pol_Vel, Gref_pol_eVel, Z, eZ, FeH, eFeH]
       else:
          Rest_frame_coo = [Sky_Coo, Gref_car_Coo, Gref_car_Vel, Gref_car_eVel, Gref_pol_Coo, Gref_pol_Vel, Gref_pol_eVel, Extra_Coo, Extra_eCoo]
    except:
       Rest_frame_coo = [Sky_Coo, Gref_car_Coo, Gref_car_Vel, Gref_car_eVel, Gref_pol_Coo, Gref_pol_Vel, Gref_pol_eVel]

    Stars_data_rf = stars_generator(Rest_frame_coo)
    
    return Stars_data_rf()

def rest_frame_generator(galaxy, read_data_options, sref_options):
            
    if (sref_options.Find_centre_sky_radius != 0) and (sref_options.Find_centre_vhel_radius !=0):
        Find_centre_coordinates_list = np.mgrid[galaxy.Sky_Coo[0]-sref_options.Find_centre_sky_radius/60. : galaxy.Sky_Coo[0]+sref_options.Find_centre_sky_radius/60. :  sref_options.Find_centre_sky_resolution/60., galaxy.Sky_Coo[1]-sref_options.Find_centre_sky_radius/60. :  galaxy.Sky_Coo[1]+sref_options.Find_centre_sky_radius/60. : sref_options.Find_centre_sky_resolution/60., galaxy.Sky_Vel-sref_options.Find_centre_vhel_radius : galaxy.Sky_Vel+sref_options.Find_centre_vhel_radius: sref_options.Find_centre_vhel_resolution]
        
        Find_centre_RA_sky_coordinates_List = Find_centre_coordinates_list[0].flatten()
        Find_centre_Dec_sky_coordinates_List = Find_centre_coordinates_list[1].flatten()
        Find_centre_Vhel_sky_coordinates_List = Find_centre_coordinates_list[2].flatten()
        Dimensions = np.shape(Find_centre_coordinates_list[0])
        
    elif (sref_options.Find_centre_sky_radius != 0) and (sref_options.Find_centre_vhel_radius == 0):
        Sky_coo_list = np.mgrid[galaxy.Sky_Coo[0]-sref_options.Find_centre_sky_radius/60. : galaxy.Sky_Coo[0]+sref_options.Find_centre_sky_radius/60. :  sref_options.Find_centre_sky_resolution/60., galaxy.Sky_Coo[1]-sref_options.Find_centre_sky_radius/60. :  galaxy.Sky_Coo[1]+sref_options.Find_centre_sky_radius/60. : sref_options.Find_centre_sky_resolution/60.]
        
        Find_centre_RA_sky_coordinates_List = Sky_coo_list[0].flatten()
        Find_centre_Dec_sky_coordinates_List = Sky_coo_list[1].flatten()
        Find_centre_Vhel_sky_coordinates_List = np.ones(len(Sky_coo_list[0].flatten()))*galaxy.Sky_Vel
        Dimensions = np.shape(Sky_coo_list[0])
    
    elif (sref_options.Find_centre_sky_radius == 0) and (sref_options.Find_centre_vhel_radius != 0):        
        Find_centre_RA_sky_coordinates_List = np.ones(len(Vhel_coo_list))*galaxy.Sky_Coo[0]
        Find_centre_Dec_sky_coordinates_List = np.ones(len(Vhel_coo_list))*galaxy.Sky_Coo[1]
        Find_centre_Vhel_sky_coordinates_List = np.arange(galaxy.Sky_Vel-sref_options.Find_centre_vhel_radius, galaxy.Sky_Vel+sref_options.Find_centre_vhel_radius, sref_options.Find_centre_vhel_resolution)
        Dimensions = 1

    else:
        Find_centre_RA_sky_coordinates_List = [galaxy.Sky_Coo[0]]
        Find_centre_Dec_sky_coordinates_List = [galaxy.Sky_Coo[1]]
        Find_centre_Vhel_sky_coordinates_List = [galaxy.Sky_Vel]
        Dimensions = 0
    
    Std_data = read_data(read_data_options)
    
    Stars_data_lis = []
    
    for Find_centre_RA_sky_coordinates, Find_centre_Dec_sky_coordinates, Find_centre_Vhel_coordinates in zip(Find_centre_RA_sky_coordinates_List, Find_centre_Dec_sky_coordinates_List, Find_centre_Vhel_sky_coordinates_List):
             
       galaxy.Sky_Coo = np.array([Find_centre_RA_sky_coordinates, Find_centre_Dec_sky_coordinates])
       galaxy.Sky_Vel = np.array([Find_centre_Vhel_coordinates])
                     
       Stars_data_lis.append(ref_frame(Std_data, galaxy, sref_options))
    
    return Stars_data_lis, Find_centre_RA_sky_coordinates_List, Find_centre_Dec_sky_coordinates_List, Find_centre_Vhel_sky_coordinates_List, Dimensions

def autoorder(Possible_Circular_Stream_Indices, Possible_NonCircular_Stream_Indices, stars_data):

    try:
       OV = stars_data.Z
    except:
       OV = stars_data.Extra_Coo
       
    Possible_Circular_Stream_Mean_OV = []
    for Circ_Index in Possible_Circular_Stream_Indices:
        Possible_Circular_Stream_Mean_OV.append(np.mean(OV[list(Circ_Index)]))

    Possible_NonCircular_Stream_Mean_OV = []
    for NonCirc_Index in Possible_NonCircular_Stream_Indices:
        Possible_NonCircular_Stream_Mean_OV.append(np.mean(OV[list(NonCirc_Index)]))        
        
    """
    Order Results By OV
    """
    Possible_Circular_Stream_Sorted_Indices = sorted(range(len(Possible_Circular_Stream_Mean_OV)), key=lambda k: Possible_Circular_Stream_Mean_OV[k])
    Possible_NonCircular_Stream_Sorted_Indices = sorted(range(len(Possible_NonCircular_Stream_Mean_OV)), key=lambda k: Possible_NonCircular_Stream_Mean_OV[k])

    Possible_Circular_Stream_Indices = [Possible_Circular_Stream_Indices[i] for i in Possible_Circular_Stream_Sorted_Indices]

    Possible_NonCircular_Stream_Indices = [Possible_NonCircular_Stream_Indices[i] for i in Possible_NonCircular_Stream_Sorted_Indices]

    return Possible_Circular_Stream_Indices, Possible_NonCircular_Stream_Indices
    
def beacon(Rest_frame_coo, beacon_options):

    """
    This routine will call OPTICS with an apropiate set of coordinates. 
    Tipically coordinates, velocities, and chemical abundances
    """

    if (beacon_options.Coordinate_system == 'pol'):
       """
       Coo are the coordinates of the stars. The angle (Theta), must be inverted, but not the radius
       """
       Coo = Rest_frame_coo.Gref_pol_Coo
       Simmet_Coo = np.array([Rest_frame_coo.Gref_pol_Coo[:,0], (Rest_frame_coo.Gref_pol_Coo[:,1]+np.pi) % (2.*np.pi)]).T
       
    if (beacon_options.Velocity_system == 'pol'):
       """
       Polar velocities are expected to be the same at both sides of
       the CM, except for the line-of-sight (z) component, that shoud be inverted.
       """
       Vel = Rest_frame_coo.Gref_pol_Vel
       Simmet_Vel = np.array([Vel[:,0], Vel[:,1], -1.*Vel[:,2]]).T

    if (beacon_options.Coordinate_system == 'car'):
       """
       If one prefers cartesian coordinates, then both have to be inverted
       """
       Coo = Rest_frame_coo.Gref_car_Coo
       Simmet_Coo = np.array([(-Rest_frame_coo.Gref_car_Coo[:,0]), -Rest_frame_coo.Gref_car_Coo[:,1]]).T

    if (beacon_options.Velocity_system == 'car'):
       """
       All cartesian components of the velocity are expected to be the
       oposite at the other side of the CM.
       """
       Vel = Rest_frame_coo.Gref_car_Vel
       Simmet_Vel = -1.*Vel    

    """
    Extra properties should be consistently equal at both sides of the CM (In case they exist)
    """
    try:
       Extra = Rest_frame_coo.Z
    except:
       pass
    try:
       Extra = Rest_frame_coo.Extra_Coo
    except:
       pass
    try:
       if Extra.ndim == 1:
          Extra = np.expand_dims(Extra, axis=1)

       Direct_Clustering_Data = np.hstack([Coo, Vel, Extra])
       Simmet_Clustering_Data = np.hstack([Simmet_Coo, Simmet_Vel, Extra])
    except:
       Direct_Clustering_Data = np.hstack([Coo, Vel])
       Simmet_Clustering_Data = np.hstack([Simmet_Coo, Simmet_Vel])

    OAL = len(Direct_Clustering_Data)   # Original Array Length

    """
    If the only coordinates being used are v_rot and v_rad then don't use the simmetric vector
    """
    if (beacon_options.Velocity_system == 'pol') & (beacon_options.Standarization_Weights[4] == 0):
      Clustering_Data = Direct_Clustering_Data
      mcs = beacon_options.Min_cluster_size
    else:
      Clustering_Data = np.vstack([Direct_Clustering_Data, Simmet_Clustering_Data])
      mcs = float(beacon_options.Min_cluster_size)*2.

    if not Clustering_Data.size:
      print 'No stars remain after cleaning, check option file!'
      quit()
    """
    Standardization
    """
    if beacon_options.Standarization_Method == 'ptp':
        Clustering_Data /= np.ptp(Direct_Clustering_Data, axis=0)
    elif beacon_options.Standarization_Method == 'std':
        Clustering_Data /= np.std(Direct_Clustering_Data, axis=0)
    elif beacon_options.Standarization_Method == 'var':
        Clustering_Data /= np.var(Direct_Clustering_Data, axis=0)

    else:
        print "Non standarization was applied!"
    
    try:
        Clustering_Data *= beacon_options.Standarization_Weights
    except:
        print "Incorrect dimension of weights. Same weight applied to all coordinates."
    
    Clustering_Data = Clustering_Data[:,np.where(beacon_options.Standarization_Weights != 0)[0]]

    """
    OPTICS call
    """
    Indices = optics_launcher(Clustering_Data, beacon_options.Smoothing, min_cluster_size = mcs, maxima_ratio = beacon_options.Maxima_ratio)
    
    Possible_Circular_Stream_Indices = []
    Possible_NonCircular_Stream_Indices = []
    for ii, Index in enumerate(Indices):
        if (sum(1 for i in Index if i >= OAL) > len(Index)/beacon_options.Epsilon) & (sum(1 for i in Index if i < OAL) > len(Index)/beacon_options.Epsilon):
            Possible_Circular_Stream = True
            Possible_Nonirc_streams = False
        elif ((Index < OAL).all()) or ((Index >= OAL).all()):
            Possible_Circular_Stream = False
            Possible_Nonirc_streams = True
        else:
            Possible_Circular_Stream = False
            Possible_Nonirc_streams = False  
        Index[Index >= OAL] -= OAL
        if Possible_Circular_Stream:
            Possible_Circular_Stream_Indices.append(set(Index))
        elif Possible_Nonirc_streams:
            Possible_NonCircular_Stream_Indices.append(set(Index))
          
    """
    Uniqueness of the groups
    """
    if beacon_options.Uniqueness_Method == 'all':    
        Possible_Circular_Stream_Indices = [list(i) for i in set(tuple(i) for i in Possible_Circular_Stream_Indices)]
        Possible_NonCircular_Stream_Indices = [list(i) for i in set(tuple(i) for i in Possible_NonCircular_Stream_Indices)]
        
    if beacon_options.Uniqueness_Method == 'any':
        from networkx.algorithms.components.connected import connected_components
        Possible_Circular_Stream_Indices = list(connected_components(lib.to_graph(Possible_Circular_Stream_Indices)))
        Possible_NonCircular_Stream_Indices = list(connected_components(lib.to_graph(Possible_NonCircular_Stream_Indices)))
    
    """
    If there are many groups, exit program
    """
    if (len(Possible_Circular_Stream_Indices) > 1000) or (len(Possible_NonCircular_Stream_Indices) > 1000):
        print 'To many groups: ',len(Possible_Circular_Stream_Indices),' aborting to prevent computer failure...'
        sys.exit(0)

    return Possible_Circular_Stream_Indices, Possible_NonCircular_Stream_Indices   


def useful_quantities(stars_data, galaxy, sref_options, Indices_list=None):
    """
    This routine calculates some useful quantities for plotting and arranging results
    """
    if Indices_list is None:
        Indices_list = [list(np.arange(len(stars_data.Gref_car_Coo)))]

    """
    Averages must be derived using linear quantities
    """
    N_Stars = []
    Velocities_stars = []
    Extra_Coo_stars = []
    Radius_stars = []
    
    Velocity_groups = []
    Extra_Coo_groups = []
    Radius_groups = []
    Angular_Momenta_groups = []

    print 'WARNING: Projected angular momentum (line-of-sight velocities) with no error!'

    for ii, Indices in enumerate(Indices_list):
        Indices = list(Indices)
        """
        Angular Momentum
        """
        L_k = lib.kinetic_pa(stars_data.Gref_car_Coo[Indices, :], stars_data.Gref_car_Vel[Indices, 2], 1./stars_data.Gref_car_eVel[Indices, 2])

        kl = lib.kinetic_pa(stars_data.Gref_car_Coo[Indices, :], stars_data.Gref_car_Vel[Indices, 2], 1./stars_data.Gref_car_eVel[Indices, 2])
        if sref_options.Use_elliptical_coordinates:
           kl[0:2] =  lib.deproject_ellip(kl[0:2], Ellipticity = 1./(1.-1./galaxy.E), PA = galaxy.PA)

        Angular_Momenta_groups.append(kl)

        try:
           Extra_Coo_stars.append([stars_data.FeH[Indices], stars_data.eFeH[Indices]])

           Avg_Std_Z = lib.avg_std(stars_data.Z[Indices], 1./(stars_data.eZ[Indices])**2)

           Extra_Coo_groups.append(lib.lin2log_metal(Avg_Std_Z[0], Avg_Std_Z[1]))
        except:
           pass
        try:
           Extra_Coo_stars.append([stars_data.Extra_Coo[Indices], stars_data.Extra_eCoo[Indices]])
           
           Extra_Coo_groups.append(lib.avg_std(stars_data.Extra_Coo[Indices], 1./(stars_data.Extra_eCoo[Indices])**2))
        except:
           pass


        N_Stars.append(len(Indices))
        Velocities_stars.append([stars_data.Gref_car_Vel[Indices],stars_data.Gref_car_eVel[Indices]])
        
        Radius_stars.append(np.sqrt(np.sum(stars_data.Gref_car_Coo[Indices, :]**2, axis = 1)))
        
        Velocity_groups.append(lib.avg_std((stars_data.Gref_car_Coo[Indices]), 1./(stars_data.Gref_car_Coo[Indices])**2))

        Radius_groups.append(np.average(np.sqrt(np.sum(stars_data.Gref_car_Coo[Indices, :]**2, axis = 1))))

    """
    CONVERTING LISTS TO ARRAYS
    """
    N_Stars = np.array(N_Stars)
    Angular_Momenta_groups = np.array(Angular_Momenta_groups)

    Radius_stars = np.array(Radius_stars)
    Velocity_groups = np.array(Velocity_groups)

    """
    These quantities may be optional
    """
    Extra_Coo_stars = np.array(Extra_Coo_stars)    
    Extra_Coo_groups =  np.array(Extra_Coo_groups)
    

    Mean_quantities = [Extra_Coo_groups, Velocity_groups, Radius_groups, Angular_Momenta_groups]
    Raw_quantities = [N_Stars, Extra_Coo_stars, Velocities_stars, Radius_stars]

    return Raw_quantities, Mean_quantities

def save_clusters(Circ_mean_quantities, File_name = 'Output/Clusters.dat'):
    Circ_Mean_Metallicity, Circ_Mean_Velocity, Circ_Mean_rDistance, Circ_KL = Circ_mean_quantities

    np.savetxt(File_name, np.array([np.arange(len(Circ_Mean_Metallicity))+1, Circ_Mean_Metallicity[:,0], Circ_Mean_Metallicity[:,1], Circ_KL[:,0], Circ_KL[:,1]]).T, fmt='%i  %.3f  %.3f  %.3f  %.3f', header='Pop   Metal   eMetal   L_angle  L_modulus')

def save_clustered_stars(stars_data, Stream_Indices, File_name = 'Output/Stars.dat'):
   Circ_clustered_stars = []
   for kk, Indices  in enumerate(Stream_Indices):
      Indices = list(Indices)

      Cluster = np.array([np.array(Indices)+1., stars_data.Sky_Coo[Indices,0], stars_data.Sky_Coo[Indices,1], stars_data.Gref_car_Vel[Indices,2], stars_data.Gref_car_eVel[Indices,2], np.squeeze(stars_data.FeH[Indices]), np.squeeze(stars_data.eFeH[Indices]), np.ones(len(stars_data.FeH[Indices]))*(kk+1)]).T

      Circ_clustered_stars.append(Cluster)
   Circ_clustered_stars = np.array([item for sublist in Circ_clustered_stars for item in sublist])

   np.savetxt(File_name, Circ_clustered_stars, fmt='%i  %.8f  %.8f  %.3f  %.3f  %.3f  %.3f  %i', header='Indices    RA    DEC    V_los    eV_los    Metal    eMetal    Pop')


def optics_launcher(Data, Smoothing, min_cluster_size = 5., maxima_ratio = 0.75):
    import optics_core as core

    RD, CD, Order = core.optics(Data,Smoothing)

    Reach_Plot = []
    Reach_Points = []

    for item in Order:
        Reach_Plot.append(RD[item])
        Reach_Points.append(Data[item])

    min_cluster_size_ratio = min_cluster_size / float(len(Data))

    rootNode = core.auto_cluster(Reach_Plot, Reach_Points, min_cluster_size_ratio, maxima_ratio)

    #get only the leaves of the tree
    leaves = core.get_leaves(rootNode, [])

    index_ordered = []
    for item in leaves:
        index_ordered.append(np.array(Order[item.start:item.end]))
    return index_ordered


"""
Andres del Pino Molina
"""
