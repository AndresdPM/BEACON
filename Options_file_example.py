"""
Input file for BEACON
"""

Galaxy_Name = "N-Body"                                   #Name of the stellar system. Will be used for naming output files.
Galaxy_Data_file = "./Catalogs/N-Body.dat"               #Path to the input catalogue.
Galaxy_WCS = False                                       #True: RA, Dec. False: x, y (One may chose false when working with N-body)

Galaxy_Coordinates = [0., 0.]                            #Center of Mass coordinates in the sky.
Galaxy_Velocity = 0.                                     #Galaxy Velocity (Galaxy_vel).
Galaxy_Distance = 136000                                 #Distance to the galaxy in parsec, it may be ignored when working with N-body
Galaxy_rh = '--'                                         #Galaxy half-light radius, just for plotting purposes. Not required.
Galaxy_rt = '--'                                         #Galaxy tidal radius, just for plotting purposes. Not required.
Galaxy_e = '--'                                          #Galaxy half-light radius, just for plotting purposes. Not required.
Galaxy_pa = '--'                                         #Galaxy half-light radius, just for plotting purposes. Not required.

Plot_maps = True                                         #Plotting options. Plot velocity maps.
Maps_plotting_radius = 2.                                #Plotting options. Defines the used limits to plot coordinates in the velocity maps.
Maps_vlos_limits = [-25, 25]                             #Plotting options. Line-of-sight velocity limits in the velocity maps. 
Draw_markers = True                                      #Plotting options. Draw markers in the velocity maps.
Markers_color = (0.1,0.1,0.1)                            #Plotting options. Color (R,G,B) used for markers in the velocity maps.
Markers_size = 5                                         #Plotting options. Used size for the markers in the velocity maps.
Markers_linewidth = 1.5                                  #Plotting options. Markers line width in the velocity maps.

Plot_angular_momenta = True                              #Plotting options. Plot angular momenta of the clusters.
Plot_extra_coo = False                                   #Plotting options. Plot Extra Coo instead of metallicity.
Auxiliar_axes_limits = [-3., 0.]                         #Plotting options. Histogram limits.
Auxiliar_axes_label = r"[Fe/H]"                          #Plotting options. Histogram axis label.

Extra_Coo_is_FeH = True                                  #Extra coordinate is [Fe/H]. This will make BEACON to convert [Fe/H] to linear Z using Solar metallicity = 0.019.
Use_elliptical_coordinates = False                       #Use ellipticity-deprojected coordinates instead of sky-projected ones.                    
Coordinate_system = "pol"                                #Clustering coordinate system for positions. "pol" = Polar, "car" = Cartesian
Velocity_system = "pol"                                  #Clustering coordinate system for velocities. "pol" = Polar, "car" = Cartesian

Clean_sample = False                                     #This must be set to False for the moment.
Admisible_range_Velocity = [-50, 50]                     #Not working for the moment.
Admisible_max_error_Velocity = 5                         #Not working for the moment.
Admisible_max_error_Extra_Coo = 0.01                     #Not working for the moment.

Smoothing = 0                                            #BEACON options: Smothing performed to the data. 0 = No smothing. More somothing help to detect more general rotation patterns.
Min_cluster_size = 9.                                    #BEACON options: Minimum Cluster Size (MCS). Minimum number of stars that can be considered as a cluster.
Maxima_ratio = 0.7                                       #BEACON options: Reachability distance ratio. 0 to 1. Values between 0.65 and 0.85 provide best results. The larger the maxima ratio is the more sensitivity BEACON has, making it prone to spureous detection
Epsilon = 4.                                             #BEACON options: 1/Epsilon is the minimum amout of stars that must be found at both sides of the center of rotation as to consider the cluster as "circular" or "both side stream"

Standarization_Method = "std"                            #BEACON options: Method used for standarization of the state vector. "std" = Standard deviation, "ptp" = dynamic range, "var" = variance
Uniqueness_Method = "any"                                #BEACON options: Uniquenes of solution. "any" = mergers clusters sharing stars into one. "all" = clusters are different if they differ in at least one star (stars could be common to more than one cluster)
Standarization_Weights = [1., 1., 1., 1., 0., 0.]        #BEACON options: Standarization weights for used coordinates. 0 - inf. Format = [Coo_1, Coo_2, v_1, v_2, v_3, Extra_coo_1, Extra_Coo_2, ...]. Coo_1, Coo_2 = r, theta. v_1, v_2, v_3 = v_r, v_theta, v_z
NStars_normed = False                                    #BEACON options: Normalize derived quantities to # stars in clusters. NOT IN USE.

Find_centre_sky_radius = 0.                              #BEACON options: Useful when the position of the center of mass of the galaxy is not known. It computes a data cube with results around the galaxy coordinates with radius "Find_centre_sky_radius" in arcmin. Use 0 when coordinates are known.
Find_centre_sky_resolution = 0.                          #BEACON options: Useful when the position of the center of mass of the galaxy is not known. Resolution of the data cube in arcmin. 
Find_centre_vhel_radius = 0                              #BEACON options: Useful when the velocity of the center of mass of the galaxy is not known. It computes a data cube with results around the galaxy velocity with radius "Find_centre_vhel_radius" in km/s. Use 0 when velocity is known.
Find_centre_vhel_resolution = 0.                         #BEACON options: Useful when the velocity of the center of mass of the ggalaxy is not known. Resolution of the data cube in km/s. 
