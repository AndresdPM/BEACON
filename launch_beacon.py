"""
This script will read and option file and launch BEACON
"""

import sys, os, warnings

def beacon_launcher(Options_file):
    import beacon as bc
    import plotting as pl
    warnings.filterwarnings("ignore")
  
    print 'Processing: %s \n'%Options_file

    print 'Reading options file...\n'
    galaxy, read_data_options, beacon_options, sref_options, plotting_options = bc.parameters_generator(Filename=Options_file)
    
    print 'Reading data file...\n'
    Stars_data_lis, Centre_RA_sky_coordinates_List, Centre_Dec_sky_coordinates_List, Centre_Vhel_sky_coordinates_List, Centre_coordinates_dimensions = bc.rest_frame_generator(galaxy, read_data_options, sref_options)

    """
    Loop over all Systems of References beging
    """
    for ii, stars in enumerate(Stars_data_lis):
        print ii, 'of', len(Stars_data_lis)

        Original_n_stars = len(stars.Sky_Coo)

        if read_data_options.Clean_sample:
            print 'Pre cleanning sample...\n'
            stars = bc.clean_data(stars, read_data_options, optics_cleanning = False, sigma_clipping = 3.)

        Clean_n_stars = len(stars.Sky_Coo)

        print 'Launching BEACON...\n' 
        Circular_Stream_Indices, NonCircular_Stream_Indices = bc.beacon(stars, beacon_options)

        print 'Arranging results...\n' 
        Circular_Stream_Indices, NonCircular_Stream_Indices = bc.autoorder(Circular_Stream_Indices, NonCircular_Stream_Indices, stars)

        CS_n_stars = len([item for sublist in Circular_Stream_Indices for item in sublist])
        NCS_n_stars = len([item for sublist in NonCircular_Stream_Indices for item in sublist])

        print '--------------------'
        print 'GENERAL STATISTICS'
        print '--------------------'
        print 'Number of input stars:',  Original_n_stars
        print 'Number of cleaned sample:', Clean_n_stars
        print 'Number of circular streams:', len(Circular_Stream_Indices), 'with a total of', CS_n_stars, 'stars.'
        print 'Number of non circular streams:', len(NonCircular_Stream_Indices), 'with a total of', NCS_n_stars, 'stars.'
        print '--------------------\n'

        """
        CREATING STREAMS DATA ARRAYS
        """
        Raw_quantities, Raw_mean_quantities = bc.useful_quantities(stars, galaxy, sref_options)
        Raw_Mean_Metallicity_linear, Raw_Mean_Velocity, Raw_Mean_rDistance, Raw_KL = Raw_mean_quantities
       
        if Circular_Stream_Indices:
          
            Circ_raw_quantities, Circ_mean_quantities = bc.useful_quantities(stars, galaxy, sref_options, Indices_list = Circular_Stream_Indices)
        
            Circ_Mean_Metallicity, Circ_Mean_Velocity, Circ_Mean_rDistance, Circ_KL = Circ_mean_quantities
            Circ_N_Stars, Circ_Metallicities_linear, Circ_Velocities, Circ_radius = Circ_raw_quantities
            
        if NonCircular_Stream_Indices:
          
            NonCirc_raw_quantities, NonCirc_mean_quantities = bc.useful_quantities(stars, galaxy, sref_options, Indices_list = NonCircular_Stream_Indices)

            NonCirc_Mean_Metallicity, NonCirc_Mean_Velocity, NonCirc_Mean_rDistance, NonCirc_KL = NonCirc_mean_quantities
            NonCirc_N_Stars, NonCirc_Metallicities_linear, NonCirc_Velocities, NonCirc_radius = NonCirc_raw_quantities

        ##########
        #PLOTTING#
        ##########
    
        """
        Raw Map
        """
        if plotting_options.Plot_maps:

            pl.plot_data_maps(galaxy, stars, plotting_options, File_name = 'Output/'+galaxy.Name+'_Data.png')
            
            RAW_histo = pl.plot_streams(galaxy, stars, plotting_options, Raw_KL, File_name = 'Output/'+galaxy.Name+'_Raw_maps.png')
            if Circular_Stream_Indices:
                pl.plot_streams(galaxy, stars, plotting_options, Circ_KL, RAW_histogram = RAW_histo, Indices_list = Circular_Stream_Indices, File_name = 'Output/'+galaxy.Name+'_Circular_Streams.png')
        
            if NonCircular_Stream_Indices:
                pl.plot_streams(galaxy, stars, plotting_options, NonCirc_KL, Indices_list = NonCircular_Stream_Indices, File_name = 'Output/'+galaxy.Name+'_NonCircular_Streams.png')

        """
        Angular Momentum
        """
        if plotting_options.Plot_angular_momenta:
            if Circular_Stream_Indices:
                pl.plot_angular_momenta(galaxy, Circ_KL, Circ_Mean_Metallicity, plotting_options, File_name = 'Output/'+galaxy.Name+'_Momenta.png')

        """
        MAIN TABLE WITH DATA ABOUT GROUPS
        """
        bc.save_clusters(Circ_mean_quantities, File_name = 'Output/'+galaxy.Name+'_cs.dat')

        """
       MAIN TABLE WITH DATA ABOUT GROUPED STARS
        """
        bc.save_clustered_stars(stars, Circular_Stream_Indices, File_name = 'Output/'+galaxy.Name+'_Stars_cs.dat')

if __name__=="__main__":
  
    print '--------'
    print 'Starting '
    print '--------\n'

    Options_file = sys.argv[1]
    beacon_launcher(Options_file)
   
    print '-----'
    print 'Done! '
    print '-----\n'

"""
Andres del Pino Molina
"""
