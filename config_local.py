'''
config_local_web.py

part of the "auroramaps" package

CONFIGURATION PARAMETERS for main program aurora_web.py 

---------------------
by C. Moestl, Austrian Space Weather Office, GeoSphere Austria.
https://github.com/helioforecast/auroramaps
https://helioforecast.space

published under GNU Lesser General Public License v3.0

last update April 2025
'''




mode=0                     # select mode: 0 for real time wind from URL, 1 for local file, 2 for OMNI2 data

time_resolution = 60        # time resolution of resulting auroramaps in minutes

frame_rate=1           #output movie frame rate frame rate 20 is good for 10 minute resolution if 3 days want to be seen quickly

# --------------------------- mode 0 settings
                                  # in real time mode, start time is always now in UTC
past_hours      =  -24          # in real time mode, start time with previous hours, negative = past
future_hours    =  -1              # in real time mode, number of hours into future, 0 for one frame only


output_directory='aurora_test'            #specify output directory of frames and movies under "results/"

#----------------------------- select map types 

map_type = 'marble'   #marble, viirs or topography

#map_type = 'topography'
#map_type = 'viirs'

#1.0 erg /cm^2 /s is the threshold in flux known that approximates visible aurora starting
equatorial_boundary_flux_threshold=1.0

# set 1 for making the map, 0 for not making the map

#flux maps
global_flux_map=0         #northern polar view
europe_flux_map=1       #Europe
canada_flux_map=0      #North America
greenland_flux_map=1       #Greenland



#probability maps
global_probability_map=0
europe_probability_map=0
canada_probability_map=0
greenland_probability_map=0       #Greenland



#------------------------------- controls for computation

window_minutes=20                #window in minutes for smoothing the coupling with a running mean; ignored if time_resolution larger than window_minutes; standard=20

#calc_mode='multi'               #'multi' or 'single' processing mode for calculating the aurora image cube
#calc_mode_frame='multi'         #'multi' or 'single' processing mode for drawing and saving the aurora frames

calc_mode='single'               #'multi' or 'single' processing mode for calculating the aurora image cube
calc_mode_frame='single'         #'multi' or 'single' processing mode for drawing and saving the aurora frames


#online source file for real time mode
predstorm_url='https://helioforecast.space/static/sync/predstorm_real_1m.txt'

# --------------------------  mode 1/2 settings

#start_time = '2019-May-14 05:00'  # in mode 1/2, set start and end times manually
#end_time   = '2019-May-14 10:00'

start_time = '2024-May-10 22:00' 
end_time   = '2024-May-12 06:00'


# --------------------------  mode 1 settings
local_input_file='data/predstorm/predstorm_real_test_2023_jul_22.txt'


# --------------------------  mode 2 settings OMNI data is automatically loaded 



#local_input_file_1hour='predstorm_output/predstorm_v1_realtime_stereo_a_save_2019-05-14-21_00.txt'
#run predstorm_l5.py for producing a real time file first
#local_input_file_1hour='/Users/chris/python/predstorm/predstorm_real.txt'
#local_input_file_1min='/Users/chris/python/predstorm/predstorm_real_1m.txt'
#local_input_file_1min='predstorm_output/predstorm_real_1m_with_May14_2019_storm.txt'




