'''
input.py

part of the "auroramaps" package

CONFIGURATION PARAMETERS for main program aurora.py 

---------------------
by C. Moestl, IWF-helio group, Graz, Austria.
https://github.com/IWF-helio/auroramaps
twitter @chrisoutofspace
https://www.iwf.oeaw.ac.at/user-site/christian-moestl/

published under GNU Lesser General Public License v3.0

last update October 2019
'''


mode=0                      # select mode: 0 for real time mode, 1 for local file, 2 for OMNI2 data

time_resolution = 30         # time resolution of resulting auroramaps in minutes

output_directory='aurora_test_package_5'            #specify output directory of frames and movies under "results/"

#--------------------------------- select map types 

#flux maps
global_flux_map=False                   #northern polar view
europe_canada_flux_map=False           #2 maps in one frame for Europa and Canada/USA


#Maximum level for flux plots erg cm-2 -s-1
max_level_flux=1.5


#probability maps
global_probability_map=False
europe_probability_map=False
canada_probability_map=True
europe_canada_probability_map=False

#marble, viirs or topography
map_type = 'marble'#'viirs' #'marble'#'topography'


#valid for flux and probability maps
equatorial_boundary_flux_threshold=1.0

#------------------------------- controls for computation

window_minutes=30                #window in minutes for smoothing the coupling with a running mean; ignored if time_resolution larger than window_minutes; standard=20

calc_mode='multi'               #multi or single processing mode for calculating the aurora image cube
#calc_mode='single'


# --------------------------- mode 0 settings

                                   # in real time mode, start time is always now in UTC
past_hours      =  -15             # in real time mode, start time with previous hours, negative = past
future_hours    =  1               # in real time mode, number of hours into future, 0 for one frame only

#online source file for real time mode
predstorm_url='https://www.iwf.oeaw.ac.at/fileadmin//staff/SP/cmoestl/readtime/predstorm_real.txt'



# --------------------------  mode 1/2 settings

#start_time = '2019-May-14 05:00'  # in mode 1/2, set start and end times manually
#end_time   = '2019-May-14 10:00'

start_time = '2017-Sep-7 13:00' 
end_time   = '2017-Sep-7 14:00'



# --------------------------  mode 1 settings
local_input_file='data/predstorm/predstorm_v1_realtime_stereo_a_save_2019-05-14-21_00.txt'


# --------------------------  mode 2 settings OMNI data is automatically loaded 







#local_input_file_1hour='predstorm_output/predstorm_v1_realtime_stereo_a_save_2019-05-14-21_00.txt'
#run predstorm_l5.py for producing a real time file first
#local_input_file_1hour='/Users/chris/python/predstorm/predstorm_real.txt'
#local_input_file_1min='/Users/chris/python/predstorm/predstorm_real_1m.txt'
#local_input_file_1min='predstorm_output/predstorm_real_1m_with_May14_2019_storm.txt'




