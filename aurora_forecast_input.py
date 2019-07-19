#CONFIGURATION PARAMETERS for aurora_forecast.py 

mode=0                        # select mode: 0 for real time mode, 1 for local file, 2 for OMNI2 data

time_resolution = 60         # time resolution of resulting auroramaps in minutes

output_directory='run_boundary_test'            #specify output directory of frames and movies under "results/"



#------------------------------- controls for computation

window_minutes=20                #window in minutes for smoothing the coupling with a running mean; ignored if time_resolution larger than window_minutes; standard=20

calc_mode='multi'               #multi or single processing mode for calculating the aurora image cube
#calc_mode='single'


# --------------------------- mode 0 settings

                                   # in real time mode, start time is always now in UTC
past_hours      =  0               # in real time mode, start time with previous hours, negative = past
future_hours    =  24               # in real time mode, number of hours into future, 0 for one frame only

#online source file for real time mode
predstorm_url='https://www.iwf.oeaw.ac.at/fileadmin//staff/SP/cmoestl/readtime/predstorm_real.txt'


# --------------------------  mode 1/2 settings

#start_time = '2019-May-14 05:00'  # in mode 1/2, set start and end times manually
#end_time   = '2019-May-14 10:00'

start_time = '2017-Sep-1 13:00' 
end_time   = '2017-Sep-10 14:00'


# --------------------------  mode 2 settings
local_input_file='data/predstorm/predstorm_v1_realtime_stereo_a_save_2019-05-14-21_00.txt'



#local_input_file_1hour='predstorm_output/predstorm_v1_realtime_stereo_a_save_2019-05-14-21_00.txt'
#run predstorm_l5.py for producing a real time file first
#local_input_file_1hour='/Users/chris/python/predstorm/predstorm_real.txt'
#local_input_file_1min='/Users/chris/python/predstorm/predstorm_real_1m.txt'
#local_input_file_1min='predstorm_output/predstorm_real_1m_with_May14_2019_storm.txt'




