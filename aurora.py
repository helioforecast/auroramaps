'''
aurora.py

main program for running the OVATION PRIME 2010 model to make an aurora forecast/hindcast 
based on the PREDSTORM solar wind prediction method (available online)
or OMNI2 data for historic events

uses the python "auroramaps" package
https://github.com/IWF-helio/auroramaps

input parameters are given in the input.py file

by C. Moestl, Rachel L. Bailey, A. J. Weiss, Helio4Cast group, Graz, Austria.
Contributions by Liam Kilcommons and Diana Morosan

twitter @chrisoutofspace
https://www.iwf.oeaw.ac.at/user-site/christian-moestl/

This package uses a rewritten version of the ovationpyme aurora model 
by Liam Kilcommons https://github.com/lkilcommons/OvationPyme

published under GNU Lesser General Public License v3.0

last update October 2019

-----------------------------------------------------------------------------------
open issues: 

core:
- check with IDL code, debuggen, both hemispheres correct?
- check with direct comparison with NOAA global images nowcast
- get flux for time -> weights for both seasons? how correctly?, compare further with ovationpyme and IDL version
- add wave flux (needs good interpolation like OP13)
- add southern hemisphere maps (with coordinate conversion etc.)
- use numba or multiprocessing somewhere further for speed? 
- multiprocessing for saving the frames does not work on MacOS! 
  maybe on linux simply use the plotting function with multiprocessing.pool
  and multiple arguments (starmap) like for the data cubes

test bottlenecks: 
>> python -m cProfile -s tottime aurora_forecast.py

or use in ipython for function run times:
>> %timeit function_name()
>> %time  function_name()


plotting ideas:
- transparent to white colormap so that it looks like viirs images for direct comparison
- split land on dayside / night lights on night side 
  this should work in global_predstorm_north by loading background only once
  need to figure out how to get pixels from nightshade day/night and how to plot 
  only specific pixels (but then each background image must be updated)
- indicate moon phase with astropy
- cloud cover for local location? https://pypi.org/project/weather-api/ ? at least for locations





-----------------------------------------------------------------------------------
'''

import matplotlib
#matplotlib.use('Qt5Agg') 
#matplotlib.use('Agg') 
#matplotlib.use('GTK3Agg')

import urllib
from urllib.request import urlopen
from io import StringIO
from sunpy.time import parse_time
import numpy as np
from datetime import datetime
from dateutil import tz
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.dates import  DateFormatter
import mpl_toolkits  
import matplotlib.dates as mdates
import sys
import datetime
import skimage.transform
import scipy
import aacgmv2
import pdb
import os
import getopt
import time
import pickle
import seaborn as sns
import pandas as pd
from numba import njit
import importlib
from multiprocessing import Pool, cpu_count, Array

from pandas.plotting import register_matplotlib_converters               
register_matplotlib_converters()                                        

#import auroramaps
#importlib.reload(auroramaps) 
from auroramaps import ovation as amo
from auroramaps import util as amu

importlib.reload(amu) #reload again while debugging
importlib.reload(amo) #reload again while debugging




##################################### FUNCTIONS #######################################

def make_aurora_cube_multi(ts,ec,k):
    '''
    make the aurora flux image cuves 
    multiprocessing version 
    - ts is a single datetime object into this function
    - ec the averaged Newell coupling inout
    - k the counter for the frames
    for debugging making one frame use: >> make_aurora_cube_multi(tsm[0],ecm[0],0)   
    '''
    
    print('Frame number and time:', k, '  ',ts)
    
    #################  (2a) get fluxes
        
    mlatN, mltN, fluxNd=de.get_flux_for_time(ts,ec)
    mlatN, mltN, fluxNm=me.get_flux_for_time(ts,ec)
    fluxN=fluxNd+fluxNm #+fluxNw
    #print(ts), print(ec), print(k), print('....')
    
    ################  (2b) coordinate conversion magnetic to geographic 
    #Coordinate conversion MLT to AACGM mlon/lat to geographic coordinates
    mlonN_1D_small=aacgmv2.convert_mlt(mltN[0],ts,m2a=True)
    mlonN_1D=np.tile(mlonN_1D_small,mlatN.shape[0])
    mlatN_1D=np.squeeze(mlatN.reshape(np.size(mltN),1))
    (glatN_1D, glonN_1D, galtN) = aacgmv2.convert_latlon_arr(mlatN_1D,mlonN_1D, 100,ts, code="A2G") #**check 100 km

    ##############  (2c) interpolate to world map 
    geo_2D=np.vstack((glatN_1D,glonN_1D)).T      #stack 2 (7680,) arrays to a single 7680,2 arrays, .T is needed
    fluxN_1D=fluxN.reshape(7680,1)   #also change flux values to 1D array

    #make a world map grid in latitude 512 pixels, longitude 1024 pixel like NOAA
    wx,wy= np.mgrid[-90:90:180/512,-180:180:360/1024]
    aimg=  np.squeeze(scipy.interpolate.griddata(geo_2D, fluxN_1D, (wx, wy), method='linear',fill_value=0))
    aimg = scipy.ndimage.gaussian_filter(aimg,sigma=(5,7),mode='wrap') #wrap means wrapping at the 180 degree edge
      
    #Array variable to be used by all processes
    ovation_img_multi[512*1024*k:512*1024*(k+1)]=aimg.reshape(512*1024)
   
    




#######################################################################################
##################################### Main ############################################
#######################################################################################


#if new background images need to be saved in auroramaps/data/wmts do this, some are included!
#amu.save_gibs_earth_image('BlueMarble_NextGeneration',300)    
#amu.save_gibs_earth_image('VIIRS_CityLights_2012',300)    
#amu.save_gibs_earth_image('BlueMarble_NextGeneration',600)    
#amu.save_gibs_earth_image('VIIRS_CityLights_2012',600)    
        



################ READ INPUT OPTIONS FROM COMMAND LINE
argv = sys.argv[1:]
opts, args = getopt.getopt(argv,"hv=",["server", "real"])#, "help", "historic=", "verbose="])

server = False
if "--server" in [o for o, v in opts]:
    server = True
    print("in server mode")

if server:
    matplotlib.use('Agg') 
else:
    matplotlib.use('Qt5Agg') # figures are shown on mac

#real switch set to True different input file
real = False

if "--real" in [o for o, v in opts]:
    real = True
    print("using input_realtime.py")
else: 
    print("using input.py")    
    
if real:
    import input_realtime
    importlib.reload(input_realtime)   #make sure it reads file again
    from input_realtime import *       #gets all variables from this file

else:
    import input
    importlib.reload(input)   #make sure it reads file again
    from input import *       #gets all variables from this file



############### (0) get input data, get PREDSTORM solar wind files mode ###############


start_all=time.time()
plt.close('all')
sns.set_style('darkgrid')

#get current time as datetime object in UTC, rounded to minute and time zone aware
utcnow=amu.round_to_minute(datetime.datetime.now(tz=tz.tzutc())) 

print()
print('Making aurora forecasts with OVATION PRIME 2010')
print('with the PREDSTORM solar wind predictions (ACE/DSCOVR+STEREO-A) or OMNI2 data')
print()
print('UTC time now:')
print(utcnow.strftime(format="%Y-%m-%d %H:%M") )
print()
    
if mode==0:
    print('mode '+str(mode)+': PREDSTORM real time')
    t0     = utcnow + datetime.timedelta(hours=past_hours)
    tend   = utcnow + datetime.timedelta(hours=future_hours)
    
if mode==1:
    print('mode '+str(mode)+': PREDSTORM local file')
    t0   = parse_time(start_time).datetime   #parse start time from string to datetime
    tend = parse_time(end_time).datetime

if mode==2:
    print('mode '+str(mode)+': OMNI2 data')
    t0   = parse_time(start_time).datetime   #parse start time from string to datetime
    tend = parse_time(end_time).datetime     

tdiff=(tend-t0)     #difference between start and end time
n_hours=int(tdiff.total_seconds()/3600)  #for how many hours the model runs

#time resolution as set in input - smallest time resolution allowed is 1 second
ts = [t0 + datetime.timedelta(seconds=i) for i in range(0, int(tdiff.total_seconds())+1,int(time_resolution*60))]

print()
print('start time:',ts[0].strftime('%Y-%m-%d %H:%M UT' ))
if mode==0: print('now time:  ',utcnow.strftime('%Y-%m-%d %H:%M UT' ))
print('end time:  ',ts[-1].strftime('%Y-%m-%d %H:%M UT' ))
print('timerange          ',n_hours, 'hours')
print('time resolution    ',time_resolution,'minutes')
print('nr of movie frames ', len(ts) )
print()

#check if all needed directories are there
if os.path.isdir('results') == False: os.mkdir('results')
if os.path.isdir('results/'+output_directory) == False: os.mkdir('results/'+output_directory)
if os.path.isdir('results/'+output_directory+'/flux_global') == False: os.mkdir('results/'+output_directory+'/flux_global')
if os.path.isdir('results/'+output_directory+'/flux_europe') == False: os.mkdir('results/'+output_directory+'/flux_europe')
if os.path.isdir('results/'+output_directory+'/flux_canada') == False: os.mkdir('results/'+output_directory+'/flux_canada')
if os.path.isdir('results/'+output_directory+'/prob_global') == False: os.mkdir('results/'+output_directory+'/prob_global')
if os.path.isdir('results/'+output_directory+'/prob_europe') == False: os.mkdir('results/'+output_directory+'/prob_europe')
if os.path.isdir('results/'+output_directory+'/prob_canada') == False: os.mkdir('results/'+output_directory+'/prob_canada')

if os.path.isdir('auroramaps/data/predstorm') == False: os.mkdir('auroramaps/data/predstorm')
if os.path.isdir('auroramaps/data/omni2') == False: os.mkdir('auroramaps/data/omni2')


# get or set input files
if mode==0:    
   try: 
       inputfile='auroramaps/data/predstorm/predstorm_real.txt'
       urllib.request.urlretrieve(predstorm_url,inputfile)
       print('loaded from', predstorm_url)
   except urllib.error.URLError as e:
       print('Failed downloading ', predstorm_url,' ',e.reason)

if mode==1:     
     inputfile=local_input_file

if mode==2:     
     amu.omni_txt_generator(ts)   #make txt file from OMNI2 data in similar format as predstorm
     inputfile='auroramaps/data/predstorm/predstorm_omni.txt'
    
print('input data file:',inputfile)

print()
print('output directory: results/'+output_directory)
print('------------------------------------------------------')












########################## (1) Initialize OVATION ########################################


########## load ovation for different types of aurora
'''
atype - str, ['diff','mono','wave','ions']
         type of aurora for which to load regression coeffients
jtype - int or str
            1:"electron energy flux",
            2:"ion energy flux",
            3:"electron number flux",
            4:"ion number flux",
            5:"electron average energy",
            6:"ion average energy"
'''
jtype = 'electron energy flux'
print('Initialize OVATION')
start = time.time()
de = amo.FluxEstimator('diff', jtype)
me = amo.FluxEstimator('mono', jtype)
#we = amo.FluxEstimator('wave', jtype)
end = time.time()
print('done, took ',np.round(end - start,2),' seconds.')
print('OVATION uses',jtype,'with diffuse + monoenergetic aurora')
print('fluxes for southern hemisphere are currently not calculated.')
print()


###################### load input solar wind
l1wind=amu.load_predstorm_wind(inputfile)
print('Solar wind data loaded from PREDSTORM input file.')
swav=amu.calc_avg_solarwind_predstorm(ts,l1wind)  # calculate average solar wind for Newell coupling 
window=int(window_minutes/time_resolution)	#when time resolution < averaging window in minutes, do moving averages
coup_cycle=4421 #average coupling for solar cycle (see e.g. Newell et al. 2010)

####################### plot current coupling, the driver behind the ovation model
fig, ax = plt.subplots(figsize=[10, 5],dpi=100)
plt.plot_date(l1wind.time,np.ones(np.size(l1wind.time)),'k--',label='cycle average',linewidth=0.7)
plt.plot_date(l1wind.time,l1wind.ec/coup_cycle,'k-', label='input solar wind')   
plt.plot_date(swav.time,swav.ec/coup_cycle,'r-',label='4-hour weighted averages',markersize=2)   
plt.title('Newell coupling Nc')
plt.ylabel(r'Nc / 4421 $\mathrm{[(km/s)^{4/3}\/nT^{2/3}]}$')
ax.xaxis.set_major_formatter( DateFormatter('%Y-%b-%d %Hh') )
plt.xticks(rotation=45)


#plot moving averages and make them if time resolution high enough
if window > 0:
   print('running mean used for Nc with time window +/- ',window_minutes,'minutes')
   ec_run_mean=swav.ec
   for i in np.arange(window,np.size(ts)-window): ec_run_mean[i]=np.mean(swav.ec[i-window:i+window])  
   plt.plot_date(swav.time,ec_run_mean/coup_cycle,'b--',label='smoothed weighted averages' )   
   #replace Ec average by running mean
   swav.ec=ec_run_mean
    
ax.set_xlim([swav.time[0]-4/24,swav.time[-1]])
plt.legend()

plt.tight_layout()
fig.savefig('results/'+output_directory+'/run_newell_coupling.png',dpi=150,facecolor=fig.get_facecolor()) #save plot in outputdirectory
















###################### (2) RUN OVATION FOR EACH TIME STEP ##############################



##################### (2a) get flux data cubes
print()
print('------------------------------------------------------')
print('Now make the aurora flux data cubes including all timesteps.')
print()

#depending on switch in input.py:

if calc_mode == 'single':
    print('Using single processing')
    print()
    start = time.time()
    print('clock run time start ...')
    ovation_img=amo.make_aurora_cube(ts,swav.ec,de,me) #ovation_img single processing - use ts and ec to make aurora world map; get back 512*1024 image cube
    print('... end run time clock:  ',np.round(time.time() - start,2),' sec total, per frame: ',np.round((time.time() - start)/np.size(ts),2) )


if calc_mode == 'multi':  
    print('Using multiprocessing, nr of cores',cpu_count())
    print()

    start = time.time()     #time measure indicator start
    print('clock run time start ...')

    #these variables will be used for pool.starmap
    tsm=ts; ecm=swav.ec; km=np.arange(np.size(ts))
    
    oshape = (512, 1024,np.size(ts))       #define shape of final array
    arrsize=oshape[0]*oshape[1]*oshape[2]  #define size of 1D array
    
    ovation_img_multi = Array('d', arrsize) #make 1D array to be used by the processes simultaneously for the ovation map
    p = Pool()                                              #make multiprocessing Pool object  
    res=p.starmap(make_aurora_cube_multi, zip(tsm,ecm,km))  #goes through all ts times, needs Ec and counter too
    p.close()
    p.join()
    oim=np.frombuffer(ovation_img_multi.get_obj())          #get calculated array values into 1D array

    #make final array 512*1024*size(ts) out of 1D array that is used in make_aurora_cube_multi
    ovation_img=amu.reshape_ovation_img_multi(np.zeros(oshape),oim,oshape) 
 
    print('... end run time clock:  ',np.round(time.time() - start,2),' sec total, per frame: ',np.round((time.time() - start)/np.size(ts),2) )


print()
print('------------------------------------------------------')

'''
#for testing conversion and image making - note that the image is upside down with north at bottom
plt.close('all')
plt.figure(1)
plt.imshow(aimg)
sys.exit()
'''


##################### (2b) get lower equatorial boundary 

print('Now make equatorial boundary')
start = time.time()
print('clock run time start ... for total number of frames:', np.size(ts))

#dictionary for putting together latitude, longitude, boundary, time and flux map data for saving
eb = {'lat':np.linspace(-90,90,512), 'long':np.linspace(-180,180,1024),       \
      'boundary':np.zeros([np.size(ts),np.size(np.linspace(-180,180,1024))]), \
      'smooth':np.zeros([np.size(ts),np.size(np.linspace(-180,180,1024))]),   \
      'time': ts}#,'flux_map': ovation_img}
     
#make the equatorial boundary
eb['boundary']=amu.make_equatorial_boundary(ovation_img,eb['boundary'],np.size(ts),eb['lat'],equatorial_boundary_flux_threshold) 
#the result is eb['data'] as function of longitude variable eb['long'] 
ebwin=15 #size of filter window, apply filter
eb['smooth']=amu.smooth_boundary(ts,eb['boundary'],ebwin)
   
print('... end run time clock:  ',np.round(time.time() - start,2),' sec total, per frame: ',np.round((time.time() - start)/np.size(ts),2))
print('total number of frames:', np.size(ts))
print()
print('------------------------------------------------------')


#saving boundary data
eb_save_filename='results/'+output_directory+'/run_data_boundary.p' #run_data_flux_map_boundary.p'
pickle.dump(eb,open(eb_save_filename, 'wb' ))
print('Boundary data saved as '+eb_save_filename)
print()
   
        
'''for debugging equatorial boundary
plt.figure(13)
plt.plot(eb['long'],eb['boundary'][0,:],'-r')   
plt.plot(eb['long'],eb['smooth'][0,:],'-b')   
'''






#####################################  (3) PLOTS and MOVIES  #########################################


############################ (3a) Make global aurora plots 
#maybe make faster with multiprocessing pool - does not work on MacOS but should on Linux
#try multiprocessing with initializer?
#https://docs.python.org/dev/library/multiprocessing.html#multiprocessing.pool.Pool

print('Make all movie frames')  
start = time.time()


if calc_mode_frame == 'single':  
    print('with single processing')  

    ############ load background image
    map_img=amu.load_high_res_background(map_type)

    #flux maps
    if global_flux_map > 0:
      amu.plot_ovation_single(ovation_img, ts, output_directory, eb, map_type, map_img, 'global', 'flux',utcnow,swav.ec)

    if europe_flux_map > 0:
      amu.plot_ovation_single(ovation_img, ts, output_directory, eb, map_type, map_img, 'europe', 'flux',utcnow,swav.ec)
    
    if canada_flux_map > 0:
      amu.plot_ovation_single(ovation_img, ts, output_directory, eb, map_type, map_img, 'canada', 'flux',utcnow,swav.ec)

    ########### same for probability maps

    #first convert flux to probability
    ovation_img_prob=amu.flux_to_probability(ovation_img)

    if global_probability_map > 0:
      amu.plot_ovation_single(ovation_img_prob, ts, output_directory, eb, map_type, map_img, 'global', 'prob',utcnow,swav.ec)

    if europe_probability_map > 0:
      amu.plot_ovation_single(ovation_img_prob, ts, output_directory, eb, map_type, map_img, 'europe', 'prob',utcnow,swav.ec)

    if canada_probability_map > 0:
      amu.plot_ovation_single(ovation_img_prob, ts, output_directory, eb, map_type, map_img, 'canada', 'prob',utcnow,swav.ec)




if calc_mode_frame == 'multi':  
    print('with multiprocessing')  

    ############ load background image
    map_img=amu.load_high_res_background(map_type)

    #flux maps
    if global_flux_map > 0:
      amu.plot_ovation_multi(ovation_img, ts, output_directory, eb, map_type, map_img, 'global', 'flux',utcnow,swav.ec)

    if europe_flux_map > 0:
      amu.plot_ovation_multi(ovation_img, ts, output_directory, eb, map_type, map_img, 'europe', 'flux',utcnow,swav.ec)
    
    if canada_flux_map > 0:
      amu.plot_ovation_multi(ovation_img, ts, output_directory, eb, map_type, map_img, 'canada', 'flux',utcnow,swav.ec)

    ########### same for probability maps

    #first convert flux to probability
    ovation_img_prob=amu.flux_to_probability(ovation_img)

    if global_probability_map > 0:
      amu.plot_ovation_multi(ovation_img_prob, ts, output_directory, eb, map_type, map_img, 'global', 'prob',utcnow,swav.ec)

    if europe_probability_map > 0:
      amu.plot_ovation_multi(ovation_img_prob, ts, output_directory, eb, map_type, map_img, 'europe', 'prob',utcnow,swav.ec)

    if canada_probability_map > 0:
      amu.plot_ovation_multi(ovation_img_prob, ts, output_directory, eb, map_type, map_img, 'canada', 'prob',utcnow,swav.ec)




######################################

number_of_maps=sum([global_flux_map, europe_flux_map, canada_flux_map, global_probability_map,europe_probability_map, canada_probability_map])
end = time.time()
print(number_of_maps, ' different maps were produced.')
print('All movie frames took ',np.round(end - start,2),'sec, per frame',np.round((end - start)/np.size(ts)*1/number_of_maps,2),' sec.')


################################# (3b) make movies 
print()
print('Make mp4 and gif movies')
print()
print('For all results see: results/'+output_directory)

#frame rate is set in input.py

if global_flux_map > 0:
  os.system('ffmpeg -r '+str(frame_rate)+' -i results/'+output_directory+'/flux_global/aurora_%05d.jpg -b:v 5000k -r '+str(frame_rate)+' results/'+output_directory+'/flux_global.mp4 -y -loglevel quiet')
  os.system('ffmpeg -r '+str(frame_rate)+' -i results/'+output_directory+'/flux_global/aurora_%05d.jpg -b:v 5000k -r '+str(frame_rate)+' results/'+output_directory+'/flux_global.gif -y -loglevel quiet')
  ########## convert mp4 to gif and makes smaller
  os.system('ffmpeg -i results/'+output_directory+'/flux_global.mp4  -vf scale=1000:-1 results/'+output_directory+'/flux_global_small.gif  -y -loglevel quiet ')

if europe_flux_map > 0:
  os.system('ffmpeg -r '+str(frame_rate)+' -i results/'+output_directory+'/flux_europe/aurora_%05d.jpg -b:v 5000k -r '+str(frame_rate)+' results/'+output_directory+'/flux_europe.mp4 -y -loglevel quiet')
  os.system('ffmpeg -r '+str(frame_rate)+' -i results/'+output_directory+'/flux_europe/aurora_%05d.jpg -b:v 5000k -r '+str(frame_rate)+' results/'+output_directory+'/flux_europe.gif -y -loglevel quiet')
  os.system('ffmpeg -i results/'+output_directory+'/flux_europe.mp4  -vf scale=1000:-1 results/'+output_directory+'/flux_europe_small.gif  -y -loglevel quiet ')

if canada_flux_map > 0:
  os.system('ffmpeg -r '+str(frame_rate)+' -i results/'+output_directory+'/flux_canada/aurora_%05d.jpg -b:v 5000k -r '+str(frame_rate)+' results/'+output_directory+'/flux_canada.mp4 -y -loglevel quiet')
  os.system('ffmpeg -r '+str(frame_rate)+' -i results/'+output_directory+'/flux_canada/aurora_%05d.jpg -b:v 5000k -r '+str(frame_rate)+' results/'+output_directory+'/flux_canada.gif -y -loglevel quiet')
  os.system('ffmpeg -i results/'+output_directory+'/flux_canada.mp4  -vf scale=1000:-1 results/'+output_directory+'/flux_canada_small.gif  -y -loglevel quiet ')


if global_probability_map > 0:
  os.system('ffmpeg -r '+str(frame_rate)+' -i results/'+output_directory+'/prob_global/aurora_%05d.jpg -b:v 5000k -r '+str(frame_rate)+' results/'+output_directory+'/prob_global.mp4 -y -loglevel quiet')
  os.system('ffmpeg -r '+str(frame_rate)+' -i results/'+output_directory+'/prob_global/aurora_%05d.jpg -b:v 5000k -r '+str(frame_rate)+' results/'+output_directory+'/prob_global.gif -y -loglevel quiet')
  os.system('ffmpeg -i results/'+output_directory+'/prob_global.mp4  -vf scale=1000:-1 results/'+output_directory+'/prob_global_small.gif  -y -loglevel quiet ')

if europe_probability_map > 0:
  os.system('ffmpeg -r '+str(frame_rate)+' -i results/'+output_directory+'/prob_europe/aurora_%05d.jpg -b:v 5000k -r '+str(frame_rate)+' results/'+output_directory+'/prob_europe.mp4 -y -loglevel quiet')
  os.system('ffmpeg -r '+str(frame_rate)+' -i results/'+output_directory+'/prob_europe/aurora_%05d.jpg -b:v 5000k -r '+str(frame_rate)+' results/'+output_directory+'/prob_europe.gif -y -loglevel quiet')
  os.system('ffmpeg -i results/'+output_directory+'/prob_europe.mp4  -vf scale=1000:-1 results/'+output_directory+'/prob_europe_small.gif  -y -loglevel quiet ')

if canada_probability_map > 0:
  os.system('ffmpeg -r '+str(frame_rate)+' -i results/'+output_directory+'/prob_canada/aurora_%05d.jpg -b:v 5000k -r '+str(frame_rate)+' results/'+output_directory+'/prob_canada.mp4 -y -loglevel quiet')
  os.system('ffmpeg -r '+str(frame_rate)+' -i results/'+output_directory+'/prob_canada/aurora_%05d.jpg -b:v 5000k -r '+str(frame_rate)+' results/'+output_directory+'/prob_canada.gif -y -loglevel quiet')
  os.system('ffmpeg -i results/'+output_directory+'/prob_canada.mp4  -vf scale=1000:-1 results/'+output_directory+'/prob_canada_small.gif  -y -loglevel quiet ')



print()
print('Run time for everything:  ',np.round((time.time() - start_all)/60,2),' min; per frame: ',np.round((time.time() - start_all)/np.size(ts)*1/number_of_maps,2),'sec' )
print()


plt.style.use('seaborn-whitegrid') #for white ticks and labels, resetting the dark background for plotting the maps


##################################### END ################################################





