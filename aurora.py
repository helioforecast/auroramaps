'''
aurora.py

main program for running the OVATION PRIME 2010 model to make an aurora forecast/hindcast 
based on the PREDSTORM solar wind prediction method (available online)
or OMNI2 data for historic events

uses the python "auroramaps" package
https://github.com/IWF-helio/auroramaps

input parameters are given in the input.py file

by C. Moestl, IWF-helio group, Graz, Austria.
Contributions by Liam Kilcommons and Diana Morosan

twitter @chrisoutofspace
https://www.iwf.oeaw.ac.at/user-site/christian-moestl/

This package uses a rewritten version of the ovationpyme aurora model 
by Liam Kilcommons https://github.com/lkilcommons/OvationPyme

published under GNU Lesser General Public License v3.0

last update October 2019

-----------------------------------------------------------------------------------
TO DO: 

core:
- use probabilities as shown by Diana
- check with IDL code, debuggen, both hemispheres correct?
- get flux for time -> weights for both seasons? how correctly?, compare further with ovationpyme and IDL version
- multiprocessing with queue?

later:
- add wave flux (needs good interpolation like OP13)
- add southern hemisphere on world map (with coordinate conversion etc.)
- use numba or multiprocessing somewhere further for speed? 
- multiprocessing for saving the frames does not work on MacOS! 
  maybe on linux simply use the plotting function with multiprocessing.pool
  and multiple arguments (starmap) like for the data cubes
- 180° longitude on world map better interpolation (interpolate on mlatgrid before? or own way to do it instead of wrap?)
- auroral power on plot directly from ovation integrated (add function in amo)


plotting:
- cut viewing line in daylight? 
- Newell solar wind coupling als parameter in plot
- add colorbars for probabilites, colormap should fade into background but oval should also visible for small values
- check with direct comparison with NOAA global images nowcast
- transparent to white colormap so that it looks like viirs images for direct comparison
- split land on dayside / night lights on night side 
  this should work in global_predstorm_north by loading background only once
  need to figure out how to get pixels from nightshade day/night and how to plot 
  only specific pixels (but then each background image must be updated)
- Newell solar wind coupling als parameter in plot
- indicate moon phase with astropy
- cloud cover how? https://pypi.org/project/weather-api/ ? at least for locations

test bottlenecks: 
>> python -m cProfile -s tottime aurora_forecast.py

or use in ipython for function run times:
>> %timeit function_name()
>> %time  function_name()






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
import cartopy.crs as ccrs
import cartopy.feature as carfeat
from cartopy.feature.nightshade import Nightshade
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import mpl_toolkits  
import matplotlib.dates as mdates
import sys
import datetime
import skimage.transform
import scipy
import aacgmv2
import pdb
import os
import time
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

import input
importlib.reload(input)   #make sure it reads file again
from input import *       #gets all variables from this file




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



############### (0) get input data, get PREDSTORM solar wind files mode ###############

start_all=time.time()
plt.close('all')

#get current time as datetime object in UTC, rounded to minute
utcnow=amu.round_to_minute(datetime.datetime.utcnow()) 

print()
print('Making aurora forecasts with OVATION PRIME 2010')
print('with the PREDSTORM solar wind predictions (DSCOVR/STEREO-A) or OMNI2 data')
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
if os.path.isdir('results/'+output_directory+'/flux_europe_canada') == False: os.mkdir('results/'+output_directory+'/flux_europe_canada')
if os.path.isdir('results/'+output_directory+'/prob_europe_canada') == False: os.mkdir('results/'+output_directory+'/prob_europe_canada')
if os.path.isdir('results/'+output_directory+'/prob_global') == False: os.mkdir('results/'+output_directory+'/prob_global')

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
fig, ax = plt.subplots(figsize=[12, 5])
plt.plot_date(swav.time,np.ones(np.size(swav.time)),'k-.',label='cycle average')
plt.plot_date(l1wind.time,l1wind.ec/coup_cycle,'k-', label='input wind')   
plt.plot_date(swav.time,swav.ec/coup_cycle,'r-',label='weighted averages',markersize=2)   
plt.title('Newell coupling')

#plot moving averages and make them if time resolution high enough
if window > 0:
   print('running mean used for Ec with time window +/- ',window_minutes,'minutes')
   ec_run_mean=swav.ec
   for i in np.arange(window,np.size(ts)-window): ec_run_mean[i]=np.mean(swav.ec[i-window:i+window])  
   plt.plot_date(swav.time,ec_run_mean/coup_cycle,'b--',label='running mean of weighted averages' )   
   #replace Ec average by running mean
   swav.ec=ec_run_mean
    
ax.set_xlim([swav.time[0]-4/24,swav.time[-1]])
plt.legend()
fig.savefig('results/'+output_directory+'/coupling.png',dpi=150,facecolor=fig.get_facecolor()) #save plot in outputdirectory


















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


##################### (2b) get lower equatorial boundary and viewing line

print('Now make equatorial boundary and viewing line')
start = time.time()
print('clock run time start ... for total number of frames:', np.size(ts))

#define the latitude longitude grid again 
all_lat=np.linspace(-90,90,512)
all_long=np.linspace(-180,180,1024)   
eb=np.zeros([np.size(ts),np.size(all_long)])    #define array of equatorial boundaries eb

#make the equatorial boundary
eb=amu.make_equatorial_boundary(ovation_img,eb,np.size(ts),all_lat,1.0) 
#the result is eb as function of longitude variable all_long

ebwin=15 #size of filter window
ebs=amu.smooth_boundary(ts,eb,ebwin)
    
print('... end run time clock:  ',np.round(time.time() - start,2),' sec total, per frame: ',np.round((time.time() - start)/np.size(ts),2))
print('total number of frames:', np.size(ts))
print()
print('------------------------------------------------------')
    
'''
plt.figure(12)
#plt.plot(np.arange(0,512*6),ebi[0,:],'or')   
#plt.plot(np.arange(0,512*6),ebis[0,:],'-k')   

plt.plot(all_long,eb[0,:],'or')   
plt.plot(all_long,ebs[0,:],'-k')   
'''




#comparison with NOAA
#global_predstorm_noaa(ovation_img)

#plot the fluxes for direct comparison to ovation prime in IDL
#for k in np.arange(0,np.size(ts)):
#    oup.global_predstorm_flux(ovation_img[:,:,k],ts[k],k)

#for k in np.arange(0,np.size(ts)):
#    oup.europe_canada_predstorm(ovation_img[:,:,k],ts[k],k, 'hot')





#####################################  (3) PLOTS and MOVIES  #########################################


############################ (3a) Make global aurora plot for comparison with NOAA nowcast
#maybe make faster with multiprocessing pool - does not work on MacOS but should on Linux
#try multiprocessing with initializer?
#https://docs.python.org/dev/library/multiprocessing.html#multiprocessing.pool.Pool

print('Make all movie frames')  
print()

start = time.time()



#global flux images
if global_flux_map==True:
  amu.ovation_global_north(ovation_img,ts,'hot',max_level_flux,output_directory,all_long,ebs)

#europe and canada/USA flux images
if europe_canada_flux_map==True:
  amu.ovation_europe_canada(ovation_img,ts,'hot',max_level_flux,output_directory,all_long,ebs)

#same for probability maps
if global_probability_map==True:
  amu.ovation_probability_global_north(ovation_img,ts,'hot',100,output_directory,all_long,ebs)

if europe_canada_probability_map==True:
  amu.ovation_probability_europe_canada(ovation_img,ts,'hot',100,output_directory,all_long,ebs)

end = time.time()
print('All movie frames took ',np.round(end - start,2),'sec, per frame',np.round((end - start)/np.size(ts),2),' sec.')


################################# (3b) make movies 
print()
print('Make mp4 and gif movies')
print()
print('For all results see: results/'+output_directory)

#frame rate 20 is good for 10 minute resolution if 3 days want to be seen quickly

if global_flux_map==True:
  os.system('ffmpeg -r 20 -i results/'+output_directory+'/flux_global/aurora_%05d.jpg -b:v 5000k -r 20 results/'+output_directory+'/flux_global.mp4 -y -loglevel quiet')
  os.system('ffmpeg -r 20 -i results/'+output_directory+'/flux_global/aurora_%05d.jpg -b:v 5000k -r 20 results/'+output_directory+'/flux_global.gif -y -loglevel quiet')
  ########## convert mp4 to gif and makes smaller
  os.system('ffmpeg -i results/'+output_directory+'/flux_global.mp4  -vf scale=1000:-1 results/'+output_directory+'/flux_global_small.gif  -y -loglevel quiet ')

if europe_canada_flux_map==True:
  os.system('ffmpeg -r 20 -i results/'+output_directory+'/flux_europe_canada/aurora_%05d.jpg -b:v 5000k -r 20 results/'+output_directory+'/flux_europe_canada.mp4 -y -loglevel quiet')
  os.system('ffmpeg -r 20 -i results/'+output_directory+'/flux_europe_canada/aurora_%05d.jpg -b:v 5000k -r 20 results/'+output_directory+'/flux_europe_canada.gif -y -loglevel quiet')
  os.system('ffmpeg -i results/'+output_directory+'/flux_europe_canada.mp4  -vf scale=1000:-1 results/'+output_directory+'/flux_europe_canada_small.gif  -y -loglevel quiet ')

if global_probability_map==True:
  os.system('ffmpeg -r 20 -i results/'+output_directory+'/prob_global/aurora_%05d.jpg -b:v 5000k -r 20 results/'+output_directory+'/prob_global.mp4 -y -loglevel quiet')
  os.system('ffmpeg -r 20 -i results/'+output_directory+'/prob_global/aurora_%05d.jpg -b:v 5000k -r 20 results/'+output_directory+'/prob_global.gif -y -loglevel quiet')
  os.system('ffmpeg -i results/'+output_directory+'/prob_global.mp4  -vf scale=1000:-1 results/'+output_directory+'/prob_global_small.gif  -y -loglevel quiet ')

if europe_canada_probability_map==True:
  os.system('ffmpeg -r 20 -i results/'+output_directory+'/prob_europe_canada/aurora_%05d.jpg -b:v 5000k -r 20 results/'+output_directory+'/prob_europe_canada.mp4 -y -loglevel quiet')
  os.system('ffmpeg -r 20 -i results/'+output_directory+'/prob_europe_canada/aurora_%05d.jpg -b:v 5000k -r 20 results/'+output_directory+'/prob_europe_canada.gif -y -loglevel quiet')
  os.system('ffmpeg -i results/'+output_directory+'/prob_europe_canada.mp4  -vf scale=1000:-1 results/'+output_directory+'/prob_europe_canada_small.gif  -y -loglevel quiet ')

print()
print('Run time for everything:  ',np.round((time.time() - start_all)/60,2),' min; per frame: ',np.round((time.time() - start_all)/np.size(ts),2),'sec' )
print()


##################################### END ################################################





