'''
Plotting an aurora forecast based on the PREDSTORM solar wind prediction method

using a rewritten version of 
ovationpyme by Liam Kilcommons https://github.com/lkilcommons/OvationPyme

by C. Moestl, IWF-helio group, Graz, Austria.
twitter @chrisoutofspace
https://www.iwf.oeaw.ac.at/user-site/christian-moestl/

last update May 2019

--------------------------

test bottlenecks: 
>> python -m cProfile -s tottime aurora_forecast.py

or use in ipython
>> %prun function_name 

TO DO: 

- check errors in ovation? how to get probabilites correctly? as Nathan Case
- add colorbars for probabilites, colormap better?
- VIIRS night mode - takes too long for each frame, delete image and replace when making movie?
- split land on dayside / night lights on night side
- historic mode with OMNI2 data (1 hour)
- code optimizen, insbesondere ovation, coordinate conversion take 2 functions used
- Newell solar wind coupling als parameter in plot
- auroral power on plot directly from predstorm_real.txt or from ovationpyme
- indicate moon phase with astropy
- cloud cover how? https://pypi.org/project/weather-api/ ?
- higher time resolution than 1 hour - use 1 min real time file, adapt calc_solarwind...
- black white colormap so that it looks like viirs images for direct comparison
- add equatorial auroral boundary case et al. 2016



'''

import matplotlib
matplotlib.use('Qt5Agg') 

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
import numba as nb
import importlib
  
  
# get own modules; make sure they are 
# automatically reloaded every time when using ipython
  
#ovation model
import ovation_prime_predstorm as opp
importlib.reload(opp) 

#extra functions and for plotting
import ovation_utilities_predstorm as oup
importlib.reload(oup) 

#input parameters 
import aurora_forecast_input
importlib.reload(aurora_forecast_input)   
from aurora_forecast_input import *



##################################### FUNCTIONS #######################################







#######################################################################################
##################################### Main ############################################
#######################################################################################

plt.close('all')



print()
print('Making aurora forecasts with OVATION')
print('From the PREDSTORM solar wind predictions or OMNI2 data')
print()

if mode ==0: 
    real_time_mode=True
    historic_mode=False
if mode ==1: 
    historic_mode=True
    real_time_mode=False

    
if real_time_mode:
    print('running in real time mode using DSCOVR data')
    #take time now + 1hour forecast and take nearest hour (as PREDSTORM is 1 hour resolution)
    t0=oup.round_to_hour(datetime.datetime.utcnow()+datetime.timedelta(hours=1))
    #or force to given time
    t0_real_forced_str='2019-May-14 00:00'
    t0 = parse_time(t0_real_forced_str)
  

if historic_mode:
    print('running in historic mode using OMNI2 data')
    #parse start time from string to datetime
    t0 = parse_time(t0_historic_str)



#make array of hours starting with t0 for which auroramaps are produced
ts = [t0 + datetime.timedelta(hours=i) for i in range(0, n_hours,1)]


print()
print('Start time:',ts[0].strftime('%Y-%m-%d %H:%M UT' ))
print('End time:  ',ts[-1].strftime('%Y-%m-%d %H:%M UT' ))
print('Number of hourly time steps: ',n_hours)
print()
print()



#check if output directories are there
if os.path.isdir('results') == False: os.mkdir('results')
if os.path.isdir('results/frames_global') == False: os.mkdir('results/frames_global')
if os.path.isdir('results/frames_europe_canada') == False: os.mkdir('results/frames_europe_canada')
if os.path.isdir('results/forecast_europe_canada') == False: os.mkdir('results/forecast_europe_canada')
if os.path.isdir('results/forecast_global') == False: os.mkdir('results/forecast_global')

 



############### (0) get PREDSTORM solar wind input files, depending on mode ##############

if real_time_mode:
   ###########   get online file or local file from predstorm
   if use_predstorm_online_file:
      try: 
        urllib.request.urlretrieve(predstorm_online_url_1hour,'predstorm_real.txt')
        inputfile='predstorm_real.txt'
        print('use solar wind input from url:')
        print(predstorm_online_url_1hour)
      except urllib.error.URLError as e:
        print('Failed downloading ', predstorm_online_url_1hour,' ',e.reason)
    
   else:
     inputfile=local_input_file
     print('use solar wind input from local file: ')
     print(inputfile)

if historic_mode:
     #make txt file in similar format as predstorm
     oup.omni_txt_generator(ts)     
     inputfile='predstorm_omni.txt'
     print('txt file with solar wind parameters generated from OMNI2 data: ')
     print(inputfile)
     #print('starting 24 hours earlier than ovation start time to make solar wind averages')
    

print()

########################## (1) Initialize OVATION ########################################

print()
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
de = opp.FluxEstimator('diff', jtype)
me = opp.FluxEstimator('mono', jtype)
#we = opp.FluxEstimator('wave', jtype)
end = time.time()
print('done, took ',np.round(end - start,2),' seconds.')

print()

 
 
print('OVATION uses ',jtype, ' with diffuse and monoenergetic aurora types')
print('Fluxes for southern hemisphere are currently not calculated.')

print()


###################### (2) RUN OVATION FOR EACH FRAME TIME ##############################

print('Now make the aurora flux data cubes for all timesteps.')


#make a world map grid in latitude 512 pixels, longitude 1024 pixel like NOAA
wx,wy=np.mgrid[-90:90:180/512,-180:180:360/1024]

#this is the array with the final world maps for each timestep
ovation_img=np.zeros([512,1024,np.size(ts)])

#this is the array with the flux grid in magnetic coordinates for each timestep
flux_img=np.zeros([80,96,np.size(ts)])

start = time.time()
print('--------------------------------------------------------------')
print('Clock run time start ...')
print()

for k in np.arange(0,np.size(ts)): #go through all times
 
    
    print('Frame number and time:', k, '  ',ts[k])

    #########################################  (2a) get solar wind 
    
    startflux=time.time()
    sw=oup.calc_avg_solarwind_predstorm(ts[k],inputfile)   #make solarwind with averaging over last 4 hours 
    #print('solar wind input last 4 hours with weights: ')
    print('Bxyz =',sw.bx[0],sw.by[0],sw.bz[0],' nT   V =',int(sw.v[0]), 'km/s')
    print('Newell coupling to average: ',np.round(sw.ec[0]/4421,1), ' Ec =',int(sw.ec[0]))    #for the coupling Ec, a cycle average is 4421
    #get fluxes for northern hemisphere and sum them
    mlatN, mltN, fluxNd=de.get_flux_for_time(ts[k],inputfile, hemi='N')
    mlatN, mltN, fluxNm=me.get_flux_for_time(ts[k],inputfile, hemi='N')
    #mlatN, mltN, fluxNw=we.get_flux_for_time(ts[k],inputfile, hemi='N')     #wave flux not correct yet
    fluxN=fluxNd+fluxNm#+fluxNw
    endflux=time.time()
    print('OVATION: ', np.round(endflux-startflux,2),' sec')

    #for debugging
    #making flux images for comparison to OVATION IDL output
    #change IDL file in this function for flux comparison
    #oup.global_ovation_flux(mlatN,mltN,fluxNw,ts[0])
    #oup.global_ovation_flux(mlatN,mltN,fluxNd+fluxNm,ts[0])

 
 
    #####################################  (2b) coordinate conversion magnetic to geographic 
    #Coordinate conversion MLT to AACGM mlon/lat to geographic coordinates
    startcoo=time.time()
    #magnetic coordinates are  mltN mlonN; convert magnetic local time (MLT) to longitude
    mltN_1D=mltN.reshape(np.size(mltN),1)  
    mlatN_1D=mlatN.reshape(np.size(mltN),1)  
    
    
    #********extract from mltN_1D first 96, convert and copy to rest
    
    mlonN_1D=aacgmv2.convert_mlt(mltN_1D,ts[k],m2a=True)
    endcoo = time.time()
    print('Coordinates 1: ', np.round(endcoo-startcoo,2),' sec')
   
    #sys.exit()
   
    startcoo=time.time()
    #magnetic coordinates are now mlatN mlonN - this is very fast
    (glatN_1D, glonN_1D, galtN)= aacgmv2.convert_latlon_arr(mlatN_1D,mlonN_1D, 100,ts[k], code="A2G")
    endcoo = time.time()
    print('Coordinates 2: ', np.round(endcoo-startcoo,2),' sec')
   


    #####################################  (2c) interpolate to world map 
    #geographic coordinates are glatN, glonN, electron flux values are fluxN
    #change shapes from glatN (80, 96) glonN  (80, 96) to a single array with 7680,2
    startworld=time.time()
 
    #glatN_1D=glatN.reshape(np.size(fluxN),1)  
    #glonN_1D=glonN.reshape(np.size(fluxN),1)    
    geo_2D=np.hstack((glatN_1D,glonN_1D))      #stack 2 (7680,1) arrays to a single 7680,2 arrays
    fluxN_1D=fluxN.reshape(7680,1)   #also change flux values to 1D array

    #interpolate to world grid, and remove 1 dimension with squeeze
    aimg=np.squeeze(scipy.interpolate.griddata(geo_2D, fluxN_1D, (wx, wy), method='linear',fill_value=0))
    #filter
    aimg = scipy.ndimage.gaussian_filter(aimg,sigma=(5,7))
    ovation_img[:,:,k]=aimg
    endworld = time.time()
    print('World map: ', np.round(endworld-startworld,2),' sec')
    print()
    print('---------------------')


 

 
end = time.time()
print()
print('... end run time clock:  ',np.round(end - start,2),' sec total, per frame: ',np.round((end - start)/np.size(ts),2) )
print('--------------------------------------------------------------')    

'''
#for testing conversion and image making - note that the image is upside down with north at bottom
plt.close('all')
plt.figure(1)
plt.imshow(aimg)
sys.exit()
'''

 

'''
#####################  (2c) Convert to probabilities, smooth etc.
 
# ************* NOT CORRECT -> how to go from erg cm-2 s-1 to probability?
#Case et al. 2016
#print('Convert to probabilites, smooth with gaussian,')
#print('round numbers, trim values < 2, cut at 100 % probability')
#convert to probabilities from energy flux in erg cm−2 s−1 Case et al. 2016
 
#pimg=8*aimg
   
#pimg=10+8*aimg

#smooth out artefacts
#pimg = scipy.ndimage.gaussian_filter(pimg,sigma=(3,5))

#round probabilities?
#pimg=np.round(pimg,0)

#trim small values
#pimg[np.where(pimg <2)]=0

#cut at 100 percent probability
#pimg[np.where(pimg >100)]=100

#add calibration of probabilities compared to NOAA
#adhoc_calibration=3
#ovation_img[:,:,k]=pimg*adhoc_calibration
 
 
for testing
plt.close()
#16/9 ration for full hd output
fig = plt.figure(figsize=[10, 5]) 
ax1 = plt.subplot(1, 1, 1, projection=crs)
fig.set_facecolor('black') 
ax1.coastlines()
ax1.imshow(pimg, vmin=0, vmax=100, transform=crs, extent=mapextent, origin='lower',zorder=3,cmap=aurora_cmap())
#ax1.set_extent(fullextent)
'''








#####################################  (3) PLOTS  #########################################


############################ (3a) Make global aurora plot for comparison with NOAA nowcast

#comparison with NOAA
#global_predstorm_noaa(ovation_img)

#plot the fluxes for direct comparison to ovation prime in IDL
#for k in np.arange(0,np.size(ts)):
#    oup.global_predstorm_flux(ovation_img[:,:,k],ts[k],k)


#************ optimize frame making

#make images and movie frames for the generated aurora image cube
for k in np.arange(0,np.size(ts)):
    oup.global_predstorm_north(ovation_img[:,:,k],ts[k],k)

for k in np.arange(0,np.size(ts)):
    oup.europe_canada_predstorm(ovation_img[:,:,k],ts[k],k)
    
        
    
#make move with frames 
os.system('ffmpeg -r 10 -i results/frames_global/aurora_%05d.jpg -b:v 5000k -r 10 results/predstorm_aurora_global.mp4 -y -loglevel quiet')
os.system('ffmpeg -r 10 -i results/frames_europe_canada/aurora_%05d.jpg -b:v 5000k -r 10 results/predstorm_aurora_europe_canada.mp4 -y -loglevel quiet')

#os.system('ffmpeg -r 10 -i results/frames_global/aurora_%05d.jpg -b:v 5000k -r 10 results/predstorm_aurora_global.gif -y -loglevel quiet')
#os.system('ffmpeg -r 10 -i results/frames_europe_canada/aurora_%05d.jpg -b:v 5000k -r 10 results/predstorm_aurora_europe_canada.gif -y -loglevel quiet')

    
    

############################ (3b) zoom North America and Europe ############################################

#europe_canada_predstorm(ovation_img)


print()
print('Fin.')


