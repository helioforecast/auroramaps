'''
Plotting an aurora forecast/hindcast 
based on the PREDSTORM solar wind prediction method
or OMNI2 data for historic events

using a rewritten version of the ovationpyme aurora model 
by Liam Kilcommons https://github.com/lkilcommons/OvationPyme
and the cartopy package

by C. Moestl, IWF-helio group, Graz, Austria.
https://github.com/IWF-helio
twitter @chrisoutofspace
https://www.iwf.oeaw.ac.at/user-site/christian-moestl/

published under GNU Lesser General Public License v3.0

last update May 2019

-----------------------------------------------------------------------------------

TO DO: 

core:
- calculate average Ec first, smooth? and then use in get_flux_for_time
- sum both hemispheres as in IDL ovation and ovationpyme
- add southern hemisphere complete
- code optimize - bottlenecks with numba (where grids are calculated); 
  and use normal numpy arrays and check for optimization (functions often used, grids...)
- check bugs in ovation in comparison to IDL version
- 180° longitude on world map better interpolation (interpolate on mlatgrid before?)

new features:
- add equatorial auroral boundary Case et al. 2016
- how to get probabilites correctly? ask Nathan Case

plotting:

- add colorbars for probabilites, colormap should fade into background but oval should also visible for small values
- make nowcast better check with direct comparison with NOAA global images
- transparent to white colormap so that it looks like viirs images for direct comparison
- split land on dayside / night lights on night side 
  this should work in global_predstorm_north by loading background only once
  need to figure out how to get pixels from nightshade day/night and how to plot 
  only specific pixels (but then each background image must be updated)
- Newell solar wind coupling als parameter in plot
- auroral power on plot directly from ovation integrated (add function in opp)
- indicate moon phase with astropy
- cloud cover how? https://pypi.org/project/weather-api/ ? at least for locations


test bottlenecks: 
>> python -m cProfile -s tottime aurora_forecast.py

or use in ipython for functions
>> %timeit function_name 

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
  
  
# make sure own modules are 
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



############### (0) get input data, get PREDSTORM solar wind files mode ###############


start_all=time.time()

plt.close('all')

utcnow=oup.round_to_minute(datetime.datetime.utcnow()) #get current time as datetime object in UTC

print()
print('Making aurora forecasts with OVATION')
print('From the PREDSTORM solar wind predictions (DSCOVR/STEREO-A) or OMNI2 data')
print('UTC time now:')
print(utcnow.strftime(format="%Y-%m-%d %H:%M") )
print()
    
if mode==0:
    print('mode',mode,': PREDSTORM real time')
    t0     = utcnow - datetime.timedelta(hours=abs(past_hours))
    tend   = utcnow + datetime.timedelta(hours=future_hours)
    
if mode==1:
    print('mode',mode,': PREDSTORM local file')
    t0   = parse_time(start_time)   #parse start time from string to datetime
    tend = parse_time(end_time)     

if mode==2:
    print('mode',mode,': OMNI2 data')
    t0   = parse_time(start_time)   #parse start time from string to datetime
    tend = parse_time(end_time)     

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
if os.path.isdir('results/'+output_directory+'/frames_global') == False: os.mkdir('results/'+output_directory+'/frames_global')
if os.path.isdir('results/'+output_directory+'/frames_europe_canada') == False: os.mkdir('results/'+output_directory+'/frames_europe_canada')
if os.path.isdir('results/'+output_directory+'/forecast_europe_canada') == False: os.mkdir('results/'+output_directory+'/forecast_europe_canada')
if os.path.isdir('results/'+output_directory+'/forecast_global') == False: os.mkdir('results/'+output_directory+'/forecast_global')

if os.path.isdir('data') == False: os.mkdir('data')
if os.path.isdir('data/predstorm') == False: os.mkdir('data/predstorm')


# get or set input files
if mode==0:    
   try: 
       inputfile='data/predstorm/predstorm_real.txt'
       urllib.request.urlretrieve(predstorm_url,inputfile)
       print('loaded from', predstorm_url)
   except urllib.error.URLError as e:
       print('Failed downloading ', predstorm_url,' ',e.reason)

if mode==1:     
     inputfile=local_input_file

if mode==2:     
     oup.omni_txt_generator(ts)   #make txt file from OMNI2 data in similar format as predstorm
     inputfile='data/predstorm/predstorm_omni.txt'
    
print('input data file:',inputfile)

print()
print('output directory: results/'+output_directory)
print('------------------------')




########################## (1) Initialize OVATION ########################################


########## load ovation for different types of aurora

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
print('OVATION uses',jtype,'with diffuse + monoenergetic aurora')
print('Fluxes for southern hemisphere are currently not calculated.')
print()


###################### load input solar wind
l1wind=oup.load_predstorm_wind(inputfile)
print('Solar wind data loaded from PREDSTORM input file.')


swav=oup.calc_avg_solarwind_predstorm(ts,l1wind)  # calculate average solar wind for Newell coupling 

window=int(window_minutes/time_resolution)	#when time resolution < averaging window in minutes, do moving averages

coup_cycle=4421 #average coupling for solar cycle (see e.g. Newell et al. 2010)

#plot current coupling, the driver behind the ovation model
fig, ax = plt.subplots(figsize=[12, 5])
plt.plot_date(swav.time,np.ones(np.size(swav.time)),'k-.',label='cycle average')
plt.plot_date(l1wind.time,l1wind.ec/coup_cycle,'k-', label='input wind')   
plt.plot_date(swav.time,swav.ec/coup_cycle,'r-',label='weighted averages',markersize=2)   
plt.title('Newell coupling')

#plot moving averages and make them 
if window > 0:
   print('running mean used for Ec with time window +/- ',window_minutes,'minutes')
   ec_run_mean=swav.ec
   for i in np.arange(window,np.size(ts)-window): ec_run_mean[i]=np.mean(swav.ec[i-window:i+window])  
   plt.plot_date(swav.time,ec_run_mean/coup_cycle,'b--',label='running mean of weighted averages' )   
   #replace Ec average by running mean
   swav.ec=ec_run_mean
     
ax.set_xlim([swav.time[0]-4/24,swav.time[-1]])
plt.legend()
#save plot in outputdirectory
fig.savefig('results/'+output_directory+'/coupling.png',dpi=150,facecolor=fig.get_facecolor())








###################### (2) RUN OVATION FOR EACH FRAME TIME ##############################

print('Now make the aurora flux data cubes for all timesteps.')

#make a world map grid in latitude 512 pixels, longitude 1024 pixel like NOAA
wx,wy=np.mgrid[-90:90:180/512,-180:180:360/1024]

#this is the array with the final world maps for each timestep
ovation_img=np.zeros([512,1024,np.size(ts)])

#this is the array with the flux grid in magnetic coordinates for each timestep
#flux_img=np.zeros([80,96,np.size(ts)])

#time measure indicator start
start = time.time()
print('--------------------------------------------------------------')
print('Clock run time start ...')
print()


for k in np.arange(0,len(ts)): #go through all times
 
    print('Frame number and time:', k, '  ',ts[k])

    #########################################  (2a) get solar wind 
    startflux=time.time()
   
    sw=oup.calc_avg_solarwind_predstorm(ts[k],l1wind)   #make solarwind with averaging over last 4 hours, only for command line 
    print('Byz =',sw.by[0],sw.bz[0],' nT   V =',int(sw.v[0]), 'km/s')
    #for the Newell coupling Ec, normalized to cycle average
    print('Ec to cycle average: ',np.round(sw.ec[0]/coup_cycle,1), ' Ec =',int(sw.ec[0]))  
    
    #get fluxes for northern hemisphere and sum them **check -> do averages and sum both
    mlatN, mltN, fluxNd=de.get_flux_for_time(ts[k],swav.ec[k], hemi='N')
    mlatN, mltN, fluxNm=me.get_flux_for_time(ts[k],swav.ec[k], hemi='N')
    #mlatN, mltN, fluxNw=we.get_flux_for_time(ts[k],inputfile, hemi='N')  #wave flux not correct yet
    fluxN=fluxNd+fluxNm #+fluxNw
    endflux=time.time()
    print('OVATION: ', np.round(endflux-startflux,2),' sec')
    
    #for debugging
    #making flux images for comparison to OVATION IDL output
    #change IDL file in this function for flux comparison
    #oup.global_ovation_flux(mlatN,mltN,fluxNw,ts[0])
    #oup.global_ovation_flux(mlatN,mltN,fluxNd,ts[k])
    #sys.exit()


 
    #####################################  (2b) coordinate conversion magnetic to geographic 
    #Coordinate conversion MLT to AACGM mlon/lat to geographic coordinates
    startcoo=time.time()
    #magnetic coordinates are  mltN mlonN in 2D grids
    #so we need to convert magnetic local time (MLT) to longitude first
    #extract from mltN first 96 MLTs for 1 latitude bin convert to magnetic longitude 
    mlonN_1D_small=aacgmv2.convert_mlt(mltN[0],ts[k],m2a=True)
    #copy result rest 80 times because MLT is same for all latitudes
    mlonN_1D=np.tile(mlonN_1D_small,mlatN.shape[0])
    
    #this will make a 1D array for the latitudes compatible with the 1D array created above    
    mlatN_1D=np.squeeze(mlatN.reshape(np.size(mltN),1))
   
    #magnetic coordinates are now mlatN mlonN, convert to geographic
    (glatN_1D, glonN_1D, galtN) = aacgmv2.convert_latlon_arr(mlatN_1D,mlonN_1D, 100,ts[k], code="A2G")
    endcoo = time.time()
    print('Coordinates: ', np.round(endcoo-startcoo,2),' sec')
   

    #####################################  (2c) interpolate to world map 
    #geographic coordinates are glatN, glonN, electron flux values are fluxN
    #change shapes from glatN (80, 96) glonN  (80, 96) to a single array with 7680,2
    startworld=time.time()
 
    #glatN_1D=glatN.reshape(np.size(fluxN),1)  
    #glonN_1D=glonN.reshape(np.size(fluxN),1)    
    geo_2D=np.vstack((glatN_1D,glonN_1D)).T      #stack 2 (7680,) arrays to a single 7680,2 arrays, .T is needed
    #smooth small array first - but strange results!
    #fluxN=scipy.ndimage.gaussian_filter(fluxN,sigma=(5,5))
    fluxN_1D=fluxN.reshape(7680,1)   #also change flux values to 1D array

    #interpolate to world grid, and remove 1 dimension with squeeze
    aimg=np.squeeze(scipy.interpolate.griddata(geo_2D, fluxN_1D, (wx, wy), method='linear',fill_value=0))
    #filter large array
    aimg = scipy.ndimage.gaussian_filter(aimg,sigma=(5,7),mode='wrap') #wrap means wrapping at the 180 degree edge
    ovation_img[:,:,k]=aimg
    endworld = time.time()
    print('World map: ', np.round(endworld-startworld,2),' sec')
    

    print()
    print('---------------------')




 
end = time.time()
print()
print('... end run time clock:  ',np.round(end - start,2),' sec total, per frame: ',np.round((end - start)/np.size(ts),2) )
print('--------------------------------------------------------------')   
print() 

'''
#for testing conversion and image making - note that the image is upside down with north at bottom
plt.close('all')
plt.figure(1)
plt.imshow(aimg)
sys.exit()
'''

 

#####################################  (3) PLOTS  #########################################


############################ (3a) Make global aurora plot for comparison with NOAA nowcast

start = time.time()
#better use color map that starts with some basic green color
oup.ovation_global_north(ovation_img,ts,'hot',1.5,output_directory)
end = time.time()
print('All movie frames took ',np.round(end - start,2),'sec, per frame',np.round((end - start)/np.size(ts),2),' sec.')

#make move with frames 
print()
print('Make movies')
os.system('ffmpeg -r 25 -i results/'+output_directory+'/frames_global/aurora_%05d.jpg -b:v 5000k -r 25 results/'+output_directory+'/predstorm_aurora_global.mp4 -y -loglevel quiet')
os.system('ffmpeg -r 25 -i results/'+output_directory+'/frames_global/aurora_%05d.jpg -b:v 5000k -r 25 results/'+output_directory+'/predstorm_aurora_global.gif -y -loglevel quiet')
print()
end_all=time.time()
print('Run time for everything:  ',np.round((end_all - start_all)/60,2),' min; per frame: ',np.round((end_all - start_all)/np.size(ts),2),'sec' )

print()




































#####################################


#make images and movie frames for the generated aurora image cube ovation_img
#for k in np.arange(0,np.size(ts)):
    #oup.global_predstorm_north(ovation_img[:,:,k],ts[k],k,'magma')
    #oup.global_predstorm_north(ovation_img[:,:,k],ts[k],k,'hot')
    #oup.global_predstorm_north(ovation_img[:,:,k],ts[k],k,oup.aurora_cmap2())

#eliminate noise when plotting
#ovation_img[ovation_img < 0.1]=np.nan  

#flux plot of last timestep
#oup.global_ovation_flux(mlatN,mltN,fluxNd,ts[k])

#comparison with NOAA
#global_predstorm_noaa(ovation_img)

#plot the fluxes for direct comparison to ovation prime in IDL
#for k in np.arange(0,np.size(ts)):
#    oup.global_predstorm_flux(ovation_img[:,:,k],ts[k],k)

#for k in np.arange(0,np.size(ts)):
#    oup.europe_canada_predstorm(ovation_img[:,:,k],ts[k],k, 'hot')

#os.system('ffmpeg -r 25 -i results/frames_europe_canada/aurora_%05d.jpg -b:v 5000k -r 25 results/predstorm_aurora_europe_canada.mp4 -y -loglevel quiet')
#os.system('ffmpeg -r 10 -i results/frames_europe_canada/aurora_%05d.jpg -b:v 5000k -r 10 results/predstorm_aurora_europe_canada.gif -y -loglevel quiet')


############################ (3b) zoom North America and Europe ############################################

#europe_canada_predstorm(ovation_img)



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
pimg = scipy.ndimage.gaussian_filter(pimg,sigma=(3,5))

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



