'''
Plotting an aurora forecast based on the PREDSTORM solar wind prediction method

using a rewritten version of 
ovationpyme by Liam Kilcommons https://github.com/lkilcommons/OvationPyme

by C. Moestl, IWF-helio group, Graz, Austria.
twitter @chrisoutofspace
https://www.iwf.oeaw.ac.at/user-site/christian-moestl/

last update May 2019

--------------------------

TO DO: 

- check errors in ovation? how to get probabilites correctly?
- colorbars for probabilites
- colormap better?
- VIIRS night mode - takes too long for each frame, delete image and replace when making movie?
- split land on dayside / night lights on night side
- historic mode with OMNI2 data (1 hour)
- code optimizen, insbesondere ovation, coordinate conversion take 2 functions used
- Newell solar wind coupling als parameter in plot
- auroral power on plot directly from predstorm_real.txt
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
  
#ovation model
import ovation_prime_predstorm as opp
#during development automatic reload
importlib.reload(opp) 

#extra functions
import ovation_utilities_predstorm as oup
#during development automatic reload
importlib.reload(oup) 

#input parameters - make sure they are reloaded every time when using ipython
import aurora_forecast_input
importlib.reload(aurora_forecast_input)   
from aurora_forecast_input import *



##################################### FUNCTIONS #######################################












def global_predstorm_north(world_image,dt,counter):


 plt.close('all')
 fig = plt.figure(1,figsize=[12, 12]) 
 fig.set_facecolor('black') 

 #axis PREDSTORM + OVATION
 ax = plt.subplot(1, 1, 1, projection=ccrs.Orthographic(global_plot_longitude, global_plot_latitude))


 #ax.add_feature(carfeat.BORDERS, color='white',alpha=0.5)
 #ax.add_feature(carfeat.LAND,color='darkslategrey')
 #ax.coastlines(alpha=0.5,zorder=3)
 #ax.add_feature(land_50m, color='darkgreen') 

 #ax.add_feature(land_50m, color='darkslategrey')
 #ax.add_feature(carfeat.LAKES,color='navy')#,zorder=2,alpha=1)
 #ax.add_feature(carfeat.OCEAN)#,zorder=2,alpha=1)
 #ax.add_feature(ocean_50m,linewidth=0.5, color='navy')
 
 #ax.add_feature(carfeat.COASTLINE,color='white')#,zorder=2,alpha=0.5)
 #ax.add_feature(carfeat.RIVERS)#,zorder=2,alpha=0.8)
 #ax.add_feature(provinces_50m,alpha=0.5)#,zorder=2,alpha=0.8)

 #ax.add_feature(carfeat.COASTLINE, alpha=0.5,color='white')#,zorder=2,alpha=0.5)

 ax.background_patch.set_facecolor('k')    
 ax.coastlines('50m',color='white',alpha=0.8)
 #ax.add_feature(provinces_50m,alpha=0.5)#,zorder=2,alpha=0.8)
 gl=ax.gridlines(linestyle='--',alpha=0.5,color='white')
 gl.n_steps=100
 #ax.stock_img()#alpha=0.2)
 #ax.add_wmts(nightmap, layer)

 ax.imshow(world_image, vmin=0, vmax=100, transform=crs, extent=mapextent, origin='lower', zorder=3, alpha=0.9, cmap=oup.aurora_cmap())
 #ax.add_feature(Nightshade(dt))

 #pdb.set_trace()
   
 fig.text(0.5,0.92,'PREDSTORM aurora forecast  '+dt.strftime('%Y-%m-%d %H:%M UT'), color='white',fontsize=15, ha='center')
 fig.text(0.99,0.03,'C. Möstl / IWF-helio, Austria', color='white',fontsize=10,ha='right')

 plt.tight_layout()  

 #save as image with timestamp in filename
 plot_Nhemi_filename='forecast_global/predstorm_aurora_real_Nhemi_'+dt.strftime("%Y_%m_%d_%H%M")  +'.jpg'
 fig.savefig(plot_Nhemi_filename,dpi=150,facecolor=fig.get_facecolor())

 #save as movie frame
 framestr = '%05i' % (counter)  
 fig.savefig('frames_global/aurora_'+framestr+'.jpg',dpi=130,facecolor=fig.get_facecolor())
 #plt.show()
 print('Saved image:  ',plot_Nhemi_filename)
















def europe_canada_predstorm(world_image,dt,counter):

 #16/9 ration for full hd output
 fig = plt.figure(2,figsize=[16, 9]) 
 fig.set_facecolor('black') 
 # ax1 Europe
 ax1 = plt.subplot(1, 2, 2, projection=ccrs.Orthographic(0, 60),position=[0.51,0.05,0.48,0.9])#[left, bottom, width, height]
 # ax2 northern America
 ax2 = plt.subplot(1, 2, 1, projection=ccrs.Orthographic(-100, 60), position=[0.01,0.05,0.48,0.9])

 #define map extents
 canada_east = -65
 canada_west = -135 
 canada_north = 75
 canada_south = 20

 europe_east = 50
 europe_west = -20
 europe_north = 75
 europe_south = 30


 #plot one axis after another
 for ax in [ax1, ax2]:

    #ax.gridlines(linestyle='--',alpha=0.2,color='white')
    #ax.coastlines(alpha=0.5,zorder=3)

    #ax.add_feature(carfeat.LAND,zorder=2,alpha=1)
    #ax.add_feature(land_50m, color='darkgreen')
    #ax.add_feature(land_50m, color='darkslategrey')
    ax.add_feature(land_50m, color='dimgrey')
    
    ax.add_feature(provinces_50m,alpha=0.5)#,zorder=2,alpha=0.8)

  
    #ax.add_feature(carfeat.COASTLINE)#,zorder=2,alpha=0.5)
    #ax.add_feature(carfeat.OCEAN,color='mediumslateblue')#,zorder=2,alpha=1)
    ax.add_feature(ocean_50m, color='darkblue')

    #ax.add_feature(carfeat.RIVERS)#,zorder=2,alpha=0.8)
    ax.add_feature(carfeat.LAKES,color='darkblue')#,zorder=2,alpha=1) color navy?
    
    ax.add_feature(carfeat.BORDERS, alpha=0.5)#,zorder=2,alpha=0.5)
  
    ax.add_feature(Nightshade(dt))
  
    #ax.stock_img()
    
    #ax.add_wmts(nightmap, layer)
  
    ax.imshow(world_image, vmin=0, vmax=100, transform=crs,extent=mapextent, origin='lower', zorder=3, alpha=0.9, cmap=oup.aurora_cmap2())
    
    if ax == ax1: ax.set_extent([europe_west, europe_east, europe_south, europe_north])
    if ax == ax2: ax.set_extent([canada_west, canada_east, canada_south, canada_north])

 fig.text(0.01,0.92,'PREDSTORM aurora forecast   '+dt.strftime('%Y-%m-%d %H:%M UT' ), color='white',fontsize=15)
 fig.text(0.99,0.02,'C. Möstl / IWF-helio, Austria', color='white',fontsize=8,ha='right')

 
 #save as image with timestamp in filename
 plot_europe_canada_filename='forecast_euro_canada/predstorm_aurora_real_'+dt.strftime("%Y_%m_%d_%H%M")  +'.jpg'
 fig.savefig(plot_europe_canada_filename,dpi=150,facecolor=fig.get_facecolor())


 #save as movie frame
 framestr = '%05i' % (counter)  
 fig.savefig('frames_euro_canada/aurora_'+framestr+'.jpg',dpi=120,facecolor=fig.get_facecolor())
 #plt.show()
 print('Saved image:  ',plot_europe_canada_filename)

















##################################### Main ############################################


plt.close('all')



print()
print('Making aurora forecasts with OVATION from the PREDSTORM solar wind predictions')
print()
print()
print('Start time:',ts[0].strftime('%Y-%m-%d %H:%M UT' ))
print('End time:',ts[-1].strftime('%Y-%m-%d %H:%M UT' ))
print('Number of hourly time steps: ',n_hours)
print()
print()



########################### (0) get PREDSTORM file ######################################

###########   get online file or local file from predstorm
if use_predstorm_online_file:
    try: 
      urllib.request.urlretrieve(predstorm_online_url_1hour,'predstorm_real.txt')
      inputfile='predstorm_real.txt'
      print('use input from url:')
      print(predstorm_online_url_1hour)
    except urllib.error.URLError as e:
        print('Failed downloading ', predstorm_online_url_1hour,' ',e.reason)
    
else:
    inputfile=local_input_file
    print('use input from local file: ')
    print(inputfile)


print()
########################## (1) LOAD OVATION OBJECTS ####################################

print('Load FluxEstimator objects first.')
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
de = opp.FluxEstimator('diff', jtype)
me = opp.FluxEstimator('mono', jtype)
we = opp.FluxEstimator('wave', jtype)



###################### (2) RUN OVATION FOR EACH FRAME TIME ##############################


#make a world map grid in latitude 512 pixels, longitude 1024 pixel like NOAA
wx,wy=np.mgrid[-90:90:180/512,-180:180:360/1024]


#this is the array with the world maps for each timestep
ovation_img=np.zeros([512,1024,np.size(ts)])


start = time.time()
print('Start time for run time check ...')


for k in np.arange(0,np.size(ts)):
 
 print('--------------------------------------------------------------')
 print('time step:', k)

 #########  (2b) get solar wind and fluxes

 print()
 print('Time of PREDSTORM frame',ts[k])
 print()
 print('solar wind input from weighting procedure:')
 #make solarwind with averaging over last 4 hours
 sw=oup.calc_avg_solarwind_predstorm(ts[k],inputfile,4)
 print(sw)
 #for Ec, a cycle average is 4421
 print('Newell coupling compared to solar cycle average: ',np.round(sw['Ec']/4421,2))
 print()
 #print('get fluxes for northern hemisphere and sum them')
 mlatN, mltN, fluxNd=de.get_flux_for_time(ts[k],inputfile, hemi='N')
 mlatN, mltN, fluxNm=me.get_flux_for_time(ts[k],inputfile, hemi='N')
 mlatN, mltN, fluxNf=we.get_flux_for_time(ts[k],inputfile, hemi='N')

 #sum all fluxes
 fluxN=fluxNd+fluxNm+fluxNf

 #print('For our aurora maps we use the ',jtype, ' with diff, mono, wave')
 #print('Fluxes for southern hemisphere are currently not calculated.')


 #########  (2b) coordinate conversion magnetic to geographic 

 print()
 #print('Coordinate conversion MLT to AACGM mlon/lat to geographic coordinates:')
 #startcoo=time.time()

 #magnetic coordinates are  mltN mlonN; convert magnetic local time (MLT) to longitude
 mlonN=aacgmv2.convert_mlt(mltN,ts[k],m2a=True)

 #magnetic coordinates are now mlatN mlonN
 (glatN, glonN, galtN)= aacgmv2.convert_latlon_arr(mlatN,mlonN, 100,ts[k], code="A2G")
 #endcoo = time.time()
 #print('Coordinate conversion takes seconds: ', np.round(endcoo-startcoo,3))

 #geographic coordinates are glatN, glonN, electron flux values are fluxN


 #change shapes from glatN (80, 96) glonN  (80, 96) to a single array with 7680,2
 glatN_1D=glatN.reshape(np.size(fluxN),1)  
 glonN_1D=glonN.reshape(np.size(fluxN),1)    
 #stack 2 (7680,1) arrays to a single 7680,2 arrays
 geo_2D=np.hstack((glatN_1D,glonN_1D)) 
 #also change flux values to 1D array
 fluxN_1D=fluxN.reshape(7680,1) 

 #interpolate to world grid, and remove 1 dimension with squeeze
 #see https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html
 aimg=np.squeeze(scipy.interpolate.griddata(geo_2D, fluxN_1D, (wx, wy), method='linear',fill_value=0))


 '''
 #for testing conversion and image making - note that the image is upside down with north at bottom
 plt.close('all')
 plt.figure(1)
 plt.imshow(aimg)
 sys.exit()
 '''



 #####################  (2c) Convert to probabilities, smooth etc.
 print()
 
 
 # ************* NOT CORRECT -> how to go from erg cm-2 s-1 to probability?
 #Case et al. 2016
 #print('Convert to probabilites, smooth with gaussian,')
 #print('round numbers, trim values < 2, cut at 100 % probability')
 #convert to probabilities from energy flux in erg cm−2 s−1 Case et al. 2016

 #*******CHECK for UNITS in ovationpyme
 
 pimg=8*aimg
   
 #pimg=10+8*aimg

 #smooth out artefacts
 pimg = scipy.ndimage.gaussian_filter(pimg,sigma=(3,5))

 #round probabilities?
 #pimg=np.round(pimg,0)

 #trim small values
 pimg[np.where(pimg <2)]=0

 #cut at 100 percent probability
 pimg[np.where(pimg >100)]=100


 #add calibration of probabilities compared to NOAA
 adhoc_calibration=3
 ovation_img[:,:,k]=pimg*adhoc_calibration



 '''
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




end = time.time()

print()
print('... end time:  ',np.round(end - start,2),' sec')
print('------------------------------------------------')
#print()
#print('Calculation for 20 frames would take:', np.round((end - start)*20/60,2),' minutes')



#####################################  (3) PLOTS  #########################################
#Night map from VIIRS
#https://wiki.earthdata.nasa.gov/display/GIBS/GIBS+Available+Imagery+Products#expand-EarthatNight4Products
nightmap = 'https://map1c.vis.earthdata.nasa.gov/wmts-geo/wmts.cgi'


#define extent of the produced world maps - defined as: west east south north
mapextent=[-180,180,-90,90]   


#not so good: layer='VIIRS_Black_Marble'
#better but takes time
layer = 'VIIRS_CityLights_2012'



land_50m = carfeat.NaturalEarthFeature('physical', 'land', '50m',
                                        edgecolor='k',
                                        facecolor=carfeat.COLORS['land'])
ocean_50m = carfeat.NaturalEarthFeature('physical', 'ocean', '50m',
                                        edgecolor='k',
                                        facecolor='steelblue')#carfeat.COLORS['water'])                                        
provinces_50m = carfeat.NaturalEarthFeature('cultural',
                                             'admin_1_states_provinces_lines',
                                             '50m',
                                             facecolor='none',edgecolor='black')
crs=ccrs.PlateCarree()



############################ (3a) Make global aurora plot for comparison with NOAA nowcast

#comparison with NOAA
#global_predstorm_noaa(ovation_img)


for k in np.arange(0,np.size(ts)):
    global_predstorm_north(ovation_img[:,:,k],ts[k],k)
    europe_canada_predstorm(ovation_img[:,:,k],ts[k],k)
    
#make move with frames 
os.system('ffmpeg -r 10 -i frames_global/aurora_%05d.jpg -b:v 5000k -r 10 predstorm_aurora_global.mp4 -y -loglevel quiet')
os.system('ffmpeg -r 10 -i frames_euro_canada/aurora_%05d.jpg -b:v 5000k -r 10 predstorm_aurora_euro_canada.mp4 -y -loglevel quiet')

os.system('ffmpeg -r 10 -i frames_global/aurora_%05d.jpg -b:v 5000k -r 10 predstorm_aurora_global.gif -y -loglevel quiet')
os.system('ffmpeg -r 10 -i frames_euro_canada/aurora_%05d.jpg -b:v 5000k -r 10 predstorm_aurora_euro_canada.gif -y -loglevel quiet')

    
    

############################ (3b) zoom North America and Europe ############################################

#europe_canada_predstorm(ovation_img)


print()
print('Fin.')


