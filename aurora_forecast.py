'''
Plotting an aurora forecast based on the PREDSTORM solar wind prediction method

using a rewritten version of 
ovationpyme by Liam Kilcommons https://github.com/lkilcommons/OvationPyme

by C. Moestl, IWF-helio, Graz, Austria.
twitter @chrisoutofspace

---

TO DO: 

- make MLT longitude conversion correct (mlong dreht sich mit
- colorbars probabilites
- movie of next 5 days mp4 / gif movie machen
- dunkler den hintergrund am Tag
- historic mode with OMNI2 data
- code optimizen, insbesondere ovation
- Newell solar wind coupling als parameter in plot
- auroral power directly from predstorm_real.txt
- moon phase astropy
- cloud cover how? https://pypi.org/project/weather-api/ ?
- higher time resolution than 1 hour

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
import time

import ovation_prime_predstorm as opp
import ovation_utilities_predstorm as oup



######################### input parameters 


#initialize

#test
#t0 = parse_time("2019-Apr-16 06:00")
#inputfile='predstorm_sample/predstorm_real.txt'

#set time

#t0 = parse_time("2019-May-01 23:00")
#t0 = parse_time("2019-Apr-29 23:00")




#take time now + 1hour forecast and take nearest hour (as PREDSTORM is 1 hour resolution)
t0=oup.round_to_hour(datetime.datetime.utcnow()+datetime.timedelta(hours=1))

#t0 = parse_time("2019-May-02 04:00")


#for a real time run - run pred5 before
inputfile='/Users/chris/python/predstorm/predstorm_real.txt'

#-100 for North America, +10 for Europe
global_plot_longitude=-100

global_plot_latitude=90
##################################### FUNCTIONS #######################################













##################################### Main ############################################


#make solarwind with averaging over last 4 hours
sw=oup.calc_avg_solarwind_predstorm(t0,inputfile,4)

print()
print('Running OVATION aurora model with PREDSTORM input')
print()
print()
print('input file',inputfile)


print('time',t0)
print()

print('solar wind input from weighting procedure:')
print(sw)
#for Ec, a cycle average is 4421
print('Newell coupling compared to solar cycle average: ',np.round(sw['Ec']/4421,2))

#sw=calc_avg_solarwind(t0,'predstorm_sample/predstorm_real.txt')
#print(sw)


#t0 = parse_time("2010-Apr-5 08:00")




#electron energy flux - mono diff wave needed

print()
print('Now run own OVATION developed from OvationPyme')

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



############################################# run ovationpyme for each frame

#ts = [t0 + datetime.timedelta(hours=i) for i in range(1, 24,1)]
#print(ts)
#for k in np.arange(1,24):
#    print(k)


start = time.time()

print('Fluxes for northern hemisphere')
mlatN, mltN, fluxNd=de.get_flux_for_time(t0,inputfile, hemi='N')
mlatN, mltN, fluxNm=me.get_flux_for_time(t0,inputfile, hemi='N')
mlatN, mltN, fluxNf=we.get_flux_for_time(t0,inputfile, hemi='N')

#sum all fluxes
fluxN=fluxNd+fluxNm+fluxNf

print('for aurora maps we use the ',jtype, ' with diff, mono, wave')




############################################# coordinate conversion

print()
print('Coordinate conversion MLT to AACGM mlon/lat to geographic coordinates ')


startcoo=time.time()
#Convert to geographic coordinates
mlonN=aacgmv2.convert_mlt(mltN,t0,m2a=True)
(glatN, glonN, galtN)= aacgmv2.convert_latlon_arr(mlatN,mlonN, 100,t0, code="A2G")
endcoo = time.time()
print('Coordinate conversion takes in seconds: ', np.round(endcoo-startcoo,3))



print('done.')



print('********* To do: interpolate ovation output correctly on world map')
print()



#***************HIER WEITER  das stimmt nicht, richtig mit dem glatN glonN grid plotten
#west, east, south, north
mapextent=(np.min(glonN), np.max(glonN),np.min(glatN), np.max(glatN))

#** it would be better to interpolate the aurora maps on a full Earth grid and then do the coordinate conversions
#then you would not need to plot a specific part of the image first on mapextent


print()
print()





print()
print('image rescaling ...')

#rescale image

aimg = skimage.transform.resize(fluxN, (512,1024), anti_aliasing=True,mode='constant')



#convert to probabilities from energy flux in erg cm−2 s−1 Case et al. 2016
#pimg=10+8*aimg ???
pimg=8*aimg


#smooth out artefacts
pimg = scipy.ndimage.gaussian_filter(pimg,sigma=(11,11))


#round probabilities
pimg=np.round(pimg,0)


#trim small values
pimg[np.where(pimg <2)]=0

#cut at 100 percent probability
pimg[np.where(pimg >100)]=100


#add calibration of probabilities compared to NOAA
adhoc_calibration=1
pimg=pimg*adhoc_calibration

#plt.close()
#16/9 ration for full hd output
#fig = plt.figure(figsize=[10, 5]) 

#ax1 = plt.subplot(1, 1, 1, projection=crs)
#fig.set_facecolor('black') 
#ax1.coastlines()
#ax1.imshow(pimg, vmin=0, vmax=100, transform=crs, extent=mapextent, origin='lower',zorder=3,cmap=aurora_cmap())
#ax1.set_extent(fullextent)






end = time.time()

print('Time needed:  ',np.round(end - start,2),' sec')

print('Calculation for 120 frames would take:', np.round((end - start)*120/60,2),' minutes')










########################################## Make aurora plot global for comparison with NOAA nowcast

plt.close()
fig = plt.figure(figsize=[15, 10]) 

fig.set_facecolor('black') 

# We choose to plot in an Orthographic projection as it looks natural
# and the distortion is relatively small around the poles where
# the aurora is most likely.

# North America centered
#ax1 = plt.subplot(1, 1, 1, projection=ccrs.Orthographic(-100, 90))

# Europe centered
#ax1 = plt.subplot(1, 2, 1, projection=ccrs.Orthographic(+10, 60))


#PREDSTORM + OVATION
ax1 = plt.subplot(1, 2, 1, projection=ccrs.Orthographic(global_plot_longitude, global_plot_latitude))


#NOAA 
ax2 = plt.subplot(1, 2, 2, projection=ccrs.Orthographic(global_plot_longitude, global_plot_latitude))


#nightmap = 'https://map1c.vis.earthdata.nasa.gov/wmts-geo/wmts.cgi'
#layer = 'VIIRS_CityLights_2012'
#try: urllib.request.urlretrieve(nightmap,'wmts.cgi')
#except urllib.error.URLError as e:
#    print('Failed downloading ', nightmap,' ',e.reason)


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



#load NOAA nowcast
img, crs, extent, origin, dt = oup.aurora_now()


for ax in [ax1,ax2]:

    ax.gridlines(linestyle='--',alpha=0.5)
    #ax.coastlines(alpha=0.5,zorder=3)
    #ax.add_feature(land_50m, color='darkgreen')
    ax.add_feature(land_50m, color='darkslategrey')
    #ax.add_feature(carfeat.LAND,zorder=2,alpha=1)
    ax.add_feature(carfeat.LAKES,color='navy')#,zorder=2,alpha=1)
    #ax.add_feature(carfeat.OCEAN)#,zorder=2,alpha=1)
    ax.add_feature(ocean_50m,linewidth=0.5, color='navy')

    ax.add_feature(carfeat.BORDERS, alpha=0.5)#,zorder=2,alpha=0.5)
    #ax.add_feature(carfeat.COASTLINE)#,zorder=2,alpha=0.5)
    #ax.add_feature(carfeat.RIVERS)#,zorder=2,alpha=0.8)
    ax.add_feature(provinces_50m,alpha=0.5)#,zorder=2,alpha=0.8)
    #ax.stock_img()
    
    #ax.add_wmts(nightmap, layer)
    if ax==ax1: 
        ax.imshow(pimg*5, vmin=0, vmax=100, transform=crs, extent=mapextent, origin='lower', zorder=3, alpha=0.8, cmap=oup.aurora_cmap())
        #contour plot?
        ax.add_feature(Nightshade(t0))

    if ax==ax2: 
        ax.imshow(img*5, vmin=0, vmax=100, transform=crs, extent=extent, origin='lower', zorder=3, alpha=0.8, cmap=oup.aurora_cmap())
        ax.add_feature(Nightshade(dt))

   
    
fig.text(0.01,0.92,'PREDSTORM aurora forecast   '+t0.strftime('%Y-%m-%d %H:%M UT' )+ '                                                            NOAA forecast  '+dt.strftime('%Y-%m-%d %H:%M UT' ), color='white',fontsize=15)
fig.text(0.99,0.02,'C. Möstl / IWF-helio, Austria', color='white',fontsize=8,ha='right')


plt.tight_layout()  



fig.savefig('test/predstorm_aurora_real_Nhemi_'+t0.strftime("%Y_%m_%d_%H%M")  +'.jpg',dpi=120,facecolor=fig.get_facecolor())
plt.show()







sys.exit()












####################################### ZOOM PLOT FOR AMERICA and EUROPE ############################################









plt.close()
#16/9 ration for full hd output
fig = plt.figure(figsize=[16, 9]) 

fig.set_facecolor('black') 

# We choose to plot in an Orthographic projection as it looks natural
# and the distortion is relatively small around the poles where
# the aurora is most likely.

# ax1 Europe
ax1 = plt.subplot(1, 2, 2, projection=ccrs.Orthographic(0, 60),position=[0.51,0.05,0.48,0.9])#[left, bottom, width, height]
# class cartopy.crs.Orthographic(central_longitude=0.0, central_latitude=0.0, globe=None)[source]
# ax2 northern America
ax2 = plt.subplot(1, 2, 1, projection=ccrs.Orthographic(-100, 60), position=[0.01,0.05,0.48,0.9])

'''
# ax2 Europe
ax1 = fig.add_subplot(1, 2, 2, projection=ccrs.PlateCarree())
# class cartopy.crs.Orthographic(central_longitude=0.0, central_latitude=0.0, globe=None)[source]
# ax1 Canada
ax2 = fig.add_subplot(1, 2, 1, projection=ccrs.PlateCarree())
'''

canada_east = -65
canada_west = -135
canada_north = 75
canada_south = 20

europe_east = 50
europe_west = -20
europe_north = 75
europe_south = 30



#nightmap = 'https://map1c.vis.earthdata.nasa.gov/wmts-geo/wmts.cgi'
#layer = 'VIIRS_CityLights_2012'
#try: urllib.request.urlretrieve(nightmap,'wmts.cgi')
#except urllib.error.URLError as e:
#    print('Failed downloading ', nightmap,' ',e.reason)


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

for ax in [ax1, ax2]:


    ax.gridlines(linestyle='--',alpha=0.5)
    #ax.coastlines(alpha=0.5,zorder=3)
    #ax.add_feature(land_50m, color='darkgreen')
    ax.add_feature(land_50m, color='darkslategrey')
    #ax.add_feature(carfeat.LAND,zorder=2,alpha=1)
    ax.add_feature(carfeat.LAKES,color='navy')#,zorder=2,alpha=1)
    #ax.add_feature(carfeat.OCEAN)#,zorder=2,alpha=1)
    ax.add_feature(ocean_50m,linewidth=0.5, color='navy')

    ax.add_feature(carfeat.BORDERS, alpha=0.5)#,zorder=2,alpha=0.5)
    #ax.add_feature(carfeat.COASTLINE)#,zorder=2,alpha=0.5)
    #ax.add_feature(carfeat.RIVERS)#,zorder=2,alpha=0.8)
    ax.add_feature(provinces_50m,alpha=0.5)#,zorder=2,alpha=0.8)
    #ax.stock_img()
    
    #ax.add_wmts(nightmap, layer)
    ax.add_feature(Nightshade(t0))
    ax.imshow(pimg, vmin=0, vmax=100, transform=crs,
    extent=mapextent, origin='lower', zorder=3, alpha=0.9,
    cmap=oup.aurora_cmap2())
    if ax == ax1: ax.set_extent([europe_west, europe_east, europe_south, europe_north])
    if ax == ax2: ax.set_extent([canada_west, canada_east, canada_south, canada_north])


fig.text(0.01,0.92,'PREDSTORM aurora forecast   '+t0.strftime('%Y-%m-%d %H:%M UT' ), color='white',fontsize=15)
fig.text(0.99,0.02,'C. Möstl / IWF-helio, Austria', color='white',fontsize=8,ha='right')

########### ****** add colorbar for probabilities


#exactly full hd resolution with dpi=120 and size 16 9
fig.savefig('forecast/predstorm_aurora_real_'+t0.strftime("%Y_%m_%d_%H%M")  +'.jpg',dpi=120,facecolor=fig.get_facecolor())
plt.show()





#interpolate onto full Earth grid later
#   1024 values covering 0 to 360 degrees in the horizontal (longitude) direction  (0.32846715 degrees/value)
#   512 values covering -90 to 90 degrees in the vertical (latitude) direction  (0.3515625 degrees/value)
#   Values range from 0 (little or no probability of visible aurora) to 100 (high probability of visible aurora)
#mlonN360=mlonN*15     
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html#scipy.interpolate.griddata
#make full Earth a 1024 (lon) x 512 (lat) image
#mlon, mlat = np.meshgrid(np.arange(0,360,360/1024.0 ),np.arange(-90,90,180/512.0))
#for i in np.arange(0,79,1):
#  points=points.append(mlatN,[024/96 ]
#mx, my = np.mgrid[-180:360/1024:180j, -90:180/512:90j]  
#aimg = scipy.interpolate.griddata(points, fluxN, (mx, my), method='cubic')#linear, nearest
#resize image
#aimg = skimage.transform.resize(fluxN, (512,1024), anti_aliasing=False)
#Aurora Specification Tabular Values
# Product: Ovation Aurora Short Term Forecast    <path to tabular data>
# Product Valid At: 2019-04-15 10:15
# Product Generated At: 2019-04-15 09:45
#
# Prepared by the U.S. Dept. of Commerce, NOAA, Space Weather Prediction Center.
# Please send comments and suggestions to SWPC.Webmaster@noaa.gov
#
# Missing Data:  (n/a)
# Cadence:   5 minutes
#
# Tabular Data is on the following grid
#
#   1024 values covering 0 to 360 degrees in the horizontal (longitude) direction  (0.32846715 degrees/value)
#   512 values covering -90 to 90 degrees in the vertical (latitude) direction  (0.3515625 degrees/value)
#   Values range from 0 (little or no probability of visible aurora) to 100 (high probability of visible aurora)

