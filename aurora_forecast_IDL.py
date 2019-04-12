"""

Plotting the Aurora based on OVATION Prime outputs in IDL

C. Moestl, IWF-helio, Graz, Austria.
twitter @chrisoutofspace



"""


try:
    from urllib2 import urlopen
except ImportError:
    from urllib.request import urlopen

from io import StringIO

import numpy as np
import sys
from datetime import datetime
import cartopy.crs as ccrs
import cartopy.feature as carfeat
from cartopy.feature.nightshade import Nightshade
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import skimage.transform
import scipy
import time


def aurora_forecast():
    """
  
    Returns
    -------
    img : numpy array
        The pixels of the image in a numpy array.
    img_proj : cartopy CRS
        The rectangular coordinate system of the image.
    img_extent : tuple of floats
        The extent of the image ``(x0, y0, x1, y1)`` referenced in
        the ``img_proj`` coordinate system.
    origin : str
        The origin of the image to be passed through to matplotlib's imshow.
    dt : datetime
        Time of forecast validity.

    """




    # GitHub gist to download the example data from
    #url = ('https://gist.githubusercontent.com/belteshassar/'
    #       'c7ea9e02a3e3934a9ddc/raw/aurora-nowcast-map.txt')
    # To plot the current forecast instead, uncomment the following line
    url = 'http://services.swpc.noaa.gov/text/aurora-nowcast-map.txt'

    response_text = StringIO(urlopen(url).read().decode('utf-8'))
    img = np.loadtxt(response_text)
    # Read forecast date and time
    response_text.seek(0)
    for line in response_text:
        if line.startswith('Product Valid At:', 2):
            dt = datetime.strptime(line[-17:-1], '%Y-%m-%d %H:%M')




    img_proj = ccrs.PlateCarree()
    img_extent = (-180, 180, -90, 90)



    return img, img_proj, img_extent, 'lower', dt


def aurora_cmap():
    """Return a colormap with aurora like colors"""
    stops = {'red': [(0.00, 0.1725, 0.1725),
                     (0.50, 0.1725, 0.1725),
                     (1.00, 0.95, 0.95)],

             'green': [(0.00, 0.9294, 0.9294),
                       (0.50, 0.9294, 0.9294),
                       (1.00, 0., 0.)],

             'blue': [(0.00, 0.3843, 0.3843),
                      (0.50, 0.3843, 0.3843),
                      (1.00, 0., 0.)],

             'alpha': [(0.00, 0.0, 0.0),
                       (0.50, 1.0, 1.0),
                       (1.00, 1.0, 1.0)]}

    return LinearSegmentedColormap('aurora', stops)









######################################### main

# url = 'http://services.swpc.noaa.gov/text/aurora-nowcast-map.txt'
# 
# response_text = StringIO(urlopen(url).read().decode('utf-8'))
# img = np.loadtxt(response_text)
# # Read forecast date and time
# response_text.seek(0)
# for line in response_text:
#   if line.startswith('Product Valid At:', 2):
#      dt = datetime.strptime(line[-17:-1], '%Y-%m-%d %H:%M')
# 
# 

start = time.time()


img_proj = ccrs.PlateCarree()
crs=img_proj	
img_extent = (-180, 180, -90, 90)
extent=img_extent
origin='lower'



#**use ovation prime 2013? https://sourceforge.net/projects/ovation-prime/
data_in_file=['ov_diff_Eflux_2017_1230_2330.txt','ov_wave_Eflux_2017_1230_2330.txt','ov_mono_Eflux_2017_1230_2330.txt']

dtovation = datetime.strptime(data_in_file[0][14:18] +'-'+data_in_file[0][19:21]+ 
      '-'+data_in_file[0][21:23]+ ' '+data_in_file[0][24:26]+':'
      +data_in_file[0][26:28] , '%Y-%m-%d %H:%M')     


#ax,ay=np.meshgrid(amlat, amlt)
#full globe
allon=np.arange(0,24,0.25)
allat=np.arange(-89.5,90,0.5)
#array values
latvals=np.arange(50,90,0.5)
aimg=np.zeros([np.size(allat),np.size(allon)]) 

#------------------------------------------------
#sum all aurora contributions on one image
for p in np.arange(0,3):
 print(p)
 data_in=np.genfromtxt('ovation_output/'+data_in_file[p], skip_header=4,max_rows=7680) 
 amlat=data_in[:,1] 
 amlt=data_in[:,0] 
 aint=data_in[:,2]
 #write 1D data on 2D image array in the right coordinate bins
 k=0
 for j in np.arange(len(allon)-1):
  for i in np.arange(279,359,1):
    #sum over all files mono, diff, wave
    aimg[i][j]=aint[k]+aimg[i][j]   
    k=k+1
#----------------------------------------------


#resize image
aimg = skimage.transform.resize(aimg, (512, 1024), anti_aliasing=False)
#convert to probabilities
pimg=10+8*aimg
#trim small values
pimg[np.where(pimg <12)]=0
#cut at 100 percent probability
pimg[np.where(pimg >100)]=100
#smooth out artefacts
pimg = scipy.ndimage.gaussian_filter(pimg,sigma=(9,9))

#TO DO: make sure transition 24 MLT is correct
#*********** coordinate conversion



#get auroral power from files

#------------------------------------




plt.close('all')


#plt.figure(2)
#plt.imshow(pimg,vmin=0, vmax=100,cmap=aurora_cmap())
#plt.imshow(pimg,vmin=0, vmax=100,cmap=aurora_cmap())


plt.figure(1)
plt.plot(pimg) 


fig = plt.figure(2,figsize=[20, 10]) 

# We choose to plot in an Orthographic projection as it looks natural
# and the distortion is relatively small around the poles where
# the aurora is most likely.

# ax1 for Northern Hemisphere
ax1 = fig.add_subplot(1, 2, 1, projection=ccrs.Orthographic(0, 60))
# class cartopy.crs.Orthographic(central_longitude=0.0, central_latitude=0.0, globe=None)[source]

# ax2 for Southern Hemisphere
ax2 = fig.add_subplot(1, 2, 2, projection=ccrs.Orthographic(-100, 60))

#img, crs, extent, origin, dt = aurora_forecast()

fig.set_facecolor((0,0,0))


url = 'https://map1c.vis.earthdata.nasa.gov/wmts-geo/wmts.cgi'
layer = 'VIIRS_CityLights_2012'



for ax in [ax1, ax2]:
 ax.coastlines(zorder=3,color='grey',alpha=0.9)
 ax.add_feature(carfeat.BORDERS,zorder=3,alpha=0.9)
 #ax.stock_img()
 #takes long!!
 #ax.add_wmts(url, layer)
 ax.gridlines(alpha=0.3)
 ax.add_feature(Nightshade(dtovation),alpha=0.9)
 #ax.imshow(img, vmin=0, vmax=100, transform=crs,	extent=extent, origin=origin, zorder=2,		cmap=aurora_cmap())
 ax.imshow(pimg, vmin=0, vmax=100, transform=crs, extent=extent, origin=origin, zorder=2, cmap=aurora_cmap())
 

plt.show()
#plt.tight_layout()
fig.savefig('predstorm_aurora_pred_'+dtovation.strftime("%Y_%m_%d_%H%M")  +'.png',dpi=300)

end = time.time()
print(end - start)