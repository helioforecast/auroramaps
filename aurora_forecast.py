"""
Plotting the PREDSTORM aurora forecast

using ovationpyme by Liam Kilcommons https://github.com/lkilcommons/OvationPyme
C. Moestl, IWF-helio, Graz, Austria.
twitter @chrisoutofspace

"""


import matplotlib
matplotlib.use('Qt5Agg') 

import urllib
from urllib.request import urlopen

from io import StringIO
import ciso8601
import numpy as np
from datetime import datetime
import cartopy.crs as ccrs
import cartopy.feature as carfeat
from cartopy.feature.nightshade import Nightshade
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import sys
import datetime
import skimage.transform
import scipy



from ovationpyme import ovation_prime
from ovationpyme import ovation_utilities
from geospacepy import satplottools, special_datetime



##################################### FUNCTIONS

def aurora_now():
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

    #GitHub gist to download the example data from
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





########################################### Main







#run ovationpyme



#electron energy flux - mono diff wave needed

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

de = ovation_prime.FluxEstimator('diff', jtype)
me = ovation_prime.FluxEstimator('mono', jtype)
we = ovation_prime.FluxEstimator('wave', jtype)


t0 = ciso8601.parse_datetime("2010-04-06 01:00")
#ts = [t0 + datetime.timedelta(hours=i) for i in range(1, 24,1)]
#print(ts)

#for k in np.arange(1,24):
#    print(k)
mlatN, mlonN, fluxN=de.get_flux_for_time(t0, hemi='N')




#resize image
aimg = skimage.transform.resize(fluxN, (512, 1024), anti_aliasing=False)
#convert to probabilities
pimg=10+8*aimg
#trim small values
pimg[np.where(pimg <12)]=0
#cut at 100 percent probability
pimg[np.where(pimg >100)]=100
#smooth out artefacts
pimg = scipy.ndimage.gaussian_filter(pimg,sigma=(9,9))


img_proj = ccrs.PlateCarree()
crs=img_proj	
img_extent = (-180, 180, -90, 90)
extent=img_extent
origin='lower'

plt.imshow(pimg, vmin=0, vmax=100)#, transform=crs)
 





sys.exit()
########################################## Make aurora plot


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



img, crs, extent, origin, dt = aurora_now()

canada_east = -65
canada_west = -135
canada_north = 75
canada_south = 20

europe_east = 50
europe_west = -20
europe_north = 75
europe_south = 30



nightmap = 'https://map1c.vis.earthdata.nasa.gov/wmts-geo/wmts.cgi'
layer = 'VIIRS_CityLights_2012'

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


for ax in [ax1, ax2]:
    if ax == ax1: ax.set_extent([europe_west, europe_east, europe_south, europe_north])
    if ax == ax2: ax.set_extent([canada_west, canada_east, canada_south, canada_north])
    ax.gridlines(linestyle='--',alpha=0.5)
    #ax.coastlines(alpha=0.5,zorder=3)
    ax.add_feature(land_50m)
    #ax.add_feature(carfeat.LAND,zorder=2,alpha=1)
    ax.add_feature(carfeat.LAKES)#,zorder=2,alpha=1)
    #ax.add_feature(carfeat.OCEAN)#,zorder=2,alpha=1)
    ax.add_feature(ocean_50m,linewidth=0.5)

    ax.add_feature(carfeat.BORDERS, alpha=0.5)#,zorder=2,alpha=0.5)
    #ax.add_feature(carfeat.COASTLINE)#,zorder=2,alpha=0.5)
    ax.add_feature(carfeat.RIVERS)#,zorder=2,alpha=0.8)
    ax.add_feature(provinces_50m,alpha=0.5)#,zorder=2,alpha=0.8)
    #ax.stock_img()
  
    
    
    #ax.add_wmts(nightmap, layer)
    ax.add_feature(Nightshade(dt))
    ax.imshow(img, vmin=0, vmax=100, transform=crs,
    extent=extent, origin=origin, zorder=3, alpha=0.9,
    cmap=aurora_cmap())

fig.text(0.01,0.92,'PREDSTORM aurora forecast   '+dt.strftime('%Y-%m-%d %H:%M UT' ), color='white',fontsize=15)
fig.text(0.99,0.02,'C. Möstl / IWF-helio, Austria', color='white',fontsize=8,ha='right')

#exactly full hd resolution with dpi=120 and size 16 9
fig.savefig('forecast/predstorm_aurora_real_'+dt.strftime("%Y_%m_%d_%H%M")  +'.jpg',dpi=120,facecolor=fig.get_facecolor())
plt.show()
