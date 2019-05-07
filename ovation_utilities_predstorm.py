"""
ovation_utilities rewritten for predstorm input


modified source code taken originally from https://github.com/lkilcommons/OvationPyme

if you change this do 
>> import importlib
>> import ovation_utilities_predstorm  
>> importlib.reload(ovation_utilities_predstorm)  
to use altered functions in main program

"""
import datetime
import numpy as np
import matplotlib.dates as mdates
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import urllib
from urllib.request import urlopen
from io import StringIO
import datetime
#from datetime import datetime
import cartopy.crs as ccrs
import cartopy.feature as carfeat
from cartopy.feature.nightshade import Nightshade
import numba as nb
from numba import njit



def calc_avg_solarwind_predstorm(dt, filein,n_hours):
    """
    Calculates a weighted average of speed and magnetic field
    n_hours (4 by default) backward
    in time from the closest hour in the predstorm forecast
    rewritten from https://github.com/lkilcommons/OvationPyme ovation_utilities.py
    because there the weighting was linearly decreasing, but in Newell et al. 2010
    its 1 0.65 0.65*0.65 ...
    
    input: 
    - datetime object dt
       
    - filein filename and path of solar wind input file
    real time  file
    file='/Users/chris/python/predstorm/predstorm_real.txt'
    sample file
    file='predstorm_sample/predstorm_real.txt'
    time    matplotlib_time B[nT] Bx   By     Bz   N[ccm-3] V[km/s] Dst[nT]   Kp   aurora [GW]
    
    - n_hours  hours previous to integrate over, usually 4
    
    """
    
    l1wind = np.loadtxt(filein)
    
    # Read forecast date and time as matplotlib date number
    l1wind_time=l1wind[:,6] 

    dt_mat=mdates.date2num(dt)  
    #find index of closest time to dt
    closest_time_ind=np.argmin(abs(l1wind_time-dt_mat))
    
    #print('input time ',dt)
    #print('closest time in predstorm',mdates.num2date(l1wind_time[closest_time_ind])     )

    bx, by, bz = l1wind[:,8],l1wind[:,9],l1wind[:,10]
    v,n = l1wind[:,12],l1wind[:,11]
    ec = calc_coupling_predstorm(bx, by, bz, v)

    prev_hour_weight = 0.65    # reduce weighting by factor of wh each hour back
    #make array with weights according to Newell et al. 2010, par 25
    weights=np.ones(1)
    for k in np.arange(1,n_hours,1):
      weights = np.append(weights,weights[k-1]*prev_hour_weight) 


    #print(weights)  
    #print(closest_time_ind)
    times_for_weight_ind = np.arange(closest_time_ind, closest_time_ind-n_hours,-1)
    #print(times_for_weight_ind)
    
    #make list of average solar wind variables
    avgsw = dict()
    #print(bx[times_for_weight_ind])
    #print(v[times_for_weight_ind])

    #pdb.set_trace()
    
    avgsw['Bx'] = np.round(np.nansum(bx[times_for_weight_ind]*weights)/ np.nansum(weights),2)
    avgsw['By'] = np.round(np.nansum(by[times_for_weight_ind]*weights)/ np.nansum(weights),2)
    avgsw['Bz'] = np.round(np.nansum(bz[times_for_weight_ind]*weights)/ np.nansum(weights),2)
    avgsw['V'] = np.round(np.nansum(v[times_for_weight_ind]*weights)/ np.nansum(weights),1)
    avgsw['Ec'] = np.round(np.nansum(ec[times_for_weight_ind]*weights)/ np.nansum(weights),1)

   
    return avgsw


@njit
def calc_coupling_predstorm(Bx, By, Bz, V):
    """
    Empirical Formula for dF/dt - i.e. the Newell coupling
    e.g. paragraph 25 in Newell et al. 2010 doi:10.1029/2009JA014805
    taken from https://github.com/lkilcommons/OvationPyme
    """
    Ec = np.zeros_like(Bx)
    Ec.fill(np.nan)
    B = np.sqrt(Bx**2 + By**2 + Bz**2)
    BT = np.sqrt(By**2 + Bz**2)
    #no 0s allowed in Bz?
    bztemp = Bz
    bztemp[Bz == 0] = 0.001
    #Caculate clock angle (theta_c = t_c)
    tc = np.abs(np.arctan2(By,bztemp))
    sintc = np.sin(tc/2.)
    Ec = (V**1.33333)*(sintc**2.66667)*(BT**0.66667)
    return Ec



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
            dt = datetime.datetime.strptime(line[-17:-1], '%Y-%m-%d %H:%M')

  
    return img, dt



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

def aurora_cmap2 ():
    """Return a colormap with aurora like colors"""
    stops = {'red': [(0.00, 0.1725, 0.1725),
                     (0.50, 0.1725, 0.1725),
                     (1.00, 0.8353, 0.8353)],

             'green': [(0.00, 0.9294, 0.9294),
                       (0.50, 0.9294, 0.9294),
                       (1.00, 0.8235, 0.8235)],

             'blue': [(0.00, 0.3843, 0.3843),
                      (0.50, 0.3843, 0.3843),
                      (1.00, 0.6549, 0.6549)],

             'alpha': [(0.00, 0.0, 0.0),
                       (0.50, 1.0, 1.0),
                       (1.00, 1.0, 1.0)]}

    return LinearSegmentedColormap('aurora', stops)
    
    
    
    

def round_to_hour(dt):
    '''
    round datetime objects to nearest hour
    '''
    dt_start_of_hour = dt.replace(minute=0, second=0, microsecond=0)
    dt_half_hour = dt.replace(minute=30, second=0, microsecond=0)

    if dt >= dt_half_hour:
        # round up
        dt = dt_start_of_hour + datetime.timedelta(hours=1)
    else:
        # round down
        dt = dt_start_of_hour
    return dt    
    
    
    
 

def global_predstorm_noaa(world_image):



 fig = plt.figure(1,figsize=[15, 10]) 
 fig.set_facecolor('black') 

 #axis PREDSTORM + OVATION
 ax1 = plt.subplot(1, 2, 1, projection=ccrs.Orthographic(global_plot_longitude, global_plot_latitude))
 #axis NOAA 
 ax2 = plt.subplot(1, 2, 2, projection=ccrs.Orthographic(global_plot_longitude, global_plot_latitude))
 #load NOAA nowcast
 noaa_img, dt = oup.aurora_now()


 for ax in [ax1,ax2]:

     ax.gridlines(linestyle='--',alpha=0.5,color='white')
     #ax.coastlines(alpha=0.5,zorder=3)
     #ax.add_feature(land_50m, color='darkgreen') 

     #ax.add_feature(land_50m, color='darkslategrey')
     #ax.add_feature(carfeat.LAND,color='darkslategrey')
     #ax.add_feature(carfeat.LAKES,color='navy')#,zorder=2,alpha=1)
     #ax.add_feature(carfeat.OCEAN)#,zorder=2,alpha=1)
     #ax.add_feature(ocean_50m,linewidth=0.5, color='navy')

     #ax.add_feature(carfeat.BORDERS, alpha=0.5)#,zorder=2,alpha=0.5)
     #ax.add_feature(carfeat.COASTLINE)#,zorder=2,alpha=0.5)
     #ax.add_feature(carfeat.RIVERS)#,zorder=2,alpha=0.8)
     #ax.add_feature(provinces_50m,alpha=0.5)#,zorder=2,alpha=0.8)
     ax.stock_img()#alpha=0.2)
    
     #ax.add_wmts(nightmap, layer)
  
     #for testing with black background
     #ax.background_patch.set_facecolor('k')    
     if ax==ax1: 
        ax.imshow(world_image, vmin=0, vmax=100, transform=crs, extent=mapextent, origin='lower', zorder=3, alpha=0.9, cmap=oup.aurora_cmap())
        ax.add_feature(Nightshade(t0))

     if ax==ax2: 
        ax.imshow(noaa_img, vmin=0, vmax=100, transform=crs, extent=mapextent, origin='lower', zorder=3, alpha=0.9, cmap=oup.aurora_cmap())
        ax.add_feature(Nightshade(dt))

   
    
 fig.text(0.01,0.92,'PREDSTORM aurora forecast   '+t0.strftime('%Y-%m-%d %H:%M UT' )+ '                                                            NOAA forecast  '+dt.strftime('%Y-%m-%d %H:%M UT' ), color='white',fontsize=15)
 fig.text(0.99,0.02,'C. Möstl / IWF-helio, Austria', color='white',fontsize=8,ha='right')
 plt.tight_layout()  
 plot_Nhemi_comparison_filename='test/predstorm_aurora_real_Nhemi_'+t0.strftime("%Y_%m_%d_%H%M")  +'.jpg'
 fig.savefig(plot_Nhemi_comparison_filename,dpi=120,facecolor=fig.get_facecolor())
 #plt.show()
 print('Saved image:  ',plot_Nhemi_comparison_filename)

   