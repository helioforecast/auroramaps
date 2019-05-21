"""
ovation_utilities for predstorm


modified source code taken originally from https://github.com/lkilcommons/OvationPyme

if you change this do 
>> import importlib
>> import ovation_utilities_predstorm  
>> importlib.reload(ovation_utilities_predstorm)  
in main program to use altered functions in main program

"""
import datetime
import numpy as np
import matplotlib
import matplotlib.dates as mdates
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import urllib
from urllib.request import urlopen
from io import StringIO
import datetime
import cartopy.crs as ccrs
import cartopy.feature as carfeat
from cartopy.feature.nightshade import Nightshade
import numba as nb
from numba import njit
import os
import pickle	
import scipy
import pdb






def calc_avg_solarwind_predstorm(dt, filein):
    """
    Calculates a weighted average of speed and magnetic field
    ave_hours (4 by default) backward, for 1 minute data time resolution in filein
    
    input: 
    - datetime object dt, 1min resolution
       
    - filein filename and path of solar wind input file
    contains:
    time    matplotlib_time B[nT] Bx   By     Bz   N[ccm-3] V[km/s] Dst[nT]   Kp   aurora [GW]
    """

    ave_hours=4                #hours previous to integrate over, usually 4
    prev_hour_weight = 0.65    # reduce weighting by factor of wh each hour back

    l1wind = np.loadtxt(filein) #*******save as pickle mit matplotlib times? dauert viel zu lange, gleich datacube übergeben
    
    # Read forecast date and time as matplotlib date number
    l1wind_time=l1wind[:,6] 
    
    #get with input dt to hour start and continute with matplotlib times
    dt_mat_hour=mdates.date2num(round_to_hour_start(dt))
    #find index of closest time to dt_mat_hour
    closest_time_ind_hour=np.argmin(abs(l1wind_time-dt_mat_hour))
    
    dt_mat=mdates.date2num(dt) #convert to matplotlib time
    closest_time_ind=np.argmin(abs(l1wind_time-dt_mat))  #find index of closest time to dt

    bx, by, bz, v = l1wind[:,8],l1wind[:,9],l1wind[:,10],l1wind[:,12]
    ec = calc_coupling_predstorm(bx, by, bz, v)

    #make array for weights, current hour 1, then go down     
    #make array with weights according to Newell et al. 2010, par 25
    weights=np.ones(1)
    for k in np.arange(1,ave_hours+1,1):
         weights = np.append(weights,weights[k-1]*prev_hour_weight) 
      
    #the weight for current time comes in at the beginning of the weights array      
    #difference between hour start and actual time, in hours between 0 and 1; this 
    #is the weight for the current hour (check IDL code for averaging)***
    weights[0]=(dt_mat-dt_mat_hour)*24
   
    #now define the times for the weights to get bxyz; or extract 5 average values for bxyz v first? **
    times_for_weight_ind=np.zeros(5,dtype=int)
    times_for_weight_ind[0]=closest_time_ind
    times_for_weight_ind[1:5] = np.arange(closest_time_ind_hour, closest_time_ind_hour-ave_hours*60,-60)
   
   
    '''    
    print('input time ',dt)
    print('hour start in predstorm',mdates.num2date(l1wind_time[closest_time_ind_hour])     )
    print('all weights are  ', weights)    
    print('time indices for weights are',times_for_weight_ind)
    '''    
    

    #******* HERE NEEDS TO BE AN AVERAGING FOR the FULL HOURS LIKE OMNI2, before! 1 time only
    
    
    #make array of average solar wind variables
    avgsw=np.recarray(1,dtype=[('bx', float), ('by', float),('bz', float),('v', float),('ec', float)])
    
    #print(bx[times_for_weight_ind])
    #print(v[times_for_weight_ind])
    
    #sum over last 4 hours with weighting
    avgsw.bx = np.round(np.nansum(bx[times_for_weight_ind]*weights)/ np.nansum(weights),2)
    avgsw.by = np.round(np.nansum(by[times_for_weight_ind]*weights)/ np.nansum(weights),2)
    avgsw.bz = np.round(np.nansum(bz[times_for_weight_ind]*weights)/ np.nansum(weights),2)
    avgsw.v = np.round(np.nansum(v[times_for_weight_ind]*weights)/ np.nansum(weights),1)
    avgsw.ec = np.round(np.nansum(ec[times_for_weight_ind]*weights)/ np.nansum(weights),1)
    
        
   
    return avgsw
 


@njit
def calc_coupling_predstorm(Bx, By, Bz, V):
    """
    Empirical Formula for dF/dt - i.e. the Newell coupling
    e.g. paragraph 25 in Newell et al. 2010 doi:10.1029/2009JA014805
    slightly rewritten from https://github.com/lkilcommons/OvationPyme
    """
    Ec = np.zeros_like(Bx)
    Ec.fill(np.nan)
    B = np.sqrt(Bx**2 + By**2 + Bz**2)
    BT = np.sqrt(By**2 + Bz**2)
    #no 0s allowed in Bz
    bztemp = Bz
    bztemp[Bz == 0] = 0.001
    #Caculate clock angle (theta_c = t_c)
    tc = np.abs(np.arctan2(By,bztemp))
    sintc = np.sin(tc/2.)
    Ec = (V**1.33333)*(sintc**2.66667)*(BT**0.66667)
    return Ec






def get_wind():
  '''
  get 1 min resolution Wind data from pickle file
  '''
  



def global_predstorm_north(wic,dt,colormap_input):
 '''
 wic is a world image cube with ovation results 512x1024
 dt are the datetimes for each frame
 colormap_input is the colormap


 '''
 
  #plotting parameters
 #-100 for North America, +10 or 0 for Europe
 #global_plot_longitude=-100
 global_plot_longitude=0
 global_plot_latitude=90

 #define extent of the produced world maps - defined as: west east south north
 mapextent=[-180,180,-90,90]   
 
 
 
 #need to set alpha linearly increasing in the colormap so that it blends into background for 
 #small values
 cmap = plt.get_cmap(colormap_input)  # Choose colormap
 my_cmap = cmap(np.arange(cmap.N))  # Get the colormap colors
 my_cmap[:,-1] = np.linspace(0, 1, cmap.N)  # Set alpha
 my_cmap = ListedColormap(my_cmap) # Create new colormap

 crs=ccrs.PlateCarree()

 fig = plt.figure(1,figsize=[12, 12],dpi=80) 
 fig.set_facecolor('black') 
 fig.text(0.99,0.01,'C. Möstl / IWF-helio, Austria', color='white',fontsize=10,ha='right',va='bottom')

 
 ax = plt.subplot(1, 1, 1, projection=ccrs.Orthographic(global_plot_longitude, global_plot_latitude))

 ax.background_patch.set_facecolor('k')    
 ax.coastlines('50m',color='white',alpha=0.8)
 gl=ax.gridlines(linestyle='--',alpha=0.5,color='white') #make grid
 gl.n_steps=100   #make grid finer
 ax.stock_img()  

 
 #these are calls that create the first object to be removed from the plot with each frame
 txt=fig.text(0.5,0.92,'')
 border=ax.add_feature(Nightshade(dt[0]))  #add day night border
 img=ax.imshow(wic[:,:,0], cmap=my_cmap)
  

 for i in np.arange(0,np.size(dt)):

     print('global movie frame',i)
 
     #plt.cla()  #clear axes
     txt.set_visible(False)  #clear previous plot title
     #plot title with time
     txt=fig.text(0.5,0.92,'PREDSTORM aurora forecast  '+dt[i].strftime('%Y-%m-%d %H:%M UT'), color='white',fontsize=15, ha='center')
     img.remove()     #remove previous wic
     border.remove()  #remove previous nightshade
 
     
     border=ax.add_feature(Nightshade(dt[i]))  #add day night border
     img=ax.imshow(wic[:,:,i], vmin=0.01, vmax=5, transform=crs, extent=mapextent, origin='lower', zorder=3, alpha=0.8, cmap=my_cmap)
  
     #save as image with timestamp in filename
     #plot_Nhemi_filename='results/forecast_global/predstorm_aurora_real_Nhemi_'+dt.strftime("%Y_%m_%d_%H%M")  +'.jpg'
     #fig.savefig(plot_Nhemi_filename,dpi=150,facecolor=fig.get_facecolor())

     #save as movie frame
     framestr = '%05i' % (i)  
     fig.savefig('results/frames_global/aurora_'+framestr+'.jpg',dpi=150,facecolor=fig.get_facecolor())
     
     #print('Saved image:  ',plot_Nhemi_filename)
     
     
 #matplotlib.use('Qt5Agg')





















#version for 1 hour data
def calc_avg_solarwind_predstorm2(dt, filein):
    """
    Calculates a weighted average of speed and magnetic field
    ave_hours (4 by default) backward
    in time from the closest hour in the predstorm forecast
    rewritten from https://github.com/lkilcommons/OvationPyme ovation_utilities.py
    because there the weighting was linearly decreasing, but in Newell et al. 2010
    its 1 0.65 0.65*0.65 ...
    
    input: 
    - datetime object dt
       
    - filein filename and path of solar wind input file
    real time  file
    file='/Users/chris/python/predstorm/predstorm_real.txt'
    contains:
    time    matplotlib_time B[nT] Bx   By     Bz   N[ccm-3] V[km/s] Dst[nT]   Kp   aurora [GW]
    """
    #hours previous to integrate over, usually 4
    ave_hours=4  
    
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
    for k in np.arange(1,ave_hours,1):
      weights = np.append(weights,weights[k-1]*prev_hour_weight) 

    #print(weights)  
    #print(closest_time_ind)
    times_for_weight_ind = np.arange(closest_time_ind, closest_time_ind-ave_hours,-1)
    print(times_for_weight_ind)
    
    #make array of average solar wind variables
    avgsw=np.recarray(1,dtype=[('bx', float), ('by', float),('bz', float),('v', float),('ec', float)])
    
    #print(bx[times_for_weight_ind])
    #print(v[times_for_weight_ind])
    
    avgsw.bx = np.round(np.nansum(bx[times_for_weight_ind]*weights)/ np.nansum(weights),2)
    avgsw.by = np.round(np.nansum(by[times_for_weight_ind]*weights)/ np.nansum(weights),2)
    avgsw.bz = np.round(np.nansum(bz[times_for_weight_ind]*weights)/ np.nansum(weights),2)
    avgsw.v = np.round(np.nansum(v[times_for_weight_ind]*weights)/ np.nansum(weights),1)
    avgsw.ec = np.round(np.nansum(ec[times_for_weight_ind]*weights)/ np.nansum(weights),1)
   
    return avgsw
 
    

   


      
    
    
    
    


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
    
 

def round_to_minute(dt):
    '''
    round datetime objects to nearest minute
    '''
    dt_start_of_min = dt.replace(second=0, microsecond=0)
    dt_half_min = dt.replace(second=30, microsecond=0)

    if dt >= dt_half_min:
        # round up
        dt = dt_start_of_min + datetime.timedelta(minutes=1)
    else:
        # round down
        dt = dt_start_of_min
    return dt        
    
 
def round_to_hour_start(dt):
    '''
    round datetime objects to start of the current hour
    '''
    dt_start_of_hour = dt.replace(minute=0, second=0, microsecond=0)
    return dt_start_of_hour
        
    
    

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


  


def omni_txt_generator(dt):

   '''
   returns solar wind data in predstorm format from OMNI2 based on datetime array dt
   starting 24 hours earlier as previous averages are used later in calc_avg_solarwind_predstorm
   '''
 
   #load all omni data
   o=omni_loader()  
   
   
   #convert to matplotlib time    
   dt_mat=mdates.date2num(dt) 
       
   #starting index of dt start time in all omni data - 24 hours 
   #needed for averaging solar wind before time dt[0]
   inds=np.argmin(abs(o.time-dt_mat[0]))-24
   
   #end index for dt in omni data
   inde=inds+24+np.size(dt)
   
   #print('omni_txt_generator')
   #print(dt[0])
   #print(dt[-1])
   
   o_time=o.time[inds:inde]
   
   vartxtout=np.zeros([inde-inds,13])

    #get date in ascii
   for i in np.arange(np.size(o_time)):
       time_dummy=mdates.num2date(o_time[i])
       vartxtout[i,0]=time_dummy.year
       vartxtout[i,1]=time_dummy.month
       vartxtout[i,2]=time_dummy.day
       vartxtout[i,3]=time_dummy.hour
       vartxtout[i,4]=time_dummy.minute
       vartxtout[i,5]=time_dummy.second

   vartxtout[:,6]=o.time[inds:inde]
   vartxtout[:,7]=o.btot[inds:inde]
   vartxtout[:,8]=o.bx[inds:inde]
   vartxtout[:,9]=o.bygsm[inds:inde]
   vartxtout[:,10]=o.bzgsm[inds:inde]
   vartxtout[:,11]=o.den[inds:inde]
   vartxtout[:,12]=o.speed[inds:inde]
   #vartxtout[:,13]=dst_temerin_li
   #vartxtout[:,14]=kp_newell
   #vartxtout[:,15]=aurora_power

   #description
   #np.savetxt(filename_save, ['time     Dst [nT]     Kp     aurora [GW]   B [nT]    Bx [nT]     By [nT]     Bz [nT]    N [ccm-3]   V [km/s]    '])
   filename_save='predstorm_omni.txt'
   np.savetxt(filename_save, vartxtout,  delimiter='',fmt='%4i %2i %2i %2i %2i %2i %10.6f %5.1f %5.1f %5.1f %5.1f   %7.0i %7.0i ', \
               header='        time      matplotlib_time B[nT] Bx   By     Bz   N[ccm-3] V[km/s] ')

   #with last 3 variables
   #np.savetxt(filename_save, vartxtout, delimiter='',fmt='%4i %2i %2i %2i %2i %2i %10.6f %5.1f %5.1f %5.1f %5.1f   %7.0i %7.0i   %5.0f %5.1f %5.1f', \
   #            header='        time      matplotlib_time B[nT] Bx   By     Bz   N[ccm-3] V[km/s] Dst[nT]   Kp   aurora [GW]')



    


def omni_loader():
   '''
   downloads all omni2 data into data folder, converts to pickle and returns object
   '''

   if os.path.isdir('data') == False: os.mkdir('data')
  
   if not os.path.exists('data/omni2_all_years.dat'):
      #see http://omniweb.gsfc.nasa.gov/html/ow_data.html
      print('download OMNI2 data from')
      omni2_url='ftp://nssdcftp.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_all_years.dat'
      print(omni2_url)
      try: urllib.request.urlretrieve(omni2_url, 'data/omni2_all_years.dat')
      except urllib.error.URLError as e:
          print(' ', omni2_url,' ',e.reason)

   #if omni2 hourly data is not yet converted and saved as pickle, do it:
   if not os.path.exists('data/omni2_all_years_pickle.p'):
       #load OMNI2 dataset from .dat file with a function from dst_module.py
       o=get_omni_data()
       #contains: o. time,day,hour,btot,bx,by,bz,bygsm,bzgsm,speed,speedx,den,pdyn,dst,kp
       #save for faster loading later
       pickle.dump(o, open('data/omni2_all_years_pickle.p', 'wb') )

   else:  o=pickle.load(open('data/omni2_all_years_pickle.p', 'rb') )
   print('loaded OMNI2 data')
   return o


def get_omni_data():
    """FORMAT(2I4,I3,I5,2I3,2I4,14F6.1,F9.0,F6.1,F6.0,2F6.1,F6.3,F6.2, F9.0,F6.1,F6.0,2F6.1,F6.3,2F7.2,F6.1,I3,I4,I6,I5,F10.2,5F9.2,I3,I4,2F6.1,2I6,F5.1)
    1963   1  0 1771 99 99 999 999 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 999.9 9999999. 999.9 9999. 999.9 999.9 9.999 99.99 9999999. 999.9 9999. 999.9 999.9 9.999 999.99 999.99 999.9  7  23    -6  119 999999.99 99999.99 99999.99 99999.99 99999.99 99999.99  0   3 999.9 999.9 99999 99999 99.9

    define variables from OMNI2 dataset
    see http://omniweb.gsfc.nasa.gov/html/ow_data.html

    omni2_url='ftp://nssdcftp.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_all_years.dat'
    """

    #check how many rows exist in this file
    f=open('data/omni2_all_years.dat')
    dataset= len(f.readlines())
    #print(dataset)
    #global Variables
    spot=np.zeros(dataset) 
    btot=np.zeros(dataset) #floating points
    bx=np.zeros(dataset) #floating points
    by=np.zeros(dataset) #floating points
    bz=np.zeros(dataset) #floating points
    bzgsm=np.zeros(dataset) #floating points
    bygsm=np.zeros(dataset) #floating points

    speed=np.zeros(dataset) #floating points
    speedx=np.zeros(dataset) #floating points
    speed_phi=np.zeros(dataset) #floating points
    speed_theta=np.zeros(dataset) #floating points

    dst=np.zeros(dataset) #float
    kp=np.zeros(dataset) #float

    den=np.zeros(dataset) #float
    pdyn=np.zeros(dataset) #float
    year=np.zeros(dataset)
    day=np.zeros(dataset)
    hour=np.zeros(dataset)
    t=np.zeros(dataset) #index time
    
    
    j=0
    print('Read OMNI2 data ...')
    with open('data/omni2_all_years.dat') as f:
        for line in f:
            line = line.split() # to deal with blank 
            #print line #41 is Dst index, in nT
            dst[j]=line[40]
            kp[j]=line[38]
            
            if dst[j] == 99999: dst[j]=np.NaN
            #40 is sunspot number
            spot[j]=line[39]
            #if spot[j] == 999: spot[j]=NaN

            #25 is bulkspeed F6.0, in km/s
            speed[j]=line[24]
            if speed[j] == 9999: speed[j]=np.NaN
          
            #get speed angles F6.1
            speed_phi[j]=line[25]
            if speed_phi[j] == 999.9: speed_phi[j]=np.NaN

            speed_theta[j]=line[26]
            if speed_theta[j] == 999.9: speed_theta[j]=np.NaN
            #convert speed to GSE x see OMNI website footnote
            speedx[j] = - speed[j] * np.cos(np.radians(speed_theta[j])) * np.cos(np.radians(speed_phi[j]))



            #9 is total B  F6.1 also fill ist 999.9, in nT
            btot[j]=line[9]
            if btot[j] == 999.9: btot[j]=np.NaN

            #GSE components from 13 to 15, so 12 to 14 index, in nT
            bx[j]=line[12]
            if bx[j] == 999.9: bx[j]=np.NaN
            by[j]=line[13]
            if by[j] == 999.9: by[j]=np.NaN
            bz[j]=line[14]
            if bz[j] == 999.9: bz[j]=np.NaN
          
            #GSM
            bygsm[j]=line[15]
            if bygsm[j] == 999.9: bygsm[j]=np.NaN
          
            bzgsm[j]=line[16]
            if bzgsm[j] == 999.9: bzgsm[j]=np.NaN    
          
          
            #24 in file, index 23 proton density /ccm
            den[j]=line[23]
            if den[j] == 999.9: den[j]=np.NaN
          
            #29 in file, index 28 Pdyn, F6.2, fill values sind 99.99, in nPa
            pdyn[j]=line[28]
            if pdyn[j] == 99.99: pdyn[j]=np.NaN      
          
            year[j]=line[0]
            day[j]=line[1]
            hour[j]=line[2]
            j=j+1     
      

    #convert time to matplotlib format
    #http://docs.sunpy.org/en/latest/guide/time.html
    #http://matplotlib.org/examples/pylab_examples/date_demo2.html

    times1=np.zeros(len(year)) #datetime time
    print('convert time start')
    for index in range(0,len(year)):
        #first to datetimeobject 
        timedum=datetime.datetime(int(year[index]), 1, 1) + datetime.timedelta(day[index] - 1) +datetime.timedelta(hours=hour[index])
        #then to matlibplot dateformat:
        times1[index] = mdates.date2num(timedum)
    print('convert time done')   #for time conversion

    print('all done.')
    print(j, ' datapoints')   #for reading data from OMNI file
    
    #make structured array of data
    omni_data=np.rec.array([times1,btot,bx,by,bz,bygsm,bzgsm,speed,speedx,den,pdyn,dst,kp], \
    dtype=[('time','f8'),('btot','f8'),('bx','f8'),('by','f8'),('bz','f8'),\
    ('bygsm','f8'),('bzgsm','f8'),('speed','f8'),('speedx','f8'),('den','f8'),('pdyn','f8'),('dst','f8'),('kp','f8')])
    
    return omni_data

 
 
 
 


def global_ovation_flux(magnetic_latitude,magnetic_local_time,flux,dt):


 
 #cmap_name='hot'
 cmap_name='jet'
 
 
 res_lat=40/80
 res_lon=360/96
 
 
 #read in IDL output for comparison
 idl_file_in='ovation_output/ov_diff_Eflux_2017_1230_2330.txt'
 #idl_file_in='ovation_output/ov_mono_Eflux_2017_1230_2330.txt'
 #idl_file_in='ovation_output/ov_wave_Eflux_2017_1230_2330.txt'
 
 #idl_file_in='ovation_output/ov_diff_Eflux_2017_1129_1300.txt'
 #idl_file_in='ovation_output/ov_mono_Eflux_2017_1129_1300.txt'

 
 ovaidl=np.loadtxt(idl_file_in,max_rows=7680)  
     
 mltN_idl=ovaidl[:,0]*15*(np.pi/180)
 mlatN_idl=-ovaidl[:,1] 
 flux_idl=ovaidl[:,2] 
 
 mltN_idl=mltN_idl.reshape(np.size(flux_idl),1)  
 mlatN_idl=mlatN_idl.reshape(np.size(flux_idl),1) 
     
 flux_idl=flux_idl.reshape(np.size(flux_idl),1)
 m2D_idl=np.hstack((mltN_idl, mlatN_idl))
   
 #idl data theta, idl data radial:
 it,ir=np.meshgrid(np.arange(0,np.pi*2+1,res_lon*np.pi/180),np.arange(-90,-50,res_lat) )

 
 #interpolate onto 2D grid 
 idlimg=np.squeeze(scipy.interpolate.griddata(m2D_idl, flux_idl,  (it, ir), method='linear',fill_value=0))
 idlimg[np.where(np.isfinite(idlimg)==False)]=0 

 
 ################ ovationpyme output
 
 #first convert so that suitable for polar plots: theta= longitude, r=latitude
 #magnetic local time to radians
 mltN=magnetic_local_time*15*(np.pi/180) 
 #invert latitude which is the radial axis so that 90N is in center of plot
 mlatN=-magnetic_latitude
 
 #make 1-D array and then 2D array of coordinates for image interpolation later
 mltN_1D=mltN.reshape(np.size(flux),1)  
 mlatN_1D=mlatN.reshape(np.size(flux),1)  
 m2D=np.hstack((mltN_1D, mlatN_1D))
 #make flux a 1-D array
 flux_1D=flux.reshape(np.size(flux),1) 

 #make grid for python ovation version
 pt,pr=np.meshgrid(np.arange(mltN.min(),mltN.max(),res_lon*np.pi/180),np.arange(mlatN.min(),mlatN.max(),res_lat))
 #interpolate onto new grid
 pyimg=np.squeeze(scipy.interpolate.griddata(m2D, flux_1D, (pt, pr), method='linear',fill_value=0))
 #####################################
 
 plt.close('all')
 plt.figure(1,figsize=[16,8])
 plt.suptitle('OVATION aurora energy flux  '+dt.strftime('%Y-%m-%d %H:%M UT'),fontsize=15)
 

 ax1=plt.subplot(121,projection='polar')
 ax1.set_title('python') 
 cs=ax1.contourf(pt, pr, pyimg, cmap=cmap_name, vmin=0,vmax=2.0,levels=100,zorder=0)
 ax1.set_facecolor('black')
 ax1.set_rlim(-90,-50)
 plt.rgrids((-90,-85,-80,-75,-70,-65,-60,-55,-50),('90N','85','80','75','70','65','60','55','50'),angle=150, fontsize=12, color='white')
 plt.thetagrids((0,45,90,135,180,225,270,315),('0','3','6','9','12','15','18','21'), fontsize=12, color='black')
 ax1.set_theta_zero_location('S')
 plt.colorbar(cs,fraction=0.03, pad=0.1)
 ax1.grid(color='white', linestyle='--', linewidth=0.5,zorder=2) 


 ax2=plt.subplot(122,projection='polar')
 ax2.set_title('IDL') 
 cs=ax2.contourf(it, ir, idlimg, cmap=cmap_name, vmin=0,vmax=2.0,levels=100,zorder=0)
 ax2.set_facecolor('black')  
 ax2.set_rlim(-90,-50)
 plt.rgrids((-90,-85,-80,-75,-70,-65,-60,-55,-50),('90N','85','80','75','70','65','60','55','50'),angle=150, fontsize=12, color='white')
 plt.thetagrids((0,45,90,135,180,225,270,315),('0','3','6','9','12','15','18','21'), fontsize=12, color='black')
 ax2.set_theta_zero_location('S')
 plt.colorbar(cs,fraction=0.03, pad=0.1)
 ax2.grid(color='white', linestyle='--', linewidth=0.5,zorder=1) 

 print('Maximum flux python',np.round(np.max(flux_1D),2))
 print('Maximum flux IDL',np.round(np.max(flux_idl),2))
 
 plt.tight_layout()
 plt.savefig('flux_test.png')











def europe_canada_predstorm(world_image,dt,counter,colormap_input):


 plt.close('all')
 #Night map from VIIRS
 #https://wiki.earthdata.nasa.gov/display/GIBS/GIBS+Available+Imagery+Products#expand-EarthatNight4Products
 #nightmap = 'https://map1c.vis.earthdata.nasa.gov/wmts-geo/wmts.cgi'


 #define extent of the produced world maps - defined as: west east south north
 mapextent=[-180,180,-90,90]   


 #not so good: layer='VIIRS_Black_Marble'
 #better but takes time
 #layer = 'VIIRS_CityLights_2012'

 
 #need to set alpha linearly increasing in the colormap
 cmap = plt.get_cmap(colormap_input)  # Choose colormap
 my_cmap = cmap(np.arange(cmap.N))  # Get the colormap colors
 my_cmap[:,-1] = np.linspace(0.3, 1, cmap.N)  # Set alpha
 my_cmap = ListedColormap(my_cmap) # Create new colormap


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



 #16/9 ration for full hd output
 fig = plt.figure(2,figsize=[16, 9],dpi=100) 
 fig.set_facecolor('black') 
 plt.cla()

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
  
    ax.imshow(world_image, vmin=0.001, vmax=5, transform=crs, extent=mapextent, origin='lower', zorder=3, alpha=0.8, cmap=my_cmap)
    
    if ax == ax1: ax.set_extent([europe_west, europe_east, europe_south, europe_north])
    if ax == ax2: ax.set_extent([canada_west, canada_east, canada_south, canada_north])

 fig.text(0.01,0.92,'PREDSTORM aurora forecast   '+dt.strftime('%Y-%m-%d %H:%M UT' ), color='white',fontsize=15)
 fig.text(0.99,0.02,'C. Möstl / IWF-helio, Austria', color='white',fontsize=8,ha='right')

 
 #save as image with timestamp in filename
 plot_europe_canada_filename='results/forecast_europe_canada/predstorm_aurora_real_'+dt.strftime("%Y_%m_%d_%H%M")  +'.jpg'
 fig.savefig(plot_europe_canada_filename,dpi=120,facecolor=fig.get_facecolor())


 #save as movie frame
 framestr = '%05i' % (counter)  
 fig.savefig('results/frames_europe_canada/aurora_'+framestr+'.jpg',dpi=120,facecolor=fig.get_facecolor())
 #plt.show()
 print('Saved image:  ',plot_europe_canada_filename)




'''

def global_predstorm_north2(world_image,dt,counter,colormap_input):


 #plotting parameters
 #-100 for North America, +10 or  for Europe
 #global_plot_longitude=-100
 global_plot_longitude=0
 global_plot_latitude=90


 #Night map from VIIRS
 #https://wiki.earthdata.nasa.gov/display/GIBS/GIBS+Available+Imagery+Products#expand-EarthatNight4Products
 nightmap = 'https://map1c.vis.earthdata.nasa.gov/wmts-geo/wmts.cgi'


 #define extent of the produced world maps - defined as: west east south north
 mapextent=[-180,180,-90,90]   
 
 #need to set alpha linearly increasing in the colormap
 cmap = plt.get_cmap(colormap_input)  # Choose colormap
 my_cmap = cmap(np.arange(cmap.N))  # Get the colormap colors
 my_cmap[:,-1] = np.linspace(0.3, 1, cmap.N)  # Set alpha
 my_cmap = ListedColormap(my_cmap) # Create new colormap



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

 fig = plt.figure(1,figsize=[12, 12],dpi=80) 
 fig.set_facecolor('black') 
 plt.clf()
 plt.cla()

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
 #ax.add_feature(provinces_50m,alpha=0.8,color='white')#,zorder=2,alpha=0.8)

 #ax.add_feature(provinces_50m,alpha=0.5)#,zorder=2,alpha=0.8)


 gl=ax.gridlines(linestyle='--',alpha=0.5,color='white') #make grid
 gl.n_steps=100   #make grid finer
 
 ax.add_feature(Nightshade(dt))  #add day night border

 
 ax.stock_img()
 #ax.add_wmts(nightmap, layer)
 
  

 #ax.imshow(world_image, vmin=0, vmax=5, transform=crs, extent=mapextent, origin='lower', zorder=3, alpha=0.9, cmap=aurora_cmap())
 ax.imshow(world_image, vmin=0.01, vmax=5, transform=crs, extent=mapextent, origin='lower', zorder=3, alpha=0.8, cmap=my_cmap)


 
 fig.text(0.5,0.92,'PREDSTORM aurora forecast  '+dt.strftime('%Y-%m-%d %H:%M UT'), color='white',fontsize=15, ha='center')
 fig.text(0.99,0.01,'C. Möstl / IWF-helio, Austria', color='white',fontsize=10,ha='right',va='bottom')

 #plt.tight_layout()  

 #save as image with timestamp in filename
 plot_Nhemi_filename='results/forecast_global/predstorm_aurora_real_Nhemi_'+dt.strftime("%Y_%m_%d_%H%M")  +'.jpg'
 fig.savefig(plot_Nhemi_filename,dpi=150,facecolor=fig.get_facecolor())

 #save as movie frame
 framestr = '%05i' % (counter)  
 fig.savefig('results/frames_global/aurora_'+framestr+'.jpg',dpi=150,facecolor=fig.get_facecolor())
 #plt.show()
 print('Saved image:  ',plot_Nhemi_filename)



'''



   