"""
ovation_utilities for predstorm

includes a few functions taken from https://github.com/lkilcommons/OvationPyme

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
from numba import njit, jit
import os
import pickle	
import scipy
import pdb





#################################### MATH ##########################################################




def calc_avg_solarwind_predstorm(dt,l1wind):
    """
    Calculates a weighted average of speed and magnetic field bx, by, bz and the Newell coupling ec
    for a number of ave_hours (4 by default) back in time
    input data time resolution in the l1wind array is 1 hour
    aurora output time resolution as given by dt can be higher
    corresponds roughly to ap_inter_sol.pro in IDL ovation
    
    input: 
    - datetime object dt (single value)
    - l1wind: the solar wind input is a recarray containing time 
      (in matplotlib format), bx, by, bz, v, ec
      
    test run time after running aurora_forecast.py with 
    %timeit oup.calc_avg_solarwind_predstorm(ts[0],l1wind)         
    """
    ave_hours=4                # hours previous to integrate over, usually 4
    prev_hour_weight = 0.65    # reduce weighting by factor with each hour back
    
    
    #initiate array of averaged solar wind variables with size of dt
    avgsw=np.recarray(np.size(dt),dtype=[('time', float),('bx', float), ('by', float),('bz', float),('v', float),('ec', float)])
    
    #make array out of dt if its scalar so subscripts work
    if np.size(dt) == 1: dt=[dt] 

    avgsw.time = mdates.date2num(dt)
    
  
    #go through all dt times:
    for i in np.arange(0,np.size(dt)):

       
        dt_mat_hour = mdates.date2num(round_to_hour_start(dt[i]))  #get with input dt to hour start and continute with matplotlib times
        closest_time_ind_hour = np.argmin(abs(l1wind.time-dt_mat_hour))  #find index of closest time to dt_mat_hour
        dt_mat = mdates.date2num(dt[i])  #convert input datetimte dt to matplotlib time
        closest_time_ind = np.argmin(abs(l1wind.time-dt_mat)) #find index of closest time to dt
    
        weights = np.ones(ave_hours)  #make array with weights     
        for k in np.arange(1,ave_hours,1):  weights[k] = weights[k-1]*prev_hour_weight

        a = (dt_mat-dt_mat_hour)*24 # fraction of hour elapsed
        #weights[0]=np.sqrt(a)  #the current fraction of hour is weighted as square root (IDL ap_inter_sol.pro)
        weights[0] = a  #the current fraction of hour is weighted as linear

        times_for_weight_ind = np.arange(closest_time_ind_hour, closest_time_ind_hour-ave_hours,-1)
        
        #sum last hours with each weight and normalize    
        avgsw[i].bx = np.round(np.nansum(l1wind.bx[times_for_weight_ind]*weights)/ np.nansum(weights),2)
        avgsw[i].by = np.round(np.nansum(l1wind.by[times_for_weight_ind]*weights)/ np.nansum(weights),2)
        avgsw[i].bz = np.round(np.nansum(l1wind.bz[times_for_weight_ind]*weights)/ np.nansum(weights),2)
        avgsw[i].v  = np.round(np.nansum(l1wind.v[times_for_weight_ind]*weights)/ np.nansum(weights),1)
        avgsw[i].ec = np.round(np.nansum(l1wind.ec[times_for_weight_ind]*weights)/ np.nansum(weights),1)

    ''' #for debugging
    #print(weights)
    #print(mdates.num2date(l1wind.time[times_for_weight_ind]))
 
    print('input wind')
    print('by',l1wind.by[times_for_weight_ind])
    print('bz',l1wind.bz[times_for_weight_ind])
    print('v',np.round(l1wind.v[times_for_weight_ind]))
    print('ec',np.round(l1wind.ec[times_for_weight_ind]))
    print()
   
    print('output wind')
    print('by',avgsw.by[0])
    print('bz',avgsw.bz[0])
    print('v',avgsw.v[0])
    print('ec',avgsw.ec[0])
    print()
    print('---------------')
    '''
   
    return avgsw




@njit
def calc_coupling_newell(by, bz, v):
    '''    
    Empirical Formula for dFlux/dt - the Newell coupling
    e.g. paragraph 25 in Newell et al. 2010 doi:10.1029/2009JA014805
    IDL ovation: sol_coup.pro - contains 33 coupling functions in total
    input: needs arrays for by, bz, v    
    ''' 
    bt = np.sqrt(by**2 + bz**2)
    bztemp = bz
    bztemp[bz == 0] = 0.001
    tc = np.arctan2(by,bztemp)     #calculate clock angle (theta_c = t_c)
    neg_tc = bt*np.cos(tc)*bz < 0  #similar to IDL code sol_coup.pro
    tc[neg_tc] = tc[neg_tc] + np.pi
    sintc = np.abs(np.sin(tc/2.))
    ec = (v**1.33333)*(sintc**2.66667)*(bt**0.66667)
    
    return ec
   
  

def round_to_hour(dt):
    '''
    round datetime objects to nearest hour
    '''
    #if not a list, change to list
    if type(dt)!=list: 
         dt=[dt]
    
    for i in np.arange(0,np.size(dt)):
       dt_start_of_hour = dt[i].replace(minute=0, second=0, microsecond=0)
       dt_half_hour = dt[i].replace(minute=30, second=0, microsecond=0)

       if dt[i] >= dt_half_hour:
          # round up
           dt[i] = dt_start_of_hour + datetime.timedelta(hours=1)
       else:
          # round down
           dt[i] = dt_start_of_hour


    return dt    
 
 
def round_to_hour_start(dt):
    '''
    round datetime objects to start of the current hour
    '''
    #if not a list, change to list
    if type(dt)!=list: 
        dt=[dt]

    for i in np.arange(0,np.size(dt)):
        dt_start_of_hour = dt[i].replace(minute=0, second=0, microsecond=0)
    return dt_start_of_hour
        


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
















########################################### DATA HANDLING ######################################



def load_predstorm_wind(file_in):
    '''
    loads the predstorm input file
    and makes a recarray of matplotlib time, Bxyz, speed and the Newell coupling
    calls calc_coupling_newell()
    '''
    l1 = np.loadtxt(file_in)       # load data file
    # define array
    wind_array=np.recarray(len(l1),dtype=[('time', float),('bx', float), ('by', float),('bz', float),('v', float),('ec', float)])
    # load parameters from file into recarray
    wind_array.time, wind_array.bx, wind_array.by, wind_array.bz, wind_array.v = l1[:,6],l1[:,8],l1[:,9],l1[:,10],l1[:,12]
    # get Newell coupling, only by, bz needed 
    wind_array.ec = calc_coupling_newell(wind_array.by, wind_array.bz, wind_array.v)    
    
    return wind_array     #return array of average solar wind variables
     

def omni_txt_generator(dt):

   '''
   returns solar wind data in predstorm format from OMNI2 based on datetime array dt
   starting 24 hours earlier as previous averages are used later in calc_avg_solarwind_predstorm
   
   calls omni_loader() (which calls get_omni_data())
   '''
 
   o=omni_loader()     #load all omni data
        
   dt_mat=mdates.date2num(dt) #convert to matplotlib time    
   
          
   #starting index of dt start time in all omni data - 24 hours 
   #needed for averaging solar wind before time dt[0]
   inds=np.argmin(abs(o.time-dt_mat[0]))-24
   inde=np.argmin(abs(o.time-dt_mat[-1]))+24 #end index for dt in omni data, add +24 hours to end
    
   o_time=o.time[inds:inde]
   
   vartxtout=np.zeros([inde-inds,13])

   for i in np.arange(np.size(o_time)):  #get date in ascii   
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
   filename_save='data/predstorm/predstorm_omni.txt'
   np.savetxt(filename_save, vartxtout,  delimiter='',fmt='%4i %2i %2i %2i %2i %2i %10.6f %5.1f %5.1f %5.1f %5.1f   %7.0f %7.0f ',\
               header='        time      matplotlib_time B[nT] Bx   By     Bz   N[ccm-3] V[km/s] ')

   #with last 3 variables
   #np.savetxt(filename_save, vartxtout, delimiter='',fmt='%4i %2i %2i %2i %2i %2i %10.6f %5.1f %5.1f %5.1f %5.1f   %7.0i %7.0i   %5.0f %5.1f %5.1f', \
   #            header='        time      matplotlib_time B[nT] Bx   By     Bz   N[ccm-3] V[km/s] Dst[nT]   Kp   aurora [GW]')

    


def omni_loader():
   '''
   downloads all omni2 data into "data/omni2" folder
   converts to pickle file and returns object with data
   '''
   if os.path.isdir('data') == False: os.mkdir('data')
   if os.path.isdir('data/omni2') == False: os.mkdir('data/omni2')
  
   if not os.path.exists('data/omni2/omni2_all_years.dat'):
      #see http://omniweb.gsfc.nasa.gov/html/ow_data.html
      print('OMNI2 .dat file not in "data" directory, so download OMNI2 data from')
      omni2_url='ftp://nssdcftp.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_all_years.dat'
      print(omni2_url)
      try: urllib.request.urlretrieve(omni2_url, 'data/omni2/omni2_all_years.dat')
      except urllib.error.URLError as e:
          print(' ', omni2_url,' ',e.reason)

   #if omni2 hourly data is not yet converted and saved as pickle, do it:
   if not os.path.exists('data/omni2/omni2_all_years_pickle.p'):
       #load OMNI2 dataset from .dat file with a function from dst_module.py
       print('OMNI2 .p file not in "data" directory, so convert to pickle')
       o=get_omni_data()
       #contains: o. time,day,hour,btot,bx,by,bz,bygsm,bzgsm,speed,speedx,den,pdyn,dst,kp
       #save for faster loading later
       pickle.dump(o, open('data/omni2/omni2_all_years_pickle.p', 'wb') )

   else:  o=pickle.load(open('data/omni2/omni2_all_years_pickle.p', 'rb') )
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
    f=open('data/omni2/omni2_all_years.dat')
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
    with open('data/omni2/omni2_all_years.dat') as f:
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

 
 


















######################################## PLOTTING ########################################


def ovation_global_north(wic,dt,colormap_input,max_level, outputdir):
 '''
 plots the ovation aurora on the northern hemisphere from a polar view, makes movie frames
 wic is a world image cube with ovation results 512x1024
 dt are the datetimes for each frame
 colormap_input is the colormap
 max_level is the maximum level for the colorbar
 '''
 #plotting parameters
 #-100 for North America, +10 or 0 for Europe
 global_plot_longitude=-100
 #global_plot_longitude=0
 global_plot_latitude=90


 provinces_50m = carfeat.NaturalEarthFeature('cultural',
                                             'admin_1_states_provinces_lines',
                                             '50m',
                                             facecolor='none',edgecolor='black')



 #define extent of the produced world maps - defined as: west east south north
 mapextent=[-180,180,-90,90]   
 
 #need to set alpha linearly increasing in the colormap so that it blends into background for 
 #small values
 cmap = plt.get_cmap(colormap_input)  # Choose colormap
 my_cmap = cmap(np.arange(cmap.N))  # Get the colormap colors
 my_cmap[:,-1] = np.linspace(0, 1, cmap.N)  # Set alpha
 my_cmap = ListedColormap(my_cmap) # Create new colormap

 crs=ccrs.PlateCarree()

 fig = plt.figure(2,figsize=[12, 12],dpi=80) 
 fig.set_facecolor('black') 
 fig.text(0.99,0.01,'C. Möstl / IWF-helio, Austria', color='white',fontsize=10,ha='right',va='bottom')
 
 ax = plt.subplot(1, 1, 1, projection=ccrs.Orthographic(global_plot_longitude, global_plot_latitude))

 ax.background_patch.set_facecolor('k')    
 
 
 gl=ax.gridlines(linestyle='--',alpha=0.5,color='white') #make grid
 gl.n_steps=100   #make grid finer
 ax.stock_img()  
 
 #ax.coastlines('50m',color='black',alpha=0.5)
 ax.coastlines('50m',color='white',alpha=0.9)
 ax.add_feature(provinces_50m,alpha=0.2)    #,zorder=2,alpha=0.8)
 ax.add_feature(carfeat.BORDERS, alpha=0.2) #,zorder=2,alpha=0.5)
  

 
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
     img=ax.imshow(wic[:,:,i], vmin=0.01, vmax=max_level, transform=crs, extent=mapextent, origin='lower', zorder=3, alpha=0.8, cmap=my_cmap)
  
     #save as image with timestamp in filename
     #plot_Nhemi_filename='results/forecast_global/predstorm_aurora_real_Nhemi_'+dt.strftime("%Y_%m_%d_%H%M")  +'.jpg'
     #fig.savefig(plot_Nhemi_filename,dpi=150,facecolor=fig.get_facecolor())

     #save as movie frame
     framestr = '%05i' % (i)  
     fig.savefig('results/'+outputdir+'/frames_global/aurora_'+framestr+'.jpg',dpi=150,facecolor=fig.get_facecolor())
     
     #print('Saved image:  ',plot_Nhemi_filename)
     
     












      
    
    
    
    


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



 
 


def global_ovation_flux(magnetic_latitude,magnetic_local_time,flux,dt):

 
 #cmap_name='hot'
 cmap_name='jet'
 
 
 res_lat=40/80
 res_lon=360/96
 
 
 #read in IDL output for comparison
 idl_file_in='data/idl_output/ov_diff_Eflux_2017_1230_2330.txt'
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
 plt.savefig('results/flux_test.png')











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





   

#@njit    
#def calc_avg_solarwind_predstorm(dt_mat_hour,dt_mat, l1wind):
def calc_avg_solarwind_predstorm_1min(dt,l1wind):

    """
    Calculates a weighted average of speed and magnetic field
    ave_hours (4 by default) backward, for 1 minute data time resolution in input recarray l1wind
    
    input: 
    - datetime object dt, 1min resolution, for the times the auroramaps are produced
    - l1wind is the solar l1wind input is a recarray containing time (in matplotlib format), bx, by, bz, v, ec
    test run time with %timeit oup.calc_avg_solarwind_predstorm(ts[0],l1wind)    
    """
    ave_hours=4                #hours previous to integrate over, usually 4
    prev_hour_weight = 0.65    # reduce weighting by factor of wh each hour back
    
    dt_mat_hour=mdates.date2num(round_to_hour_start(dt))  #get with input dt to hour start and continute with matplotlib times
    closest_time_ind_hour=np.argmin(abs(l1wind.time-dt_mat_hour))     #find index of closest time to dt_mat_hour
   
    dt_mat=mdates.date2num(dt) #convert all dates to input matplotlib time
    closest_time_ind=np.argmin(abs(l1wind.time-dt_mat))  #find index of closest time to dt

    weights=np.ones(5)     #make array for weights, most recent hour is 1, then go down by 0.65, 0.65*0.65 ...
    #make array with weights according to Newell et al. 2010, par 25
    for k in np.arange(2,ave_hours+1,1):  weights[k] = weights[k-1]*prev_hour_weight
      
    #the weight for current time comes in at the beginning of the weights array      
    #difference between hour start and actual time, in hours between 0 and 1; this 
    #is the weight for the current hour (check IDL real time code for this)***
    weights[0]=(dt_mat-dt_mat_hour)*24
   
    #now define the times for the weights to sum bxyz, v and ec
    times_for_weight_ind=np.zeros(5,dtype=int)
    times_for_weight_ind[0]=closest_time_ind
    
    #time resolution 1 hour
    if l1wind.time[1]-l1wind.time[0] > 0.01:
        times_for_weight_ind[1:5] = np.arange(closest_time_ind_hour, closest_time_ind_hour-ave_hours,-1)
    else: #time resolution 1 minute	
        times_for_weight_ind[1:5] = np.arange(closest_time_ind_hour, closest_time_ind_hour-ave_hours*60,-60)
   

    #for debugging
    print('--')
    print('input time ',dt)
    print('hour start in predstorm',mdates.num2date(l1wind.time[closest_time_ind_hour])     )
    print('all weights are  ', weights)    
    print('time indices for weights are',times_for_weight_ind)
    weight_times=mdates.num2date(l1wind.time[times_for_weight_ind])
    print('times are',weight_times[0],weight_times[1], weight_times[2], weight_times[3], weight_times[4] )
    print('--')

   
    #******* HERE NEEDS TO BE AN AVERAGING FOR the FULL HOURS LIKE OMNI2, before! 1 time only?    
    #all l1wind.bx, ... need to be averaged 
    
    #make array of average solar l1wind variables
    avgsw=np.recarray(1,dtype=[('bx', float), ('by', float),('bz', float),('v', float),('ec', float)])
    
    #sum over last 4 hours with weighting
    avgsw.bx = np.round(np.nansum(l1wind.bx[times_for_weight_ind]*weights) / np.nansum(weights),2)
    avgsw.by = np.round(np.nansum(l1wind.by[times_for_weight_ind]*weights) / np.nansum(weights),2)
    avgsw.bz = np.round(np.nansum(l1wind.bz[times_for_weight_ind]*weights) / np.nansum(weights),2)
    avgsw.v  = np.round(np.nansum(l1wind.v [times_for_weight_ind]*weights) / np.nansum(weights),1)
    avgsw.ec = np.round(np.nansum(l1wind.ec[times_for_weight_ind]*weights) / np.nansum(weights),1)
            
    #return averaged solar l1wind for the input time dt
    return avgsw
 
'''



   