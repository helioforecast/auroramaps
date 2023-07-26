"""
ovation.py

part of the "auroramaps" package

The module for the ovation prime 2010 model in python

Modified source code taken originally as OvationPyme
from https://github.com/lkilcommons/OvationPyme

This module contains the main model routines for 
Ovation Prime (historically called season_epoch.pro in the IDL version)

by C. Moestl, Austrian Space Weather Office, GeoSphere Austria.
https://github.com/helioforecast/auroramaps
twitter @chrisoutofspace
https://helioforecast.space

using a rewritten version of the ovationpyme aurora model 
by Liam Kilcommons https://github.com/lkilcommons/OvationPyme

published under GNU Lesser General Public License v3.0

last update October 2019

 

"""
import os
import datetime
import numpy as np
from scipy import interpolate
from numba import njit, jit
import pdb
import sys
import pickle
import time	
import aacgmv2
import scipy





def make_aurora_cube(ts, ec,diff,mono):
  '''
  single processor version
  '''
  #this is the array with the final world maps for each timestep
  ovation_img=np.zeros([512,1024,np.size(ts)])
  #make a world map grid in latitude 512 pixels, longitude 1024 pixel like NOAA
  wx,wy=np.mgrid[-90:90:180/512,-180:180:360/1024]

  #if ts is a single not a list, change to list to make it subscriptable
  if type(ts)!=list: ts=[ts]

  for k in np.arange(0,np.size(ts)): #go through all times
 
    print('Frame number and time:', k, '  ',ts[k])
    

    #########################################  (2a) get solar wind 
    startflux=time.time()
   
    #sw=oup.calc_avg_solarwind_predstorm(ts[k],l1wind)   #make solarwind with averaging over last 4 hours, only for command line 
    #print('Byz =',sw.by[0],sw.bz[0],' nT   V =',int(sw.v[0]), 'km/s')
    #for the Newell coupling Ec, normalized to cycle average, smoothed for high time resolution
    #print('Ec to cycle average: ',np.round(swav.ec[k]/coup_cycle,1), ' <Ec> =',int(swav.ec[k]))  
    
    #get fluxes for northern hemisphere 
    mlatN, mltN, fluxNd=diff.get_flux_for_time(ts[k],ec[k])
    mlatN, mltN, fluxNm=mono.get_flux_for_time(ts[k],ec[k])
    #mlatN, mltN, fluxNw=we.get_flux_for_time(ts[k],inputfile, hemi='N')  #wave flux not correct yet
    fluxN=fluxNd+fluxNm #+fluxNw
    #endflux=time.time()
    #print('OVATION: ', np.round(endflux-startflux,2),' sec')
    
    #for debugging
    #making flux images for comparison to OVATION IDL output
    #change IDL file in this function for flux comparison
    #oup.global_ovation_flux(mlatN,mltN,fluxNw,ts[0])
    #oup.global_ovation_flux(mlatN,mltN,fluxNd,ts[k])
    #sys.exit()
  
    #if k==3: sys.exit()

 
    #####################################  (2b) coordinate conversion magnetic to geographic 
    #Coordinate conversion MLT to AACGM mlon/lat to geographic coordinates
    #startcoo=time.time()
    #magnetic coordinates are  mltN mlonN in 2D grids
    #so we need to convert magnetic local time (MLT) to longitude first
    #extract from mltN first 96 MLTs for 1 latitude bin convert to magnetic longitude 
    mlonN_1D_small=aacgmv2.convert_mlt(mltN[0],ts[k],m2a=True)
    #copy result  80 times because MLT is same for all latitudes
    mlonN_1D=np.tile(mlonN_1D_small,mlatN.shape[0])
    
    #this will make a 1D array for the latitudes compatible with the 1D array created above    
    mlatN_1D=np.squeeze(mlatN.reshape(np.size(mltN),1))
   
    #magnetic coordinates are now mlatN mlonN, convert to geographic
    (glatN_1D, glonN_1D, galtN) = aacgmv2.convert_latlon_arr(mlatN_1D,mlonN_1D, 100,ts[k], method_code="A2G")
    #endcoo = time.time()
    #print('Coordinates: ', np.round(endcoo-startcoo,2),' sec')
   

    #####################################  (2c) interpolate to world map 
    #geographic coordinates are glatN, glonN, electron flux values are fluxN
    #change shapes from glatN (80, 96) glonN  (80, 96) to a single array with 7680,2
    #startworld=time.time()
 
    #glatN_1D=glatN.reshape(np.size(fluxN),1)  
    #glonN_1D=glonN.reshape(np.size(fluxN),1)    
    geo_2D=np.vstack((glatN_1D,glonN_1D)).T      #stack 2 (7680,) arrays to a single 7680,2 arrays, .T is needed
    #smooth small array first - but strange results!
    #fluxN=scipy.ndimage.gaussian_filter(fluxN,sigma=(5,5))
    fluxN_1D=fluxN.reshape(7680,1)   #also change flux values to 1D array

    #bottleneck, maybe make own interpolation and accelerate with numba
    #https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.griddata.html

    #interpolate to world grid, and remove 1 dimension with squeeze
    aimg=np.squeeze(scipy.interpolate.griddata(geo_2D, fluxN_1D, (wx, wy), method='linear',fill_value=0))
    #filter large array
    aimg = scipy.ndimage.gaussian_filter(aimg,sigma=(5,7),mode='wrap') #wrap means wrapping at the 180 degree edge
    ovation_img[:,:,k]=aimg
    #endworld = time.time()
    #print('World map: ', np.round(endworld-startworld,2),' sec')
    
   

  return ovation_img






























class FluxEstimator(object):
    """
    A class which estimates auroral flux based on the Ovation Prime regressions,
    at arbitrary locations and times.

    Locations are in magnetic latitude and local time, and are interpolated using a B-spline
    representation.
    
    Ovation works like this: for each seasons there are different coefficients to make the flux bins
    The coefficients for all 4 seasons are loaded from txt files
    Later, each season is given a weight, those seasons with weight > 0 are evaluated
    and the fluxes are summed with the weight for each season
    
    time: 334 ms
        
    """
    def __init__(self, atype, jtype):
        """
        Initialize class with two variables:
        (1) Type of aurora for which to load regression coeffients
            atype - str, ['diff','mono','wave','ions']

        (2) Type of flux you want to estimate
            jtype - int or str
            1:"electron energy flux",
            2:"ion energy flux",
            3:"electron number flux",
            4:"ion number flux",
            5:"electron average energy",
            6:"ion average energy"
        """
        self.atype = atype #Type of aurora
        self.jtype = jtype #Type of flux 
        self.seasons = ['spring', 'summer', 'fall', 'winter']
       
        #create list of objects, one for each season
        self.seasonal_flux_estimators={}
        for i in self.seasons:
           self.seasonal_flux_estimators[i]=SeasonalFluxEstimator(i, self.atype, self.jtype)

  

    def get_flux_for_time(self, dt, ec_average):
        """
        returns grid_mlats, grid_mlts, gridflux for 
        given time, predstorm_input file, and hemisphere 
        only works for northern hemisphere at the moment
        
        dt: time of auroramap
        ec_average: Newell coupling for this timestep    
        
        calls: season_weights(), seasonal_flux_estimator.get_gridded_flux
  
        """
          
        doy = dt.timetuple().tm_yday        #get doy of current date dt which is a datetime object
        weightsN = self.season_weights(doy) #get weights for northern hemisphere
        weightsS = self.season_weights(365.-doy) #weights for southern hemisphere are simply phase-shifted by half a year
   
        #print(weightsN)
        #print(weightsS)
    
        
        #determines size of the grid - 40° latitude 50° to 90°, in 0.5° steps, and 24 hours time in 0.25h steps
        #taken from the SeasonalFluxEstimator class where this is defined
        n_mlat_bins = list(self.seasonal_flux_estimators.items())[0][1].n_mlat_bins//2 #div by 2 because combined N&S hemispheres 
        n_mlt_bins= list(self.seasonal_flux_estimators.items())[0][1].n_mlt_bins
        gridflux = np.zeros((n_mlat_bins, n_mlt_bins))         #makes the flux grid in magnetic coordinates with all values 0
        
        #this works only for N hemisphere at the moment
        #go through each season        
        for season in self.seasons:
            #if the weight for a season is greater 0, go on (otherwise the flux from this season stays 0)
            if np.logical_or(weightsN[season] > 0, weightsS[season] > 0):
                #print(season)
                flux_estimator = self.seasonal_flux_estimators[season]          #get one of the four SeasonalFluxEstimator objects, defined when calling the FluxEstimator class           
                #calculate fluxes for season and solar wind coupling ec average and for both hemispheres
                #grid_mlatsN, grid_mltsN, grid_fluxN, gridmlatsS, grid_mltsS, grid_fluxS = flux_estimator.get_gridded_flux(ec_average)
                grid_mlatsN, grid_mltsN, grid_fluxN = flux_estimator.get_gridded_flux(ec_average)                
                #add for this doy each season with weights
                gridflux = gridflux + grid_fluxN*weightsN[season][0]
                #gridflux = gridflux + grid_fluxN*float(weightsN[season])+grid_fluxS*float(weightsS[season])
        
        return grid_mlatsN, grid_mltsN, gridflux

        
     
    def season_weights(self,doy):
        """
        Determines the relative weighting of the model coefficients for the various 
        seasons for a particular day of year (doy). Nominally, weights the seasons
        based on the difference between the doy and the peak of the season (solstice/equinox)

        Returns:
        A numpy structured array with a key for each season, values between 0 and 1
        
        #time: 3.3 microsec
        """
        weight = np.zeros(1,dtype=[('winter', float), ('spring', float),('summer', float),('fall', float)])        
        
        if doy >= 79. and doy < 171:
            weight['summer'] = 1. - (171.-doy)/92.
            weight['spring'] = 1. - weight['summer']

        elif doy >= 171. and doy < 263.:
            weight['fall'] = 1. - (263.-doy)/92.
            weight['summer'] = 1. - weight['fall']

        elif doy >= 263. and doy < 354.:
            weight['winter'] = 1. - (354.-doy)/91.
            weight['fall'] = 1. - weight['winter']

        elif doy >= 354 or doy < 79:
            #For days of year > 354, subtract 365 to get negative
            #day of year values for computation
            doy0 = doy- 365. if doy >= 354 else doy
            weight['spring'] = 1. - (79.-doy0)/90.
            weight['winter'] = 1. - weight['spring']

        return weight
  
     
 
 
 
 
 
 
 
 
 
        


class SeasonalFluxEstimator(object):
    """
    A class to calculate predictions from the regression coeffecients
    which are tabulated in the auroramaps/data/premodel/{season}_{atype}_*.txt
    files.
    A function is included that loads these files into a single array that 
    is dumped in the python pickle file auroramaps/data/premodel/all_premodel_python.p

    Given a particular season, type of aurora ( one of ['diff','mono','wave','ions'])
    and type of flux, extracts all fit parameters upon initialization
    """
    
    def __init__(self, season, atype, jtype):
        """
        season - str,['winter','spring','summer','fall']
            season for which to load regression coeffients

        atype - str, ['diff','mono','wave','ions']
            type of aurora for which to load regression coeffients

        jtype - int or str
            1:"electron energy flux",
            2:"ion energy flux",
            3:"electron number flux",
            4:"ion number flux",
            5:"electron average energy",
            6:"ion average energy"
        """

        self.premodel_directory='data/premodel/'         #define premodel directory
        nmlt = 96              #number of mag local times in arrays (resolution of 15 minutes)
        nmlat = 160            #number of mag latitudes in arrays (resolution of 0.5 degree, for two hemispheres))
        nEc = 12               #number of coupling strength bins
        self.jtype, self.atype = jtype, atype
        self.n_mlt_bins, self.n_mlat_bins, self.n_Ec_bins = nmlt, nmlat, nEc


        #The mlat bins are organized like -50:-dlat:-90, 50:dlat:90
        self.mlats = np.concatenate([np.linspace(-90., -50., self.n_mlat_bins//2)[::-1],
                                     np.linspace(50., 90., self.n_mlat_bins//2)])

        self.mlts = np.linspace(0., 24., self.n_mlt_bins)

        self.fluxtypes = {1:"electron energy flux",
                          2:"ion energy flux",
                          3:"electron number flux",
                          4:"ion number flux",
                          5:"electron average energy",
                          6:"ion average energy"}
 
        #check whether premodel data has been loaded from the txts into one pickle file, do if not
        if os.path.isfile(self.premodel_directory+'/all_premodel_python.p') == False:
            self.make_premodel_pickle()          
      
        #load pickle file       
        ov=pickle.load(open(self.premodel_directory+'all_premodel_python.p', 'rb' ))
        #print('loaded pickle')
        #print(atype)
        #print(season)    

        #get current position in array for this type and season
        ov_ind=np.where(np.logical_and(ov['season']==season, ov['type']==atype))  
        
        #extract parameters at correct position; squeeze to remove first dimension
        self.b1a = np.squeeze(ov['b1a'][ov_ind])
        self.b2a = np.squeeze(ov['b2a'][ov_ind])
                
        #Determine if suffix is there
        file_suffix = '_n' if (jtype in [3, 4] or 'number flux' in jtype) else ''
        # if yes replace b1a b2a with b1an b2an
        if file_suffix == '_n':
            self.b1a = np.squeeze(ov['b1an'][ov_ind])
            self.b2a = np.squeeze(ov['b2an'][ov_ind])
        
        self.b1p = np.squeeze(ov['b1p'][ov_ind])
        self.b2p = np.squeeze(ov['b2p'][ov_ind])
        self.prob= np.squeeze(ov['prob'][ov_ind])



  
  
  
   
    def get_gridded_flux(self, Ec):
        """
        Return the flux interpolated onto arbitrary locations
        in mlats and mlts for a given Newell coupling (dF)

        Average the fluxes for northern and southern hemisphere
        and use them for both hemispheres (this is what standard
        ovation prime does always I think, so I've made it default)
        The original code says that this result is appropriate for 
        the northern hemisphere, and to use 365 - actual doy to
        get a combined result appropriate for the southern hemisphere
      
        - calls: make_flux_fast
         
         
        check time in main program 
        de=opp.FluxEstimator('diff', 'electron energy flux')   
        %timeit a,b,c=de.seasonal_flux_estimators['spring'].get_gridded_flux(swav.ec[0])  
        
        """

        #make grid and flux array for both hemispheres
        #north latitudes from +50 to +90 MLT from 0 to 24
        mlatgridN, mltgridN = np.meshgrid(self.mlats[self.n_mlat_bins//2:], self.mlts, indexing='ij')
        #south latitudes from -50 to -90, MLT from 0 to 24
        mlatgridS, mltgridS = np.meshgrid(self.mlats[:self.n_mlat_bins//2], self.mlts, indexing='ij')
         
   
        fluxgridN = np.zeros((self.n_mlat_bins//2, self.n_mlt_bins)) #// means results is integer
        fluxgridS = np.zeros((self.n_mlat_bins//2, self.n_mlt_bins)) #// means results is integer

        #numba optimized version for electron energy flux; key for using numba is to not define arrays but pass all of them to the function
        fluxgridN,fluxgridS = make_flux_fast(fluxgridN, fluxgridS, self.n_mlt_bins, self.n_mlat_bins, Ec,self.b1a,self.b2a,self.b1p,self.b2p,self.prob,self.n_Ec_bins)
 
        #wedge interpolation
        fluxgridN = self.interp_wedge(mlatgridN, mltgridN, fluxgridN)
  
        #return mlatgridN, mltgridN, fluxgridN
        return mlatgridN, mltgridN, (fluxgridN+fluxgridS)/2.
      
        
     
      
       
        
    def interp_wedge(self, mlatgridN, mltgridN, fluxgridN):
        """
        Interpolates across the wedge shaped data gap around 50 magnetic latitude and 23-4 MLT.
        Interpolation is performed individually across each magnetic latitude ring,
        only missing flux values are filled with the using the interpolant
        """
        #Constants copied from IDL code
        x_mlt_min=-1.0   #minimum MLT for interpolation [hours]
        x_mlt_max=4.0    #maximum MLT for interpolation [hours]
        x_mlat_min=50.0  #minimum MLAT for interpolation [degrees]
        x_mlat_max=75.0  #maximum MLAT for interpolation [degrees] --change if desired (LMK increased this from 67->75)

        valid_interp_mlat_bins = np.logical_and(mlatgridN[:, 0]>=x_mlat_min, mlatgridN[:, 0]<=x_mlat_max)
       
        #go through all rings of constant latitude
        for i_mlat_bin in np.flatnonzero(valid_interp_mlat_bins):
            
            #get array for this latitude and this flux
            this_mlat = mlatgridN[i_mlat_bin, 0]
            this_flux = fluxgridN[i_mlat_bin, :]
            this_mlt = mltgridN[i_mlat_bin, :]
 
            #Change from 0-24 MLT to -12 to 12 MLT, so that there is no discontinuity at midnight when we interpolate
            this_mlt[this_mlt>12.] = this_mlt[this_mlt>12.]-24.
            valid_interp_mlt_bins = np.logical_and(this_mlt >= x_mlt_min, this_mlt <= x_mlt_max)     #select grid points from -1 to 4 h MLT
           
            #bins where flux is missing
            mlt_bins_missing_flux = np.logical_not(this_flux > 0.)
            #bins where flux is missing and which are in the valid range of the wedge for interpolation
            interp_bins_missing_flux = np.logical_and(valid_interp_mlt_bins, mlt_bins_missing_flux)
            
 
            if np.count_nonzero(interp_bins_missing_flux) > 0: #go through all bins with missing data for flux
                #Bins right next to missing wedge probably have bad statistics, so don't include them
                interp_bins_missing_flux_inds = np.flatnonzero(interp_bins_missing_flux)
                
                #edge function see below in numba section
                interp_bins_missing_flux = edge(interp_bins_missing_flux_inds,interp_bins_missing_flux)
                
                interp_source_bins = np.flatnonzero(np.logical_not(interp_bins_missing_flux))
                #interpolate all missing data points - this is the computational bottleneck
                flux_interp = interpolate.interp1d(this_mlt[interp_source_bins], this_flux[interp_source_bins], kind='linear')
                fluxgridN[i_mlat_bin, interp_bins_missing_flux] = flux_interp(this_mlt[interp_bins_missing_flux])

        return fluxgridN    
        
        
   
      
      
      

   
    def make_premodel_pickle(self):        
            """
            takes the premodel txt files and makes a pickle file for much faster initialization
            """
            print('load all premodel files and put into pickle file')
            
            self.seasons = ['spring', 'summer', 'fall', 'winter']
            self.types   = ['diff','mono','wave','ions']

            nx=self.n_mlt_bins
            ny=self.n_mlat_bins
            
            #define structured array for ovation data from premodel files    season type b1a b2a b1an b2an b1p b2p prob          
            #4*4 =16 because of 4 types and 4 seasons in total
            self.ov = np.zeros(16,dtype=[('season', 'U6'),('type', 'U4'), ('b1a', 'f8',(nx,ny)),('b2a', 'f8',(nx,ny)),
                                    ('b1an', 'f8',(nx,ny)),('b2an', 'f8',(nx,ny)),('b1p', 'f8',(nx,ny)),
                                    ('b2p', 'f8',(nx,ny)),('prob', 'f8',(nx, ny, self.n_dF_bins)) ] )

            counter=0
            
            
            #go through all seasons and types and load all files
            for season_counter in self.seasons:
                for flux_counter in self.types:
                
                    self.ov['season'][counter]=season_counter        
                    self.ov['type'][counter]=flux_counter        
                    print(self.ov['season'][counter])
                    print(self.ov['type'][counter])
                              
                    afile = self.premodel_directory+'{0}_{1}.txt'.format(season_counter, flux_counter)
                    adata=np.loadtxt(afile,skiprows=1)
                    b1a, b2a = np.zeros((self.n_mlt_bins, self.n_mlat_bins)), np.zeros((self.n_mlt_bins, self.n_mlat_bins))

                    #column 0 refers to mlt bins, column 1 to mlat bins
                    mlt_bin_inds, mlat_bin_inds = adata[:, 0].astype(int), adata[:, 1].astype(int)
                    #data is in column 2 and 3
                    b1a[mlt_bin_inds, mlat_bin_inds] = adata[:, 2]
                    b2a[mlt_bin_inds, mlat_bin_inds] = adata[:, 3]
     
                    anfile = self.premodel_directory+'{0}_{1}_n.txt'.format(season_counter, flux_counter)
                    andata=np.loadtxt(anfile,skiprows=1)
                    b1an, b2an = np.zeros((self.n_mlt_bins, self.n_mlat_bins)), np.zeros((self.n_mlt_bins, self.n_mlat_bins))
                    b1an[mlt_bin_inds, mlat_bin_inds] = andata[:, 2]
                    b2an[mlt_bin_inds, mlat_bin_inds] = andata[:, 3]
                     
                    print(afile)
                    print(anfile)
                    self.ov['b1a'][counter]=b1a
                    self.ov['b2a'][counter]=b2a
                    self.ov['b1an'][counter]=b1an
                    self.ov['b2an'][counter]=b2an

                    #if one of three type where prob data files are available:
                    if flux_counter in ['diff', 'mono', 'wave']:
           
                        pfile = self.premodel_directory+'{0}_prob_b_{1}.txt'.format(season_counter, flux_counter)       
                        #load first part of this file into 2 arrays   
                        #pdata has 2 columns, b1, b2 for first 15360 rows
                        pdata_b = np.loadtxt(pfile, skiprows=1,max_rows=self.n_mlt_bins*self.n_mlat_bins) 
                        print(pfile)
                        #write into b1p array, shape is 96 * 160; format is similar to the afile
                        #load pfile
                        b1p, b2p = np.zeros((self.n_mlt_bins, self.n_mlat_bins)), np.zeros((self.n_mlt_bins, self.n_mlat_bins))
                        b1p[mlt_bin_inds, mlat_bin_inds]=pdata_b[:, 0]
                        b2p[mlt_bin_inds, mlat_bin_inds]=pdata_b[:, 1]

                        #load 2nd part of the file 
                        #184320 rows in 1 column = 160 latitude bins *96 localtime bins *12 ndF
                        pdata_p = np.loadtxt(pfile, skiprows=self.n_mlt_bins*self.n_mlat_bins+1) 
                        
                        #in the file the probability is stored with coupling strength bin
                        #varying fastest (this is Fortran indexing order)
                        #reshape from 184320, to 15360,12
                        pdata_p2 = pdata_p.reshape((-1, self.n_dF_bins), order='F')

                        #self.prob has shape  (96, 160, 12)
                        prob = np.zeros((self.n_mlt_bins,self.n_mlat_bins, self.n_dF_bins))
                        for idF in range(self.n_dF_bins):
                            prob[mlt_bin_inds, mlat_bin_inds, idF]=pdata_p2[:, idF]

                        #write to ov array
                        self.ov['b1p'][counter]=b1p
                        self.ov['b2p'][counter]=b2p
                        self.ov['prob'][counter]=prob

                    print(counter)
                    counter=counter+1
                    print()                    

            pickle.dump(self.ov,open(self.premodel_directory+'all_premodel_python.p', 'wb' ))
            print('done.')     
        
           
           
           
           
           
           
           
           
           
           
           
           
           
               
#################### NUMBA optimization in get_flux_for_time for electron energy flux

		

@njit
def make_flux_fast(fluxgridN,fluxgridS, n_mlt_bins, n_mlat_bins, Ec, b1a, b2a, b1p, b2p, prob, n_Ec_bins):  
    '''
    get electron energy flux for each grid point optimized with numba for southern hemisphere
    '''
    #go through each grid point
    for i_mlt in np.arange(n_mlt_bins):                #all magnetic local times 96
          for j_mlat in np.arange(n_mlat_bins//2):     #all magnetic latitudes 80

              #northern hemisphere 
              p = prob_estimate(Ec, b1p[i_mlt,j_mlat+n_mlat_bins//2],b2p[i_mlt,j_mlat+n_mlat_bins//2],prob[i_mlt,j_mlat+n_mlat_bins//2],n_Ec_bins) 
              flux=(b1a[i_mlt, j_mlat+n_mlat_bins//2]+b2a[i_mlt, j_mlat+n_mlat_bins//2]*Ec)*p
              fluxgridN[j_mlat, i_mlt] = correct_flux(flux)

              #southern hemisphere
              p = prob_estimate(Ec, b1p[i_mlt,j_mlat],b2p[i_mlt,j_mlat],prob[i_mlt,j_mlat],n_Ec_bins) 
              flux=(b1a[i_mlt, j_mlat]+b2a[i_mlt, j_mlat]*Ec)*p
              fluxgridS[j_mlat, i_mlt] = correct_flux(flux)
              
    return fluxgridN,fluxgridS


@njit  #ok mit ovationpyme
def correct_flux(flux):
    #for electron energy flux only - ***add others
    if flux < 0.: flux = 0.
    if flux > 10.:
        flux = 0.5
    elif flux > 5.:
        flux = 5.             
    return flux

     
@njit      
def prob_estimate(Ec, b1,b2,prob,n_Ec_bins): 

    p1 = b1 + b2*Ec 

    if p1 < 0.: p1 = 0.
    if p1 > 1.: p1 = 1.

    if b1 == 0. and b2 == 0.:
        i_Ecbin = which_Ec_bin(Ec,n_Ec_bins)
        p1 = prob[i_Ecbin]
        if p1 == 0.:
           i_Ecbin_1 = i_Ecbin - 1 if i_Ecbin > 0 else i_Ecbin+2 
           i_Ecbin_2 = i_Ecbin + 1 if i_Ecbin < n_Ec_bins-1 else i_Ecbin-2 
           p1 = (prob[i_Ecbin_1] + prob[i_Ecbin_2])/2.
    return p1
        
@njit  #***check in IDL version and ovationpyme
def which_Ec_bin(Ec,n_Ec_bins):
        i_Ecbin = np.round(Ec/(4421./8.))
        if i_Ecbin < 0: i_Ecbin = 0
        if i_Ecbin > 0: i_Ecbin = n_Ec_bins-1    
        return int(i_Ecbin)



@njit
def edge(interp_bins_missing_flux_inds,interp_bins_missing_flux):
     for edge_offset in np.arange(1, 7):
       lower_edge_ind = interp_bins_missing_flux_inds[0]-edge_offset
       upper_edge_ind = np.mod(interp_bins_missing_flux_inds[-1]+edge_offset, len(interp_bins_missing_flux))       
       interp_bins_missing_flux[lower_edge_ind] = interp_bins_missing_flux[interp_bins_missing_flux_inds[0]]
       interp_bins_missing_flux[upper_edge_ind] = interp_bins_missing_flux[interp_bins_missing_flux_inds[-1]]
     return interp_bins_missing_flux
		






