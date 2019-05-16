"""
ovation_prime_predstorm.py

The module for the ovation prime model in python

Modified source code taken originally as OvationPyme
from https://github.com/lkilcommons/OvationPyme

This module contains the main model routines for 
Ovation Prime (historically called season_epoch.pro in the IDL version)


"""
import os
import datetime
import numpy as np
from scipy import interpolate
from numba import njit, jit, jitclass
import pdb
import sys

import ovation_utilities_predstorm as oup


#reload for debugging
#import importlib
#importlib.reload(oup) 



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

  

    def get_flux_for_time(self, dt, fileinput, hemi):
        """
        returns grid_mlats, grid_mlts, gridflux for 
        given time, predstorm_input file, and hemisphere 
        only works for northern hemisphere at the moment
        """
    
        #calculates the average solar wind of the last 4 hours with weighting as
        #described in Newell et al. 2010 JGR Seasonal variations, paragraph 25     
        #time dt, file with predstorm output; only avgsw.ec is used later
        avgsw = oup.calc_avg_solarwind_predstorm(dt,fileinput)
      
        doy = dt.timetuple().tm_yday        #get doy of current date dt which is a datetime object
        weightsN = self.season_weights(doy) #get weights for northern hemisphere
        
        #determines size of the grid - 40° latitude 50° to 90°, in 0.5° steps, and 24 hours time in 0.25h steps
        #taken from the SeasonalFluxEstimator class where this is defined
        n_mlat_bins = list(self.seasonal_flux_estimators.items())[0][1].n_mlat_bins//2 #div by 2 because combined N&S hemispheres 
        n_mlt_bins= list(self.seasonal_flux_estimators.items())[0][1].n_mlt_bins
        gridflux = np.zeros((n_mlat_bins, n_mlt_bins))         #makes the flux grid in magnetic coordinates with all values 0

        #go through each season        
        for season in self.seasons:
            #if the weight for a season is greater 0, go on (otherwise the flux from this season stays 0)
            if weightsN[season] > 0:
                flux_estimator = self.seasonal_flux_estimators[season]          #get one of the four SeasonalFluxEstimator objects, defined when calling the FluxEstimator class           
                #calculate fluxes for season and solar wind coupling avgsw.ec and for both hemispheres
                #grid_mlatsN, grid_mltsN, grid_fluxN, gridmlatsS, grid_mltsS, grid_fluxS = flux_estimator.get_gridded_flux(avgsw.ec)
                grid_mlatsN, grid_mltsN, grid_fluxN = flux_estimator.get_gridded_flux(avgsw.ec[0])                
                #add for this doy each season with weights
                gridflux = gridflux+ grid_fluxN*float(weightsN[season])
        
        return grid_mlatsN, grid_mltsN, gridflux



        #older code
        #get season weights for current time dt
        #if hemi=='N':
        #    weightsN = self.season_weights(doy)
        #    weightsS = self.season_weights(365.-doy)
        #elif hemi=='S':        
            ######### ****************************** correct or reverse?
        #    weightsN = self.season_weights(doy)
        #    weightsS = self.season_weights(365.-doy)
            #This is already handled in the organization of the
            #ptype datafiles (re Betsey Mitchell)- LMK, 6-15-2017
            #weightsN = self.season_weights(365.-doy)
            #weightsS = self.season_weights(doy)
        #other older code:    
        #if hemi == 'S':
        # doy = 365.-doy #Use opposite season coefficients to get southern hemisphere results
        #print('North: season weights for doy: ',doy,'  ', weightsN) 
        #print('South: season weights for doy: ',doy,'  ', weightsS) 
        
        #******check summarizing both hemispheres and then dividing flux by 0.5
        #gridflux = gridflux+ grid_fluxN*float(weightsN[season])+grid_fluxS*float(weightsS[season])
        
        #old:
        #if hemi == 'S':
        #    grid_mlats = -1.*grid_mlats #by default returns positive latitudes


        
     
    def season_weights(self,doy):
        """
        Determines the relative weighting of the model coefficients for the various 
        seasons for a particular day of year (doy). Nominally, weights the seasons
        based on the difference between the doy and the peak of the season (solstice/equinox)

        Returns:
        A numpy structured array with a key for each season, values between 0 and 1
        """
        weight = np.zeros(1,dtype=[('winter', float), ('spring', float),('summer', float),('fall', float)])        
        
        #*******************check
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
    which are tabulated in the data/premodel/{season}_{atype}_*.txt
    files.

    Given a particular season, type of aurora ( one of ['diff','mono','wave','ions'])
    and type of flux, returns
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

        premodel_directory='data/premodel/'         #define premodel directory
        nmlt = 96              #number of mag local times in arrays (resolution of 15 minutes)
        nmlat = 160            #number of mag latitudes in arrays (resolution of 0.5 degree, for two hemispheres))
        ndF = 12               #number of coupling strength bins
        self.jtype, self.atype = jtype, atype

        self.n_mlt_bins, self.n_mlat_bins, self.n_dF_bins = nmlt, nmlat, ndF

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
   
        #Determine file names
        file_suffix = '_n' if (jtype in [3, 4] or 'number flux' in jtype) else ''
        self.afile = premodel_directory+'{0}_{1}{2}.txt'.format(season, atype, file_suffix)
        self.pfile = premodel_directory+'{0}_prob_b_{1}.txt'.format(season, atype)
        
        #load file for flux coefficients
        adata=np.loadtxt(self.afile,skiprows=1)

        self.b1a, self.b2a = np.zeros((nmlt, nmlat)), np.zeros((nmlt, nmlat))
        mlt_bin_inds, mlat_bin_inds = adata[:, 0].astype(int), adata[:, 1].astype(int)
        self.b1a[mlt_bin_inds, mlat_bin_inds] = adata[:, 2]
        self.b2a[mlt_bin_inds, mlat_bin_inds] = adata[:, 3]

        self.b1p, self.b2p = np.zeros((nmlt, nmlat)), np.zeros((nmlt, nmlat))
        self.prob = np.zeros((nmlt, nmlat, ndF))
        #pdata has 2 columns, b1, b2 for first 15361 rows
        #pdata has nmlat*nmlt rows (one for each positional bin)
        
        #if one of three type where prob data files are available:
        if atype in ['diff', 'mono', 'wave']:
                
            #load file into 2 arrays   
            pdata_b = np.loadtxt(self.pfile,skiprows=1,max_rows=nmlt*nmlat) 
            pdata_p = np.loadtxt(self.pfile, skiprows=nmlt*nmlat+1)#max_rows=nmlt*nmlat) 
            
            #in the file the probability is stored with coupling strength bin
            #varying fastest (this is Fortran indexing order)
            pdata_p_column_dFbin = pdata_p.reshape((-1, ndF), order='F')

            #mlt is first dimension
            self.b1p[mlt_bin_inds, mlat_bin_inds]=pdata_b[:, 0]
            self.b2p[mlt_bin_inds, mlat_bin_inds]=pdata_b[:, 1]
            for idF in range(ndF):
                self.prob[mlt_bin_inds, mlat_bin_inds, idF]=pdata_p_column_dFbin[:, idF]

        
        

    def which_dF_bin(self, dF):
        """
        Given a coupling strength value, finds the bin it falls into
        """
        #4421 is is a solar cycle average of the Ec coupling parameter, e.g. Newell et al. 2009
        i_dFbin = np.round(dF/(4421./8.))
        if i_dFbin < 0: i_dFbin = 0
        if i_dFbin > 0: i_dFbin = self.n_dF_bins-1
    
        return int(i_dFbin)
        
     
        

    def prob_estimate(self, dF, i_mlt_bin, i_mlat_bin): ########### ****** OPTIMIZE
        """
        Estimate probability of <something> by using tabulated
        linear regression coefficients ( from prob_b files )
        WRT coupling strength dF (which are different for each position bin)

        If p doesn't come out sensible by the initial regression,
        (i.e both regression coefficients are zero)
        then tries loading from the probability array. If the value
        in the array is zero, then estimates a value using adjacent
        coupling strength bins in the probability array
        """

        #Look up the regression coefficients
        b1, b2 = self.b1p[i_mlt_bin, i_mlat_bin], self.b2p[i_mlt_bin, i_mlat_bin]

        p1 = b1 + b2*dF #What is this the probability of?
        
        #pdb.set_trace()
        
        #check if p1 is in correct range 0 to 1, otherwise correct    
        if p1 < 0.: p1 = 0.
        if p1 > 1.: p1 = 1.

        if b1 == 0. and b2 == 0.:
            i_dFbin = self.which_dF_bin(dF)
            #Get the tabulated probability
            p1 = self.prob[i_mlt_bin, i_mlat_bin, i_dFbin]

            if p1 == 0.:
                #If no tabulated probability we must estimate by interpolating
                #between adjacent coupling strength bins
                i_dFbin_1 = i_dFbin - 1 if i_dFbin > 0 else i_dFbin+2 #one dF bin before by preference, two after in extremis
                i_dFbin_2 = i_dFbin + 1 if i_dFbin < self.n_dF_bins-1 else i_dFbin-2 #one dF bin after by preference, two before in extremis
                p1 = (self.prob[i_mlt_bin, i_mlat_bin, i_dFbin_1] + self.prob[i_mlt_bin, i_mlat_bin, i_dFbin_2])/2.

        return p1
        
   
        
            
    def correct_flux(self, flux):
        """
        A series of magical (unexplained, unknown) corrections to flux given 
        a particular type of flux
        """
        fluxtype = self.jtype

        if flux < 0.: flux = 0.

        if self.atype is not 'ions':
            #Electron Energy Flux
            if fluxtype in [1, self.fluxtypes[1]]:
                if flux > 10.:
                    flux = 0.5
                elif flux > 5.:
                    flux = 5.
            #Electron Number Flux
            if fluxtype in [3, self.fluxtypes[3]]:
                if flux > 2.0e9:
                    flux = 1.0e9
                elif flux > 2.0e10:
                    flux = 0.
        else:
            #Ion Energy Flux
            if fluxtype in [2, self.fluxtypes[2]]:
                if flux > 2.:
                    flux = 2.
                elif flux > 4.:
                    flux = 0.25
            #Ion Number Flux
            if fluxtype in [4, self.fluxtypes[4]]:
                if flux > 1.0e8:
                    flux = 1.0e8
                elif flux > 5.0e8:
                    flux = 0.
        return flux



        

    def estimate_auroral_flux(self, dF, i_mlt_bin, i_mlat_bin): ########### ****** OPTIMIZE
        """
        estimates the flux using the regression coeffecients in the 'a' files
        """
        p = self.prob_estimate(dF, i_mlt_bin, i_mlat_bin)
        #print(p, b1, b2, dF)
        flux = (self.b1a[i_mlt_bin, i_mlat_bin]+self.b2a[i_mlt_bin, i_mlat_bin]*dF)*p
        return self.correct_flux(flux)
        
        
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
        #x_mlat_max=67.0
        x_mlat_max=75.0  #maximum MLAT for interpolation [degrees] --change if desired (LMK increased this from 67->75)

        valid_interp_mlat_bins = np.logical_and(mlatgridN[:, 0]>=x_mlat_min, mlatgridN[:, 0]<=x_mlat_max)
        #inwedge = np.zeros(fluxgridN.shape, dtype=bool) #Store where we did interpolation
       
        #go through all rings of constant latitude
        for i_mlat_bin in np.flatnonzero(valid_interp_mlat_bins):
            
            #get array for this latitude and this flux
            this_mlat = mlatgridN[i_mlat_bin, 0]
            this_flux = fluxgridN[i_mlat_bin, :]
            this_mlt = mltgridN[i_mlat_bin, :]

            #Change from 0-24 MLT to -12 to 12 MLT, so that there is no discontinuity at midnight when we interpolate
            this_mlt[this_mlt>12.] = this_mlt[this_mlt>12.]-24.
        
            #select grid points from -1 to 4 h MLT
            valid_interp_mlt_bins = np.logical_and(this_mlt >= x_mlt_min, this_mlt <= x_mlt_max)
            
    
            #bins where flux is missing
            mlt_bins_missing_flux = np.logical_not(this_flux > 0.)
            #bins where flux is missing and which are in the valid range of the wedge for interpolation
            interp_bins_missing_flux = np.logical_and(valid_interp_mlt_bins, mlt_bins_missing_flux)
            
            #inwedge[i_mlat_bin, :] = interp_bins_missing_flux     
            
            if np.count_nonzero(interp_bins_missing_flux) > 0: #go through all bins with missing data for flux
                #Bins right next to missing wedge probably have bad statistics, so don't include them
                interp_bins_missing_flux_inds = np.flatnonzero(interp_bins_missing_flux)
                for edge_offset in np.arange(1, 7):
                        lower_edge_ind = interp_bins_missing_flux_inds[0]-edge_offset
                        upper_edge_ind = np.mod(interp_bins_missing_flux_inds[-1]+edge_offset, len(interp_bins_missing_flux))       
                        interp_bins_missing_flux[lower_edge_ind] = interp_bins_missing_flux[interp_bins_missing_flux_inds[0]]
                        interp_bins_missing_flux[upper_edge_ind] = interp_bins_missing_flux[interp_bins_missing_flux_inds[-1]]
                
                interp_source_bins = np.flatnonzero(np.logical_not(interp_bins_missing_flux))
                 
                #interpolate all missing data points 
                flux_interp = interpolate.interp1d(this_mlt[interp_source_bins], this_flux[interp_source_bins], kind='linear')
                fluxgridN[i_mlat_bin, interp_bins_missing_flux] = flux_interp(this_mlt[interp_bins_missing_flux])

        return fluxgridN      
      


    


   
   
    def get_gridded_flux(self, dF, combined_N_and_S=False, interp_N=True):
        """
        Return the flux interpolated onto arbitrary locations
        in mlats and mlts for a given Newell coupling (dF)

        combined_N_and_S, bool, optional
            Average the fluxes for northern and southern hemisphere
            and use them for both hemispheres (this is what standard
            ovation prime does always I think, so I've made it default)
            The original code says that this result is appropriate for 
            the northern hemisphere, and to use 365 - actual doy to
            get a combined result appropriate for the southern hemisphere

        interp_N, bool, optional
            Interpolate flux linearly for each latitude ring in the wedge
            of low coverage in northern hemisphere dawn/midnight region
        """

        #make arrays for northern hemisphere
        fluxgridN = np.zeros((self.n_mlat_bins//2, self.n_mlt_bins)) #// means results is integer
        #Make grid coordinates
        mlatgridN, mltgridN = np.meshgrid(self.mlats[self.n_mlat_bins//2:], self.mlts, indexing='ij')

        #make arrays for southern hemisphere
        #fluxgridS = np.zeros((self.n_mlat_bins//2, self.n_mlt_bins))
        #Make grid coordinates
        #mlatgridS, mltgridS = np.meshgrid(self.mlats[:self.n_mlat_bins//2], self.mlts, indexing='ij')
        
        #go through each grid point of both hemispheres to get the flux with estimate_auroral_flux
        for i_mlt in np.arange(self.n_mlt_bins):         #all magnetic local times
            for j_mlat in np.arange(self.n_mlat_bins//2):     #all magnetic latitudes
                #The mlat bins for the northern hemisphere start at 80, southern at 0
                fluxgridN[j_mlat, i_mlt] = self.estimate_auroral_flux(dF, i_mlt, self.n_mlat_bins//2+j_mlat)
                #for southern hemisphere do
                #fluxgridS[j_mlat, i_mlt] = self.estimate_auroral_flux(dF, i_mlt, j_mlat)

        
        #******check why mlts < 0 happen here
        #interpolate wedge
        if interp_N: fluxgridN = self.interp_wedge(mlatgridN, mltgridN, fluxgridN)
    
        return mlatgridN, mltgridN, fluxgridN

        #if not combined_N_and_S:
        #    return mlatgridN, mltgridN, fluxgridN, mlatgridS, mltgridS, fluxgridS
        #else:
        #    return mlatgridN, mltgridN, (fluxgridN+fluxgridS)/2.
      
      
      
      
# functions to optimize with numba      
      
      
      
      
      
      
      
      
      
      
def which_dF_bin2(self, dF):
        """
        Given a coupling strength value, finds the bin it falls into
        """
        #4421 is is a solar cycle average of the Ec coupling parameter, e.g. Newell et al. 2009
        i_dFbin = np.round(dF/(4421./8.))
        if i_dFbin < 0: i_dFbin = 0
        if i_dFbin > 0: i_dFbin = self.n_dF_bins-1
    
        return int(i_dFbin)  
      
      
      
@njit      
def interp_wedge2(mlatgridN, mltgridN, fluxgridN):
        """
        Interpolates across the wedge shaped data gap around 50 magnetic latitude and 23-4 MLT.
        Interpolation is performed individually across each magnetic latitude ring,
        only missing flux values are filled with the using the interpolant
        """
        #Constants copied from IDL code
        x_mlt_min=-1.0   #minimum MLT for interpolation [hours]
        x_mlt_max=4.0    #maximum MLT for interpolation [hours]
        x_mlat_min=50.0  #minimum MLAT for interpolation [degrees]
        #x_mlat_max=67.0
        x_mlat_max=75.0  #maximum MLAT for interpolation [degrees] --change if desired (LMK increased this from 67->75)

        valid_interp_mlat_bins = np.logical_and(mlatgridN[:, 0]>=x_mlat_min, mlatgridN[:, 0]<=x_mlat_max)
        #inwedge = np.zeros(fluxgridN.shape, dtype=bool) #Store where we did interpolation
        
        pdb.set_trace()
        #go through all rings of constant latitude
        for i_mlat_bin in np.flatnonzero(valid_interp_mlat_bins):
            
            #get array for this latitude and this flux
            this_mlat = mlatgridN[i_mlat_bin, 0]
            this_flux = fluxgridN[i_mlat_bin, :]
            this_mlt = mltgridN[i_mlat_bin, :]

            #Change from 0-24 MLT to -12 to 12 MLT, so that there is no discontinuity at midnight when we interpolate
            this_mlt[this_mlt>12.] = this_mlt[this_mlt>12.]-24.
        
            #select grid points from -1 to 4 h MLT
            valid_interp_mlt_bins = np.logical_and(this_mlt >= x_mlt_min, this_mlt <= x_mlt_max)
            
    
            #bins where flux is missing
            mlt_bins_missing_flux = np.logical_not(this_flux > 0.)
            #bins where flux is missing and which are in the valid range of the wedge for interpolation
            interp_bins_missing_flux = np.logical_and(valid_interp_mlt_bins, mlt_bins_missing_flux)
            
            #inwedge[i_mlat_bin, :] = interp_bins_missing_flux     
            
            if np.count_nonzero(interp_bins_missing_flux) > 0: #go through all bins with missing data for flux
                #Bins right next to missing wedge probably have bad statistics, so don't include them
                interp_bins_missing_flux_inds = np.flatnonzero(interp_bins_missing_flux)
                for edge_offset in np.arange(1, 7):
                        lower_edge_ind = interp_bins_missing_flux_inds[0]-edge_offset
                        upper_edge_ind = np.mod(interp_bins_missing_flux_inds[-1]+edge_offset, len(interp_bins_missing_flux))       
                        interp_bins_missing_flux[lower_edge_ind] = interp_bins_missing_flux[interp_bins_missing_flux_inds[0]]
                        interp_bins_missing_flux[upper_edge_ind] = interp_bins_missing_flux[interp_bins_missing_flux_inds[-1]]
                
                interp_source_bins = np.flatnonzero(np.logical_not(interp_bins_missing_flux))
                 
                #interpolate all missing data points 
                flux_interp = interpolate.interp1d(this_mlt[interp_source_bins], this_flux[interp_source_bins], kind='linear')
                fluxgridN[i_mlat_bin, interp_bins_missing_flux] = flux_interp(this_mlt[interp_bins_missing_flux])

        return fluxgridN      
      