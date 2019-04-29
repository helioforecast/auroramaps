"""
This module contains the main model routines for 
Ovation Prime (historically called season_epoch.pro)


modified source code taken originally from https://github.com/lkilcommons/OvationPyme
"""
import os
import datetime
import numpy as np
from scipy import interpolate

import ovation_utilities_predstorm

#previously needed in ovationpyme but not anymore?
#import geospacepy
#from geospacepy import special_datetime,astrodynamics2,satplottools
#import aacgmv2 #available on pip
#import apexpy

#Determine where this module's source file is located
#to determine where to look for the tables
src_file_dir = os.path.dirname(os.path.realpath(__file__))
ovation_datadir = os.path.join(src_file_dir,'data')







class FluxEstimator(object):
    """
    A class which estimates auroral flux
    based on the Ovation Prime regressions,
    at arbitrary locations and times.

    Locations are in magnetic latitude and local
    time, and are interpolated using a B-spline
    representation
    """
    def __init__(self, atype, jtype, seasonal_estimators=None, start_dt=None, end_dt=None):
        """
        doy - int
            day of year

        atype - str, ['diff','mono','wave','ions']
            type of aurora for which to load regression coeffients

        jtype - int or str
            1:"electron energy flux",
            2:"ion energy flux",
            3:"electron number flux",
            4:"ion number flux",
            5:"electron average energy",
            6:"ion average energy"

            Type of flux you want to estimate

        seasonal_estimators - dict, optional
            A dictionary of SeasonalFluxEstimators for seasons
            'spring','fall','summer','winter', if you
            don't want to create them
            (for efficiency across multi-day calls)
        """
        self.atype = atype #Type of aurora
        self.jtype = jtype #Type of flux
        if start_dt is not None and end_dt is not None:
            self.oi = geospacepy.omnireader.omni_interval(start_dt-datetime.timedelta(days=1),
                                                            end_dt+datetime.timedelta(days=1),
                                                            '5min', silent=True) #Give 1 day +- buffer because we need avg SW
            #Pre-create an omni interval (for speed if you are estimating auroral flux across many days)
        else:
            self.oi = None #omni_interval objects will be created on-the-fly (slow, but fine for single calls to get_flux)

        seasons = ['spring', 'summer', 'fall', 'winter']

        if seasonal_estimators is None:
            #Make a seasonal estimator for each season with nonzero weight
            self.seasonal_flux_estimators = {season:SeasonalFluxEstimator(season, atype, jtype) for season in seasons}
        else:
            #Ensure the passed seasonal estimators are approriate for this atype and jtype
            for season, estimator in seasonal_flux_estimators.items():
                jtype_atype_ok = jtype_atype_ok and (self.jtype == estimator.jtype and self.atype == estimator.atype)
            if not jtype_atype_ok:
                raise RuntimeError('Auroral and flux type of SeasonalFluxEstimators do not match {0} and {1}!'.format(self.atype, self.jtype))

    def get_flux_for_time(self, dt, fileinput, hemi,  return_dF=False):
        """
        doy must be single value
        mlats and mlts can be arbitary shape, but both must be same shape
        """
        doy = dt.timetuple().tm_yday
        #if hemi == 'S':
        #   doy = 365.-doy #Use opposite season coefficients to get southern hemisphere results

        if hemi=='N':
            weightsN = self.season_weights(doy)
            weightsS = self.season_weights(365.-doy)
        elif hemi=='S':
            weightsN = self.season_weights(doy)
            weightsS = self.season_weights(365.-doy)
            #This is already handled in the organization of the
            #ptype datafiles (re Betsey Mitchell)- LMK, 6-15-2017
            #weightsN = self.season_weights(365.-doy)
            #weightsS = self.season_weights(doy)

        avgsw = ovation_utilities_predstorm.calc_avg_solarwind_predstorm(dt,fileinput,4)
        #avgsw = ovation_predstorm.calc_avg_solarwind(dt, oi=self.oi)

        dF = avgsw['Ec']
        seasonal_flux = {}
        n_mlat_bins = list(self.seasonal_flux_estimators.items())[0][1].n_mlat_bins//2 #div by 2 because combined N&S hemispheres
        n_mlt_bins = list(self.seasonal_flux_estimators.items())[0][1].n_mlt_bins
        gridflux = np.zeros((n_mlat_bins, n_mlt_bins))
        for season in weightsN:
            flux_estimator = self.seasonal_flux_estimators[season]
            grid_mlats, grid_mlts, grid_fluxN, gridmlatsS, grid_mltsS, grid_fluxS = flux_estimator.get_gridded_flux(dF)
            gridflux += grid_fluxN*weightsN[season]+grid_fluxS*weightsS[season]

        #Because we added together both hemispheres we now must divide everything by two
        #to get the proper opposite season for other hemisphere weighting
        gridflux *= .5

        if hemi == 'S':
            grid_mlats = -1.*grid_mlats #by default returns positive latitudes

        if not return_dF:
            return grid_mlats, grid_mlts, gridflux
        else:
            return grid_mlats, grid_mlts, gridflux, dF

    def season_weights(self, doy):
        """
        Determines the relative weighting of the
        model coeffecients for the various seasons for a particular
        day of year (doy). Nominally, weights the seasons
        based on the difference between the doy and the peak
        of the season (solstice/equinox)

        Returns:
            a dictionary with a key for each season.
            Each value in the dicionary is a float between 0 and 1
        """
        weight = {'winter':0., 'spring':0., 'summer':0., 'fall':0.}
        winter_w, spring_w, summer_w, fall_w = 0., 0., 0., 0.

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






















class LatLocaltimeInterpolator(object):
    def __init__(self, mlat_grid, mlt_grid, var):
        self.mlat_orig = mlat_grid
        self.mlt_orig = mlt_grid
        self.zvar = var
        n_north, n_south = np.count_nonzero(self.mlat_orig>0.), np.count_nonzero(self.mlat_orig<0.)
        if n_south == 0.:
            self.hemisphere = 'N'
        elif n_north == 0.:
            self.hemisphere = 'S'
        else:
            raise ValueError('Latitude grid contains northern (N={0}) and southern (N={1}) values.'.format(n_north,n_south)+\
                                                ' Can only interpolate one hemisphere at a time.')

    def interpolate(self, new_mlat_grid, new_mlt_grid ,method='nearest'):
        """
        Rectangularize and Interpolate (using Linear 2D interpolation)
        """
        X0, Y0 = satplottools.latlt2cart(self.mlat_orig.flatten(), self.mlt_orig.flatten(),self.hemisphere)
        X, Y = satplottools.latlt2cart(new_mlat_grid.flatten(), new_mlt_grid.flatten(),self.hemisphere)
        interpd_zvar = interpolate.griddata((X0,Y0), self.zvar.flatten(), (X,Y), method=method, fill_value=0.)
        return interpd_zvar.reshape(new_mlat_grid.shape)


class BinCorrector(object):
    """
    We've found that often there are strange outlier bins that show up in
    OvationPyme results. This attempts to identify them by computing a numerical
    derivative around each ring of constant latitude.
    """
    def __init__(self, mlat_grid, mlt_grid):
        self.mlat_grid = mlat_grid
        self.mlats = self.mlat_grid[:, 0].flatten()
        self.mlt_grid = mlt_grid
        self.mlts = self.mlt_grid[0, :].flatten()
        self.dy_thresh = None

    def fix(self, y_grid, min_mlat=49, max_mlat=75, label=''):
        """
        Compute derivatives and attempt to identify bad bins
        Assumes mlat varies along the first dimension of the gridded location
        arrays
        """
        debug=False
        plot=False
        bad_bins = np.zeros_like(y_grid, dtype=bool)
        y_grid_corr = y_grid.copy()
        if self.dy_thresh is None:
                self.dy_thresh = 3.*np.nanstd(np.diff(y_grid.flatten()))
        wraparound = lambda x, nwrap: np.concatenate([x[-1*(nwrap+1):-1], x, x[:nwrap]])

        for i_mlat, mlat in enumerate(self.mlats):
            if not(np.abs(mlat)>=min_mlat and np.abs(mlat)<=max_mlat):
                if debug:
                    print('MLAT ring at {0} mlat is not between'.format(mlat)
                              +' {0} and {1}'.format(min_mlat, max_mlat)
                              +' skipping')
                continue
            mlts_nowrap = self.mlt_grid[i_mlat, :].copy()
            mlts_nowrap[mlts_nowrap<0] += 24
            mlts_nowrap[-1] = 23.9
            y = y_grid[i_mlat, :]
            #Wrap around first and last nwarp indicies in MLT
            #this prevents out of bounds errors in the spline/derviative
            nwrap = 4 # Pchip is cubic so order+1
            mlts = wraparound(mlts_nowrap, nwrap)
            mlts[:nwrap] -= 24. #to keep mlt in increasing order
            mlts[-1*nwrap:] += 24.
            y = wraparound(y, nwrap)
            #y_i = interpolate.PchipInterpolator(mlts, y)
            dy = np.diff(np.concatenate([y[:1], y])) # compute 1st derivative of spline
            i_dy = interpolate.interp1d(mlts,dy,kind='nearest')
            mlt_mask = np.ones_like(mlts,dtype=bool)
            for i_mlt,mlt in enumerate(mlts_nowrap.flatten()):
                if np.abs(i_dy(mlt)) > self.dy_thresh:
                    bad_bins[i_mlat, i_mlt] = True
                    mlt_mask[i_mlt+nwrap] = False

            y_corr_i = interpolate.PchipInterpolator(mlts[mlt_mask], y[mlt_mask])
            y_corr = y_corr_i(mlts)
            y_grid_corr[i_mlat, :] = y_corr_i(mlts_nowrap)
            if plot:
                self.plot_single_spline(mlat, mlts, y, dy, mlt_mask, y_corr, label=label)

        return y_grid_corr

    def plot_single_spline(self, mlat, mlts, y, dy, mlt_mask, y_corr, label=''):
        import matplotlib.pyplot as plt
        f = plt.figure(figsize=(8, 6))
        ax = f.add_subplot(111)
        ax.plot(mlts, y, 'bo', label='data')
        ax.plot(mlts, dy, 'r-', label='Deriv')
        bad_bins = np.logical_not(mlt_mask)
        ax.plot(mlts, y_corr, 'g.', label='After Correction')
        ax.plot(mlts[bad_bins], y[bad_bins], 'rx',
                label='Bad@dy>{0:.1f}'.format(self.dy_thresh))
        ax.set_title('Spline fit (mlat={0:.1f})'.format(mlat))
        ax.set_xlabel('MLT')
        ax.legend()
        if not os.path.exists('/tmp/ovationpyme'):
                os.makedirs('/tmp/ovationpyme')
        f.savefig('/tmp/ovationpyme/ovationpyme_spline_{0}_{1}.png'.format(label, np.floor(mlat*10)))
        plt.close(f)

    def __call__(self, y):
        """

        """
        return self.fix(y)


class ConductanceEstimator(object):
    """
    Implements the 'Robinson Formula'
    for estimating Pedersen and Hall height integrated
    conductivity (conducance) from
    average electron energy and
    total electron energy flux
    (assumes a Maxwellian electron energy distribution)
    """
    def __init__(self, start_dt, end_dt, fluxtypes=['diff', 'mono', 'wave']):
        #Use diffuse aurora only
        self.numflux_estimator = {}
        self.eavg_estimator = {}
        for fluxtype in fluxtypes:
            self.numflux_estimator[fluxtype] = FluxEstimator(fluxtype, 'electron number flux', start_dt=start_dt, end_dt=end_dt)
            #self.energyflux_estimator = FluxEstimator('diff', 'electron energy flux', start_dt=start_dt, end_dt=end_dt)
            self.eavg_estimator[fluxtype] = FluxEstimator(fluxtype, 'electron average energy', start_dt=start_dt, end_dt=end_dt)

        #Need hourly omni data for F10.7
        self.oi = geospacepy.omnireader.omni_interval(start_dt, end_dt, 'hourly', silent=True) 
        self.omjd = special_datetime.datetimearr2jd(self.oi['Epoch'])
        self.omf107 = self.oi['F10_INDEX']

    #def remove_bad_bins(self, numflux, eavg, max_replace=10):
    #    flat_nflux = numflux.flatten()
    #    flat_eavg = eavg.flatten()
    #    for i in range(max_replace):
    #        i_max = np.nanargmax(flat_nflux)
    #        flat_nflux[i_max]=np.nan
    #        flat_eavg[i_max]=np.nan

    #    #Calculate replacement values
    #    rep_nflux = np.nanmax(flat_nflux)
    #    i_replace = np.nanargmin(flat_nflux-rep_nflux)
    #    rep_eavg = flat_eavg[i_replace]
    #    #Find and replace any bind simultaneously about the rep values
    #    bad_bins = np.logical_or(flat_nflux>rep_nflux, flat_eavg>rep_eavg)
    #    flat_nflux[bad_bins] = rep_nflux
    #    flat_eavg[bad_bins] = rep_eavg
    #    fixed_numflux = flat_nflux.reshape(numflux.shape)
    #    fixed_eavg = flat_eavg.reshape(eavg.shape)
    #    print("Set %d bad bins to nflux %.2e eavg %.2f" % (np.count_nonzero(bad_bins),
    #                                                                                                            rep_nflux,rep_eavg))
    #    return fixed_numflux,fixed_eavg

    def get_closest_f107(self, dt):
        """
        Finds closest F10.7 value from hourly omni data to match with datetime
        Brekke and Moen describe using the daily F10.7 in the parameterizaton, 
        so I just do the mean for all of the 1 hour values for the day.
        """
        jd = special_datetime.datetime2jd(dt)
        imatch = np.floor(self.omjd.flatten())==np.floor(jd)
        return np.nanmean(self.omf107[imatch])

    def robinson_formula(self, numflux_grid, eavg_grid):
        #From E. Cousins IDL code
        #Implement the Robinson formula
        #Assume all of the particles come in at the average energy??

        energyflux_grid = numflux_grid*1.6022e-9*eavg_grid #keV to ergs, * #/(cm^2 s)
        #energyflux_grid *= 1.6022e-9
        sigp_auroral = 40.*eavg_grid/(16+eavg_grid**2) * np.sqrt(energyflux_grid)
        sigh_auroral = 0.45*eavg_grid**0.85*sigp_auroral
        return sigp_auroral, sigh_auroral

    def get_conductance(self, dt, hemi='N', solar=True, auroral=True,  background_p=None, background_h=None,
                        conductance_fluxtypes=['diff'], interp_bad_bins=True, return_dF=False):
        """
        Compute total conductance using Robinson formula and emperical solar conductance model
        """
        print("Getting conductance with solar {0}, aurora {1}, fluxtypes {2}, background_ped: {3}, background_hall {4}".format(solar,
                auroral, conductance_fluxtypes, background_p, background_h))

        all_sigp_auroral, all_sigh_auroral = [], []
        #Create a bin interpolation corrector
        for fluxtype in conductance_fluxtypes:
            mlat_grid, mlt_grid, numflux_grid, dF = self.numflux_estimator[fluxtype].get_flux_for_time(dt, hemi=hemi, return_dF=True)
            #mlat_grid, mlt_grid, energyflux_grid = self.energyflux_estimator.get_flux_for_time(dt, hemi=hemi)
            mlat_grid, mlt_grid, eavg_grid = self.eavg_estimator[fluxtype].get_flux_for_time(dt, hemi=hemi)

            if interp_bad_bins:
                #Clean up any extremely large bins
                fixer = BinCorrector(mlat_grid, mlt_grid)

                #Fix numflux
                fixer.dy_thresh = 1.0e8
                numflux_grid = fixer.fix(numflux_grid, label='nflux_{0}'.format(fluxtype))

                #Fix avg energy
                fixer.dy_thresh = .3
                eavg_grid = fixer.fix(eavg_grid, label='eavg_{0}'.format(fluxtype))

                #zero out lowest latitude numflux row because it makes no sense
                #has some kind of artefact at post midnight
                bad = np.abs(mlat_grid) < 52.0 
                numflux_grid[bad] = 0.

            #raise RuntimeError('Debug stop!')

            this_sigp_auroral, this_sigh_auroral = self.robinson_formula(numflux_grid, eavg_grid)
            all_sigp_auroral.append(this_sigp_auroral)
            all_sigh_auroral.append(this_sigh_auroral)

        sigp_solar, sigh_solar =  self.solar_conductance(dt, mlat_grid, mlt_grid)
        total_sigp_sqrd = np.zeros_like(sigp_solar)
        total_sigh_sqrd = np.zeros_like(sigh_solar)

        if solar:
            total_sigp_sqrd += sigp_solar**2
            total_sigh_sqrd += sigh_solar**2

        if auroral:
            #Sum up all contributions (sqrt of summed squares)
            for sigp_auroral, sigh_auroral in zip(all_sigp_auroral, all_sigh_auroral):
                    total_sigp_sqrd += sigp_auroral**2
                    total_sigh_sqrd += sigh_auroral**2

            #sigp_auroral *= 1.5
            #sigh_auroral *= 1.5

        #Now take square root to get hall and pedersen conductance
        if solar or auroral:
            sigp = np.sqrt(total_sigp_sqrd)
            sigh = np.sqrt(total_sigh_sqrd)
        else:
            #No conductance except flat background
            sigp = total_sigp_sqrd
            sigh = total_sigh_sqrd

        if background_h is not None and background_p is not None:
            #Cousins et. al. 2015, nightside artificial background of 4 Siemens
            #Ellen found this to be the background nightside conductance level which
            #best optimizes the SuperDARN ElePot AMIE ability to predict AMPERE deltaB data, and
            #the AMPERE MagPot AMIE ability to predict SuperDARN LOS V
            sigp[sigp<background_p]=background_p
            sigh[sigh<background_h]=background_h

        if return_dF:
            return mlat_grid, mlt_grid, sigp, sigh, dF
        else:
            return mlat_grid, mlt_grid, sigp, sigh

    def solar_conductance(self, dt, mlats, mlts):
        """
        Estimate the solar conductance using methods from:
            Cousins, E. D. P., T. Matsuo, and A. D. Richmond (2015), Mapping
            high-latitude ionospheric electrodynamics with SuperDARN and AMPERE

            --which cites--

            Asgeir Brekke, Joran Moen, Observations of high latitude ionospheric conductances

            Maybe is not good for SZA for southern hemisphere? Don't know
            Going to use absolute value of latitude because that's what's done
            in Cousins IDL code.
        """
        #Find the closest hourly f107 value
        #to the current time to specifiy the conductance
        f107 = self.get_closest_f107(dt)
        #print "F10.7 = {0}".format(f107)

        #Convert from magnetic to geocentric using the AACGMv2 python library
        flatmlats, flatmlts = mlats.flatten(), mlts.flatten()
        #flatmlons = (flatmlts-zero_lon_mlt)/12*180.
        flatmlons = aacgmv2.convert_mlt(flatmlts, dt, m2a=True)
        glats, glons = aacgmv2.convert(flatmlats, flatmlons, 110.*np.ones_like(flatmlats),
                                                                        date=dt, a2g=True, geocentric=False)
        szas = astrodynamics2.solar_zenith_angle(dt, glats, glons)
        szas_rad = szas/180.*np.pi

        sigp, sigh = np.zeros_like(glats), np.zeros_like(glats)

        cos65 = np.cos(65/180.*np.pi)
        sigp65  = .5*(f107*cos65)**(2./3)
        sigh65  = 1.8*np.sqrt(f107)*cos65
        sigp100 = sigp65-0.22*(100.-65.)

        in_band = szas <= 65.
        #print "%d/%d Zenith Angles < 65" % (np.count_nonzero(in_band), len(in_band))
        sigp[in_band] = .5*(f107*np.cos(szas_rad[in_band]))**(2./3)
        sigh[in_band] = 1.8*np.sqrt(f107)*np.cos(szas_rad[in_band])

        in_band = np.logical_and(szas >= 65., szas < 100.)
        #print "%d/%d Zenith Angles > 65 and < 100" % (np.count_nonzero(in_band), len(in_band))
        sigp[in_band] = sigp65-.22*(szas[in_band]-65.)
        sigh[in_band] = sigh65-.27*(szas[in_band]-65.)

        in_band = szas > 100.
        #print "%d/%d Zenith Angles > 100" % (np.count_nonzero(in_band), len(in_band))
        sigp[in_band] = sigp100-.13*(szas[in_band]-100.)
        sigh[in_band] = sigh65-.27*(szas[in_band]-65.)

        sigp[sigp<.4] = .4
        sigh[sigh<.8] = .8

        #correct for inverse relationship with magnetic field from AMIE code
        #(conductance_models.f90)
        theta = np.radians(90.-glats)
        bbp = np.sqrt(1. - 0.99524*np.sin(theta)**2)*(1. + 0.3*np.cos(theta)**2)
        bbh = np.sqrt(1. - 0.01504*(1.-np.cos(theta)) - 0.97986*np.sin(theta)**2)*(1.0+0.5*np.cos(theta)**2)
        sigp = sigp*1.134/bbp
        sigh = sigh*1.285/bbh

        sigp_unflat = sigp.reshape(mlats.shape)
        sigh_unflat = sigh.reshape(mlats.shape)

        return sigp_unflat, sigh_unflat



class SeasonalFluxEstimator(object):
    """
    A class to hold and caculate predictions from the regression coeffecients
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

        nmlt = 96                           #number of mag local times in arrays (resolution of 15 minutes)
        nmlat = 160                         #number of mag latitudes in arrays (resolution of 1/4 of a degree (.25))
        ndF = 12                                                        #number of coupling strength bins
        self.jtype, self.atype = jtype, atype

        self.n_mlt_bins, self.n_mlat_bins, self.n_dF_bins = nmlt, nmlat, ndF

        #The mlat bins are orgainized like -50:-dlat:-90, 50:dlat:90
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
        self.afile = os.path.join(ovation_datadir, 'premodel/{0}_{1}{2}.txt'.format(season, atype, file_suffix))
        self.pfile = os.path.join(ovation_datadir, 'premodel/{0}_prob_b_{1}.txt'.format(season, atype))
        #Defualt values of header (don't know why need yet)
        # b1 = 0.
        # b2 = 0.
        # yend = 1900
        # dend = 1
        # y0 = 1900
        # d0 = 1
        # files_done = 0
        # sf0 = 0
        self.valid_atypes = ['diff', 'mono', 'wave', 'ions']

        with open(self.afile, 'r') as f:
            aheader = f.readline() # y0,d0,yend,dend,files_done,sf0
            #print "Read Auroral Flux Coefficient File %s,\n Header: %s" % (self.afile,aheader)
            # Don't know if it will read from where f pointer is after reading header line
            adata = np.genfromtxt(f, max_rows=nmlat*nmlt)
            #print "First line was %s" % (str(adata[0,:]))

        self.b1a, self.b2a = np.zeros((nmlt, nmlat)), np.zeros((nmlt, nmlat))
        self.b1a.fill(np.nan)
        self.b2a.fill(np.nan)
        mlt_bin_inds, mlat_bin_inds = adata[:, 0].astype(int), adata[:, 1].astype(int)
        self.b1a[mlt_bin_inds, mlat_bin_inds] = adata[:, 2]
        self.b2a[mlt_bin_inds, mlat_bin_inds] = adata[:, 3]

        self.b1p, self.b2p = np.zeros((nmlt, nmlat)), np.zeros((nmlt, nmlat))
        self.prob = np.zeros((nmlt, nmlat, ndF))
        self.b1p.fill(np.nan)
        self.b2p.fill(np.nan)
        self.prob.fill(np.nan)
        #pdata has 2 columns, b1, b2 for first 15361 rows
        #pdata has nmlat*nmlt rows (one for each positional bin)

        if atype in ['diff', 'mono', 'wave']:
            with open(self.pfile, 'r') as f:
                pheader = f.readline() #y0,d0,yend,dend,files_done,sf0
                # Don't know if it will read from where f pointer is after reading header line
                pdata_b = np.genfromtxt(f, max_rows=nmlt*nmlat) # 2 columns, b1 and b2
                #print "Shape of b1p,b2p should be nmlt*nmlat=%d, is %s" % (nmlt*nmlat,len(pdata_b[:,0]))
                pdata_p = np.genfromtxt(f, max_rows=nmlt*nmlat*ndF) # 1 column, pval

            #in the file the probability is stored with coupling strength bin
            #varying fastest (this is Fortran indexing order)
            pdata_p_column_dFbin = pdata_p.reshape((-1, ndF), order='F')

            #I don't know why this is not used for atype 'ions'
            #mlt is first dimension
            self.b1p[mlt_bin_inds, mlat_bin_inds]=pdata_b[:, 0]
            self.b2p[mlt_bin_inds, mlat_bin_inds]=pdata_b[:, 1]
            for idF in range(ndF):
                self.prob[mlt_bin_inds, mlat_bin_inds, idF]=pdata_p_column_dFbin[:, idF]

            #if season=='spring' and atype=='diff' and jtype=='electron energy flux':
                #print(self.b1p[22:26, 138:142])
                #print(self.prob[22:26, 138:142, 5])

        #IDL original read
        #readf,20,i,j,b1,b2,rF
        #;;   b1a_all(atype, iseason,i,j) = b1
        #;;   b2a_all(atype, iseason,i,j) = b2
        #adata has 5 columns, mlt bin number, mlat bin number, b1, b2, rF
        #adata has nmlat*nmlt rows (one for each positional bin)

    def which_dF_bin(self, dF):
        """
        Given a coupling strength value, finds the bin it falls into
        """
        dFave = 4421. #Magic numbers!
        dFstep = dFave/8.
        i_dFbin = np.round(dF/dFstep)
        #Range check 0 <= i_dFbin <= n_dF_bins-1
        if i_dFbin < 0 or i_dFbin > self.n_dF_bins-1: 
            i_dFbin = 0 if i_dFbin < 0 else self.n_dF_bins-1
        return int(i_dFbin)

    def prob_estimate(self, dF, i_mlt_bin, i_mlat_bin):
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

        p = b1 + b2*dF #What is this the probability of?

        #range check 0<=p<=1
        if p < 0. or p > 1.:
            p = 1. if p > 1. else 0.

        if b1 == 0. and b2 == 0.:
            i_dFbin = self.which_dF_bin(dF)
            #Get the tabulated probability
            p = self.prob[i_mlt_bin, i_mlat_bin, i_dFbin]

            if p == 0.:
                #If no tabulated probability we must estimate by interpolating
                #between adjacent coupling strength bins
                i_dFbin_1 = i_dFbin - 1 if i_dFbin > 0 else i_dFbin+2 #one dF bin before by preference, two after in extremis
                i_dFbin_2 = i_dFbin + 1 if i_dFbin < self.n_dF_bins-1 else i_dFbin-2 #one dF bin after by preference, two before in extremis
                p = (self.prob[i_mlt_bin, i_mlat_bin, i_dFbin_1] + self.prob[i_mlt_bin, i_mlat_bin, i_dFbin_2])/2.

        return p

    def estimate_auroral_flux(self, dF, i_mlt_bin, i_mlat_bin):
        """
        Does what it says on the tin,
        estimates the flux using the regression coeffecients in the 'a' files
        """
        b1, b2 = self.b1a[i_mlt_bin, i_mlat_bin], self.b2a[i_mlt_bin, i_mlat_bin]
        p = self.prob_estimate(dF, i_mlt_bin, i_mlat_bin)
        #print(p, b1, b2, dF)
        flux = (b1+b2*dF)*p
        return self.correct_flux(flux)
            
    def correct_flux(self, flux):
        """
        A series of magical (unexplained, unknown) corrections to flux given a particular
        type of flux
        """
        fluxtype = self.jtype

        if flux < 0.:
            flux = 0.

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

    def get_gridded_flux(self, dF, combined_N_and_S=False, interp_N=True):
        """
        Return the flux interpolated onto arbitary locations
        in mlats and mlts

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

        fluxgridN = np.zeros((self.n_mlat_bins//2, self.n_mlt_bins))
        fluxgridN.fill(np.nan)
        #Make grid coordinates
        mlatgridN, mltgridN = np.meshgrid(self.mlats[self.n_mlat_bins//2:], self.mlts, indexing='ij')

        fluxgridS = np.zeros((self.n_mlat_bins//2, self.n_mlt_bins))
        fluxgridS.fill(np.nan)
        #Make grid coordinates
        mlatgridS, mltgridS = np.meshgrid(self.mlats[:self.n_mlat_bins//2], self.mlts, indexing='ij')
        #print(self.mlats[:self.n_mlat_bins//2])

        for i_mlt in range(self.n_mlt_bins):
            for j_mlat in range(self.n_mlat_bins//2):
                #The mlat bins are orgainized like -50:-dlat:-90,50:dlat:90
                fluxgridN[j_mlat, i_mlt] = self.estimate_auroral_flux(dF, i_mlt, self.n_mlat_bins//2+j_mlat)
                fluxgridS[j_mlat, i_mlt] = self.estimate_auroral_flux(dF, i_mlt, j_mlat)

        if interp_N:
            fluxgridN, inwedge = self.interp_wedge(mlatgridN, mltgridN, fluxgridN)
            self.inwedge = inwedge

        if not combined_N_and_S:
            return mlatgridN, mltgridN, fluxgridN, mlatgridS, mltgridS, fluxgridS
        else:
            return mlatgridN, mltgridN, (fluxgridN+fluxgridS)/2.

    def interp_wedge(self, mlatgridN, mltgridN, fluxgridN):
        """
        Interpolates across the wedge shaped data gap
        around 50 magnetic latitude and 23-4 MLT.
        Interpolation is performed individually
        across each magnetic latitude ring,
        only missing flux values are filled with the
        using the interpolant
        """
        #Constants copied verbatim from IDL code
        x_mlt_min=-1.0   #minimum MLT for interpolation [hours] --change if desired
        x_mlt_max=4.0    #maximum MLT for interpolation [hours] --change if desired
        x_mlat_min=49.0  #minimum MLAT for interpolation [degrees]
        #x_mlat_max=67.0
        x_mlat_max=75.0  #maximum MLAT for interpolation [degrees] --change if desired (LMK increased this from 67->75)

        valid_interp_mlat_bins = np.logical_and(mlatgridN[:, 0]>=x_mlat_min, mlatgridN[:, 0]<=x_mlat_max).flatten()
        inwedge = np.zeros(fluxgridN.shape, dtype=bool) #Store where we did interpolation

        for i_mlat_bin in np.flatnonzero(valid_interp_mlat_bins).tolist():
            #Technically any row in the MLT grid would do, but for consistancy use the i_mlat_bin-th one
            this_mlat = mlatgridN[i_mlat_bin, 0]
            this_mlt = mltgridN[i_mlat_bin, :]
            this_flux = fluxgridN[i_mlat_bin, :]

            #Change from 0-24 MLT to -12 to 12 MLT, so that there is no discontiunity at midnight
            #when we interpolate
            this_mlt[this_mlt>12.] = this_mlt[this_mlt>12.]-24.

            valid_interp_mlt_bins = np.logical_and(this_mlt>=x_mlt_min, this_mlt<=x_mlt_max).flatten()
            mlt_bins_missing_flux = np.logical_not(this_flux>0.).flatten()

            interp_bins_missing_flux = np.logical_and(valid_interp_mlt_bins, mlt_bins_missing_flux)

            inwedge[i_mlat_bin, :] = interp_bins_missing_flux

            if np.count_nonzero(interp_bins_missing_flux) > 0:
                #Bins right next to missing wedge probably have bad statistics, so
                #don't include them
                interp_bins_missing_flux_inds = np.flatnonzero(interp_bins_missing_flux)
                nedge=6
                for edge_offset in range(1, nedge+1):
                        lower_edge_ind = interp_bins_missing_flux_inds[0]-edge_offset
                        upper_edge_ind = np.mod(interp_bins_missing_flux_inds[-1]+edge_offset, len(interp_bins_missing_flux))
                        interp_bins_missing_flux[lower_edge_ind] = interp_bins_missing_flux[interp_bins_missing_flux_inds[0]]
                        interp_bins_missing_flux[upper_edge_ind] = interp_bins_missing_flux[interp_bins_missing_flux_inds[-1]]

                interp_source_bins = np.flatnonzero(np.logical_not(interp_bins_missing_flux))

                #flux_interp = interpolate.PchipInterpolator(this_mlt[interp_source_bins], this_flux[interp_source_bins])
                flux_interp = interpolate.interp1d(this_mlt[interp_source_bins], this_flux[interp_source_bins], kind='linear')
                fluxgridN[i_mlat_bin, interp_bins_missing_flux] = flux_interp(this_mlt[interp_bins_missing_flux])

                #print fluxgridN[i_mlat_bin,interp_bins_missing_flux]
                #print "For latitude %.1f, replaced %d flux bins between MLT %.1f and %.1f with interpolated flux..." % (this_mlat,
                #       np.count_nonzero(interp_bins_missing_flux),np.nanmin(this_mlt[interp_bins_missing_flux]),
                #       np.nanmax(this_mlt[interp_bins_missing_flux]))

        #notwedge = np.logical_not(inwedge)
        #crazy = np.logical_and(inwedge,fluxgridN>np.nanpercentile(fluxgridN[notwedge],90.))
        #fluxgridN[crazy]=np.nanpercentile(fluxgridN[notwedge],90.)

        return fluxgridN, inwedge

