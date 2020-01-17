# auroramaps

This is an open-source version of the OVATION Prime 2010 (OP10) aurora model in python.

Current status (January 2020): **Work in progress!** Debugging is ongoing for OP10, and updating to OP13 is underway.

by C. Moestl, Rachel L. Bailey, IWF-helio group, Graz, Austria. https://www.iwf.oeaw.ac.at/en/user-site/christian-moestl/  
Contributions by  Diana E. Morosan and Liam Kilcommons.

If you want to use parts of this code for generating results for peer-reviewed scientific publications, please contact me per email (christian.moestl@oeaw.ac.at) or via twitter @chrisoutofspace.

The solar wind input comes either from OMNI2 historic data or from the PREDSTORM L1 solar wind forecast, which is based on data from the DSCOVR or ACE and STEREO-A spacecraft (https://github.com/IWF-helio/PREDSTORM). The results are plotted with cartopy on different types of world maps. 
OVATION has been largely rewritten based on https://github.com/lkilcommons/OvationPyme and open source versions thankfully made available by NOAA and the UK MetOffice.

## Installation

I use a python 3.7.6 miniconda installation with standard packages (numpy, numba, scipy, matplotlib) as well as cartopy (0.17.0), aacgmv2 (2.6.0), and heliosat (0.3.1). It is recommended to make a conda environment to run the code. For full installation instructions with such an environment see: 
https://github.com/helioforecast/helio4cast

After installation, do in a directory of your choice:

    git clone https://github.com/IWF-helio/auroramaps

*auroramaps* is currently not distributed via PyPi.

## Usage

Run this file - the input variables are controlled in input.py and can be changed there:

    python aurora.py

For matplotlib backend 'Agg' use this:
    
    python aurora.py --server

For a real time version use this (note the solar wind input file needs to be specified in *input_realtime.py*):

    python aurora.py --real

    
This produces a folder in the results directory named as given in *input.py* or *input_realtime.py* that contains aurora movies (gif, mp4) and frames for the event, and a plot on the Newell coupling. For the animation, it is assumed that ffmpeg (https://ffmpeg.org/download.html) is available on the command line.



## Documentation


3 modes are available, selected in input.py

0 - uses the real time PREDSTORM solar wind prediction  
1 - uses a local file of the PREDSTORM solar wind output  
2 - for using any interval in the OMNI2 data since 1963

In mode 2, the OMNI2 data are downloaded automatically from https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_all_years.dat on the first time, and then converted to a pickle for faster processing in new runs. If you want to update the OMNI2 data, just delete both files in "auroramaps/data/omni2/" to force a new download and conversion.

There are 3 types of maps available - a global image of the northern hemisphere and high-resolution maps of Europe and North-America. 
For both maps, flux and probability images can be made, andt here are 3 different background images available (blue marble, VIIRS night band and a topography image).
Control all this with input.py. 

These are samples with blue marble background for the flux map (global polar view) and the probability map (North America and Europe):

![Sample image](samples/global_flux_sample.jpg)
![Sample image](samples/canada_prob_sample.jpg)
![Sample image](samples/europe_prob_sample.jpg)

