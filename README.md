AURORAMAPS
==========

This is an open-source version of the OVATION Prime 2010 (OP10) aurora model in python.

Current status (October 2019): ongoing debugging.

The solar wind input comes either from OMNI2 historic data or from the PREDSTORM L1 solar wind forecast, which is based on data from the DSCOVR and STEREO-A spacecraft. The results are plotted with cartopy on different types of world maps.

by C. Moestl, IWF-helio group, Graz, Austria. https://www.iwf.oeaw.ac.at/en/user-site/christian-moestl/

Contributions by Liam Kilcommons and Diana Morosan

OVATION has been largely rewritten based on https://github.com/lkilcommons/OvationPyme

Python packages used on top of a python 3.7 anaconda installation: cartopy 0.17.0, sunpy 1.0.3, scikit-learn 0.20.3, aacgmv2 2.5.1. For their installation, see below.

Note that for installing cartopy you need an anaconda installation.

Run time optimization with the numba and multiprocessing packages is implemented.




Installation
------------

For the dependent packages:

    pip install cartopy
    pip install aacgmv2
    conda install sunpy
    conda install scikit-learn

when you have them:

    git clone https://github.com/IWF-helio/auroramaps
    cd auroramaps
    python setup.py install


Usage
-----

Run this file - the input variables are controlled in input.py and can be changed there:

    python aurora.py


This produces a folder in the results directory named as given in input.py that contains aurora movies (gif, mp4) and frames for the event, and a plot on the Newell coupling. For the animation, it is assumed that ffmpeg (https://ffmpeg.org/download.html) is available on the command line.



Documentation
-------------

3 modes are available, selected in input.py

0 - uses the real time PREDSTORM solar wind prediction
1 - uses a local file of the PREDSTORM solar wind output
2 - for using any interval in the OMNI2 data since 1963


In mode 2, the OMNI2 data are downloaded automatically from https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_all_years.dat on the first time, and then converted to a pickle for faster processing in new runs. If you want to update the OMNI2 data, just delete both files in "auroramaps/data/omni2/" to force a new download and conversion.



![Sample image](https://raw.githubusercontent.com/cmoestl/auroramaps/samples/sample_polar_north.png)





