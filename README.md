AURORAMAPS
==========

This is an open-source version of the OVATION Prime 2010 (OP10) aurora model in python.

ONGOING DEBUGGING!

The solar wind input comes either from OMNI2 historic data or from the PREDSTORM L1 solar wind forecast, which is based on data from the DSCOVR and STEREO-A spacecraft. The results are plotted with cartopy on different types of world maps.

Status: work in progress, October 2019, C. Moestl/IWF-helio 
https://www.iwf.oeaw.ac.at/en/user-site/christian-moestl/

OVATION has been largely rewritten based on https://github.com/lkilcommons/OvationPyme

Run time optimization with the numba and multiprocessing packages is implemented.


Installation
------------

    git clone https://github.com/IWF-helio/auroramaps
    cd auroramaps
    python setup.py install


Demo
----

Run this file, controlled by input.py:

    python aurora.py







