#! /usr/bin/env python
# 2014.07.05 S. Rodney
__author__ = 'rodney'

import os
from matplotlib import pylab as pl
import numpy as np
import stardust

_datadir = os.path.abspath( '.' )
_RA, _DEC = 53.178264, -27.801989

_STACKFILE = 'bush_JH_stack-e00_sub_masked.fits'
_SNANADATFILE = 'HST_CANDELS2_bush.dat'

def dophot( ra=_RA, dec=_DEC, datadir=_datadir, stackfile=_STACKFILE,
             snanastyle=True, verbose=False ):
    """ Measure aperture photometry for SN Camille from diff images.
    """
    from phot import getmags
    getmags( ra, dec, datadir, stackfile, snanastyle, verbose )



