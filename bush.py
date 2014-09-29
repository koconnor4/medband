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


def sncosmo_circlefig( simIa=None, simCC=None,
                       simIapkl='bush_SncosmoSim_Ia.pkl',
                       simCCpkl='bush_SncosmoSim_CC.pkl',
                       z_range=[1.16,2.36], nsim=1000,
                       verbose=True, clobber=False ):
    """  Construct a color-color circle figure for SN Colfax, with observed
     photometry included.

    :param simIa:
    :param simCC:
    :param simIapkl:
    :param simCCpkl:
    :param z_range:
    :param nsim:
    :param verbose:
    :param clobber:
    :return:
    """

    import medband_classtest
    import os
    import cPickle
    from matplotlib import pyplot as pl
    # from matplotlib import patheffects as pe
    import numpy as np
    from pytools import plotsetup
    fig = plotsetup.fullpaperfig( 1, figsize=[8,4])

    pl.ioff()
    mjdpk = 55797.
    mjdmedband = 55804.

    t0_range = [ mjdmedband-mjdpk-3,mjdmedband-mjdpk+3 ]
    t0 = mjdmedband-mjdpk

    if simIa is not None :
        pass
    elif os.path.isfile( simIapkl ) and not clobber>1 :
        if verbose: print("Loading Ia simulation from pickle : %s"%simIapkl)
        fin = open( simIapkl, 'rb' )
        simIa = cPickle.load( fin )
        fin.close()
    else :
        if verbose: print("Running a new Ia simulation, then saving to pickle : %s"%simIapkl)
        simIa = medband_classtest.SncosmoSim( 'Ia' , z_range=z_range, t0_range=t0_range, nsim=nsim )
        fout = open( simIapkl, 'wb' )
        cPickle.dump( simIa, fout, protocol=-1 )
        fout.close()

    if simCC is not None :
        pass
    elif os.path.isfile( simCCpkl ) and not clobber>1 :
        if verbose: print("Loading CC simulation from pickle : %s"%simCCpkl)
        fin = open( simCCpkl, 'rb' )
        simCC = cPickle.load(fin)
        fin.close()
    else :
        if verbose: print("Running a new CC simulation, then saving to pickle : %s"%simCCpkl)
        simCC = medband_classtest.SncosmoSim( 'CC' , z_range=z_range, t0_range=t0_range, nsim=nsim )
        fout = open( simCCpkl, 'wb' )
        cPickle.dump( simCC, fout, protocol=-1 )
        fout.close()

    import getredshift

    # classify.plotcontours( simIa, simCC, plotstyle='points' )
    fig = pl.gcf()
    ax1 = fig.add_subplot(1,2,1)
    mkcirclepoints(zrange=z_range, t0=t0, colorselect=[1,0], coloredaxislabels=False, marker='o' )
    ax2 = fig.add_subplot(1,2,2, sharex=ax1)
    mkcirclepoints(zrange=z_range, t0=t0, colorselect=[0,2], coloredaxislabels=False, marker='o' )

    mag4 = {  # from psf model based on the multi-epoch stack
              'f125w':25.9432,
              'f127m':25.9555,
              'f139m':25.5784,
              'f140w':26.2892,
              'f153m':25.3580,
              'f160w':25.9067,
              }

    mag6 = { # from psf model based on the multi-epoch stack
             'f125w':25.9464,
             'f127m':25.9610,
             'f139m':25.6361,
             'f140w':26.2892,
             'f153m':25.3578,
             'f160w':25.9097,
             }

    magerr = { # from drop method
         'f125w':0.106712,
         'f127m':0.184979,
         'f139m':0.214549,
         'f140w':0.1162	 ,
         'f153m':0.160024,
         'f160w':0.138843,
         }
    mag = mag6

    f25 = {}
    ferr25 = {}
    for k in mag.keys() :
        f25[k] = 10**(-0.4*(mag[k]-25.))
        ferr25[k] = magerr[k] * f25[k] / 2.5*np.log10(np.e)

    deltamag = { 'f127m':mag['f127m']-mag['f125w'],
                 'f139m':mag['f139m']-mag['f140w'],
                 'f153m':mag['f153m']-mag['f160w'],
                 }
    deltamagerr = { 'f127m':np.sqrt(magerr['f127m']**2+magerr['f125w']**2),
                    'f139m':np.sqrt(magerr['f139m']**2+magerr['f140w']**2),
                    'f153m':np.sqrt(magerr['f153m']**2+magerr['f160w']**2),
                    }

    ax1.errorbar( deltamag['f139m'], deltamag['f127m'],
                  deltamagerr['f127m'], deltamagerr['f139m'],
                  marker='D', ms=10, elinewidth=2, capsize=0, color='darkorange' )
    ax2.errorbar( deltamag['f127m'], deltamag['f153m'],
                  deltamagerr['f153m'], deltamagerr['f127m'],
                  marker='D', ms=10, elinewidth=2, capsize=0, color='darkorange' )

    ax1.set_xlim( -0.35, 0.3 )
    ax1.set_ylim( -0.7, 0.22 )
    ax2.set_ylim( -0.4, 0.3 )

    ax2.text(  deltamag['f139m']+0.05, deltamag['f153m']-0.05, 'GND12Bus' ,
               ha='left',va='top', color='darkorange',fontsize=15,)
               # path_effects=[pe.withStroke( linewidth=3,foreground='k')] )
    # pl.legend( loc='upper right', numpoints=2, handlelength=0.3)
    pl.draw()
    pl.ion()
    return( simIa, simCC )


