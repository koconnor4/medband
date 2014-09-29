#! /usr/bin/env python
# S.Rodney
# 2014.05.04
from circlefig import plotcontours

_SNANADATFILE = 'HST_CANDELS2_colfax.dat'

def mkTemplates( clobber=False, verbose=False):
    """ Make med-band templates from sndrizpipe epoch 00 broad-band images.
    """
    from hstsntools import imageops, filters

    f125to127 = filters.computeScaling( 'F125W','F127M')
    print("F127M ~ F125W * %.3f"%f125to127)
    f160to153 = filters.computeScaling( 'F160W','F153M')
    print("F153M ~ F160W * %.3f"%f160to153)

    f125and160to140 = filters.computeScaling2to1( 'F125W','F160W','F140W')
    print("F140W ~ (F125W+F160W) * %.3f"%f125and160to140)

    f140to139 = filters.computeScaling( 'F140W','F139M')
    print("F139M ~ F140W * %.3f"%f140to139)

    # make epoch 00 med-band templates
    imageops.imscaleflux( "colfax.e00/colfax_f125w_e00_reg_drz_sci.fits",
                          f125to127,
                          outfile='colfax.e00/colfax_f127m_e00_reg_drz_sci.fits',
                          clobber=clobber, verbose=verbose)
    imageops.imscaleflux( "colfax.e00/colfax_f160w_e00_reg_drz_sci.fits",
                          f160to153,
                          outfile='colfax.e00/colfax_f153m_e00_reg_drz_sci.fits',
                          clobber=clobber, verbose=verbose)

    imageops.imsum( "colfax.e00/colfax_f125w_e00_reg_drz_sci.fits",
                    "colfax.e00/colfax_f160w_e00_reg_drz_sci.fits",
                    "colfax.e00/colfax_f125w+f160w_e00_reg_drz_sci.fits",
                    clobber=clobber, verbose=verbose )

    imageops.imscaleflux( "colfax.e00/colfax_f125w+f160w_e00_reg_drz_sci.fits",
                          f125and160to140,
                          outfile='colfax.e00/colfax_f140w_e00_reg_drz_sci.fits',
                          clobber=clobber, verbose=verbose)
    imageops.imscaleflux( "colfax.e00/colfax_f140w_e00_reg_drz_sci.fits",
                          f140to139,
                          outfile='colfax.e00/colfax_f139m_e00_reg_drz_sci.fits',
                          clobber=clobber, verbose=verbose)

def mkSubs( clobber=False, verbose=False ):
    """ Make med-band diff images.
    """
    from hstsntools import imageops

    for band in ['f127m','f139m','f153m']:
        template='colfax.e00/colfax_%s_e00_reg_drz_sci.fits'%band
        snimage='colfax.e03/colfax_%s_e03_reg_drz_sci.fits'%band
        diffimage='colfax.e03/colfax_%s_e03-e00_sub_sci.fits'%band

        imageops.imsubtract( template, snimage, diffimage,
                             clobber=clobber, verbose=verbose )

    template='colfax.e00/colfax_f140w_e00_reg_drz_sci.fits'
    for epoch in [3,4,5]:
        snimage='colfax.e%02i/colfax_f140w_e%02i_reg_drz_sci.fits'%(epoch,epoch)
        diffimage='colfax.e%02i/colfax_f140w_e%02i-e00_sub_sci.fits'%(epoch,epoch)
        imageops.imsubtract( template, snimage, diffimage,
                             clobber=clobber, verbose=verbose )



def pseudocolor_vs_z( simgridIa=None, clobber=0 ):
    """  Plot medium-broad pseudo-colors from a grid simulation as a function
     of redshift.

    :param simgridIa: SNANA Simulation Table, or None to make/load it anew
    :param clobber:  passed to snanasim.dosimGrid() to re-run the SNANA sims
    :return: new or existing SNANA sim table (a stardust.SimTable object)
    """
    import stardust
    import snanasim
    import mkplots
    from matplotlib import pyplot as pl
    sn = stardust.SuperNova('HST_CANDELS2_colfax.dat')

    if simgridIa is None :
        if clobber :
            simgridIa = snanasim.dosimGrid(sn,ngridz=20, clobber=clobber, x1range=[-2,2], crange=[-0.2,0.5], trestrange=[-5,5] )
        else :
            simgridIa = stardust.SimTable( 'sim_colfax_medbandGrid_Ia')
    mkplots.pseudocolor_vs_z( simgridIa, medbands='OPQ', broadbands='JNH')
    pl.suptitle('MED BAND GRID SIM FOR SN COLFAX @ z=2.1+- 0.2', fontsize=20)
    return( simgridIa )


def circlefig( sn=None, simdataMC=None, simIaGrid=None, fast=True, clobber=False ):
    import mkplots
    if sn is None :
        sn = _SNANADATFILE
    sn,simdataMC,simIaGrid = mkplots.circlefig(
        sn, simdataMC, simIaGrid,
        plotstyle='points', showGridline=True,
        color1='O-J',color2='P-N',color3='Q-H',
        fast=fast, clobber=clobber )
    return( sn,simdataMC,simIaGrid  )



def circlefigPoints( ):
    """  Add measured medband photometry to an existing circle diagram in the
    current figure.
    :return:
    """
    import numpy as np
    from matplotlib import pyplot as pl

    OJ, OJerr =  25.502 - 26.254, np.sqrt( 0.199**2 + 0.125**2) # F127M - F125W
    PN, PNerr =  24.963 - 25.078, np.sqrt( 0.082**2 + 0.052**2) # F139M - F140W
    QH, QHerr =  25.105 - 25.147, np.sqrt( 0.116**2 + 0.113**2) # F153M - F160W

    fig = pl.gcf()
    ax1 = fig.add_subplot( 2,2,1 )
    pt1 = ax1.errorbar( OJ, PN, PNerr, OJerr, marker='D', ms=12, capsize=0, elinewidth=1, color='k' )

    ax2 = fig.add_subplot( 2,2,2 )
    pt2 = ax2.errorbar( OJ, QH, QHerr, OJerr, marker='D', ms=12, capsize=0, elinewidth=1, color='k' )

    ax3 = fig.add_subplot( 2,2,3 )
    pt3 = ax3.errorbar( OJ, QH, QHerr, OJerr, marker='D', ms=12, capsize=0, elinewidth=1, color='k' )
    return( pt1, pt2, pt3 )


def mkcirclefigGrid( simdataGrid=None, clobber=0,
                     snanadatfile=_SNANADATFILE ):
    """  Plot the results of a SNANA Grid simulation as a circle diagram.
    """
    import snanasim
    import mkplots
    import stardust
    import numpy as np
    from matplotlib import pyplot as pl

    sn = stardust.SuperNova(snanadatfile)
    if simdataGrid is None :
        if clobber :
            simdataGrid = snanasim.dosimGrid(sn,ngridz=20, clobber=clobber, x1range=[-2,2], crange=[-0.2,0.5], trestrange=[-5,5] )
        else :
            simdataGrid = stardust.SimTable( 'sim_colfax_medbandGrid_Ia')

    mjdmedband = sn.MJD[ np.where( (sn.FLT=='7') | (sn.FLT=='8') | (sn.FLT=='P')) ]
    if len(mjdmedband)>0 :
        mjdobs = np.median( mjdmedband )
    else :
        mjdobs = sn.pkmjd

    fig = pl.gcf()
    ax1 = fig.add_subplot( 2,2,1 )
    mkplots.gridsim_circleplot( simdataGrid,  'O-J', 'P-N' )
    ax2 = fig.add_subplot( 2,2,2 )
    mkplots.gridsim_circleplot( simdataGrid, 'Q-H', 'P-N' )
    ax2 = fig.add_subplot( 2,2,3 )
    mkplots.gridsim_circleplot( simdataGrid, 'O-J', 'Q-H' )

    pl.draw()
    return simdataGrid


def sncosmo_circlefig( simIa=None, simCC=None,
                          simIapkl='colfax_SncosmoSim_Ia.pkl',
                          simCCpkl='colfax_SncosmoSim_CC.pkl',
                          z_range=[1.9,2.35], nsim=1000,
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
    mjdpk = 56078
    mjdmedband = 56084

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

    plotcontours( simIa, simCC, plotstyle='points' )
    fig = pl.gcf()
    ax1 = fig.add_subplot(1,2,1)
    mkcirclepoints(zrange=z_range, t0=t0, colorselect=[1,0], coloredaxislabels=False, marker='o' )
    ax2 = fig.add_subplot(1,2,2, sharex=ax1)
    mkcirclepoints(zrange=z_range, t0=t0, colorselect=[1,2], coloredaxislabels=False, marker='o' )



    mag4 = {  # using the multi-epoch stack for the psf model
        'f127m':25.6728,
        'f125w':26.1096,
        'f139m':25.5964,
        'f140w':25.7223,
        'f153m':25.7903,
        'f160w':25.9248,
        }

    magerr = { # using the drop-and-recover psf fitting err method
        'f125w':0.126760,
        'f127m':0.102688,
        'f139m':0.132692,
        'f140w':0.043144,
        'f153m':0.170073,
        'f160w':0.117535,
        }

    mag2 = { # Using steve's 0.3"-aperture photometry
    'f127m':25.668,
    'f125w':26.172,
    'f139m':25.601,
    'f140w':25.747,
    'f153m':25.769,
    'f160w':25.958,
    }

    magerr2 = { # Using steve's 0.3"-aperture photometry
        'f127m':0.095,
        'f125w':0.118,
        'f139m':0.096,
        'f140w':0.063,
        'f153m':0.125,
        'f160w':0.129,
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
    ax2.errorbar( deltamag['f139m'], deltamag['f153m'],
                  deltamagerr['f153m'], deltamagerr['f139m'],
                  marker='D', ms=10, elinewidth=2, capsize=0, color='darkorange' )

    ax1.set_xlim( -0.35, 0.3 )
    ax1.set_ylim( -0.7, 0.22 )
    ax2.set_ylim( -0.4, 0.3 )

    ax2.text(  deltamag['f139m']+0.05, deltamag['f153m']-0.05, 'GND12Col' ,
               ha='left',va='top', color='darkorange',fontsize=15,)
               # path_effects=[pe.withStroke( linewidth=3,foreground='k')] )
    # pl.legend( loc='upper right', numpoints=2, handlelength=0.3)
    pl.draw()
    pl.ion()
    return( simIa, simCC )