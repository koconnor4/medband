#! /usr/bin/env python
# 2014.06.29 S. Rodney

import os
from matplotlib import pylab as pl
import numpy as np
import stardust
from circlefig import plotcontours

_datadir = os.path.abspath( '.' )
_RA, _DEC = 189.31989, 62.278179

def getPhot( ra=_RA, dec=_DEC, datadir=_datadir,
             snanastyle=True, verbose=False ):
    """ Measure aperture photometry for SN Stone from diff images.
    """
    import glob
    from hstphot import hstapphot

    topdir = os.path.abspath( '.' )
    os.chdir( datadir )

    stackfile = 'stone_f160w_e05-e00_sub_masked.fits'
    xc,yc = hstapphot.getxycenter( stackfile, ra, dec, radec=True, verbose=1)
    rac,decc = hstapphot.xy2radec( stackfile, xc, yc )

    if verbose :
        print( "Recentered to x,y= %.1f %.1f on %s"%(xc,yc,stackfile))
        print( "User Input RA,Dec = %.6f  %.6f"%(ra,dec) )
        print( "Recentered RA,Dec = %.6f  %.6f"%(rac,decc) )

    if verbose and snanastyle :
        print("VARLIST:  MJD  FLT FIELD  FLUXCAL  FLUXCALERR   MAG    MAGERR")
    elif verbose :
        print( " MJD  FILT  APER  MAG  MAGERR   FLUX FLUXERR   FSKY  FSKYERR  " )

    magdict = {}
    sublist = glob.glob( '*sub_masked.fits')
    for subfile in sublist :
        camera = hstapphot.getcamera(subfile)
        if camera == 'ACS-WFC' : aperture = 0.5
        elif camera == 'WFC3-UVIS' : aperture = 0.15
        elif camera == 'WFC3-IR' : aperture = 0.3
        else : continue
        xc,yc = hstapphot.radec2xy( subfile, rac, decc )
        mag,magerr,magstring = hstapphot.getmag(
            subfile, xc, yc, aperture, system='Vega', verbose=verbose>2, snanastyle=snanastyle )
        key = subfile.split('_')[1] +'.' + subfile.split('_')[2]
        magdict[ key ] = [mag,magerr]
        print(magstring[0])

    os.chdir( topdir )

    return( magdict )


def dosim( nsim=2000, clobber=0, verbose=False ):
    import snanasim
    sn = stardust.SuperNova('HST_CANDELS2_stone.dat')
    simname = snanasim.dosimMC( sn, Nsim=nsim, clobber=clobber, verbose=verbose )
    simdata = snanasim.readSimDataMC(simname=simname)
    return( simdata )

def mkfig( simdata=None, linelevels = [ 0, 0.82 ], plotstyle='contourf',
           Nbins=80, nsim=2000, clobber=0, verbose=1  ):
    if simdata is None :
        simdata = dosim( nsim=nsim, clobber=clobber, verbose=verbose )

    pl.clf()
    sn = stardust.SuperNova('HST_CANDELS2_stone.dat')
    simIa, simIbc, simII = simdata

    mjdmedband = sn.MJD[ np.where( (sn.FLT=='P') | (sn.FLT=='Q') ) ]
    mjdobs = np.median( mjdmedband )

    print('Binning up MC sim for color-color diagram...')
    stardust.simplot.plotColorColor( simIbc, 'Q-H','P-N', binrange=[[-0.4,0.5],[-0.4,0.5]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
    stardust.simplot.plotColorColor( simII,  'Q-H','P-N', binrange=[[-0.4,0.5],[-0.4,0.5]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
    stardust.simplot.plotColorColor( simIa,  'Q-H','P-N', binrange=[[-0.4,0.5],[-0.4,0.5]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )

    plotPhot()
    pl.draw()
    return(simdata)

def getzpoints( sn, simIa, dz=0.1 ):
    """ Bin up the simulated Ia's by redshift and determine the point in
    color-color space ovvupied by the median member of each bin.
    """
    if 'FLT' not in simIa.__dict__ :
        simIa.getLightCurves()
    zbins = np.arange( sn.z - sn.zerr, sn.z + sn.zerr + dz, dz )
    izbinarg = np.digitize( simIa.z, zbins )
    iQ = np.where( simIa.FLT=='Q' )
    iH = np.where( simIa.FLT=='H' )
    iP = np.where( simIa.FLT=='P' )
    iN = np.where( simIa.FLT=='N' )
    Q,H,P,N = simIa.MAG[iQ],simIa.MAG[iH],simIa.MAG[iP],simIa.MAG[iN]
    QHpts = np.array([ (Q-H)[izbinarg==iz] for iz in range(len(zbins))])
    PNpts = np.array([ (P-N)[izbinarg==iz] for iz in range(len(zbins))])
    return( zbins, QHpts, PNpts )

def plotzpoints(sn,simIa,dz=0.1):
    """  plot the
    :param sn:
    :param simIa:
    :param dz:
    :return:
    """
    from matplotlib import cm
    from mpltools import color
    from pytools import plotsetup

    if 'FLT' not in simIa.__dict__ :
        simIa.getLightCurves()
    zbins = np.arange( sn.z - sn.zerr, sn.z + sn.zerr + dz, dz )
    izbinarg = np.digitize( simIa.z, zbins, right=True )-1
    iQ = np.where( simIa.FLT=='Q' )
    iH = np.where( simIa.FLT=='H' )
    iP = np.where( simIa.FLT=='P' )
    iN = np.where( simIa.FLT=='N' )
    Q,H,P,N = simIa.MAG[iQ],simIa.MAG[iH],simIa.MAG[iP],simIa.MAG[iN]
    QH = Q-H
    PN = P-N

    plotsetup.fullpaperfig( 1, [5,5] )
    pl.clf()
    ax1 = pl.gca()
    # ax1.plot(QH,PN, ls=' ', marker='o',color='k', alpha=0.3 )
    color.cycle_cmap( len(zbins), cmap=cm.gist_rainbow_r, ax=ax1)
    for iz in range(len(zbins)-1) :
        z0, z1 = zbins[iz],zbins[iz+1]
        QHz = QH[izbinarg==iz]
        PNz = PN[izbinarg==iz]
        if dz>=0.1 : zmid = round(np.mean([z0,z1]),1)
        elif dz>=0.01 : zmid = round(np.mean([z0,z1]),2)
        ax1.plot(QHz,PNz, marker='o', alpha=0.3, ls=' ',mew=0,label='%.1f-%.1f'%(z0,z1) )
    ax1.legend(loc='upper left', ncol=2, numpoints=1, frameon=False,
               handletextpad=0.3, handlelength=0.2 )
    ax1.set_xlabel('F153M-F160W')
    ax1.set_ylabel('F139M-F140W')
    plotPhot(label=True)
    ax1.set_xlim(-0.35,0.53)
    ax1.set_ylim(-0.35,0.53)
    fig = pl.gcf()
    fig.subplots_adjust(bottom=0.12,left=0.16,right=0.95,top=0.95)
    pl.draw()


def plotPhot(label=False):
    """ Plot the SN Stone Med Band Photometry on a color-color diagram
    """
    # OBS: 56475.01   Q   stone   16.990    1.424      24.180    0.091
    # OBS: 56474.53   H   stone   19.852    0.989      23.999    0.054

    # F153M - F160W = Q - H
    QH = 24.180 - 23.999
    QHerr = np.sqrt( 0.091**2 + 0.054**2)

    # OBS: 56474.48   P   stone   15.884    1.444      24.268    0.099
    # OBS: 56474.54   N   stone   16.704    1.012      24.212    0.066
    # F139M - F140W = P - N
    PN = 24.268 - 24.212
    PNerr = np.sqrt( 0.099**2 + 0.066**2)

    print("F153M-F160W = %.2f +- %.2f"%(QH,QHerr))
    print("F139M-F140W = %.2f +- %.2f"%(PN,PNerr))
    ax = pl.gca()
    ax.errorbar( QH, PN, PNerr, QHerr, color='k', marker='o', ms=15, capsize=0 )
    if label :
        ax.text( 0.95,0.95,'GND13Sto', transform=ax.transAxes, ha='right',va='top',fontsize='large')



def sncosmo_circlefig( simIa=None, simCC=None,
                       simIapkl='stone_SncosmoSim_Ia.pkl',
                       simCCpkl='stone_SncosmoSim_CC.pkl',
                       z_range=[1.6,2.0], nsim=1000,
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
    mjdpk = 56482.
    mjdmedband = 56475.

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

    magerr = { # from drop method
        'f139m':0.182164 ,
        'f140w':0.0920638,
        'f153m':0.185014 ,
        'f160w':0.0731440,
        }

    mag = { # from idl5 psf fitting
        'f139m':25.4198,
        'f140w':25.3622,
        'f153m':25.5056,
        'f160w':25.2278,
        }


    f25 = {}
    ferr25 = {}
    for k in mag.keys() :
        f25[k] = 10**(-0.4*(mag[k]-25.))
        ferr25[k] = magerr[k] * f25[k] / 2.5*np.log10(np.e)

    deltamag = {
                 'f139m':mag['f139m']-mag['f140w'],
                 'f153m':mag['f153m']-mag['f160w'],
                 }
    deltamagerr = {
                    'f139m':np.sqrt(magerr['f139m']**2+magerr['f140w']**2),
                    'f153m':np.sqrt(magerr['f153m']**2+magerr['f160w']**2),
                    }

    # ax1.errorbar( deltamag['f139m'], deltamag['f127m'],
    #              deltamagerr['f127m'], deltamagerr['f139m'],
    #              marker='D', ms=10, elinewidth=2, capsize=0, color='darkorange' )
    ax2.errorbar( deltamag['f139m'], deltamag['f153m'],
                  deltamagerr['f153m'], deltamagerr['f139m'],
                  marker='D', ms=10, elinewidth=2, capsize=0, color='darkorange' )

    ax1.set_xlim( -0.35, 0.3 )
    ax1.set_ylim( -0.7, 0.22 )
    ax2.set_ylim( -0.4, 0.3 )

    ax2.text(  deltamag['f139m']+0.05, deltamag['f153m']-0.05, 'GND13Sto' ,
               ha='left',va='top', color='darkorange',fontsize=15,)
               # path_effects=[pe.withStroke( linewidth=3,foreground='k')] )
    # pl.legend( loc='upper right', numpoints=2, handlelength=0.3)
    pl.draw()
    pl.ion()
    return( simIa, simCC )
