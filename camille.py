#! /usr/bin/env python
# 2014.07.05 S. Rodney

import os
from matplotlib import pylab as pl
import numpy as np
import stardust

_datadir = os.path.abspath( '.' )
_RA, _DEC = 189.2807,  62.174123

_STACKFILE = 'camille_JNH_stack-e00_sub_masked.fits'
_SNANADATFILE = 'HST_CANDELS2_camille.dat'

def getPhot( ra=_RA, dec=_DEC, datadir=_datadir, stackfile=_STACKFILE,
             snanastyle=True, verbose=False ):
    """ Measure aperture photometry for SN Camille from diff images.
    """
    import glob
    from hstphot import hstapphot

    topdir = os.path.abspath( '.' )
    os.chdir( datadir )

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
    sublist = glob.glob( '*e??-e??_sub_masked.fits')
    for subfile in sublist :
        camera = hstapphot.getcamera(subfile)
        if camera == 'ACS-WFC' : aperture = 0.5
        elif camera == 'WFC3-UVIS' : aperture = 0.15
        elif camera == 'WFC3-IR' : aperture = 0.3
        else : continue
        xc,yc = hstapphot.radec2xy( subfile, rac, decc )
        mag,magerr,magstring = hstapphot.dophot(
            subfile, xc, yc, aperture, system='AB', skyannarcsec=[6.0,10.0], verbose=verbose>2, snanastyle=snanastyle )
        key = subfile.split('_')[1] +'.' + subfile.split('_')[2]
        magdict[ key ] = [mag,magerr]
        print(magstring[0])

    os.chdir( topdir )

    return( magdict )


def dosim( nsim=2000, clobber=0, verbose=False,
           snanadatfile=_SNANADATFILE, smear=False ):
    """ Run the MC simulation.
    Set smear=False to get no magnitude smearing.  For Ia this means no range
    of x1 or c, and for both Ia and CC it means no intrinsic mag smearing.
    """
    import medbandfig
    sn = stardust.SuperNova( snanadatfile )
    sn.REDSHIFT_FINAL = '1.2 +- 0.2'
    simname = medbandfig.doMedBandSim( sn, Nsim=nsim, smear=smear, clobber=clobber, verbose=verbose )
    simdata = medbandfig.readSimData(simname=simname)
    return( simdata )

def mkfig( simdata=None, linelevels = [ 0, 0.82 ], plotstyle='contourf',
           Nbins=80, nsim=2000, clobber=0, verbose=1,
           snanadatfile=_SNANADATFILE ):
    if simdata is None :
        simdata = dosim( nsim=nsim, clobber=clobber, verbose=verbose )

    pl.clf()
    sn = stardust.SuperNova(snanadatfile)
    simIa, simIbc, simII = simdata

    mjdmedband = sn.MJD[ np.where( (sn.FLT=='7') | (sn.FLT=='8') | (sn.FLT=='P')) ]
    mjdobs = np.median( mjdmedband )


    print('Binning up MC sim for color-color diagram...')
    ax1 = pl.subplot( 2,2,1 )
    stardust.simplot.plotColorColor( simIa,  '7-I','P-N', binrange=[[-0.2,0.8],[-0.5,0.5]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
    stardust.simplot.plotColorColor( simIbc, '7-I','P-N', binrange=[[-0.2,0.8],[-0.5,0.5]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
    stardust.simplot.plotColorColor( simII,  '7-I','P-N', binrange=[[-0.2,0.8],[-0.5,0.5]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
    ax2 = pl.subplot( 2,2,2 )
    stardust.simplot.plotColorColor( simIa,   '8-I','P-N', binrange=[[-0.6,0.4],[-0.5,0.5]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
    stardust.simplot.plotColorColor( simIbc,  '8-I','P-N', binrange=[[-0.6,0.4],[-0.5,0.5]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
    stardust.simplot.plotColorColor( simII,   '8-I','P-N', binrange=[[-0.6,0.4],[-0.5,0.5]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
    ax3 = pl.subplot( 2,2,3 )
    stardust.simplot.plotColorColor( simIa,   '7-I','8-I', binrange=[[-0.2,0.8],[-0.6,0.4]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
    stardust.simplot.plotColorColor( simIbc,  '7-I','8-I', binrange=[[-0.2,0.8],[-0.6,0.4]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
    stardust.simplot.plotColorColor( simII,   '7-I','8-I', binrange=[[-0.2,0.8],[-0.6,0.4]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )

    pl.draw()
    return(simdata)


def mkfig3d( simdata, dz=0.05, snanadatfile=_SNANADATFILE, ):
    from mpl_toolkits.mplot3d.axes3d import Axes3D
    from matplotlib import cm
    import mpltools
    from pytools import plotsetup
    fig = plotsetup.fullpaperfig( 1, [10,10] )

    sn = stardust.SuperNova(snanadatfile)
    simIa, simIbc, simII = simdata
    ax1 = fig.add_axes( [0.12,0.12,0.85,0.85], projection='3d')

    c7I,c8I, cPN = get_colors( simII )
    ax1.plot3D(c7I,c8I,cPN, marker='o', color='k',
               alpha=1, ls=' ', mew=0,label='Type II' )

    c7I,c8I, cPN = get_colors( simIbc )
    ax1.plot3D(c7I,c8I,cPN, marker='o', color='sienna',
               alpha=1, ls=' ', mew=0,label='Type Ibc' )

    zbins, c7I,c8I, cPN = get_colors_binned_by_z( sn, simIa, dz=dz )
    mpltools.color.cycle_cmap( len(zbins), cmap=cm.gist_rainbow_r, ax=ax1)
    for iz in range(len(zbins)-1) :
        z0, z1 = zbins[iz],zbins[iz+1]
        #if dz>=0.1 : zmid = round(np.mean([z0,z1]),1)
        #elif dz>=0.01 : zmid = round(np.mean([z0,z1]),2)
        ax1.plot3D(c7I[iz],c8I[iz],cPN[iz], marker='o', alpha=1, ls=' ',
                 mew=0,label='%.2f-%.2f'%(z0,z1) )

    ax1.set_xlim(-0.1,0.5)
    ax1.set_ylim(-0.4,0.2)
    ax1.set_zlim(-0.3,0.3)


    ax1.legend(loc='upper left', ncol=2, numpoints=1, frameon=False,
               handletextpad=0.3, handlelength=0.2 )
    ax1.set_xlabel('F763M-F814W')
    ax1.set_ylabel('F845M-F140W')
    ax1.set_zlabel('F139M-F140W')


    #ax1.plot3D( c7I[0], c8I[0], zs=cPN[0], zdir='z' )



def get_colors_binned_by_z( sn, simIa, dz=0.1 ):
    """ Bin up the simulated Ia's by redshift and determine the point in
    color-color space occupied by the median member of each bin.
    """
    if 'FLT' not in simIa.__dict__ :
        simIa.getLightCurves()
    zbins = np.arange( sn.z - sn.zerr, sn.z + sn.zerr + dz, dz )
    izbinarg = np.digitize( simIa.z, zbins )
    i7 = np.where( simIa.FLT=='7' )
    i8 = np.where( simIa.FLT=='8' )
    iI = np.where( simIa.FLT=='I' )
    iP = np.where( simIa.FLT=='P' )
    iN = np.where( simIa.FLT=='N' )
    m7,m8,mI,mP,mN = simIa.MAG[i7],simIa.MAG[i8],simIa.MAG[iI],simIa.MAG[iP],simIa.MAG[iN]
    c7Ipts = np.array([ (m7-mI)[izbinarg==iz] for iz in range(len(zbins))])
    c8Ipts = np.array([ (m8-mI)[izbinarg==iz] for iz in range(len(zbins))])
    cPNpts = np.array([ (mP-mN)[izbinarg==iz] for iz in range(len(zbins))])
    return( zbins, c7Ipts, c8Ipts, cPNpts )


def get_colors( simtable, dz=0.1 ):
    """ Bin up the simulated Ia's by redshift and determine the point in
    color-color space occupied by the median member of each bin.
    """
    if 'FLT' not in simtable.__dict__ :
        simtable.getLightCurves()
    i7 = np.where( simtable.FLT=='7' )
    i8 = np.where( simtable.FLT=='8' )
    iI = np.where( simtable.FLT=='I' )
    iP = np.where( simtable.FLT=='P' )
    iN = np.where( simtable.FLT=='N' )
    m7,m8,mI,mP,mN = simtable.MAG[i7],simtable.MAG[i8],simtable.MAG[iI],simtable.MAG[iP],simtable.MAG[iN]
    c7Ipts = m7-mI
    c8Ipts = m8-mI
    cPNpts = mP-mN
    return( c7Ipts, c8Ipts, cPNpts )






def plotzpoints(sn,simIa,dz=0.05):
    """  plot the
    :param sn:
    :param simIa:
    :param dz:
    :return:
    """
    from matplotlib import cm
    import mpltools
    from pytools import plotsetup

    if 'FLT' not in simIa.__dict__ :
        simIa.getLightCurves()

    zbins, c7Ipts, c8Ipts, cPNpts = get_colors_binned_by_z(sn,simIa,dz)
    plotsetup.fullpaperfig( 1, [5,5] )
    pl.clf()
    ax1 = pl.gca()
    # ax1.plot(QH,PN, ls=' ', marker='o',color='k', alpha=0.3 )
    mpltools.color.cycle_cmap( len(zbins), cmap=cm.gist_rainbow_r, ax=ax1)
    for iz in range(len(zbins)-1) :
        z0, z1 = zbins[iz],zbins[iz+1]
        #if dz>=0.1 : zmid = round(np.mean([z0,z1]),1)
        #elif dz>=0.01 : zmid = round(np.mean([z0,z1]),2)
        ax1.plot(c7Ipts[iz],cPNpts[iz], marker='o', alpha=0.3, ls=' ',
                 mew=0,label='%.2f-%.2f'%(z0,z1) )

    ax1.legend(loc='upper left', ncol=2, numpoints=1, frameon=False,
               handletextpad=0.3, handlelength=0.2 )
    ax1.set_xlabel('F763M-F814W')
    ax1.set_ylabel('F139M-F140W')
    # plotPhot(label=True)
    ax1.set_xlim(-0.8,0.8)
    ax1.set_ylim(-0.8,0.8)
    fig = pl.gcf()
    fig.subplots_adjust(bottom=0.12,left=0.16,right=0.95,top=0.95)
    pl.draw()


def plotPhot(label=False):
    """ Plot the SN Camille Med Band Photometry on a color-color diagram
    """
    m7 = 25.813 ; dm7 = 0.142
    m8 = 25.430 ; dm8 = 0.161
    mI = 25.359 ; dmI = 0.090
    mP = 25.593 ; dmP = 0.117
    mN = 25.142 ; dmN = 0.052

    # F763M - F814W = 7 - I
    c7I = m7-mI
    c7Ierr = np.sqrt( dm7**2 + dmI**2)

    # F845M - F814W = 8 - I
    c8I = m8-mI
    c8Ierr = np.sqrt( dm8**2 + dmI**2)

    # F139M - F140W = P - N
    cPN = mP-mN
    cPNerr = np.sqrt( dmP**2 + dmN**2)

    print("F763M-F814W = %.2f +- %.2f"%(c7I,c7Ierr))
    print("F845M-F814W = %.2f +- %.2f"%(c8I,c8Ierr))
    print("F139M-F140W = %.2f +- %.2f"%(cPN,cPNerr))
    fig = pl.gcf()
    ax1 = fig.add_subplot(2,2,1)
    ax1.errorbar( c7I, cPN, cPNerr, c7Ierr, color='k', marker='o', ms=15, capsize=0 )
    ax2 = fig.add_subplot(2,2,2)
    ax2.errorbar( c8I, cPN, cPNerr, c8Ierr, color='k', marker='o', ms=15, capsize=0 )
    ax3 = fig.add_subplot(2,2,3)
    ax3.errorbar( c7I, c8I, c8Ierr, c7Ierr, color='k', marker='o', ms=15, capsize=0 )

    if label :
        ax2.text( 0.95,0.95,'GND13Cam', transform=ax2.transAxes, ha='right',va='top',fontsize='large')


def plotbands(z = 1.2):
    import stardust
    from medbandfig import w763,f763,w845,f845,w139,f139
    w1a, f1a = stardust.snsed.getsed( sedfile='/usr/local/SNDATA_ROOT/snsed/Hsiao07.dat', day=0 )

    w1az = w1a * (1+z)
    f1az = f1a / f1a.max() / 2.
    ax18 = pl.gca() # subplot(3,2,1)
    ax18.plot(w1az, f1az, ls='-', lw=0.7, color='0.5', label='_nolegend_')
    ax18.plot(w763, f763, ls='-', color='DarkOrchid',label='F763M')
    ax18.plot(w845, f845, ls='-',color='Teal', label='F845M')
    ax18.plot(w139, f139, ls='-',color='Maroon', label='F139M')
    #ax18.fill_between( w1az, np.where(f763>0.01,f1az,0), color='DarkOrchid', alpha=0.3 )
    #ax18.fill_between( w1az, f1az, where=((w1az>13500) & (w1az<14150)), color='teal', alpha=0.3 )
    #ax18.fill_between( w1az, f1az, where=((w1az>15000) & (w1az<15700)), color='Maroon', alpha=0.3 )
    ax18.text(0.95,0.4, 'SNIa\n@ z=%.1f'%(z), color='k',ha='right',va='bottom',fontweight='bold', transform=ax18.transAxes, fontsize='large')
    ax18.set_xlim( 6500, 16000 )
    pl.setp(ax18.get_xticklabels(), visible=False)
    pl.setp(ax18.get_yticklabels(), visible=False)
    ax18.text( 7630, 0.65, 'F763M', ha='right', va='center', color='DarkOrchid', fontweight='bold')
    ax18.text( 8450, 0.65, 'F845M', ha='center', va='center', color='Teal', fontweight='bold')
    ax18.text( 13900, 0.65, 'F139M', ha='left', va='center', color='Maroon', fontweight='bold')
