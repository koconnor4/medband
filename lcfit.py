__author__ = 'rodney'

from astropy import cosmology
from astropy.io import ascii
import numpy as np
from scipy.stats.stats import chisqprob

import sncosmo
# from sncosmost import hstbandpasses,ccsnmodels

cosmo=cosmology.FlatLambdaCDM( name="WMAP9", H0=70.0, Om0=0.3 )
dm = lambda z : cosmo.distmod( z ).value
MBmodel = -19.223 # AB mags ( = -19.12 Vega )

def fitbush( ):
    sn, fit, res = dofit( datfile='DATFILES/HST_CANDELS4_bush.sncosmo.dat',  z=1.76, dz=0.53, t0=55803.1, dt0=25.0 )
    snid = 'GSD11Bus'


def dofitIa( snid='colfax' ):
    """  fit the given SN with a SALT2 Ia model
    :param snid:
    :return:
    """
    if snid=='bush' :
        sn, fit, res = dofit( datfile='DATFILES/HST_CANDELS4_bush.sncosmo.dat',  z=1.76, dz=0.53, t0=55803.1, dt0=25.0 )
    elif snid=='colfax' :
        sn, fit, res = dofit( datfile='DATFILES/HST_CANDELS4_colfax.sncosmo.dat', z=2.13, dz=0.3, t0=56075, dt0=20.0 )
    elif snid=='stone' :
        sn, fit, res = dofit( datfile='DATFILES/HST_CANDELS4_stone.sncosmo.dat',  z=1.9, dz=0.3, t0=56467.5, dt0=15.0 )
    else :
        return(-1)
    return( sn, fit, res)

def lcplot( snid='colfax', yunits='flux', snfitres=None ):
    from matplotlib import pyplot as pl, ticker
    from pytools import plotsetup, colorpalette as cp
    import time

    start = time.time()

    if snfitres is not None :
        sn,fit,res = snfitres
    else :
        sn,fit,res = dofitIa( snid )

    if snid=='bush' :
        snid = 'GSD11Bus'
        ymax=0.8
        labelcorner='upper right'
    elif snid=='colfax' :
        snid = 'GND12Col'
        ymax=0.7
        labelcorner='upper right'
    elif snid=='stone' :
        snid = 'GND13Sto'
        ymax=1.2
        labelcorner='upper right'
    else :
        return(-1)

    # optbands = ['f350lp','f606w','f814w','f850lp']
    allbands = np.unique( sn['filter'] )
    broadbands = [ band for band in allbands if band in ['f105w','f125w','f140w','f160w'] ]

    nax = len(broadbands)
    medbanddict = {'f105w':'f098m','f125w':'f127m', 'f140w':'f139m', 'f160w':'f153m'}
    bands = sn['filter']
    markers = ['^','o','s','d']
    colors = [cp.purple, cp.bluegray,cp.darkgreen,cp.red]


    mjd = sn['mjd']
    mag = sn['mag']
    magerr = sn['magerr']
    flux = sn['flux'] * 10**(-0.4*(sn['zpt']-25))
    fluxerr = sn['fluxerr'] * 10**(-0.4*(sn['zpt']-25))

    plotsetup.fullpaperfig( 2, [8,3] )
    pl.clf()
    fig = pl.gcf()

    mjdmin = mjd.min()-10
    mjdmax = mjd.max()+25
    mjdmod = np.arange( mjdmin, mjdmax, 1 )

    z = res['parameters'][ res['param_names'].index('z') ]
    mjdpeak = res['parameters'][ res['param_names'].index('t0') ]

    trestmin = (mjdmin - mjdpeak)/(1+z)
    trestmax = (mjdmax - mjdpeak)/(1+z)

    iax = 0
    for bb, marker, color in zip( broadbands, markers, colors ) :
        t1 = time.time()
        iax += 1
        ibb = np.where( bands == bb )
        mb = medbanddict[bb]
        imb = np.where( bands == mb )

        if iax==1 :
            ax = fig.add_subplot( 1, nax, iax )
            ax1 = ax
            axtop = ax.twiny()
            ax1.text( -0.3, 1.22, snid, color='k', fontsize='large', fontweight='heavy', ha='left',va='top', transform=ax.transAxes )

        elif iax==2 :
            ax = fig.add_subplot( 1, nax, iax, sharex=ax1, sharey=ax1 )
            axtop = ax.twiny()
            axtop.set_xlabel( 'rest frame time (days from peak)')
            ax.set_xlabel( 'observer frame time (MJD)')
            pl.setp( ax.get_yticklabels(), visible=False )

        elif iax==nax :
            ax = fig.add_subplot( 1, nax, iax, sharex=ax1, sharey=ax1 )
            axtop = ax.twiny()
            ax.yaxis.set_ticks_position('right')
            ax.yaxis.set_ticks_position('both')

        else :
            ax = fig.add_subplot( 1, nax, iax, sharex=ax1, sharey=ax1 )
            axtop = ax.twiny()
            pl.setp( ax.get_yticklabels(), visible=False )


        if yunits=='mag' :
            yval = mag
            yerr = magerr
            ymod_bb = fit.bandmag( bb, 'ab', mjdmod )
            ymod_mb = fit.bandmag( mb, 'ab', mjdmod )
            if iax==1 :
                ax.set_ylabel( 'AB mag' )
            # ax.set_ylim( 30.6, 25.2 )
        else :
            yval = flux
            yerr = fluxerr
            ymod_bb = fit.bandflux( bb, mjdmod , zp=25, zpsys='ab' )
            ymod_mb = fit.bandflux( mb, mjdmod , zp=25, zpsys='ab' )
            if iax==1 :
                ax.set_ylabel( 'flux (zp$_{AB}$=25)' )
            ax.set_ylim( -0.09, ymax )

        ax.plot( mjdmod, ymod_bb, marker=None, color='0.5', ls='-'  , zorder=-10 )
        ax.errorbar( mjd[ibb], yval[ibb], yerr[ibb], ms=8, marker='o', color=color, mfc=color, mec=color, capsize=0, ls=' ' , zorder=10 )

        if mb in allbands :
            ax.plot( mjdmod, ymod_mb, marker=None, color='0.5', ls='--' , zorder=-10 )
            ax.errorbar( mjd[imb], yval[imb], yerr[imb], ms=9, marker='o', color='k', mfc='w',   mec='k', capsize=0, ls=' ' , alpha=0.5, zorder=20 )

        if labelcorner=='lower right' :
            ax.text( 0.92,0.42, '%s'%(bb.upper()) , color=color, fontsize='medium', fontweight='heavy', ha='right',va='top', transform=ax.transAxes )
            if mb in allbands :
                ax.text( 0.92,0.32, '%s'%(mb.upper()) , color='k', fontsize='medium', fontweight='heavy', ha='right',va='top', transform=ax.transAxes )
        elif labelcorner=='upper left' :
            ax.text( 0.08,0.92, '%s'%(bb.upper()) , color=color, fontsize='medium', fontweight='heavy', ha='left',va='top', transform=ax.transAxes )
            if mb in allbands :
                ax.text( 0.08,0.82, '%s'%(mb.upper()) , color='k', fontsize='medium', fontweight='heavy', ha='left',va='top', transform=ax.transAxes )
        if labelcorner=='upper right' :
            ax.text( 0.92,0.92, '%s'%(bb.upper()) , color=color, fontsize='medium', fontweight='heavy', ha='right',va='top', transform=ax.transAxes )
            if mb in allbands :
                ax.text( 0.92,0.82, '%s'%(mb.upper()) , color='k', fontsize='medium', fontweight='heavy', ha='right',va='top', transform=ax.transAxes )


        ax.set_xlim( mjdmin, mjdmax )
        axtop.set_xlim(  trestmin, trestmax )
        ax.xaxis.set_major_locator( ticker.MultipleLocator( 100 ) )
        ax.xaxis.set_minor_locator( ticker.MultipleLocator( 20 ) )
        axtop.xaxis.set_major_locator( ticker.MultipleLocator( 20 ) )
        axtop.xaxis.set_minor_locator( ticker.MultipleLocator( 5 ) )
        ax.yaxis.set_major_locator( ticker.MultipleLocator( 0.2 ) )
        ax.yaxis.set_minor_locator( ticker.MultipleLocator( 0.1 ) )

        for tick in axtop.get_xaxis().get_major_ticks():
            tick.set_pad(5.)
            tick.label1 = tick._get_text1()

        t2 = time.time()
        print( "%s  %.1f %.1f "%( bb.upper(), t2-start, t2-t1))


    fig.subplots_adjust( left=0.1, bottom=0.18, right=0.95, top=0.82, wspace=0., hspace=0.15 )
    pl.draw()
    t3 = time.time()
    print( "All  %.1f"%( t3-start) )


def dofit( datfile='HST_CANDELS4_bush.sncosmo.dat', z=1.76, dz=0.53,
           t0=55803.1, dt0=25.0, model='Ia', debug=False) :

    # TODO : read in the redshift, etc from the header.

    # read in the obs data
    sn = ascii.read( datfile, format='commented_header', header_start=-1, data_start=0 )


    if model == 'Ia' :
        # define SALT2 models and set initial guesses for z and t0
        # salt2ex = sncosmo.Model( source='salt2-extended')
        salt2ex = sncosmo.Model( source='salt2')
        salt2ex.source.set_peakmag( 0., 'bessellb', 'ab' )
        x0_AB0 = salt2ex.get('x0')
        salt2ex.set( z=z, t0=t0, x1=0.1, c=-0.2 )
        # salt2ex.set( z=1.33, t0=56814.6, hostebv=0.05, hostr_v=3.1 )

        # Do a bounded fit :
        # salt2res, salt2fit = sncosmo.fit_lc( sn, salt2, ['z','t0','x0','x1','c'], bounds={'z':(1.28,1.37),'t0':(56804,56824)} )
        res, fit = sncosmo.fit_lc( sn, salt2ex, ['z','t0','x0','x1','c'],
                                                 bounds={'z':(z-dz,z+dz),'t0':(t0-dt0,t0+dt0),
                                                         'x1':(-5.,5.), 'c':(-0.5,3.0) })
        x0 = fit.get( 'x0' )
        z = fit.get( 'z' )
        mB = -2.5*np.log10(  x0 / x0_AB0 )
        distmod = mB - MBmodel
        deltamuLCDM = distmod - dm(z)
        print( "mB = %.2f"%mB )
        print( "dist.mod. = %.2f"%distmod)
        print( "Delta.mu_LCDM = %.2f"%deltamuLCDM)

        chi2 = res.chisq
        ndof = res.ndof
        pval = chisqprob( chi2, ndof )

        print( "chi2/dof= %.3f"% (chi2/float(ndof) ) )
        print( "p-value = %.3f"% pval )
        return( sn, fit, res )


    elif model.lower() in ['cc','ib','ic','ii','ibc','iip','iin'] :
        # remove the blue filters from the sn data
        bandlist = sn['filter'].data
        igood = np.array( [ band.startswith('f1') for band in bandlist ] )
        sn = sn.copy()[igood]

        # define a host-galaxy dust model
        dust = sncosmo.CCM89Dust( )
        version = '1.0'

        if model.lower()=='cc' : classlist = ['Ib','Ic','IIP','IIn']
        elif model.lower()=='ii' : classlist = ['IIP','IIn']
        elif model.lower()=='ibc' : classlist = ['Ib','Ic']
        else : classlist = [model]

        # find the best-fit from each CC sub-class
        chi2list, reslist, fitlist  = [],[],[]
        for snclass in classlist :
            for tempnum in range( 1, 10 ):
                Av = 0.2
                modname = snclass.lower() + '.%02i'%tempnum
                modkey = ( sncosmo.Source, modname, version )
                if modkey not in sncosmo.registry._loaders : continue
                ccmodel = sncosmo.Model( source=modname, effects=[dust],
                                         effect_names=['host'], effect_frames=['rest'])
                ccmodel.set( z=z, t0=t0, hostr_v=3.1, hostebv=Av/3.1 )
                # Do a bounded fit :
                res, fit  = sncosmo.fit_lc(
                    sn, ccmodel, ['z','t0','amplitude','hostebv' ], debug=debug,
                    bounds={'z':(z-dz,z+dz),'t0':(t0-dt0,t0+dt0),
                            'hostebv':(0.0,1.0) } )

                chi2 = res.chisq
                ndof = res.ndof
                pval = chisqprob( chi2, ndof )

                print( "%s  chi2/dof= %.3f  p=%.3f"%(modname, chi2/float(ndof), pval  ) )
                chi2list.append( chi2/float(ndof) )
                reslist.append( res )
                fitlist.append( fit )
        ichi2min = np.argmin( chi2list )
        res, fit = reslist[ichi2min], fitlist[ichi2min]
        return( sn, fit, res )




def sncosmoplot( sn, fit, res ):
    from pytools import plotsetup
    plotsetup.presfig()
    sncosmo.plot_lc( sn, model=fit, errors=res.errors )