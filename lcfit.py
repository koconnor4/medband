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
        ymax=1.2
        labelcorner='upper right'

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
           t0=55803.1, dt0=25.0, x1=None, c=None,
           model='Ia', noUV=True, debug=False) :

    # TODO : read in the redshift, etc from the header.

    # read in the obs data
    sn = ascii.read( datfile, format='commented_header', header_start=-1, data_start=0 )


    if model == 'Ia' :
        # define SALT2 models and set initial guesses for z and t0
        if noUV :
            salt2ex = sncosmo.Model( source='salt2')
        else :
            salt2ex = sncosmo.Model( source='salt2-extended')
        salt2ex.source.set_peakmag( 0., 'bessellb', 'ab' )
        x0_AB0 = salt2ex.get('x0')
        salt2ex.set( z=z, t0=t0, x1=0.1, c=-0.2 )
        # salt2ex.set( z=1.33, t0=56814.6, hostebv=0.05, hostr_v=3.1 )

        # Do a bounded fit :
        # salt2res, salt2fit = sncosmo.fit_lc( sn, salt2, ['z','t0','x0','x1','c'], bounds={'z':(1.28,1.37),'t0':(56804,56824)} )
        varlist = varlist = ['z','t0','x0']
        bounds={  'z':(z-dz,z+dz), 't0':(t0-dt0,t0+dt0) }
        if x1 is not None:
            salt2ex.set( x1=x1 )
            bounds['x1'] = (x1-1e-6,x1+1e-6)
            varlist.append( 'x1' )
        else :
            bounds['x1'] = (-5,5)
            varlist.append( 'x1' )
        if c is not None:
            salt2ex.set( c=c )
        else :
            bounds['c'] = (-0.5,3.0)
            varlist.append( 'c' )

        res, fit = sncosmo.fit_lc( sn, salt2ex, varlist, bounds )

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

        if ndof>0:
            print( "chi2/dof= %.3f"% (chi2/float(ndof) ) )
            print( "p-value = %.3f"% pval )
        else :
            print( "chi2/dof= %.3f/%i"%( chi2, ndof) )
            print( "p-value = %.3f"% pval )

        print( "z = %.3f"% fit.get('z') )
        print( "t0 = %.3f"% fit.get('t0') )
        print( "x0 = %.3e"% fit.get('x0') )
        print( "x1 = %.3f"% fit.get('x1') )
        print( "c = %.3f"% fit.get('c') )

    elif model.lower() in ['cc','ib','ic','ii','ibc','iip','iin']:
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
    else : # 'nugent-sn91bg'
        # remove the blue filters from the sn data
        bandlist = sn['filter'].data
        igood = np.array( [ band.startswith('f1') for band in bandlist ] )
        sn = sn.copy()[igood]

        # define a host-galaxy dust model
        dust = sncosmo.CCM89Dust( )
        version = '1.0'

        Av = 0.2
        altmodel = sncosmo.Model( source=model, effects=[dust],
                                 effect_names=['host'], effect_frames=['rest'])
        altmodel.set( z=z, t0=t0, hostr_v=3.1, hostebv=Av/3.1 )
        # Do a bounded fit :
        res, fit  = sncosmo.fit_lc(
            sn, altmodel, ['z','t0','amplitude','hostebv' ], debug=debug,
            bounds={'z':(z-dz,z+dz),'t0':(t0-dt0,t0+dt0),
                    'hostebv':(0.0,1.0) } )

        chi2 = res.chisq
        ndof = res.ndof
        pval = chisqprob( chi2, ndof )

        print( "%s  chi2/dof= %.3f  p=%.3f"%(model, chi2/float(ndof), pval  ) )

    return( sn, fit, res )




def sncosmoplot( sn, fit, res ):
    from pytools import plotsetup
    plotsetup.presfig()
    sncosmo.plot_lc( sn, model=fit, errors=res.errors )

def dofitSCP0401( datfile='HST_SCP_0401.sncosmo.dat', z=1.713,
                  t0=53080.0, dt0=50.0 ) :

    # TODO : read in the redshift, etc from the header.

    # read in the obs data
    sn = ascii.read( datfile, format='commented_header', header_start=-1, data_start=0 )

    # define SALT2 models and set initial guesses for z and t0
    salt2ex = sncosmo.Model( source='salt2-extended')
    salt2ex.source.set_peakmag( 0., 'bessellb', 'ab' )
    x0_AB0 = salt2ex.get('x0')
    x0_from_mB = lambda m : x0_AB0 * 10**(-0.4*(m) )
    salt2ex.set( z=1.713, t0=53090.0, x0=x0_from_mB(26.14), x1=0.2, c=-0.1 )
    # salt2ex.set( z=1.33, t0=56814.6, hostebv=0.05, hostr_v=3.1 )

    # Do a bounded fit :
    #res, fit = sncosmo.fit_lc( sn, salt2ex, ['z','t0','x0','x1','c'],
    #                           bounds={'z':(1.712,1.714),'t0':(t0-dt0,t0+dt0),
    #                                   'x1':(-5.,5.), 'c':(-0.5,3.0) })
    res, fit = sncosmo.fit_lc( sn, salt2ex, ['z','t0','x0'],
                               bounds={'z':(1.712,1.714),'t0':(t0-dt0,t0+dt0)})

    x0 = fit.get( 'x0' )
    mB = -2.5*np.log10(  x0 / x0_AB0 )
    distmod = mB - -19.19 #  MBmodel from Rubin et al 2013
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


def plot_light_curve_fit_from_SNANA( snid='colfax', plotmags=False ) :
    """ make a plot showing the light curve fit, extracted from a SNANA
    output .TEXT file generated by D.Scolnic"""
    from pytools import plotsetup, colorpalette as cp
    from astropy.io import ascii
    import numpy as np
    from matplotlib import pyplot as pl
    from scipy.interpolate import interp1d
    from matplotlib import ticker
    import os
    import sys
    from hstphot import hstzpt

    alpha2filter = { 'H':'f160w','N':'f140w','J':'f125w','Y':'f105w',
                     'Z':'f850l','I':'f814w','V':'f606w',
                     'L':'f098m','O':'f127m','P':'f139m','Q':'f153m',
                     }
    colordict = { 'H':'k', 'J':cp.darkgold, 'Y':cp.cadetblue,
                  'N':cp.darkgreen, 'V':cp.coral }

    fluxdat = ascii.read( 'CANDELs-%s.LCPLOT.TEXT'%snid )
    name = fluxdat['col1']
    mjd = fluxdat['col2']
    tobs = fluxdat['col3']
    fluxcal = fluxdat['col4']
    fluxcalerr = fluxdat['col5']
    obsflag = fluxdat['col6']
    bandletter = fluxdat['col7']

    flux    = fluxcal * 0.1
    fluxerr = fluxcalerr * 0.1
    # mag = np.where( flux>0, -2.5*np.log10( flux ) + 25, 35 )
    # magerr = np.where( flux>0, 1.0857 * fluxerr / flux, 0.2 )

    # medbanddict = {'f105w':'f098m','f125w':'f127m', 'f140w':'f139m', 'f160w':'f153m'}
    colors = [cp.purple, cp.bluegray,cp.darkgreen,cp.red]


    plotsetup.fullpaperfig( [8,3] )
    pl.clf()
    fig = pl.gcf()

    mjdmin = mjd.min()-10
    mjdmax = mjd.max()+25
    mjdmod = np.arange( mjdmin, mjdmax, 1 )

    if snid=='colfax' :
        ymax=0.7
        ymaxMB=0.3
        mjdmin= 56025
        mjdmax = 56220
        z, mjdpeak = 2.1, 56080.
        snlabel='GND12Col'
    elif snid=='stone' :
        ymax=1.2
        ymaxMB=0.35
        mjdmin= 56430
        mjdmax = 56600
        z,mjdpeak = 1.8, 56485.
        snlabel='GND13Sto'

    trestmin = (mjdmin - mjdpeak)/(1+z)
    trestmax = (mjdmax - mjdpeak)/(1+z)

    ax1w = pl.subplot2grid( [4,3], [0,0], rowspan=3 )
    ax2w = pl.subplot2grid( [4,3], [0,1], rowspan=3, sharex=ax1w, sharey=ax1w )
    ax3w = pl.subplot2grid( [4,3], [0,2], rowspan=3, sharex=ax1w, sharey=ax1w  )
    ax1m = pl.subplot2grid( [4,3], [3,0], sharex=ax1w )
    ax2m = pl.subplot2grid( [4,3], [3,1], sharex=ax1w, sharey=ax1m )
    ax3m = pl.subplot2grid( [4,3], [3,2], sharex=ax1w, sharey=ax1m )

    iax = 0
    for bbl, mbl, color in zip( ['Y','J','N','H'], ['L','O','P','Q'], colors ) :
        if bbl not in bandletter : continue
        iax += 1
        filternamew = alpha2filter[ bbl ].upper()
        filternamem = alpha2filter[ mbl ].upper()

        if iax==1 :
            axw = ax1w
            axm = ax1m
            axw.text( 0.05, 0.92, snlabel, color='k', fontsize='large', fontweight='heavy', ha='left',va='top', transform=axw.transAxes )
            # axw.text( -0.28, 1.2, snlabel, backgroundcolor='w', color='k', fontsize='large', fontweight='heavy', ha='left',va='top', transform=axw.transAxes, zorder=1000 )
            axw.set_ylabel( 'flux (zp$_{AB}$=25)' )
            axm.set_ylabel( '$\Delta$f$_{25}$' )

        elif iax==2 :
            axw = ax2w
            axm = ax2m
        else :
            axw = ax3w
            axm = ax3m
            axw.yaxis.set_ticks_position('right')
            axw.yaxis.set_ticks_position('both')
            axm.yaxis.set_ticks_position('right')
            axm.yaxis.set_ticks_position('both')

        axtop = axw.twiny()
        axtop.set_xlim(  trestmin, trestmax )
        if iax==2 :
            axtop.set_xlabel( 'rest frame time (days from peak)')
            axm.set_xlabel( 'observer frame time (MJD)')
            pl.setp( axw.get_yticklabels(), visible=False )
            pl.setp( axm.get_yticklabels(), visible=False )

        pl.setp( axw.get_xticklabels(), visible=False )

        axw.text( 0.95,0.92, filternamew, ha='right', va='top', color=color, transform=axw.transAxes )

        iobs = np.where( (bandletter==bbl) & (obsflag>0) )[0]
        imod = np.where( (bandletter==bbl) & (obsflag==0) )[0]

        # import pdb; pdb.set_trace()
        axw.errorbar( mjd[iobs], flux[iobs], fluxerr[iobs], color=color, marker='D', ls=' ', zorder=-100 )
        axw.plot( mjd[imod], flux[imod], color=color, marker=' ', ls='-', zorder=-100 )
        axw.fill_between( mjd[imod], flux[imod]+fluxerr[imod], flux[imod]-fluxerr[imod],
                          color=color, alpha=0.3, zorder=-1000 )

        if mbl not in bandletter : continue
        axm.text( 0.95,0.5, filternamem, ha='right', va='center', backgroundcolor='w', color=color, transform=axm.transAxes )

        iobsmb = np.where( (bandletter==mbl) & (obsflag>0) )[0]
        imodmb = np.where( (bandletter==mbl) & (obsflag==0) )[0]

        mbfluxinterp = interp1d( mjd[imodmb], flux[imodmb], fill_value=0, bounds_error=False )
        mbfluxdiff = flux[iobsmb] - mbfluxinterp( mjd[iobsmb] )

        axm.errorbar( mjd[iobsmb], mbfluxdiff, fluxerr[iobsmb], color=color, marker='o', mfc='w', ls=' ' )
        axm.fill_between( mjd[imodmb], fluxerr[imodmb], -fluxerr[imodmb],
                          color=color, alpha=0.3, zorder=-1000 )
        axm.axhline( 0, color=color, lw=0.8, ls='--')

    ax1w.set_xlim( mjdmin, mjdmax )
    ax1w.set_ylim( -0.08, ymax )
    ax1m.set_ylim( -ymaxMB, ymaxMB )

    ax1w.xaxis.set_major_locator( ticker.MultipleLocator( 100 ) )
    ax1w.xaxis.set_minor_locator( ticker.MultipleLocator( 25 ) )
    ax1w.yaxis.set_major_locator( ticker.MultipleLocator( 0.2 ) )
    ax1w.yaxis.set_minor_locator( ticker.MultipleLocator( 0.1 ) )
    ax1m.xaxis.set_major_locator( ticker.MultipleLocator( 100 ) )
    ax1m.xaxis.set_minor_locator( ticker.MultipleLocator( 25 ) )
    ax1m.yaxis.set_major_locator( ticker.MultipleLocator( 0.2 ) )
    ax1m.yaxis.set_minor_locator( ticker.MultipleLocator( 0.1 ) )


    fig.subplots_adjust( left=0.1, bottom=0.18, right=0.95, top=0.82, wspace=0., hspace=0 )

    pl.draw()


def mk_bush_lc_plot( snfitres=None ):
    from astropy.io import ascii
    import numpy as np
    from scipy.interpolate import interp1d
    from matplotlib import pyplot as pl, ticker
    from pytools import plotsetup, colorpalette as cp
    import time

    start = time.time()

    if snfitres is not None :
        sn,fit,res = snfitres
    else :
        sn,fit,res = dofit('/Users/rodney/Dropbox/MEDBAND/DATFILES/HST_CANDELS4_bushALL.sncosmo.dat',
                           z=1.15, dz=0.001, t0=55803.1, dt0=25.0,
                           model='s11-2006fo', noUV=True )
    snid = 'GSD11Bus'
    ymax=0.8

    # optbands = ['f350lp','f606w','f814w','f850lp']
    allbands = np.unique( sn['filter'] )
    broadbands = [ band for band in allbands if band in ['f814w','f105w','f125w','f140w','f160w'] ]

    nax = 5
    medbanddict = {'f105w':'f098m','f125w':'f127m', 'f140w':'f139m', 'f160w':'f153m'}
    bands = sn['filter']

    mjd = sn['mjd']
    flux = sn['flux'] * 10**(-0.4*(sn['zpt']-25))
    fluxerr = sn['fluxerr'] * 10**(-0.4*(sn['zpt']-25))

    alpha2filter = { 'H':'f160w','N':'f140w','J':'f125w','Y':'f105w',
                     'Z':'f850l','I':'f814w','V':'f606w',
                     'L':'f098m','O':'f127m','P':'f139m','Q':'f153m',
                     }

    plotsetup.fullpaperfig( [8,3] )
    pl.clf()
    fig = pl.gcf()

    #mjdmin = mjd.min()-10
    #mjdmax = mjd.max()+25

    ymax=0.7
    ymaxMB=0.39
    mjdmin= 55701
    mjdmax = 55995
    mjdmod = np.arange( mjdmin, mjdmax, 1 )
    z, mjdpeak = 1.15, 55802.
    snlabel='GSD11Bus'

    trestmin = (mjdmin - mjdpeak)/(1+z)
    trestmax = (mjdmax - mjdpeak)/(1+z)

    ax1w = pl.subplot2grid( [4,4], [0,0], rowspan=3 )
    ax2w = pl.subplot2grid( [4,4], [0,1], rowspan=3, sharex=ax1w, sharey=ax1w )
    ax3w = pl.subplot2grid( [4,4], [0,2], rowspan=3, sharex=ax1w, sharey=ax1w  )
    ax4w = pl.subplot2grid( [4,4], [0,3], rowspan=3, sharex=ax1w, sharey=ax1w  )
    #ax5w = pl.subplot2grid( [4,5], [0,4], rowspan=3, sharex=ax1w, sharey=ax1w  )
    ax1m = pl.subplot2grid( [4,4], [3,0], sharex=ax1w )
    ax2m = pl.subplot2grid( [4,4], [3,1], sharex=ax1w, sharey=ax1m )
    ax3m = pl.subplot2grid( [4,4], [3,2], sharex=ax1w, sharey=ax1m )
    ax4m = pl.subplot2grid( [4,4], [3,3], sharex=ax1w, sharey=ax1m )
    #ax5m = pl.subplot2grid( [4,5], [3,4], sharex=ax1w, sharey=ax1m )

    iax = 0
    for bb, mb in zip( ['f105w','f125w','f140w','f160w'], ['f098m','f127m','f139m','f153m'] ) :
        iax += 1
        ibb = np.where( bands == bb )
        imb = np.where( bands == mb )

        yval = flux
        yerr = fluxerr
        ymod_bb = fit.bandflux( bb, mjdmod , zp=25, zpsys='ab' )

        if iax==1 :
            axw = ax1w
            axm = ax1m
            axw.text( 0.05, 0.92, snlabel, color='k', fontsize='large', fontweight='heavy', ha='left',va='top', transform=axw.transAxes )
            # axw.text( -0.28, 1.2, snlabel, backgroundcolor='w', color='k', fontsize='large', fontweight='heavy', ha='left',va='top', transform=axw.transAxes, zorder=1000 )
            axw.set_ylabel( 'flux (zp$_{AB}$=25)' )
            axm.set_ylabel( '$\Delta$f$_{25}$' )

        elif iax==2 :
            axw = ax2w
            axm = ax2m
        elif iax==3 :
            axw = ax3w
            axm = ax3m
        elif iax==4 :
            axw = ax4w
            axm = ax4m
            axw.yaxis.set_ticks_position('right')
            axw.yaxis.set_ticks_position('both')
            axm.yaxis.set_ticks_position('right')
            axm.yaxis.set_ticks_position('both')

        axtop = axw.twiny()
        axtop.set_xlim(  trestmin, trestmax )
        if iax in [2,3] :
            pl.setp( axw.get_yticklabels(), visible=False )
            pl.setp( axm.get_yticklabels(), visible=False )
        if iax == 2 :
            axtop.set_xlabel( 'rest frame time (days from peak)')
            axm.set_xlabel( 'observer frame time (MJD)')

        pl.setp( axw.get_xticklabels(), visible=False )

        axw.text( 0.95,0.92, bb.upper(), ha='right', va='top', color='k', transform=axw.transAxes )

        color='k'
        axw.plot( mjdmod, ymod_bb, marker=None, color='0.5', ls='-'  , zorder=-10 )
        axw.errorbar( mjd[ibb], yval[ibb], yerr[ibb], ms=8, marker='o', color=color, mfc=color, mec=color, capsize=0, ls=' ' , zorder=10 )

        if not mb : continue
        axm.text( 0.95,0.5, mb.upper(), ha='right', va='center', backgroundcolor='w', color=color, transform=axm.transAxes )

        mbfluxdiff = yval[imb] - fit.bandflux( mb, mjd[imb], zp=25, zpsys='ab' )

        axm.errorbar( mjd[imb], mbfluxdiff, fluxerr[imb], color='k', marker='o', mfc='w', ls=' ' )
        # axm.fill_between( mjdmod, fluxerr[imodmb], -fluxerr[imodmb],
        #                  color=color, alpha=0.3, zorder=-1000 )
        axm.axhline( 0, color='k', lw=0.8, ls='--')

    ax1w.set_xlim( mjdmin, mjdmax )
    ax1w.set_ylim( -0.08, ymax )
    ax1m.set_ylim( -ymaxMB, ymaxMB )

    ax1w.xaxis.set_major_locator( ticker.MultipleLocator( 100 ) )
    ax1w.xaxis.set_minor_locator( ticker.MultipleLocator( 25 ) )
    ax1w.yaxis.set_major_locator( ticker.MultipleLocator( 0.2 ) )
    ax1w.yaxis.set_minor_locator( ticker.MultipleLocator( 0.1 ) )
    ax1m.xaxis.set_major_locator( ticker.MultipleLocator( 100 ) )
    ax1m.xaxis.set_minor_locator( ticker.MultipleLocator( 25 ) )
    ax1m.yaxis.set_major_locator( ticker.MultipleLocator( 0.2 ) )
    ax1m.yaxis.set_minor_locator( ticker.MultipleLocator( 0.1 ) )


    fig.subplots_adjust( left=0.1, bottom=0.18, right=0.95, top=0.82, wspace=0., hspace=0 )

    pl.draw()
