__author__ = 'rodney'

import sncosmo
from matplotlib import pyplot as pl
import numpy as np
from scipy import interpolate as scint

def mkdemofig():
    """ Construct a demonstration figure showing how the medium band filters
    pick out features in the Type Ia SN spectrum and how pseudo-colors can
    provide redshift information.  (This is Figure 1 in the paper)
    """
    from matplotlib import ticker
    from pytools import plotsetup
    plotsetup.fullpaperfig(figsize=[8,4])
    pl.clf()

    ax1 = pl.subplot( 3, 2, 1 )
    mksedplot( 1.8, color='darkblue' )
    ax2 = pl.subplot( 3, 2, 3, sharex=ax1, sharey=ax1 )
    mksedplot( 2.0, color='green')
    ax3 = pl.subplot( 3, 2, 5, sharex=ax1, sharey=ax1 )
    mksedplot( 2.2, color='darkred' )

    pl.setp( ax1.get_xticklabels(), visible=False )
    pl.setp( ax2.get_xticklabels(), visible=False )

    ax3.set_xlim( 1.12, 1.69 )
    ax3.set_ylim( 0.0, 0.7 )
    ax3.set_yticklabels([])
    ax3.set_xlabel('wavelength ($\mu$m)')
    ax2.set_ylabel('SN Flux or Filter Transmission\n (arbitrary units)')

    ax1.text(1.27,0.6,'F127M',color='darkmagenta',fontsize=9, ha='center',va='center')
    ax1.text(1.39,0.6,'F139M',color='teal',fontsize=9, ha='center',va='center')
    ax1.text(1.53,0.6,'F153M',color='darkorange',fontsize=9, ha='center',va='center')

    ax1.text(0.95,0.45,'z=1.8', color='darkblue', transform=ax1.transAxes,ha='right',va='top')
    ax2.text(0.95,0.45,'z=2.0', color='green', transform=ax2.transAxes,ha='right',va='top')
    ax3.text(0.95,0.55,'z=2.2', color='darkred', transform=ax3.transAxes,ha='right',va='top')

    ax4 = pl.subplot( 2, 2, 2 )
    mkcirclepoints(zrange=[1.8,2.2], colorselect=[1,0], source='hsiao',
                   marker='o', ls=' ', mew=0.1, mec='w' )
    ax4.text( 0.19, -0.13, 'z=1.8', color='darkblue', ha='right', va='center', rotation=0)
    ax4.text( -0.06, 0.02, 'z=2.0', color='green', ha='center', va='center', rotation=0)
    ax4.text( 0.15, -0.34, 'z=2.2', color='darkred', ha='left', va='center', rotation=0)

    pl.setp( ax4.get_xticklabels(), visible=False )

    ax4.yaxis.set_ticks_position( 'right')
    ax4.yaxis.set_ticks_position( 'both')
    ax4.yaxis.set_label_position( 'right')

    ax5 = pl.subplot( 2, 2, 4, sharex=ax4 )
    mkcirclepoints(zrange=[1.8,2.2], colorselect=[1,2], source='hsiao',
                   marker='o', ls=' ', mew=0.1, mec='w' )

    ax5.yaxis.set_ticks_position( 'right')
    ax5.yaxis.set_ticks_position( 'both')
    ax5.yaxis.set_label_position( 'right')

    ax5.text( 0.2, 0.02, 'z=1.8', color='darkblue', ha='center', va='center', rotation=0)
    ax5.text( -0.07,-0.09, 'z=2.0', color='green', ha='center', va='center', rotation=0)
    ax5.text( 0.11, 0.13, 'z=2.2', color='darkred', ha='center', va='center', rotation=0)

    ax5.set_xlim( -0.2, 0.25 )
    ax4.set_ylim( -0.39, 0.19 )
    ax5.set_ylim( -0.15, 0.29 )
    pl.setp( ax5.yaxis.get_label(), rotation=-90 )
    pl.setp( ax4.yaxis.get_label(), rotation=-90 )
    ax4.yaxis.labelpad= 20
    ax5.yaxis.labelpad= 20

    ax1.xaxis.set_major_locator( ticker.MultipleLocator( 0.1 ) )
    ax1.xaxis.set_minor_locator( ticker.MultipleLocator( 0.05 ) )
    ax5.xaxis.set_major_locator( ticker.MultipleLocator( 0.1 ) )
    ax5.xaxis.set_minor_locator( ticker.MultipleLocator( 0.05 ) )
    ax4.yaxis.set_major_locator( ticker.MultipleLocator( 0.1 ) )
    ax4.yaxis.set_minor_locator( ticker.MultipleLocator( 0.05 ) )
    ax5.yaxis.set_major_locator( ticker.MultipleLocator( 0.1 ) )
    ax5.yaxis.set_minor_locator( ticker.MultipleLocator( 0.05 ) )

    fig = pl.gcf()
    fig.subplots_adjust( left=0.08, right=0.88, bottom=0.15, top=0.97, hspace=0, wspace=0.08 )

def mksedplot( z=1.8, color='k' ):
    """ make a set of plots showing the SN Ia SED (Hsiao template)
    at various redshifts, with bandpasses overlaid.
    """
    snIa = sncosmo.Model( source='hsiao' )
    snIa.set( z=z, t0=0 )
    snwave = np.arange( 6000., 20000., 10. )
    snflux = snIa.flux( 0, snwave )
    snwave = snwave / 10000.
    snflux = 0.5 * snflux / snflux.max()

    pl.plot( snwave, snflux, color=color, ls='-')

    f105w = sncosmo.get_bandpass( 'f105w' )
    f098m = sncosmo.get_bandpass( 'f098m' )
    f127m = sncosmo.get_bandpass( 'f127m' )
    f139m = sncosmo.get_bandpass( 'f139m' )
    f153m = sncosmo.get_bandpass( 'f153m' )

    wf127m = f127m.wave / 10000.
    wf139m = f139m.wave / 10000.
    wf153m = f153m.wave / 10000.

    pl.plot( wf127m, f127m.trans, color='darkmagenta', ls='-')
    pl.plot( wf139m, f139m.trans, color='teal', ls='-')
    pl.plot( wf153m, f153m.trans, color='darkorange', ls='-')

    intf127m = scint.interp1d( wf127m, f127m.trans, bounds_error=False, fill_value=0 )
    overlap = np.min( [snflux, intf127m(snwave)], axis=0 )
    pl.fill_between( snwave, np.zeros(len(snwave)), overlap, color='darkmagenta' )

    intf139m = scint.interp1d( wf139m, f139m.trans, bounds_error=False, fill_value=0 )
    overlap = np.min( [snflux, intf139m(snwave)], axis=0 )
    pl.fill_between( snwave, np.zeros(len(snwave)), overlap, color='teal' )

    intf153m = scint.interp1d( wf153m, f153m.trans, bounds_error=False, fill_value=0 )
    overlap = np.min( [snflux, intf153m(snwave)], axis=0 )
    pl.fill_between( snwave, np.zeros(len(snwave)), overlap, color='darkorange' )




def mkcirclepoints(zrange=[1.8,2.2], t0=0, colorselect=[0,1],
                   source='hsiao', **plotargs ):
    """ Make a string of colored points in pseudo-color-color space,
    showing where a fiducial Type Ia SN lands as a function of redshift.
    """
    import mpltools
    from matplotlib import cm
    from copy import copy
    labelcolors = ['darkmagenta','teal','darkorange']
    plotargs10 = copy(plotargs)
    plotargs10.update( mec='w', marker='D', zorder=1000, ms=10 )

    zsteps = np.arange( zrange[0], zrange[1], 0.01 )
    mpltools.color.cycle_cmap( length=len(zsteps), cmap=cm.jet )
    for z in np.arange( zrange[0], zrange[1], 0.01 ) :
        sn = sncosmo.Model( source=source )
        sn.set( z=z, t0=t0 )
        f127m = sn.bandmag('f127m', 'ab', 0.)
        f139m = sn.bandmag('f139m', 'ab', 0.)
        f153m = sn.bandmag('f153m', 'ab', 0.)
        f125w = sn.bandmag('f125w', 'ab', 0.)
        f140w = sn.bandmag('f140w', 'ab', 0.)
        f160w = sn.bandmag('f160w', 'ab', 0.)
        med = np.array([ f127m,f139m,f153m] )
        broad = np.array([ f125w,f140w,f160w] )
        pseudocolor = med-broad
        colorx = pseudocolor[colorselect[0]]
        colory = pseudocolor[colorselect[1]]

        if (int(np.round(z*100))%10) == 0 :
            pl.plot( colorx, colory, **plotargs10 )
        else :
            pl.plot( colorx, colory, **plotargs )

    colorlabel = ['F127M-F125W','F139M-F140W','F153M-F160W']
    ax = pl.gca()
    ax.set_xlabel(colorlabel[colorselect[0]], color=labelcolors[colorselect[0]])
    ax.set_ylabel(colorlabel[colorselect[1]], color=labelcolors[colorselect[1]])


