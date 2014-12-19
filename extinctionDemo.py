__author__ = 'rodney'


def mkExtinctionDemoFigSmall( z=2.0 ):
    """ make a set of plots showing the SN Ia SED (Hsiao template)
    at three extinction values, with bandpasses overlaid.
    :return:
    """
    import numpy as np
    import sncosmo
    # from sncosmost import hstbandpasses, ccsnmodels
    from matplotlib import rc
    rc('text',usetex=True)
    rc('text.latex', preamble='\usepackage[usenames]{xcolor}')
    from matplotlib import pyplot as pl
    from matplotlib import ticker
    from pytools import plotsetup
    from scipy import interpolate as scint


    fig = plotsetup.fullpaperfig( 1, [8,3] )


    # load the O'Donnell 1994 dust model
    dust = sncosmo.OD94Dust()
    snIa = sncosmo.Model( source='hsiao', effects=[dust],
                          effect_names=['host'], effect_frames=['rest'])

    ax1 = pl.gca()


    f127m = sncosmo.get_bandpass( 'f127m' )
    f139m = sncosmo.get_bandpass( 'f139m' )
    f153m = sncosmo.get_bandpass( 'f153m' )

    f125w = sncosmo.get_bandpass( 'f125w' )
    f140w = sncosmo.get_bandpass( 'f140w' )
    f160w = sncosmo.get_bandpass( 'f160w' )

    wf127m = f127m.wave / 10000.
    wf139m = f139m.wave / 10000.
    wf153m = f153m.wave / 10000.

    wf125w = f125w.wave / 10000.
    wf140w = f140w.wave / 10000.
    wf160w = f160w.wave / 10000.

    # ax2 = ax1.twinx()
    ax2 = ax1
    ax2.plot( wf127m, f127m.trans, color='darkmagenta', ls='-', lw=2)
    ax2.plot( wf153m, f153m.trans, color='darkorange', ls='-', lw=2)

    ax2.plot( wf125w, f125w.trans, color='darkmagenta', ls='--', lw=2)
    ax2.plot( wf160w, f160w.trans, color='darkorange', ls='--', lw=2)

    intf127m = scint.interp1d( wf127m, f127m.trans, bounds_error=False, fill_value=0 )
    intf153m = scint.interp1d( wf153m, f153m.trans, bounds_error=False, fill_value=0 )

    colorlist1, colorlist2 = [], []
    for Av,ls, alpha in zip([2,1,0],[':','--','-'],[0.1,0.3,0.5]):
        snIa.set( z=z, t0=0, hostr_v=3.1, hostebv=Av/3.1 )
        colorlist1.append( snIa.bandmag( 'f127m','ab', 0) - snIa.bandmag( 'f125w','ab', 0) )
        colorlist2.append( snIa.bandmag( 'f153m','ab', 0) - snIa.bandmag( 'f160w','ab', 0) )

        snwave = np.arange( 6000., 20000., 10. )
        snflux = snIa.flux( 0, snwave )
        snwave = snwave / 10000.
        snflux = 0.12 * snflux / snflux[400]
        ax1.plot( snwave, snflux, color='k', ls=ls, lw=1, label='%.1f'%Av )
        overlap127 = np.min( [snflux, intf127m(snwave)], axis=0 )
        ax2.fill_between( snwave, np.zeros(len(snwave)), overlap127, color='darkmagenta', alpha=alpha )
        overlap153 = np.min( [snflux, intf153m(snwave)], axis=0 )
        pl.fill_between( snwave, np.zeros(len(snwave)), overlap153, color='darkorange', alpha=alpha )


    ax1.legend(loc='upper left', bbox_to_anchor=(0.0,0.9),frameon=False,fontsize=11 )
    ax1.text( 0.08, 0.88, 'A$_V$', transform=ax1.transAxes, ha='left',va='bottom',fontsize=11 )
    ax1.text( 0.13, 0.88, '$\Delta$m$_{127}$', color='darkmagenta',transform=ax1.transAxes, ha='left',va='bottom',fontsize=11 )
    ax1.text( 0.23, 0.88, '$\Delta$m$_{153}$', color='darkorange',transform=ax1.transAxes, ha='left',va='bottom',fontsize=11 )

    ax1.text( 0.14, 0.78, '%.3f'%colorlist1[0], color='darkmagenta',transform=ax1.transAxes, ha='left',va='bottom',fontsize=11 )
    ax1.text( 0.23, 0.78, '%.3f'%colorlist2[0], color='darkorange',transform=ax1.transAxes, ha='left',va='bottom',fontsize=11 )

    ax1.text( 0.14, 0.68, '%.3f'%colorlist1[1], color='darkmagenta',transform=ax1.transAxes, ha='left',va='bottom',fontsize=11 )
    ax1.text( 0.23, 0.68, '%.3f'%colorlist2[1], color='darkorange',transform=ax1.transAxes, ha='left',va='bottom',fontsize=11 )

    ax1.text( 0.14, 0.58, '%.3f'%colorlist1[2], color='darkmagenta',transform=ax1.transAxes, ha='left',va='bottom',fontsize=11 )
    ax1.text( 0.23, 0.58, '%.3f'%colorlist2[2], color='darkorange',transform=ax1.transAxes, ha='left',va='bottom',fontsize=11 )

    #           title=+
    #                 '\\textcolor{DarlMagenta}{W}' +
    #                 '\\textcolor{F153M-F160W')#, handlelength=0.5, numpoints=3)
    # ax1.text( 0.15,0.95,,ha='left',va='bottom')

    ax1.yaxis.set_major_locator( ticker.MultipleLocator( 0.1 ) )
    ax1.yaxis.set_minor_locator( ticker.MultipleLocator( 0.05 ) )
    ax1.xaxis.set_major_locator( ticker.MultipleLocator( 0.2 ) )
    ax1.xaxis.set_minor_locator( ticker.MultipleLocator( 0.1 ) )
    ax1.set_xlabel('wavelength ($\mu$m)')
    ax1.set_ylabel('SN Flux or Filter Transmission\n (arbitrary units)')

    ax1.set_xlim(0.6,2.0)
    ax1.set_ylim( 0.0, 0.7 )
    ax1.set_yticklabels([])

    ax1.text(1.27,0.6,'F127M,F125W',color='darkmagenta',fontsize=9, ha='center',va='center')
    ax1.text(1.53,0.6,'F153M,F160W',color='darkorange',fontsize=9, ha='center',va='center')

    fig.subplots_adjust( left=0.12, right=0.95, bottom=0.18, top=0.92 )


def mkExtinctionDemoFig():
    """  Make a demo figure that shows how the pseudo-colors
    do not change much even under extreme host galaxy extinction
    :return:
    """
    from matplotlib import pyplot as pl
    from matplotlib import ticker
    from pytools import plotsetup
    plotsetup.fullpaperfig(figsize=[8,4])
    pl.clf()

    ax1 = pl.subplot( 3, 2, 1 )
    mksedplot( z=2.0, Av=0.0, color='darkblue' )
    ax2 = pl.subplot( 3, 2, 3, sharex=ax1, sharey=ax1 )
    mksedplot( z=2.0, Av=1.0, color='green')
    ax3 = pl.subplot( 3, 2, 5, sharex=ax1, sharey=ax1 )
    mksedplot( z=2.0, Av=2.0, color='darkred' )

    pl.setp( ax1.get_xticklabels(), visible=False )
    pl.setp( ax2.get_xticklabels(), visible=False )

    ax3.set_xlim( 1.12, 1.69 )
    ax3.set_ylim( 0.0, 0.7 )
    ax3.set_yticklabels([])
    #ax3.set_xticks([1.1,1.3,1.5,1.7])
    ax3.set_xlabel('wavelength ($\mu$m)')
    ax2.set_ylabel('SN Flux or Filter Transmission\n (arbitrary units)')

    ax1.text(1.27,0.6,'F127M',color='darkmagenta',fontsize=9, ha='center',va='center')
    ax1.text(1.39,0.6,'F139M',color='teal',fontsize=9, ha='center',va='center')
    ax1.text(1.53,0.6,'F153M',color='darkorange',fontsize=9, ha='center',va='center')

    ax1.text(0.95,0.45,'A$_V$=0.0', color='darkblue', transform=ax1.transAxes,ha='right',va='top')
    ax2.text(0.95,0.65,'A$_V$=1.0', color='green', transform=ax2.transAxes,ha='right',va='top')
    ax3.text(0.95,0.75,'A$_V$=2.0', color='darkred', transform=ax3.transAxes,ha='right',va='top')

    ax4 = pl.subplot( 2, 2, 2 )
    mkcirclepoints( colorselect=[1,0], source='hsiao', marker='o', ls=' ', mew=0 )
    pl.setp( ax4.get_xticklabels(), visible=False )

    ax4.yaxis.set_ticks_position( 'right')
    ax4.yaxis.set_ticks_position( 'both')
    ax4.yaxis.set_label_position( 'right')

    ax5 = pl.subplot( 2, 2, 4, sharex=ax4 )
    mkcirclepoints( colorselect=[1,2], source='hsiao', marker='o', ls=' ', mew=0 )

    ax5.yaxis.set_ticks_position( 'right')
    ax5.yaxis.set_ticks_position( 'both')
    ax5.yaxis.set_label_position( 'right')

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

def mksedplot( z=1.8, Av=0, color='k' ):
    """ make a set of plots showing the SN Ia SED (Hsiao template)
    at various redshifts, with bandpasses overlaid.
    :return:
    """
    import numpy as np
    import sncosmo
    from sncosmohst import hstbandpasses, ccsnmodels
    from matplotlib import pyplot as pl
    from matplotlib import ticker
    from pytools import plotsetup
    from scipy import interpolate as scint

    # load the O'Donnell 1994 dust model
    dust = sncosmo.OD94Dust()
    snIa = sncosmo.Model( source='hsiao', effects=[dust],
                          effect_names=['host'], effect_frames=['rest'])
    snIa.set( z=z, t0=0, hostr_v=3.1, hostebv=Av/3.1 )
    snwave = np.arange( 6000., 20000., 10. )
    snflux = snIa.flux( 0, snwave )
    snwave = snwave / 10000.
    snflux = 0.5 * snflux / snflux.max()

    pl.plot( snwave, snflux, color=color, ls='-')

    f105w = sncosmo.get_bandpass( 'wfc3f105w' )
    f098m = sncosmo.get_bandpass( 'wfc3f098m' )
    f127m = sncosmo.get_bandpass( 'wfc3f127m' )
    f139m = sncosmo.get_bandpass( 'wfc3f139m' )
    f153m = sncosmo.get_bandpass( 'wfc3f153m' )

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


def mkcirclepoints(z = 2.0, avrange=[0.0,2.0], colorselect=[0,1],
                   source='hsiao', **plotargs ):
    import numpy as np
    from matplotlib import pyplot as pl
    from mpltools import color
    from matplotlib import cm
    import sncosmo
    from sncosmohst import hstbandpasses, ccsnmodels

    # load the O'Donnell 1994 dust model
    dust = sncosmo.OD94Dust()
    sn = sncosmo.Model( source='hsiao', effects=[dust],
                        effect_names=['host'], effect_frames=['rest'])

    avsteps = np.arange( avrange[0], avrange[1], 0.05 )
    color.cycle_cmap( length=len(avsteps), cmap=cm.jet )

    f127m, f139m, f153m = [], [], []
    f125w, f140w, f160w = [], [], []
    for av in avsteps :
        sn.set( z=z, t0=0, hostr_v=3.1, hostebv=av/3.1 )
        f127m.append( sn.bandmag('wfc3f127m', 'ab', 0.) )
        f139m.append( sn.bandmag('wfc3f139m', 'ab', 0.) )
        f153m.append( sn.bandmag('wfc3f153m', 'ab', 0.) )
        f125w.append( sn.bandmag('wfc3f125w', 'ab', 0.) )
        f140w.append( sn.bandmag('wfc3f140w', 'ab', 0.) )
        f160w.append( sn.bandmag('wfc3f160w', 'ab', 0.) )
    med = np.array([ f127m,f139m,f153m] )
    broad = np.array([ f125w,f140w,f160w] )
    color = med-broad
    pl.plot( color[colorselect[0]], color[colorselect[1]], **plotargs )

    colorlabel = ['F127M-F125W','F139M-F140W','F153M-F160W']
    ax = pl.gca()
    ax.set_xlabel(colorlabel[colorselect[0]])
    ax.set_ylabel(colorlabel[colorselect[1]])

    labelcolors = ['darkmagenta','teal','darkorange']

    colorlabel = ['F127M-F125W','F139M-F140W','F153M-F160W']
    ax = pl.gca()
    ax.set_xlabel(colorlabel[colorselect[0]], color=labelcolors[colorselect[0]])
    ax.set_ylabel(colorlabel[colorselect[1]], color=labelcolors[colorselect[1]])