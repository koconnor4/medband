__author__ = 'rodney'


#def plot_snana_light_curve_fit(
#    """ make a figure showing the SNANA light curve fit """





def lcplot( sn='colfax', datadir='/Users/rodney/MEDBAND/LIGHTCURVEFIG',
            yunits='flux', snfitres=None ):
    import time
    import os
    import sys
    from pytools import plotsetup, colorpalette as cp
    from astropy.io import ascii
    import numpy as np
    from matplotlib import pyplot as pl
    from matplotlib import ticker

    plotsetup.fullpaperfig( figsize=[8,3] )
    pl.clf()
    fig = pl.gcf()

    thisfile = sys.argv[0]
    if 'ipython' in thisfile :
        thisfile = __file__
    thisdir = os.path.abspath( os.path.dirname(thisfile))
    datadir = os.path.abspath( datadir )


    alpha2filter = { 'H':'f160w','N':'f140w','J':'f125w','Y':'f105w',
                     'L':'f098m','O':'f127m','P':'f139m','Q':'f153m',
                     'Z':'f850l','I':'f814w','V':'f606w'    }
    filter2alpha = {     'f160w':'H','f140w':'N','f125w':'J','f105w':'Y',
                         'f098m':'L','f127m':'O','f139m':'P','f153m':'Q',
                         'f850l':'Z','f814w':'I','f606w':'V' }

    colordict = { 'f160w':cp.black, 'f125w':cp.darkgold, 'f105w':cp.cadetblue,
                  'f105w':cp.black, 'f127m':cp.darkgold, 'f139m':cp.cadetblue, 'f153m':cp.beige,
                  'f814w':cp.darkgreen, 'f606w':cp.coral }

    fluxdat = ascii.read( os.path.join( datadir, 'HSThighz-%s.LCPLOT.TEXT'%sn) )

    mjd = fluxdat['col2']
    tobs = fluxdat['col3']
    fluxcal = fluxdat['col4']
    fluxcalerr = fluxdat['col5']
    dataflag  = fluxdat['col6']
    snanaband = fluxdat['col7']
    obsfilterlist = [ alpha2filter[b] for b in snanaband ]

    wideIRfilters = [ fn.lower() for fn in np.unique(obsfilterlist) if fn.lower() in ['f105w','f125w','f140w','f160w'] ]
    nax = len(wideIRfilters)
    medbanddict = {'f105w':'f098m','f125w':'f127m', 'f140w':'f139m', 'f160w':'f153m'}

    #fluxab      = fluxcal * 10**(-0.4*(27.5-zptab))
    #fluxerrab   = fluxcalerr * 10**(-0.4*(27.5-zptab))

    flux    = fluxcal * 0.1  # flux with zp_AB=25
    fluxerr = fluxcalerr * 0.1  # flux with zp_AB=25

    start = time.time()

    if sn=='bush' :
        snid = 'GSD11Bus'
        z = 1.76
        ymax=0.8
        labelcorner='upper right'
    elif sn=='colfax' :
        snid = 'GND12Col'
        z = 2.07
        ymax=0.7
        labelcorner='upper right'
    elif sn=='stone' :
        snid = 'GND13Sto'
        z = 1.8
        ymax=0.95
        labelcorner='upper right'
    else :
        ymax=1.2
        labelcorner='upper right'

    # markers = ['^','o','s','d']
    # colors = [cp.purple, cp.bluegray,cp.darkgreen,cp.red]

    mjdmin = mjd.min()-10
    mjdmax = mjd.max()+25

    mjdpeak = mjd[ flux.argmax() ]

    trestmin = (mjdmin - mjdpeak)/(1+z)
    trestmax = (mjdmax - mjdpeak)/(1+z)

    iax = 0

    color = cp.darkbluegray
    for widefiltname in  wideIRfilters :
        t1 = time.time()
        iax += 1

        medfiltname = medbanddict[widefiltname]

        widesnanaband = filter2alpha[ widefiltname ]
        medsnanaband = filter2alpha[ medfiltname ]

        iwide = np.where( snanaband == widefiltname )
        imed = np.where( snanaband == medfiltname )

        # filterletter = filter2alpha[widefiltname]
        # color = colordict[widefiltname]

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


        y,yerr = flux, fluxerr
        if iax==1 :
            ax.set_ylabel( 'flux (zp$_{AB}$=25)' )
        ax.set_ylim( -0.09, ymax )


        # Plot the wide band data and model
        iobs = np.where( (snanaband==widesnanaband) & (dataflag!=0) )[0]
        imod = np.where( (snanaband==widesnanaband) & (dataflag==0) )[0]
        ax.errorbar( mjd[iobs], y[iobs], yerr[iobs], color=color, marker='D', ls=' ', zorder=-100 )

        ax.plot( mjd[imod], y[imod], color=color, marker=' ', ls='-', zorder=-100 )
        ax.fill_between( mjd[imod], y[imod]+yerr[imod], y[imod]-yerr[imod],
                         color=color, alpha=0.3, zorder=-1000 )



        # ax.plot( mjdmod, ymod_bb, marker=None, color='0.5', ls='-'  , zorder=-10 )
        # ax.errorbar( mjd[iwide], yval[iwide], yerr[iwide], ms=8, marker='o', color=color, mfc=color, mec=color, capsize=0, ls=' ' , zorder=10 )

        # Plot the med band data and model :
        # import pdb; pdb.set_trace()
        if medfiltname in np.unique( obsfilterlist ):
            iobsmed = np.where( (snanaband==medsnanaband) & (dataflag!=0) )[0]
            imodmed = np.where( (snanaband==medsnanaband) & (dataflag==0) )[0]
            ax.errorbar( mjd[iobsmed], y[iobsmed], yerr[iobsmed], color='k', marker='o', mfc='w', ls=' ', alpha=0.5, zorder=200)

            ax.plot( mjd[imodmed], y[imodmed], color='0.5', marker=' ', ls='-', zorder=-100 )
            ax.fill_between( mjd[imodmed], y[imodmed]+yerr[imodmed], y[imodmed]-yerr[imodmed],
                             color='0.5', alpha=0.3, zorder=-1000 )


        if labelcorner=='lower right' :
            ax.text( 0.92,0.42, '%s'%(widefiltname.upper()) , color=color, fontsize='medium', fontweight='heavy', ha='right',va='top', transform=ax.transAxes )
            if medfiltname in obsfilterlist :
                ax.text( 0.92,0.32, '%s'%(medfiltname.upper()) , color='k', fontsize='medium', fontweight='heavy', ha='right',va='top', transform=ax.transAxes )
        elif labelcorner=='upper left' :
            ax.text( 0.08,0.92, '%s'%(widefiltname.upper()) , color=color, fontsize='medium', fontweight='heavy', ha='left',va='top', transform=ax.transAxes )
            if medfiltname in obsfilterlist :
                ax.text( 0.08,0.82, '%s'%(medfiltname.upper()) , color='k', fontsize='medium', fontweight='heavy', ha='left',va='top', transform=ax.transAxes )
        if labelcorner=='upper right' :
            ax.text( 0.92,0.92, '%s'%(widefiltname.upper()) , color=color, fontsize='medium', fontweight='heavy', ha='right',va='top', transform=ax.transAxes )
            if medfiltname in obsfilterlist :
                ax.text( 0.92,0.82, '%s'%(medfiltname.upper()) , color='k', fontsize='medium', fontweight='heavy', ha='right',va='top', transform=ax.transAxes )


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
        print( "%s  %.1f %.1f "%( widefiltname.upper(), t2-start, t2-t1))


    fig.subplots_adjust( left=0.1, bottom=0.18, right=0.95, top=0.82, wspace=0., hspace=0.15 )
    pl.draw()
    t3 = time.time()
    print( "All  %.1f"%( t3-start) )

