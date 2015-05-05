"""
S.Rodney
2014.12.22
making a figure showing how the SN classifications
change as med band data are added.
"""


def readdata(datfile = 'lightcurve_subset_classtest.dat'):
    """  read in the summary data file
    :param datfile:
    :return:
    """
    from astropy.io import ascii
    import os
    import sys
    import exceptions

    if not os.path.isfile( datfile ):
        thisfile = sys.argv[0]
        if 'ipython' in thisfile :
            thisfile = __file__
        thisdir = os.path.dirname( thisfile )
        datfile = os.path.join( thisdir, datfile )
    if not os.path.isfile( datfile ):
        raise exceptions.RuntimeError("No such file %s"%datfile)


    dat = ascii.read( datfile, format='commented_header', header_start=-1, data_start=0)
    return( dat )


def plotdata( dat, hostprior=True) :
    """ make a single-panel plot showing how P(Ia) changes with additions
     of new data.
    :param dat:
    :return:
    """
    from matplotlib import pyplot as pl
    import numpy as np
    from pytools import colorpalette as cp

    for sn,marker,color in zip( ['stone','colfax','bush'],
                                ['d','s','o'],
                                [cp.teal, cp.darkred, cp.black ] ):
        if sn=='bush' : continue
        # with med band data
        imb = np.where( (dat['x']>0) & (dat['x']<4) & (dat['mb']>0) & (dat['sn']==sn) & (dat['hostprior']==hostprior) )[0]

        # no med band data
        inmb = np.where( (dat['x']>0) & (dat['x']<4) & (dat['mb']<1) & (dat['sn'] == sn) & (dat['hostprior']==hostprior))[0]

        pl.plot( dat['x'][imb], dat['PIa'][imb],   color=color, marker=marker, mec=color, mfc=color, ms=10, ls='-', alpha=0.5 )
        pl.plot( dat['x'][inmb], dat['PIa'][inmb], color=color, marker=marker, mec=color, mfc='w', ms=10, ls='--', alpha=0.5 )

    ax = pl.gca()
    ax.set_xlim( -1,5)
    ax.set_ylim( 0.0, 1.1 )
    ax.set_xticks( [0,1,2,3,4] )
    ax.set_xticklabels( [] )
    ax.set_ylabel('P(Ia$|${\\bf D})')


def mkfigure2():
    """ Two-panel figure
    :return:
    """
    from matplotlib import pyplot as pl
    import numpy as np
    from pytools import plotsetup, colorpalette as cp

    plotsetup.halfpaperfig( figsize=[4,7])
    fig = pl.gcf()
    fig.clf()
    fig.subplots_adjust( bottom=0.05, top=0.92, left=0.15, right=0.97, hspace=0 )

    dat = readdata()
    ax1 = fig.add_subplot( 2, 1, 1 )
    plotdata( dat, hostprior=False )

    #ax1.set_xlim(0.5,4.5)
    #ax1.set_ylim(0.02,1.18)
    ax1.set_xticks( [1,2,3,4] )
    ax1.set_xticklabels( [] )
    ax1.set_ylabel('P(Ia$|${\\bf D})')


    ax1.text( 0.4, 1.32, '{\\bf D} =', color='0.5', ha='center', va='top')#, size='small')
    # ax1.text( 0, 1.25, 'Med. Band\nPseudo-colors\nOnly', color='0.5', ha='center', va='top', size='small')
    ax1.text( 1, 1.35, 'F160W\nLight Curve', color='0.5', ha='center', va='top')#, size='small')
    ax1.text( 2, 1.35, 'F125W\nLight Curve', color='0.5', ha='center', va='top')#, size='small')
    ax1.text( 3, 1.35, 'All Wide\nBand IR', color='0.5', ha='center', va='top')#, size='small')
    # ax1.text( 4, 1.25, 'All Wide\nBand\nIR+ACS', color='0.5', ha='center', va='top', size='small')

    ax1.text( 2.5, 0.95, 'GND13Sto', color=cp.teal, rotation=18, ha='center', va='center')
    ax1.text( 2.6, 0.53, 'GND11Col', color=cp.darkred, rotation=20, ha='center', va='center')
    ax1.text( 0.85, 0.1, 'No host $z$ prior.', color='k', ha='right',va='bottom',transform=ax1.transAxes)

    ax1.text( 1.8, 1.0, 'With med. band\npseudo-colors', color='0.7', ha='right',va='bottom')
    ax1.plot( [1.5,2.0], [0.98,0.84], color='0.7', marker=' ', ls='-', lw=0.7 )
    ax1.text( 1.3, 0.85, 'Without', color='0.7', ha='right',va='bottom')
    ax1.plot( [1.15,2.0], [0.83,0.65], color='0.7', marker=' ', ls='-', lw=0.7 )


    ax2 = fig.add_subplot( 2, 1, 2, sharex=ax1, sharey=ax1 )
    plotdata( dat, hostprior=True )

    ax2.set_xlim(0.5,3.5)
    ax2.set_ylim(0.01,1.18)
    ax2.set_xticks( [1,2,3] )
    ax2.set_xticklabels( ['1','2','3'] )
    ax2.set_ylabel('P(Ia$|${\\bf D})')

    ax2.text( 0.85, 0.1, 'With host $z$ prior.', color='k', ha='right',va='bottom',transform=ax2.transAxes)

    pl.draw()


def mkfigure1():
    """  single-panel figure
    :return:
    """
    from matplotlib import pyplot as pl
    import numpy as np
    from pytools import plotsetup, colorpalette as cp

    plotsetup.halfpaperfig( figsize=[4,4])
    fig = pl.gcf()
    fig.clf()
    fig.subplots_adjust( bottom=0.09, top=0.85, left=0.15, right=0.97, hspace=0 )

    dat = readdata()
    ax1 = fig.add_subplot( 1, 1, 1 )
    plotdata( dat, hostprior=True )

    #ax1.set_xlim(0.5,4.5)
    #ax1.set_ylim(0.02,1.18)
    ax1.set_xticks( [1,2,3,4] )
    ax1.set_xticklabels( [] )
    ax1.set_ylabel('P(Ia$|${\\bf D})')

    ax1.text( 0.25, 1.32, '{\\bf D} =', color='0.5', ha='center', va='top')#, size='small')
    # ax1.text( 0, 1.25, 'Med. Band\nPseudo-colors\nOnly', color='0.5', ha='center', va='top', size='small')
    ax1.text( 1, 1.35, 'F160W\nLight Curve', color='0.5', ha='center', va='top')#, size='small')
    ax1.text( 2, 1.35, 'F125W\nLight Curve', color='0.5', ha='center', va='top')#, size='small')
    ax1.text( 3, 1.35, 'All Wide\nBand IR', color='0.5', ha='center', va='top')#, size='small')
    # ax1.text( 4, 1.25, 'All Wide\nBand\nIR+ACS', color='0.5', ha='center', va='top', size='small')

    ax1.text( 1.45, 0.47, 'GND13Sto', color=cp.teal, rotation=46, ha='center', va='center')
    ax1.text( 1.65, 0.28, 'GND11Col', color=cp.darkred, rotation=48, ha='center', va='center')

    ax1.text( 1.4, 1.0, 'With med. band\npseudo-colors', color='0.7', ha='right',va='bottom')
    ax1.plot( [1.5,2.0], [1.06,1.01], color='0.7', marker=' ', ls='-', lw=0.7 )
    ax1.text( 0.85, 0.65, 'Without', color='0.7', ha='right',va='center')
    ax1.plot( [0.88,1.97], [0.65,0.6], color='0.7', marker=' ', ls='-', lw=0.7 )

    ax1.set_xlim(0.01,3.5)
    ax1.set_ylim(0.0,1.18)
    ax1.set_xticks( [1,2,3] )
    ax1.set_xticklabels( ['1','2','3'] )
    ax1.set_ylabel('P(Ia$|${\\bf D})')

    # ax2.text( 0.85, 0.1, 'With host $z$ prior.', color='k', ha='right',va='bottom',transform=ax2.transAxes)

    pl.draw()
