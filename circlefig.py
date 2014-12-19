from matplotlib import pyplot as pl
import numpy as np
import sncosmo
# from sncosmost import hstbandpasses, ccsnmodels
from medband_classtest import SncosmoSim, scumsum

__author__ = 'rodney'


_COLOR1 = '#882255' # maroon
_COLOR2 = '#44AA99' # teal
_COLOR3 = '#DDCC77' # beige

class supernova( object ):
    def __init__( self, name ):
        self.name = name
        self.t0err = 3

    @property
    def t0_range( self ):
        return( [ self.mjdmedband-self.mjdpk-self.t0err,self.mjdmedband-self.mjdpk+self.t0err ] )

    @property
    def t0( self ):
        return( self.mjdmedband-self.mjdpk )

DEMO = supernova('demo')
DEMO.mjdpk = 0.0
DEMO.mjdmedband = 0.0
DEMO.z_range = [1.8, 2.2]

JWST2 = supernova('jwst2')
JWST2.mjdpk = 0.0
JWST2.mjdmedband = 0.0
#JWST.z_range = [2.8, 3.2]
JWST2.z_range = [1.6, 2.4]

JWST3 = supernova('jwst3')
JWST3.mjdpk = 0.0
JWST3.mjdmedband = 0.0
#JWST.z_range = [2.8, 3.2]
JWST3.z_range = [2.6, 3.4]


STONE = supernova('stone')
STONE.mjdpk = 56482.
STONE.mjdmedband = 56475.
STONE.z_range = [1.6, 2.2]
STONE.mag = { # from idl5 psf fitting
            'f139m':25.375,
            'f140w':25.344, # 25.349, # <- from fit.
            #'f140w':25.344,
            'f153m':25.464,
            'f160w':25.233, # 25.354, # <- from fit.
            #'f160w':25.233,
            }
STONE.magerr = { # from drop method
            'f139m':0.248,
            'f140w':0.155, # 0.01, # <- from fit.
            'f153m':0.188,
            'f160w':0.08, # 0.01, # <- from fit.
            }

COLFAX = supernova('colfax')
COLFAX.mjdpk = 56078.
COLFAX.mjdmedband = 56084.
COLFAX.z_range = [1.9, 2.4]
COLFAX.mag = {  # final photometry 2014.10.21
                'f127m':25.677,
                'f139m':25.628,
                'f153m':25.791,
                'f125w':26.106,
                'f140w':25.734,
                'f160w':25.931,
                }
COLFAX.magerr = { # using the drop-and-recover psf fitting err method
                  'f127m':0.124,
                  'f139m':0.157,
                  'f153m':0.190,
                  'f125w':0.163,
                  'f140w':0.100,
                  'f160w':0.160,
                  }

BUSH = supernova('bush')
BUSH.mjdpk = 55797.
BUSH.mjdmedband = 55803.5
BUSH.z_range = [ 0.77, 1.34 ]  # host photoz 68% conf. region
# BUSH.z_range = [ 0.26, 2.21 ]
BUSH.mag = { # from psf phot  (entered 2010.10.20)
             'f098m':26.474,
             'f127m':25.961,
             'f139m':25.569,
             'f153m':25.463,
             'f105w':26.068,
             'f125w':25.953,
             'f160w':25.891,
             'f140w':25.937,
             }
BUSH.magerr = { # from drop-method psf errors  (entered 2010.10.20)
                'f098m':0.175,
                'f127m':0.241,
                'f139m':0.215,
                'f153m':0.211,
                'f105w':0.132,
                'f125w':0.131,
                'f160w':0.148,
                'f140w':0.144,
                }

def medband_matching_filter( medband ):
    if medband=='f098m': return('f105w')
    if medband=='f127m': return('f125w')
    if medband=='f139m': return('f140w')
    if medband=='f153m': return('f160w')
    if medband=='f140m': return('f150w')
    if medband=='f162m': return('f150w')
    if medband=='f182m': return('f200w')
    if medband=='f210m': return('f200w')
    if medband=='f250m': return('f277w')
    if medband=='f300m': return('f277w')
    if medband=='f335m': return('f356w')
    if medband=='f360m': return('f356w')
    if medband=='f410m': return('f444w')
    if medband=='f430m': return('f444w')
    if medband=='f460m': return('f444w')
    if medband=='f480m': return('f444w')



def sncosmo_sim( snroot='demo', filterset='hst',
                 z_range=[1.8,2.2], t0_range=[-3,3],
                 nsim=1000,
                 verbose=True, clobber=False ):
    """  Construct a color-color circle figure for SN Colfax, with observed
     photometry included.

    :param simIa:
    :param simCC:
    :param simIapkl:
    :param simIIpkl:
    :param z_range:
    :param nsim:
    :param verbose:
    :param clobber:
    :return:
    """
    import medband_classtest
    import os
    import cPickle

    simIapkl='%s_SncosmoSim_Ia.pkl'%snroot
    simIIpkl='%s_SncosmoSim_II.pkl'%snroot
    simIbcpkl='%s_SncosmoSim_Ibc.pkl'%snroot

    if os.path.isfile( simIapkl ) and not clobber>1 :
        if verbose: print("Loading Ia simulation from pickle : %s"%simIapkl)
        fin = open( simIapkl, 'rb' )
        simIa = cPickle.load( fin )
        fin.close()
    else :
        if verbose: print("Running a new Ia simulation, then saving to pickle : %s"%simIapkl)
        simIa = medband_classtest.SncosmoSim( 'Ia' , z_range=z_range, t0_range=t0_range, nsim=nsim, filterset=filterset )
        fout = open( simIapkl, 'wb' )
        cPickle.dump( simIa, fout, protocol=-1 )
        fout.close()

    if os.path.isfile( simIIpkl ) and not clobber>1 :
        if verbose: print("Loading II simulation from pickle : %s"%simIIpkl)
        fin = open( simIIpkl, 'rb' )
        simII = cPickle.load(fin)
        fin.close()
    else :
        if verbose: print("Running a new II simulation, then saving to pickle : %s"%simIIpkl)
        simII = medband_classtest.SncosmoSim( 'II' , z_range=z_range, t0_range=t0_range, nsim=nsim, filterset=filterset )
        fout = open( simIIpkl, 'wb' )
        cPickle.dump( simII, fout, protocol=-1 )
        fout.close()

    if os.path.isfile( simIbcpkl ) and not clobber>1 :
        if verbose: print("Loading Ibc simulation from pickle : %s"%simIbcpkl)
        fin = open( simIbcpkl, 'rb' )
        simIbc = cPickle.load(fin)
        fin.close()
    else :
        if verbose: print("Running a new Ibc simulation, then saving to pickle : %s"%simIbcpkl)
        simIbc = medband_classtest.SncosmoSim( 'Ibc' , z_range=z_range, t0_range=t0_range, nsim=nsim, filterset=filterset )
        fout = open( simIbcpkl, 'wb' )
        cPickle.dump( simIbc, fout, protocol=-1 )
        fout.close()


    return( simIa, simII, simIbc )



def _plot_colorcolor_singlesim( snsim, medbandx, medbandy,
                            plotstyle='points', nbins=None, **plotargs ):
    """ plot the med band pseudo-colors in a color-color diagram
    :param snsim:
    :return:
    """
    import numpy as np
    from matplotlib import pyplot as pl
    from matplotlib import ticker

    widebandx = medband_matching_filter(medbandx)
    widebandy = medband_matching_filter(medbandy)

    mag = np.array([ -2.5*np.log10( snlc['flux'] ) + snlc['zp'] for snlc in snsim.lightcurves ])

    flt = np.array( [snlc['band'] for snlc in snsim.lightcurves] )
    imx = np.where( flt==medbandx)
    imy = np.where( flt==medbandy )
    iwx = np.where( flt==widebandx)
    iwy = np.where( flt==widebandy)

    ax = pl.gca()
    if plotstyle=='points':
        plotargfinal = {'marker':'o', 'alpha':0.3, 'color':'darkorange', 'ls':' '}
        plotargfinal.update( **plotargs )
        ax.plot( mag[imx]-mag[iwx], mag[imy]-mag[iwy], **plotargfinal )
    elif plotstyle.startswith('contour') or plotstyle=='gradient':
        xarray = mag[imy]-mag[iwy]
        nsim = len(xarray)
        if nbins is None : nbins = int( np.sqrt( nsim  ) )
        plotargfinal = {'levels':[0.0,0.68,0.95],'colors':['r','g','b'],'ls':'-',
                        'alpha':0.5, 'extend':'neither'}
        plotargfinal.update( **plotargs )

        # Plot filled contours, showing  the full extent of the population,
        # and contour lines containing 68% of the population.
        # First, bin the points into a 2-d histogram:
        # (Note that we reverse the x-y order here to get the binned arrays
        #  plotted in the correct direction )
        count,y,x = np.histogram2d( mag[imy]-mag[iwy],mag[imx]-mag[iwx],
                                    bins=nbins, range=[[-1,1],[-1,1]] )

        # Renormalize relative to the sum of all SNe in this class :
        count /= count.sum()

        # Now set up an array 'cabove' such that  the cell value in cabove[i,j]
        # is equal to the sum of all cells that have a value higher than c[i,j]
        cabove = scumsum( count )

        # solid lines give probability contours at specified levels
        # (defaults to 0.68 for "1-sigma contours")
        ax.contour( x[:-1], y[:-1], cabove, **plotargfinal )
        if plotstyle=='contourf' :
            ax.contourf( x[:-1], y[:-1], cabove, **plotargfinal )

    ax.set_xlabel( '%s - %s'%( medbandx.upper(), widebandx.upper()))
    ax.set_ylabel( '%s - %s'%( medbandy.upper(), widebandy.upper()))

    ax.xaxis.set_major_locator( ticker.MultipleLocator( 0.1 ) )
    ax.xaxis.set_minor_locator( ticker.MultipleLocator( 0.05 ) )
    ax.yaxis.set_major_locator( ticker.MultipleLocator( 0.1 ) )
    ax.yaxis.set_minor_locator( ticker.MultipleLocator( 0.05 ) )
    return( ax )


def plotcontours( sim1, sim2, sim3=None, medbandx='f127m', medbandy='f139m', nbins=None, **plotargs ):
    """ Make a circle diagram, i.e. a med-wide band pseudo-color-color plot,
    showing both Type Ia and CC simulations over the given redshift range.
    :param snsim:
    :return:
    """
    plotargs1 = { 'levels':[0.,0.68,0.95], 'colors':[_COLOR1], 'alpha':0.3 }
    plotargs1.update( **plotargs )

    plotargs2 = { 'levels':[0.,0.68,0.95], 'colors':[_COLOR2], 'alpha':0.3 }
    plotargs2.update( **plotargs )

    plotargs3 = { 'levels':[0.,0.68,0.95], 'colors':[_COLOR3], 'alpha':0.3 }
    plotargs3.update( **plotargs )

    ax = _plot_colorcolor_singlesim(sim1, medbandx, medbandy, nbins=nbins, plotstyle='contourf', **plotargs1 )
    ax = _plot_colorcolor_singlesim(sim2, medbandx, medbandy, nbins=nbins, plotstyle='contourf', **plotargs2 )
    if sim3 is not None :
        ax = _plot_colorcolor_singlesim(sim3, medbandx, medbandy, nbins=nbins, plotstyle='contourf', **plotargs3 )
    return( ax )

def plotgradient( sim1, sim2, sim3=None, medbandx='f127m', medbandy='f139m', nbins=None, **plotargs ):
    """ Make a circle diagram, i.e. a med-wide band pseudo-color-color plot,
    showing both Type Ia and CC simulations over the given redshift range.
    :param snsim:
    :return:
    """
    plotargs1 = { 'levels':[0.,0.68,0.95], 'colors':[_COLOR1], 'alpha':0.3 }
    plotargs1.update( **plotargs )

    plotargs2 = { 'levels':[0.,0.68,0.95], 'colors':[_COLOR2], 'alpha':0.3 }
    plotargs2.update( **plotargs )

    plotargs3 = { 'levels':[0.,0.68,0.95], 'colors':[_COLOR3], 'alpha':0.3 }
    plotargs3.update( **plotargs )

    ax = _plot_colorcolor_singlesim(sim1, medbandx, medbandy, nbins=nbins, plotstyle='gradient', **plotargs1 )
    ax = _plot_colorcolor_singlesim(sim2, medbandx, medbandy, nbins=nbins, plotstyle='gradient', **plotargs2 )
    if sim3 is not None :
        ax = _plot_colorcolor_singlesim(sim3, medbandx, medbandy, nbins=nbins, plotstyle='gradient', **plotargs3 )
    return( ax )



def plotpoints( sim1, sim2, sim3=None, medbandx='f127m', medbandy='f139m', **plotargs ):
    """ Make a circle diagram, i.e. a med-wide band pseudo-color-color plot,
    showing both Type Ia and CC simulations over the given redshift range.
    :param snsim:
    :return:
    """
    ax = _plot_colorcolor_singlesim(sim2, medbandx, medbandy, plotstyle='points', marker='o', ls=' ', color=_COLOR2, **plotargs )
    if sim3 is not None :
        ax = _plot_colorcolor_singlesim(sim3, medbandx, medbandy, plotstyle='points', marker='o', ls=' ', color=_COLOR3, **plotargs )
    ax = _plot_colorcolor_singlesim(sim1, medbandx, medbandy, plotstyle='points', marker='o', ls=' ', color=_COLOR1, **plotargs )

    return( ax )



def plot_redshift_circle(z_range=[1.8,2.2], t0=0,
                         medbandx='f139m', medbandy='f127m',
                         source='salt2', coloredaxislabels=True, **plotargs ):
    """  plot a color-color circle showing how the position of a SNIa in
    color-color space changes with redshift, using markers color-coded by z.
    :param z_range:
    :param t0:
    :param medbandx:
    :param medbandy:
    :param source:
    :param coloredaxislabels:
    :param plotargs:
    :return:
    """
    from mpltools import color
    from matplotlib import cm

    zsteps = np.arange( z_range[0], z_range[1], 0.01 )
    color.cycle_cmap( length=len(zsteps), cmap=cm.jet )
    for z in zsteps :
        sn = sncosmo.Model( source=source )
        sn.set( z=z, t0=t0 )
        widebandx = medband_matching_filter(medbandx)
        widebandy = medband_matching_filter(medbandy)

        colorx = sn.bandmag(medbandx, 'ab', t0) - sn.bandmag(widebandx, 'ab', t0)
        colory = sn.bandmag(medbandy, 'ab', t0) - sn.bandmag(widebandy, 'ab', t0)
        pl.plot( colorx, colory, **plotargs )


    ax = pl.gca()
    ax.set_xlabel( '%s - %s'%( medbandx.upper(), widebandx.upper()))
    ax.set_ylabel( '%s - %s'%( medbandy.upper(), widebandy.upper()))
    return( ax )

def singlecircle( sn=None, sim1=None, sim2=None, sim3=None,
                  medbandx='f139m', medbandy='f153m',
                  contours=True, redshiftcircle=True,
                  clobber=False, filterset='hst', nsim=2000, **plotargs ) :
    """  make a single color-color circle diagram from sncosmo monte carlo sims.
    :param sn:
    :param sim1:
    :param sim2:
    :param contours:
    :param dotlines:
    :param clobber:
    :return:
    """
    import numpy as np
    from matplotlib import pyplot as pl

    if sn == 'stone' :
        sndat = STONE
    elif sn == 'colfax' :
        sndat = COLFAX
    elif sn == 'bush' :
        sndat = BUSH
    elif sn == 'jwst2' :
        sndat = JWST2
    elif sn == 'jwst3' :
        sndat = JWST3
    elif sn == 'demo' or sn=='none' or sn is None :
        sndat = DEMO

    if sim1 is None :
        sim1, sim2, sim3 = sncosmo_sim( snroot = sndat.name,
                                        z_range= sndat.z_range,
                                        t0_range=sndat.t0_range, nsim=nsim,
                                        clobber=clobber, filterset=filterset )

    if contours :
        plotcontours( sim1, sim2, sim3, medbandx=medbandx, medbandy=medbandy, **plotargs )
    else :
        plotpoints( sim1, sim2, sim3, medbandx=medbandx, medbandy=medbandy, **plotargs )

    if redshiftcircle :
        plot_redshift_circle(z_range=sndat.z_range, t0=sndat.t0,
                             medbandx=medbandx, medbandy=medbandy,
                             marker='o' )

    widebandx = medband_matching_filter(medbandx)
    widebandy = medband_matching_filter(medbandy)

    if sndat.name not in ['demo','none','jwst','jwst2','jwst3']:
        snmag= sndat.mag
        snmagerr = sndat.magerr
        deltasnmagx = snmag[medbandx] - snmag[widebandx]
        deltasnmagy = snmag[medbandy] - snmag[widebandy]

        deltasnmagerrx = np.sqrt(snmagerr[medbandx]**2+snmagerr[widebandx]**2)
        deltasnmagerry = np.sqrt(snmagerr[medbandy]**2+snmagerr[widebandy]**2)

        ax = pl.gca()
        ax.errorbar( deltasnmagx, deltasnmagy,
                     deltasnmagerrx, deltasnmagerry,
                     marker='D', ms=10, elinewidth=2, capsize=0,
                     color='k' )

    pl.draw()
    if sim3 is not None :
        return sim1, sim2, sim3
    return sim1,sim2



def doublecircle( sn='stone', sim1=None, sim2=None, sim3=None,
                  contours=True, redshiftcircle=True, clobber=False,
                  circle1bands=['f139m','f127m'], circle2bands=['f139m','f153m'],
                  filterset='hst', nsim=2000, **plotargs ):
    """  Two circle diagrams in side-by-side plots
    :param sn:
    :param sim1:
    :param sim2:
    :param contours:
    :param redshiftcircle:
    :param clobber:
    :param plotargs:
    :return:
    """
    from matplotlib import pyplot as pl
    fig = pl.gcf()
    fig.clf()
    ax1 = pl.subplot( 1, 2, 1 )
    simlist  = singlecircle( sn=sn, sim1=sim1, sim2=sim2, sim3=sim3,
                             medbandx=circle1bands[0], medbandy=circle1bands[1],
                             contours=contours, redshiftcircle=redshiftcircle,
                             clobber=clobber, filterset=filterset, nsim=nsim, **plotargs )
    sim1,sim2,sim3 = simlist

    ax2 = pl.subplot( 1, 2, 2 )
    singlecircle( sn=sn, sim1=sim1, sim2=sim2, sim3=sim3,
                  medbandx=circle2bands[0], medbandy=circle2bands[1],
                  contours=contours, redshiftcircle=redshiftcircle,
                  clobber=clobber, filterset=filterset, **plotargs )
    ax2.yaxis.set_ticks_position('right')
    ax2.yaxis.set_ticks_position('both')
    ax2.yaxis.set_label_position('right')
    ax2.yaxis.get_label().set_rotation(-90)

    return sim1, sim2, sim3


def sixcircles( sn='bush', sim1=None, sim2=None, sim3=None,
                  contours=True, redshiftcircle=True, clobber=False,
                  **plotargs ):
    """  Two circle diagrams in side-by-side plots
    :param sn:
    :param sim1:
    :param sim2:
    :param contours:
    :param redshiftcircle:
    :param clobber:
    :param plotargs:
    :return:
    """
    from matplotlib import pyplot as pl
    fig = pl.gcf()
    fig.clf()
    ax1 = pl.subplot( 2, 3, 1 )
    simlist  = singlecircle( sn=sn, sim1=sim1, sim2=sim2, sim3=sim3,
                             medbandx='f098m', medbandy='f127m',
                             contours=contours, redshiftcircle=redshiftcircle,
                             clobber=clobber, **plotargs )
    sim1,sim2,sim3 = simlist

    for iax,mbpair in zip( [2,3,4,5,6], [['f098m','f139m'], ['f098m','f153m'],
                                         ['f127m','f139m'], ['f127m','f153m'],
                                         ['f139m','f153m'] ] ) :
        ax = pl.subplot( 2, 3, iax )
        singlecircle( sn=sn, sim1=sim1, sim2=sim2, sim3=sim3,
                      medbandx=mbpair[0], medbandy=mbpair[1],
                      contours=contours, redshiftcircle=redshiftcircle,
                      clobber=False, **plotargs )
        #ax.yaxis.set_ticks_position('right')
        #ax.yaxis.set_ticks_position('both')
        #ax.yaxis.set_label_position('right')
        #ax.yaxis.get_label().set_rotation(-90)

    return sim1, sim2, sim3


def stonefig( simIa=None, simIbc=None, simII=None, contours=True, redshiftcircle=True,
              clobber=False, **plotargs ):
    from pytools import plotsetup
    fig = plotsetup.halfpaperfig(3,[4,4])

    simIa, simIbc, simII = singlecircle( 'stone', simIa, simIbc, simII,
                                         clobber=clobber,
                                         medbandx='f139m', medbandy='f153m',
                                         contours=contours, redshiftcircle=redshiftcircle,
                                         **plotargs )
    fig = pl.gcf()
    ax = pl.gca()
    fig.subplots_adjust( left=0.18, bottom=0.12, right=0.95, top=0.95, wspace=0.1 )

    ax.set_xlim( -0.35, 0.35 )
    ax.set_ylim( -0.3, 0.45 )

    ax.text(  0.95, 0.95, 'GND13Sto' , transform=ax.transAxes,
               ha='right',va='top', color='k',fontsize=15,)

    pl.draw()
    return simIa, simIbc, simII


def colfaxfig( simIa=None, simIbc=None, simII=None, contours=True, redshiftcircle=True,
               clobber=False, **plotargs ):
    from pytools import plotsetup
    from matplotlib import pyplot as pl
    fig = plotsetup.fullpaperfig([8,4])

    simIa, simIbc, simII = doublecircle( 'colfax', simIa, simIbc, simII,
                                         clobber=clobber,
                                         circle1bands=['f139m','f127m'],
                                         circle2bands=['f139m','f153m'],
                                         contours=contours, redshiftcircle=redshiftcircle,
                                         **plotargs )
    fig = pl.gcf()
    ax1,ax2 = fig.axes
    fig.subplots_adjust( left=0.1, bottom=0.18, right=0.89, top=0.95, wspace=0.1 )

    ax1.set_xlim( -0.5, 0.4 )
    ax1.set_ylim( -0.7, 0.4 )
    ax2.set_xlim( -0.5, 0.4 )
    ax2.set_ylim( -0.4, 0.4 )
    ax2.text(  0.95, 0.95, 'GND12Col' , transform=ax2.transAxes,
               ha='right',va='top', color='k',fontsize=15,)
    ax2.yaxis.labelpad=20

    pl.draw()

    return simIa, simIbc, simII

def bushfig( simIa=None, simIbc=None, simII=None, contours=True, redshiftcircle=False,
             clobber=False, **plotargs ):
    from pytools import plotsetup
    from matplotlib import pyplot as pl, ticker

    fig = plotsetup.fullpaperfig( [8,4] )

    simIa, simIbc, simII = doublecircle( 'bush', simIa, simIbc, simII, clobber=clobber,
                                         circle1bands=['f098m','f153m'], circle2bands=['f127m','f139m'],
                                         contours=contours, redshiftcircle=redshiftcircle,
                                         **plotargs )
    fig = pl.gcf()
    ax1, ax2 = fig.axes
    fig.subplots_adjust( left=0.1, bottom=0.18, right=0.89, top=0.95, wspace=0.1 )

    ax1.set_xlim( -0.6, 0.95 )
    ax1.set_ylim( -0.9, 0.5 )
    ax2.set_xlim( -0.6, 0.6 )
    ax2.set_ylim( -0.8, 0.6 )
    ax2.yaxis.labelpad=20

    for ax in [ax1,ax2]:
        ax.xaxis.set_major_locator( ticker.MultipleLocator( 0.2 ) )
        ax.xaxis.set_minor_locator( ticker.MultipleLocator( 0.1 ) )
        ax.yaxis.set_major_locator( ticker.MultipleLocator( 0.2 ) )
        ax.yaxis.set_minor_locator( ticker.MultipleLocator( 0.1 ) )

    #ax.text( 0.1, 0.1, 'SN Ia', color=color1,
    #          ha='left',va='top')#,transform=ax1.transAxes,fontsize=14)
    #ax.text( 0., -0.05, 'CC SN', color=color2,
    #          ha='center',va='top')#,transform=ax1.transAxes,fontsize=14)
    ax2.text(  0.95, 0.95, 'GSD11Bus' , transform=ax2.transAxes,
               ha='right',va='top', color='k',fontsize=15,)
    pl.draw()

    return simIa, simIbc, simII


def classdemofig( sim1=None, sim2=None, sim3=None,
                  contours=True,
                  clobber=False, **plotargs ):
    """  Make a circle figure demonstrating the segregation of SN types
    in pseudo-color space for classification.

    :param sim1:
    :param sim2:
    :param sim3:
    :param contours:
    :param clobber:
    :param plotargs:
    :return:
    """
    from pytools import plotsetup
    fig = plotsetup.fullpaperfig( figsize=[8,4] )

    sim1, sim2, sim3  = doublecircle( None, sim1, sim2, sim3, clobber=clobber,
                                      contours=contours, redshiftcircle=False,
                                      **plotargs )
    ax1 = fig.add_subplot( 1, 2, 1 )
    ax2 = fig.add_subplot( 1, 2, 2, sharex=ax1 )
    if contours :
        ax1.text( -0.25, -0.3, 'SN Ia', color=_COLOR1, ha='right',va='top', fontsize=14)
        ax1.text( 0.0, -0.08, 'SN II', color=_COLOR2, ha='center',va='center')
        ax1.text( -0.2, 0.25, 'SN Ib/c', color=_COLOR3, ha='right',va='top')
    else :
        ax1.text( -0.25, -0.3, 'SN Ia', color=_COLOR1, ha='right',va='top', fontsize=14)
        ax1.text( -0.25, -0.08, 'SN II', color=_COLOR2, ha='right',va='center')
        ax1.text( -0.25, 0.25, 'SN Ib/c', color=_COLOR3, ha='right',va='top')

    fig.subplots_adjust( left=0.12, bottom=0.14, right=0.88, top=0.95, wspace=0.1 )

    ax1.set_xlim( -0.45, 0.4 )
    ax1.set_ylim( -0.45, 0.3 )
    ax2.set_ylim( -0.25, 0.3 )


    return( sim1, sim2, sim3 )


def jwstfig_z3( simIa=None, simIbc=None, simII=None,
                contours=True,
                redshiftcircle=True, clobber=False, nsim=2000,
                **plotargs ):
    from pytools import plotsetup
    from matplotlib import pyplot as pl
    fig = plotsetup.fullpaperfig([8,4])

    simIa, simIbc, simII = doublecircle( 'jwst3', simIa, simIbc, simII,
                                         clobber=clobber,
                                         circle1bands=['f182m','f300m'],
                                         circle2bands=['f210m','f300m'],
                                         contours=contours, redshiftcircle=redshiftcircle,
                                         filterset='jwst_z3', nsim=nsim, **plotargs)
    ax1,ax2 = fig.axes
    fig.subplots_adjust( left=0.1, bottom=0.18, right=0.89, top=0.95, wspace=0.1 )
    ax2.yaxis.labelpad=20
    pl.draw()
    return simIa, simIbc, simII


def jwstfig_z2( simIa=None, simIbc=None, simII=None,
                contours=True,
                redshiftcircle=True, clobber=False, nsim=2000,
                **plotargs ):
    from pytools import plotsetup
    from matplotlib import pyplot as pl
    fig = plotsetup.fullpaperfig([8,4])

    simIa, simIbc, simII = doublecircle( 'jwst2', simIa, simIbc, simII,
                                         clobber=clobber,
                                         circle1bands=['f182m','f210m'],
                                         circle2bands=['f210m','f300m'],
                                         contours=contours, redshiftcircle=redshiftcircle,
                                         filterset='jwst_z3', nsim=nsim, **plotargs)
    ax1,ax2 = fig.axes
    fig.subplots_adjust( left=0.1, bottom=0.18, right=0.89, top=0.95, wspace=0.1 )
    ax2.yaxis.labelpad=20
    pl.draw()
    return simIa, simIbc, simII
