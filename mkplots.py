__author__ = 'rodney'
import mpltools
from matplotlib import cm
import numpy as np
from matplotlib import pyplot as pl
from __init__ import _ALPHA2FILTER


def pseudocolor_vs_z( simgridIa, medbands='OPQ', broadbands='JNH' ):
    """  Plot the med-broad pseudo colors vs z.
    :param simgridIa:
    :param medbands:
    :param broadbands:
    :param age:
    :return:
    """
    # TODO : allow user to input colors, not separated bands
    import mpltools
    cmap = cm.jet
    z = simgridIa.z

    nrows = len(broadbands)
    for medband, broadband, irow in zip(medbands,broadbands,range(nrows)) :
        medfilt = _ALPHA2FILTER[medband]
        broadfilt = _ALPHA2FILTER[broadband]
        iB = simgridIa.BANDS.index( broadband )
        iM = simgridIa.BANDS.index( medband )

        irv=0

        ax1 = pl.subplot(nrows,3,1+irow*3)
        ax2 = pl.subplot(nrows,3,2+irow*3,sharex=ax1, sharey=ax1)
        ax3 = pl.subplot(nrows,3,3+irow*3,sharex=ax1, sharey=ax1)

        mpltools.color.cycle_cmap( len(simgridIa.LUMIPAR),  cmap=cmap, ax=ax1)
        mpltools.color.cycle_cmap( len(simgridIa.COLORPAR), cmap=cmap, ax=ax2)

        iage = np.argmin(np.abs(simgridIa.TREST-0))
        ic = np.argmin(np.abs(simgridIa.COLORPAR-0))
        for x1 in simgridIa.LUMIPAR :
            ix1 = np.argmin(np.abs(simgridIa.LUMIPAR-x1))

            B = simgridIa.LCMATRIX[ ix1, irv, ic, :, iB, iage ]
            M = simgridIa.LCMATRIX[ ix1, irv, ic, :, iM, iage ]
            ax1.plot( z, M-B,  marker='o', ls=' ' )

        ix1 = np.argmin(np.abs(simgridIa.LUMIPAR-0))
        iage = np.argmin(np.abs(simgridIa.TREST-0))
        for c in simgridIa.COLORPAR :
            ic = np.argmin(np.abs(simgridIa.COLORPAR-c))
            B = simgridIa.LCMATRIX[ ix1, irv, ic, :, iB, iage ]
            M = simgridIa.LCMATRIX[ ix1, irv, ic, :, iM, iage ]
            ax2.plot( z, M-B, marker='o', ls=' ' )

        ix1 = np.argmin(np.abs(simgridIa.LUMIPAR-0))
        ic = np.argmin(np.abs(simgridIa.COLORPAR-0))
        for age in simgridIa.TREST :
            iage = np.argmin(np.abs(simgridIa.TREST-age))
            B = simgridIa.LCMATRIX[ ix1, irv, ic, :, iB, iage ]
            M = simgridIa.LCMATRIX[ ix1, irv, ic, :, iM, iage ]
            ax3.plot( z, M-B, marker='o', ls=' ' )

        # TODO : use http://tonysyu.github.io/mpltools/auto_examples/color/plot_color_mapper.html

        if irow==0 :
            ax1.text( 0.15, 0.05, 'x1 : %.1f'%simgridIa.LUMIPAR[0], color=cmap(0.0) , transform=ax1.transAxes)
            ax1.text( 0.4, 0.05, '.. %.1f'%simgridIa.LUMIPAR[-1], color=cmap(1.0) , transform=ax1.transAxes)
            ax1.text( 0.05, 0.95, 'c=0,age=0', color='k', transform=ax1.transAxes, ha='left',va='top')
            ax2.text( 0.15, 0.05, 'c : %.1f' %simgridIa.COLORPAR[0], color=cmap(0.0), transform=ax2.transAxes )
            ax2.text( 0.4, 0.05, '.. %.1f'%simgridIa.COLORPAR[-1], color=cmap(1.0), transform=ax2.transAxes )
            ax2.text( 0.05, 0.95, 'x1=0,age=0', color='k', transform=ax2.transAxes, ha='left',va='top')
            ax3.text( 0.15, 0.05, 'age : %.1f'%simgridIa.TREST[0], color=cmap(0.0), transform=ax3.transAxes )
            ax3.text( 0.4, 0.05, '.. %.1f'% simgridIa.TREST[-1], color=cmap(1.0), transform=ax3.transAxes )
            ax3.text( 0.05, 0.95, 'x1=0,c=0', color='k', transform=ax3.transAxes, ha='left',va='top')
        if irow == nrows-1 :
            ax1.set_xlabel('redshift')
            ax2.set_xlabel('redshift')
            ax3.set_xlabel('redshift')
        ax1.set_ylabel('%s-%s'%(medfilt, broadfilt) )
        pl.setp( ax2.yaxis.get_ticklabels(), visible=False )
    return


def circlefig( sn, simdataMC=None, simIaGrid=None,
               color1='O-J', color2='P-N', color3=None,
               showMCpoints=True, showGridline=False,
               snanadatfile=None, plotstyle='points',
               fast=False, clobber=0, verbose=1, ):
    """  Plot the results of a SNANA monte carlo simulation as a circle diagram.

    :param simdataMC:
    :param linelevels:
    :param plotstyle:
    :param Nbins:
    :param nsim:
    :param clobber:
    :param verbose:
    :param snanadatfile:
    :return:
    """
    import snanasim
    import stardust
    import numpy as np
    from matplotlib import pyplot as pl

    # Read in the SN, classify it, compute marginalized posterior distributions
    # for the light curve parameters  (x1,c,tpk,Av,etc)
    assert sn is not None, "Must provide a light curve .dat file!"
    if isinstance( sn, str ) :
        sn = stardust.SuperNova(sn)
    if 'ClassSim' not in sn.__dict__ :
        print("Running doGridClassify for %s with medium bands excluded."%sn.name)
        sn.doGridClassify( bands=sn.bands.translate(None,'78LOPQ'),
                           useLuminosityPrior=True,
                           nlogz=(fast and 10) or 15,
                           ncolorpar=(fast and 5) or 15,
                           ncolorlaw=1,
                           nlumipar=(fast and 6) or 15,
                           npkmjd=(fast and 10) or 20,
                           omitTemplateII='IIL', clobber=clobber )
    if 'x1_maxprob' not in sn.ClassSim.Ia.__dict__ :
        # TODO : plotting should not be needed, once the marginalization steps are extracted from plotClassStatsGrid
        sn.plotClassStatsGrid()

    # Define the date at which the med bands were observed
    mjdmedband = sn.MJD[ np.where( (sn.FLT=='7') | (sn.FLT=='8') | (sn.FLT=='L')
                         | (sn.FLT=='O') | (sn.FLT=='P') | (sn.FLT=='Q')) ]
    if len(mjdmedband)>0 :
        mjdobs = np.median( mjdmedband )
    else :
        mjdobs = sn.pkmjd

    # Best-fit SNIa LC parameter values, as determined without med-band data
    z, zerr  = sn.ClassSim.Ia.z_maxprob, sn.ClassSim.Ia.z_maxprob_err,
    x1, x1err  = sn.ClassSim.Ia.x1_maxprob, sn.ClassSim.Ia.x1_maxprob_err,
    c, cerr= sn.ClassSim.Ia.c_maxprob, sn.ClassSim.Ia.c_maxprob_err
    mjdpk, mjdpkerr =  sn.pkmjd+sn.ClassSim.Ia.tpk_maxprob, sn.ClassSim.Ia.tpk_maxprob_err
    x1range = [ max(-3.0,x1-3*x1err), min(3.0,x1+3*x1err) ]
    crange = [ max(-0.5,c-3*cerr), min(1.5,c+3*cerr) ]
    # mjdpkrange = [ mjdpk-3*mjdpkerr, mjdpk+3*mjdpkerr ]
    age, ageerr = (mjdobs-mjdpk)/(1+z), (mjdpkerr)/(1+z)
    agerange = [ age-3*ageerr, age+3*ageerr ]

    if fast :
        Nbins=30
        Nsim=500
        linelevels = [ 0, 0.82 ]
    else :
        Nbins=80
        Nsim=2000
        linelevels = [ 0, 0.82 ]

    if showMCpoints :
        # Run or read the med band Monte Carlo simulation :
        if simdataMC is None :
            simdataMC = snanasim.dosimMC( sn,  Nsim=Nsim, bands='XI78YJNHLOPQ',
                                       clobber=clobber, verbose=verbose )
            simIaMC, simIbcMC, simIIMC = simdataMC
        elif isinstance( simdataMC, str ):
            simIaMC = stardust.SimTable(simdataMC+'_Ia')
            simIbcMC = stardust.SimTable(simdataMC+'_Ibc')
            simIIMC = stardust.SimTable(simdataMC+'_II')
        else :
            simIaMC, simIbcMC, simIIMC = simdataMC

    if showGridline :
        # Run or read in the grid simulation
        if simIaGrid is None :
            simIaGrid = snanasim.dosimGrid( sn, ngridz=50, bands='X7I8LYOJPNQH',
                                             x1range=x1range, crange=crange,
                                             trestrange = agerange,
                                             clobber=clobber, verbose=verbose )
        elif isinstance( simIaGrid, str ) :
            simIaGrid = stardust.SimTable(simIaGrid)


    if verbose and plotstyle.startswith('contour'):
        print('Binning up MC sim for color-color diagram...')
    pl.clf()
    if color3 is not None :
        ax1 = pl.subplot( 2,2,1 )
    if showMCpoints :
        stardust.simplot.plotColorColor( simIaMC, color1,color2, binrange=[[-0.2,0.8],[-0.5,0.5]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
        stardust.simplot.plotColorColor( simIbcMC, color1,color2, binrange=[[-0.2,0.8],[-0.5,0.5]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
        stardust.simplot.plotColorColor( simIIMC,  color1,color2, binrange=[[-0.2,0.8],[-0.5,0.5]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
    if showGridline :
        gridsim_circleplot( simIaGrid, color1, color2, x1=x1, c=c, age=age )

    if color3 is not None :
        ax2 = pl.subplot( 2,2,2 )
        if showMCpoints :
            stardust.simplot.plotColorColor( simIaMC,  color2, color3, binrange=[[-0.6,0.4],[-0.5,0.5]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
            stardust.simplot.plotColorColor( simIbcMC, color2, color3, binrange=[[-0.6,0.4],[-0.5,0.5]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
            stardust.simplot.plotColorColor( simIIMC,  color2, color3, binrange=[[-0.6,0.4],[-0.5,0.5]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
        if showGridline :
            gridsim_circleplot( simIaGrid, color2, color3, x1=x1, c=c, age=age )

        ax3 = pl.subplot( 2,2,3 )
        if showMCpoints :
            stardust.simplot.plotColorColor( simIaMC,  color1,color3, binrange=[[-0.2,0.8],[-0.6,0.4]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
            stardust.simplot.plotColorColor( simIbcMC, color1, color3, binrange=[[-0.2,0.8],[-0.6,0.4]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
            stardust.simplot.plotColorColor( simIIMC,  color1, color3, binrange=[[-0.2,0.8],[-0.6,0.4]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
        if showGridline :
            gridsim_circleplot( simIaGrid, color1,color3, x1=x1, c=c, age=age )

    pl.draw()
    return (sn,simdataMC,simIaGrid)



def gridsim_circleplot( simgridIa, color1, color2, x1=0, c=0, age=0 ):
    """ Plot a SNANA grid simulation into a color-color circle diagram.
    :param simdatGrid:
    :param color1:
    :param color2:
    :return:
    """

    medband1,broadband1 = color1.split('-')
    # medfilt1 = ALPHA2FILTER[medband1]
    # broadfilt1 = ALPHA2FILTER[broadband1]
    iB1 = simgridIa.BANDS.index( broadband1 )
    iM1 = simgridIa.BANDS.index( medband1 )

    medband2,broadband2 = color2.split('-')
    # medfilt2 = ALPHA2FILTER[medband2]
    # broadfilt2 = ALPHA2FILTER[broadband2]
    iB2 = simgridIa.BANDS.index( broadband2 )
    iM2 = simgridIa.BANDS.index( medband2 )

    irv=0
    cmap = cm.jet
    z = simgridIa.z

    iage = np.argmin(np.abs(simgridIa.TREST-age))
    ic = np.argmin(np.abs(simgridIa.COLORPAR-c))
    ix1 = np.argmin(np.abs(simgridIa.LUMIPAR-x1))
    B1 = simgridIa.LCMATRIX[ ix1, irv, ic, :, iB1, iage ]
    M1 = simgridIa.LCMATRIX[ ix1, irv, ic, :, iM1, iage ]

    B2 = simgridIa.LCMATRIX[ ix1, irv, ic, :, iB2, iage ]
    M2 = simgridIa.LCMATRIX[ ix1, irv, ic, :, iM2, iage ]

    mpltools.color.cycle_cmap( len(z),  cmap=cmap )
    for iz in range(len(z)) :
        pl.plot( M1[iz:iz+2]-B1[iz:iz+2], M2[iz:iz+2]-B2[iz:iz+2], marker=' ', ls='-', lw=2 )

    return
