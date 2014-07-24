__author__ = 'rodney'
import mpltools
from matplotlib import cm
import numpy as np
from matplotlib import pyplot as pl
from __init__ import _ALPHA2FILTER


def plotbands(z = 1.2):
    """ plot a SNIa spectrum at the given z, overlaid with medium bands
    :param z:
    :return:
    """
    import stardust
    from demofig import w763,f763,w845,f845,w139,f139
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


def pseudocolor_vs_z( simgridIa, medbands='OPQ', broadbands='JNH', age=0 ):
    """  Plot the med-broad pseudo colors vs z.
    :param simgridIa:
    :param medbands:
    :param broadbands:
    :param age:
    :return:
    """
    # TODO : allow user to input colors, not separated bands
    cmap = cm.jet

    nrows = len(broadbands)
    for medband, broadband, irow in zip(medbands,broadbands,range(nrows)) :
        medfiltname = _ALPHA2FILTER[medband]
        broadfiltname = _ALPHA2FILTER[broadband]
        iB = simgridIa.BANDS.index( broadband )
        iM = simgridIa.BANDS.index( medband )

        irv=0
        iage = np.argmin(np.abs(simgridIa.TREST-age))

        ax1 = pl.subplot(nrows,2,1)
        ax2 = pl.subplot(nrows,2,2,sharex=ax1, sharey=ax1)

        mpltools.color.cycle_cmap( len(simgridIa.LUMIPAR),  cmap=cmap, ax=ax1)
        mpltools.color.cycle_cmap( len(simgridIa.COLORPAR), cmap=cmap, ax=ax2)

        ic = np.argmin(np.abs(simgridIa.COLORPAR-0))
        for x1 in simgridIa.LUMIPAR :
            ix1 = np.argmin(np.abs(simgridIa.LUMIPAR-x1))

            B = simgridIa.LCMATRIX[ ix1, irv, ic, :, iB, iage ]
            M = simgridIa.LCMATRIX[ ix1, irv, ic, :, iM, iage ]
            ax1.plot( z, M-B,  marker='o', ls=' ' )

        ax1.text( 0.15, 0.05, 'x1=%.1f'%simgridIa.LUMIPAR[0], color=cmap(0.0) , transform=ax1.transAxes)
        ax1.text( 0.35, 0.05, '.. %.1f'%simgridIa.LUMIPAR[len(simgridIa.LUMIPAR)/2], color=cmap(0.5) , transform=ax1.transAxes)
        ax1.text( 0.55, 0.05, '.. %.1f'%simgridIa.LUMIPAR[-1], color=cmap(1.0) , transform=ax1.transAxes)
        ax2.text( 0.15, 0.05, 'c=%.1f' %simgridIa.COLORPAR[0], color=cmap(0.0), transform=ax2.transAxes )
        ax2.text( 0.35, 0.05, '.. %.1f'%simgridIa.COLORPAR[len(simgridIa.LUMIPAR)/2], color=cmap(0.5), transform=ax2.transAxes )
        ax2.text( 0.55, 0.05, '.. %.1f'%simgridIa.COLORPAR[-1], color=cmap(1.0), transform=ax2.transAxes )

        ix1 = np.argmin(np.abs(simgridIa.LUMIPAR-0))
        for c in simgridIa.COLORPAR :
            ic = np.argmin(np.abs(simgridIa.COLORPAR-c))

            B = simgridIa.LCMATRIX[ ix1, irv, ic, :, iB, iage ]
            M = simgridIa.LCMATRIX[ ix1, irv, ic, :, iM, iage ]

            z = simgridIa.z
            ax2.plot( z, M-B, marker='o', ls=' ' )
        ax1.set_xlabel('redshift')
        ax2.set_xlabel('redshift')
        ax1.set_ylabel('%s-%s'%(medfiltname, broadfiltname) )
        pl.setp( ax2.yaxis.get_ticklabels(), visible=False )

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
    import gridsim
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
            simIaGrid = gridsim.dosimGrid( sn, ngridz=50, bands='X7I8LYOJPNQH',
                                             x1range=x1range, crange=crange,
                                             trestrange = agerange,
                                             clobber=clobber, verbose=verbose )
        elif isinstance( simIaGrid, str ) :
            simIaGrid = stardust.SimTable(simIaGrid)


    if verbose and plotstyle.startswith('contour'):
        print('Binning up MC sim for color-color diagram...')
    pl.clf()
    if color3 is not None : ax1 = pl.subplot( 2,2,1 )
    if showMCpoints :
        stardust.simplot.plotColorColor( simIaMC, color1,color2, binrange=[[-0.2,0.8],[-0.5,0.5]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
        stardust.simplot.plotColorColor( simIbcMC, color1,color2, binrange=[[-0.2,0.8],[-0.5,0.5]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
        stardust.simplot.plotColorColor( simIIMC,  color1,color2, binrange=[[-0.2,0.8],[-0.5,0.5]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
    if showGridline :
        gridsim.plotgridCircle( simIaGrid, color1, color2, x1=x1, c=c, age=age )

    if color3 is not None :
        ax2 = pl.subplot( 2,2,2 )
        if showMCpoints :
            stardust.simplot.plotColorColor( simIaMC,  color2, color3, binrange=[[-0.6,0.4],[-0.5,0.5]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
            stardust.simplot.plotColorColor( simIbcMC, color2, color3, binrange=[[-0.6,0.4],[-0.5,0.5]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
            stardust.simplot.plotColorColor( simIIMC,  color2, color3, binrange=[[-0.6,0.4],[-0.5,0.5]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
        if showGridline :
            gridsim.plotgridCircle( simIaGrid, color2, color3, x1=x1, c=c, age=age )

        ax3 = pl.subplot( 2,2,3 )
        if showMCpoints :
            stardust.simplot.plotColorColor( simIaMC,  color1,color3, binrange=[[-0.2,0.8],[-0.6,0.4]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
            stardust.simplot.plotColorColor( simIbcMC, color1, color3, binrange=[[-0.2,0.8],[-0.6,0.4]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
            stardust.simplot.plotColorColor( simIIMC,  color1, color3, binrange=[[-0.2,0.8],[-0.6,0.4]], mjdrange=[mjdobs,mjdobs], tsample=1, plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
        if showGridline :
            gridsim.plotgridCircle( simIaGrid, color1,color3, x1=x1, c=c, age=age )

    pl.draw()
    return (sn,simdataMC,simIaMC,simIaGrid)



