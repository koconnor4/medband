# 2013.02.12
# S.Rodney
"""
Figure to demonstrate the use of medium band filters 
for efficient SN classification and redshift estimation
"""
import numpy as np
from matplotlib import pyplot as pl
import os
import stardust

sndataroot = os.environ['SNDATA_ROOT']
topdir = os.path.abspath( '.' )

_medbanddat = 'ffsn_exampleIa_medband.dat'

# Read in the filter transmission curves
try: 
    os.chdir( sndataroot+'/filters/HST/HST_CANDELS')

    w350, f350 = np.loadtxt( 'WFC3_UVIS_F350LP.dat', unpack=True )
    w763, f763 = np.loadtxt( 'WFC3_UVIS_F763M.dat', unpack=True )
    w845, f845 = np.loadtxt( 'WFC3_UVIS_F845M.dat', unpack=True )
    w814, f814 = np.loadtxt( 'ACS_WFC_F814W.dat', unpack=True )
    w127, f127 = np.loadtxt( 'WFC3_IR_F127M.dat', unpack=True )
    w125, f125 = np.loadtxt( 'WFC3_IR_F125W.dat', unpack=True )
    w160, f160 = np.loadtxt( 'WFC3_IR_F160W.dat', unpack=True )
    w153, f153 = np.loadtxt( 'WFC3_IR_F153M.dat', unpack=True )
    w139, f139 = np.loadtxt( 'WFC3_IR_F139M.dat', unpack=True )
    w140, f140 = np.loadtxt( 'WFC3_IR_F140W.dat', unpack=True )

    #os.chdir( sndataroot+'/snsed') 
finally:
    os.chdir(topdir)





def dosimMC( sn, Nsim=500, bands='XI78YJNHLOPQ',
             smear=False, clobber=0, verbose=1 ) :
    """ set up and run a SNANA monte carlo simulation to show the 
    medium band technique.
    smear : when false, don't include any heterogeneity for the Ia population
        i.e. only simulate x1=0, c=0.  And don't add additional mag-smearing
        to either the Ia or CC models.
    """
    if 'ClassSim' not in sn.__dict__ :
        print("Running doGridClassify for %s with medium bands excluded."%sn.name)
        sn.doGridClassify( bands=sn.bands.translate(None,'78LOPQ'),
                           useLuminosityPrior=True, clobber=clobber,
                           # nlogz=15, ncolorpar=15, ncolorlaw=1, nlumipar=15, npkmjd=20,
                           nlogz=10, ncolorpar=5, ncolorlaw=1, nlumipar=6, npkmjd=10,
                           omitTemplateII='IIL' )

    if 'x1_maxprob' not in sn.ClassSim.Ia.__dict__ :
        # TODO : this should not be needed, once the marginalization is  extracted from plotClassStatsGrid
        sn.plotClassStatsGrid()

    mjdmedband = sn.MJD[ np.where( (sn.FLT=='7') | (sn.FLT=='8') |
                                   (sn.FLT=='L') | (sn.FLT=='O') |
                                   (sn.FLT=='P') | (sn.FLT=='Q') ) ]
    if len(mjdmedband)>0 :
        mjdobs = np.median( mjdmedband )
    else :
        mjdobs = sn.pkmjd

    # First, make a .simlib file
    simname = 'sim_'+sn.nickname+'_medbandMC'
    stardust.simulate.mksimlib( simname+'.simlib', survey='HST', field='default',
                                bands=bands, etimes=[10000]*len(bands), mjd0=mjdobs,
                                cadence=1, Nepoch=1, clobber=clobber, perfect=True )


    # -------------------------------------
    #  Type Ia simulation
    # -------------------------------------
    # define the parameter ranges allowed by the broad-band LC fit
    x1range = [ max(-3.0,sn.ClassSim.Ia.x1_maxprob-3*sn.ClassSim.Ia.x1_maxprob_err),
                min(3.0,sn.ClassSim.Ia.x1_maxprob+3*sn.ClassSim.Ia.x1_maxprob_err) ]
    crange = [ max(-0.5,sn.ClassSim.Ia.c_maxprob-3*sn.ClassSim.Ia.c_maxprob_err),
               min(1.5,sn.ClassSim.Ia.c_maxprob+3*sn.ClassSim.Ia.c_maxprob_err) ]
    mjdpkrange = [ sn.pkmjd+sn.ClassSim.Ia.tpk_maxprob-3*sn.ClassSim.Ia.tpk_maxprob_err,
                   sn.pkmjd+sn.ClassSim.Ia.tpk_maxprob+3*sn.ClassSim.Ia.tpk_maxprob_err ]
    # we allow a generous range of possible redshifts, up to dz=0.2 or 10-sigma
    zrange = [ sn.ClassSim.Ia.z_maxprob-min(0.2,10*sn.ClassSim.Ia.z_maxprob_err),
               sn.ClassSim.Ia.z_maxprob+min(0.2,10*sn.ClassSim.Ia.z_maxprob_err) ]
    # make the .input file for the Ia sim
    stardust.simulate.mkinput( simname+'_Ia', inputfile=simname+'_Ia.input',
                               simlibfile=simname+'.simlib',GENFILTERS=bands,
                               GENRANGE_SALT2x1=x1range,GENRANGE_SALT2c=crange,
                               GENRANGE_PEAKMJD=mjdpkrange,
                               GENPERFECT=2, GENRANGE_REDSHIFT=zrange, smear=smear,
                               simtype='Ia', ratemodel='flat', Nsim=Nsim, clobber=clobber )

    # -------------------------------------
    #  Type Ibc simulation
    # -------------------------------------
    # define the parameter ranges allowed by the broad-band LC fit
    avrange = [ max(0,sn.ClassSim.Ibc.Av_maxprob-3*sn.ClassSim.Ibc.Av_maxprob_err),
                sn.ClassSim.Ibc.Av_maxprob+3*sn.ClassSim.Ibc.Av_maxprob_err ]
    mjdpkrange = [ sn.pkmjd+sn.ClassSim.Ibc.tpk_maxprob-3*sn.ClassSim.Ibc.tpk_maxprob_err,
                   sn.pkmjd+sn.ClassSim.Ibc.tpk_maxprob+3*sn.ClassSim.Ibc.tpk_maxprob_err ]
    # we allow a generous range of possible redshifts
    zrange = [ sn.ClassSim.Ibc.z_maxprob-min(0.2,10*sn.ClassSim.Ibc.z_maxprob_err),
               sn.ClassSim.Ibc.z_maxprob+min(0.2,10*sn.ClassSim.Ibc.z_maxprob_err) ]
    stardust.simulate.mkinput( simname+'_Ibc', inputfile=simname+'_Ibc.input',
                               simlibfile=simname+'.simlib',GENFILTERS=bands,
                               GENPERFECT=2, GENRANGE_REDSHIFT=zrange,
                               GENRANGE_PEAKMJD=mjdpkrange, GENRANGE_AV=avrange,
                               smear=smear,
                               simtype='Ibc', ratemodel='flat', Nsim=Nsim, clobber=clobber )

    # -------------------------------------
    #  Type II simulation
    # -------------------------------------
    # define the parameter ranges allowed by the broad-band LC fit
    avrange = [ max(0,sn.ClassSim.II.Av_maxprob-3*sn.ClassSim.II.Av_maxprob_err),
                sn.ClassSim.II.Av_maxprob+3*sn.ClassSim.II.Av_maxprob_err ]
    mjdpkrange = [ sn.pkmjd+sn.ClassSim.II.tpk_maxprob-3*sn.ClassSim.II.tpk_maxprob_err,
                   sn.pkmjd+sn.ClassSim.II.tpk_maxprob+3*sn.ClassSim.II.tpk_maxprob_err ]
    # we allow a generous range of possible redshifts
    zrange = [ sn.ClassSim.II.z_maxprob-min(0.2,10*sn.ClassSim.II.z_maxprob_err),
               sn.ClassSim.II.z_maxprob+min(0.2,10*sn.ClassSim.II.z_maxprob_err) ]
    stardust.simulate.mkinput( simname+'_II', inputfile=simname+'_II.input',
                               simlibfile=simname+'.simlib',GENFILTERS=bands,
                               GENPERFECT=2, GENRANGE_REDSHIFT=zrange, GENRANGE_AV=avrange,
                               GENRANGE_PEAKMJD=mjdpkrange, smear=smear,
                               simtype='II', rate='flat', Nsim=Nsim, clobber=clobber )

    stardust.simulate.dosim('%s_Ia.input'%simname, perfect=True, verbose=verbose )
    stardust.simulate.dosim('%s_Ibc.input'%simname, perfect=True, verbose=verbose )
    stardust.simulate.dosim('%s_II.input'%simname, perfect=True, verbose=verbose )

    return readSimData( simname )

def readSimData(simname='simFF_medband'):
    simIa = stardust.SimTable( simname+'_Ia')
    simIbc = stardust.SimTable( simname+'_Ibc')
    simII = stardust.SimTable( simname+'_II')
    return [simIa, simIbc, simII]


#def circlefig( simdata, , linelevels = [ 0, 0.82 ], plotstyle='contourf', Nbins=80, showsn=False ):

def mkMedBandFig( simdata, linelevels = [ 0, 0.82 ], plotstyle='contourf', Nbins=80, showsn=False ):
    """ construct the medium-band demo figure for the FrontierSN proposal: 
    filter bandpasses on the left for three redshifts, two color-color plots on the 
    right with SN observation points overlaid """

    w1a, f1a = stardust.snsed.getsed( sedfile='/usr/local/SNDATA_ROOT/snsed/Hsiao07.dat', day=0 )

    # first the band-pass plots on the left
    z = 1.8
    w1az = w1a * (1+z)
    f1az = f1a / f1a.max() / 2.
    ax18 = pl.subplot(3,2,1)
    ax18.plot(w1az, f1az, ls='-', lw=0.7, color='0.5', label='_nolegend_')
    ax18.plot(w127, f127, ls='-', color='DarkOrchid',label='F127M')
    ax18.plot(w139, f139, ls='-',color='Teal', label='F139M')
    ax18.plot(w153, f153, ls='-',color='Maroon', label='F153M')
    ax18.fill_between( w1az, f1az, where=((w1az>12400) & (w1az<13120)), color='DarkOrchid', alpha=0.3 ) 
    ax18.fill_between( w1az, f1az, where=((w1az>13500) & (w1az<14150)), color='teal', alpha=0.3 )
    ax18.fill_between( w1az, f1az, where=((w1az>15000) & (w1az<15700)), color='Maroon', alpha=0.3 )
    ax18.text(0.95,0.4, 'SNIa\n@ z=%.1f'%(z), color='k',ha='right',va='bottom',fontweight='bold', transform=ax18.transAxes, fontsize='large')
    pl.setp(ax18.get_xticklabels(), visible=False)
    pl.setp(ax18.get_yticklabels(), visible=False)
    ax18.text( 12700, 0.65, 'F127M', ha='right', va='center', color='DarkOrchid', fontweight='bold')
    ax18.text( 13900, 0.65, 'F139M', ha='center', va='center', color='Teal', fontweight='bold')
    ax18.text( 15300, 0.65, 'F153M', ha='left', va='center', color='Maroon', fontweight='bold')

    z = 2.0
    w1az = w1a * (1+z)
    f1az = f1a / f1a.max() / 2.
    ax20 = pl.subplot(3,2,3, sharex=ax18)
    ax20.plot(w127, f127, ls='-', color='DarkOrchid',label='F127M')
    ax20.plot(w139, f139, ls='-',color='Teal', label='F139M')
    ax20.plot(w153, f153, ls='-',color='Maroon', label='F153M')
    ax20.plot(w1az, f1az, ls='-', lw=0.7, color='0.5', label='_nolegend_')
    ax20.fill_between( w1az, f1az, where=((w1az>12400) & (w1az<13120)), color='DarkOrchid', alpha=0.3 )
    ax20.fill_between( w1az, f1az, where=((w1az>13500) & (w1az<14150)), color='Teal', alpha=0.3 )
    ax20.fill_between( w1az, f1az, where=((w1az>15000) & (w1az<15700)), color='Maroon', alpha=0.3 )
    ax20.text(0.95,0.4, 'z=%.1f'%(z), color='k',ha='right',va='bottom',fontweight='bold', transform=ax20.transAxes, fontsize='large')
    pl.setp(ax20.get_xticklabels(), visible=False)
    pl.setp(ax20.get_yticklabels(), visible=False)

    z = 2.2
    w1az = w1a * (1+z)
    f1az = f1a / f1a.max() / 2.
    ax22 = pl.subplot(3,2,5, sharex=ax18)
    ax22.plot(w127, f127, ls='-', color='DarkOrchid',label='F127M')
    ax22.plot(w139, f139, ls='-',color='Teal', label='F139M')
    ax22.plot(w153, f153, ls='-',color='Maroon', label='F153M')
    ax22.plot(w1az, f1az, ls='-', lw=0.7, color='0.5', label='_nolegend_')
    ax22.fill_between( w1az, f1az, where=((w1az>12400) & (w1az<13120)), color='DarkOrchid', alpha=0.3 )
    ax22.fill_between( w1az, f1az, where=((w1az>13500) & (w1az<14150)), color='Teal', alpha=0.3 )
    ax22.fill_between( w1az, f1az, where=((w1az>15000) & (w1az<15700)), color='Maroon', alpha=0.3 )
    ax22.text(0.95,0.4, 'z=%.1f'%(z), color='k',ha='right',va='bottom',fontweight='bold', transform=ax22.transAxes, fontsize='large')
    pl.setp(ax22.get_yticklabels(), visible=False)
    ax22.set_xlabel('Observed Wavelength [\AA]')

    ax18.set_xlim( 9000, 19900 )

    # ------------------------------------------------------------
    # Now the color-color plots on the right

    if type(simdata)==str : simdata = readSimData( simdata ) 
    simIa,simIbc,simII = simdata

    if showsn : 
        sn = stardust.SuperNova(_medbanddat)
        H,dH = sn.MAG[ np.where(sn.FLT=='H')[0][0] ], sn.MAGERR[ np.where(sn.FLT=='H')[0][0] ]
        J,dJ = sn.MAG[ np.where(sn.FLT=='J')[0][0] ], sn.MAGERR[ np.where(sn.FLT=='J')[0][0] ]
        N,dN = sn.MAG[ np.where(sn.FLT=='N')[0][0] ], sn.MAGERR[ np.where(sn.FLT=='N')[0][0] ]
        Y,dY = sn.MAG[ np.where(sn.FLT=='Y')[0][0] ], sn.MAGERR[ np.where(sn.FLT=='Y')[0][0] ]
        L,dL = sn.MAG[ np.where(sn.FLT=='L')[0][0] ], sn.MAGERR[ np.where(sn.FLT=='L')[0][0] ]
        O,dO = sn.MAG[ np.where(sn.FLT=='O')[0][0] ], sn.MAGERR[ np.where(sn.FLT=='O')[0][0] ]
        P,dP = sn.MAG[ np.where(sn.FLT=='P')[0][0] ], sn.MAGERR[ np.where(sn.FLT=='P')[0][0] ]
        Q,dQ = sn.MAG[ np.where(sn.FLT=='Q')[0][0] ], sn.MAGERR[ np.where(sn.FLT=='Q')[0][0] ]

    ax1 = pl.subplot(2,2,2)
    stardust.simplot.plotColorColor( simII,  'Q-H','O-J', plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
    stardust.simplot.plotColorColor( simIbc, 'Q-H','O-J', plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
    stardust.simplot.plotColorColor( simIa,  'Q-H','O-J', plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
    if showsn : 
        ax1.errorbar( Q-H, O-J, np.sqrt(dO**2+dJ**2), np.sqrt(dQ**2+dH**2), color='k',mfc='w',mec='k',mew=2,elinewidth=2,marker='D'  )

    ax1.xaxis.set_ticks_position('top')
    ax1.xaxis.set_ticks_position('both')
    ax1.xaxis.set_label_position('top')
    ax1.yaxis.set_ticks_position('right')
    ax1.yaxis.set_ticks_position('both')
    ax1.yaxis.set_label_position('right')

    ax2 = pl.subplot(2,2,4, sharex=ax1)
    stardust.simplot.plotColorColor( simII,  'Q-H','P-N', plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
    stardust.simplot.plotColorColor( simIbc, 'Q-H','P-N', plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
    stardust.simplot.plotColorColor( simIa,  'Q-H','P-N', plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False )
    if showsn : 
        ax2.errorbar( Q-H, P-N, np.sqrt(dP**2+dN**2), np.sqrt(dQ**2+dH**2), color='k',mfc='w',mec='k',mew=2,elinewidth=2,marker='D'  )
    ax2.yaxis.set_ticks_position('right')
    ax2.yaxis.set_ticks_position('both')
    ax2.yaxis.set_label_position('right')
    ax2.text( 0.2, -0.18, 'Ia', ha='left', va='top', color='DarkRed', fontsize='large', fontweight='bold' )
    ax2.text( -0.2, 0.25, 'Ib/c', ha='left', va='top', color='DarkGreen', fontsize='large', fontweight='bold' )
    ax2.text( -0.05, 0.21, 'II', ha='left', va='bottom', color='DarkBlue', fontsize='large', fontweight='bold' )
    ax2.set_ylabel(ax2.get_ylabel(),rotation=-90, color='teal' )
    ax2.set_xlabel(ax2.get_xlabel(), color='maroon' )
    ax1.set_ylabel(ax1.get_ylabel(),rotation=-90,color='DarkOrchid' )
    ax1.set_xlabel(ax1.get_xlabel(), color='maroon' )

    fig = pl.gcf()
    fig.subplots_adjust( left=0.03, right=0.90, top=0.90, bottom=0.12, hspace=0, wspace=0.08 )

    ax1.set_xlim(-0.28,0.38)
    ax1.set_ylim(-0.54,0.2)
    ax2.set_ylim(-0.3,0.32)
    ax1.text( 0.25, -0.19, 'z=1.8', color='k',ha='left',va='bottom',fontweight='bold', fontsize='large')
    ax1.text( -0.15, -0.32, 'z=2.0', color='k',ha='right',va='bottom',fontweight='bold', fontsize='large')
    ax1.text( 0.29, -0.51, 'z=2.2', color='k',ha='right',va='bottom',fontweight='bold', fontsize='large')
    ax1.plot( [0.1,0.24],[-0.2,-0.175], ls='-', marker=' ', color='k', lw=0.8 )
    ax1.plot( [-0.14,-0.07],[-0.28,-0.185], ls='-', marker=' ', color='k', lw=0.8 )

    ax2.text( 0.2, 0.21,  '1.8', color='k',ha='left',va='bottom',fontweight='bold', fontsize='large')
    ax2.text( -0.18, -0.27, '2.0', color='k',ha='right',va='bottom',fontweight='bold', fontsize='large')
    ax2.text( 0.3, 0.01,  '2.2', color='k',ha='left',va='bottom',fontweight='bold', fontsize='large')


def mkCirclePlots( simdata, linelevels = [ 0, 0.82 ], plotstyle='contourf',
                   mjd=55500, Nbins=80 ):
    """ construct circle diagrams for the redshift range of interest. """
    if type(simdata)==str : simdata = readSimData( simdata ) 
    simIa,simIbc,simII = simdata

    ax1 = pl.subplot(2,2,1)
    stardust.simplot.plotColorColor( simII,  'L-Y','O-J', plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False, mjdrange=[mjd,mjd] )
    stardust.simplot.plotColorColor( simIbc, 'L-Y','O-J', plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False, mjdrange=[mjd,mjd] )
    stardust.simplot.plotColorColor( simIa,  'L-Y','O-J', plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False, mjdrange=[mjd,mjd] )

    ax1.xaxis.set_ticks_position('top')
    ax1.xaxis.set_ticks_position('both')
    ax1.xaxis.set_label_position('top')
    ax1.yaxis.set_ticks_position('left')
    ax1.yaxis.set_ticks_position('both')
    ax1.yaxis.set_label_position('left')

    ax2 = pl.subplot(2,2,2)
    stardust.simplot.plotColorColor( simII,  'P-N','O-J', plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False, mjdrange=[mjd,mjd] )
    stardust.simplot.plotColorColor( simIbc, 'P-N','O-J', plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False, mjdrange=[mjd,mjd] )
    stardust.simplot.plotColorColor( simIa,  'P-N','O-J', plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False, mjdrange=[mjd,mjd] )
    ax2.xaxis.set_ticks_position('top')
    ax2.xaxis.set_ticks_position('both')
    ax2.xaxis.set_label_position('top')
    ax2.yaxis.set_ticks_position('right')
    ax2.yaxis.set_ticks_position('both')
    ax2.yaxis.set_label_position('right')

    ax3 = pl.subplot(2,2,3)
    stardust.simplot.plotColorColor( simII,  'L-Y','Q-H', plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False, mjdrange=[mjd,mjd] )
    stardust.simplot.plotColorColor( simIbc, 'L-Y','Q-H', plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False, mjdrange=[mjd,mjd] )
    stardust.simplot.plotColorColor( simIa,  'L-Y','Q-H', plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False, mjdrange=[mjd,mjd] )
    ax3.yaxis.set_ticks_position('left')
    ax3.yaxis.set_ticks_position('both')
    ax3.yaxis.set_label_position('left')

    ax4 = pl.subplot(2,2,4, sharex=ax1)
    stardust.simplot.plotColorColor( simII,  'P-N','Q-H', plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False, mjdrange=[mjd,mjd]  )
    stardust.simplot.plotColorColor( simIbc, 'P-N','Q-H', plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False, mjdrange=[mjd,mjd]  )
    stardust.simplot.plotColorColor( simIa,  'P-N','Q-H', plotstyle=plotstyle, linelevels=linelevels, Nbins=Nbins, sidehist=False, mjdrange=[mjd,mjd]  )
    ax4.yaxis.set_ticks_position('right')
    ax4.yaxis.set_ticks_position('both')
    ax4.yaxis.set_label_position('right')

    if False : 
        ax4.text( 0.2, -0.18, 'Ia', ha='left', va='top', color='DarkRed', fontsize='large', fontweight='bold' )
        ax4.text( -0.2, 0.25, 'Ib/c', ha='left', va='top', color='DarkGreen', fontsize='large', fontweight='bold' )
        ax4.text( -0.05, 0.21, 'II', ha='left', va='bottom', color='DarkBlue', fontsize='large', fontweight='bold' )

    ax2.set_ylabel(ax2.get_ylabel(),rotation=-90, color='teal' )
    ax2.set_xlabel(ax2.get_xlabel(), color='maroon' )
    ax1.set_ylabel(ax1.get_ylabel(),rotation=-90,color='DarkOrchid' )
    ax1.set_xlabel(ax1.get_xlabel(), color='maroon' )

    fig = pl.gcf()
    fig.subplots_adjust( left=0.12, right=0.88, top=0.88, bottom=0.12, hspace=0, wspace=0.00 )

    if False : 
        ax1.set_xlim(-0.28,0.38)
        ax1.set_ylim(-0.54,0.2)
        ax2.set_ylim(-0.3,0.32)
        ax1.text( 0.25, -0.19, 'z=1.8', color='k',ha='left',va='bottom',fontweight='bold', fontsize='large')
        ax1.text( -0.15, -0.32, 'z=2.0', color='k',ha='right',va='bottom',fontweight='bold', fontsize='large')
        ax1.text( 0.29, -0.51, 'z=2.2', color='k',ha='right',va='bottom',fontweight='bold', fontsize='large')
        ax1.plot( [0.1,0.24],[-0.2,-0.175], ls='-', marker=' ', color='k', lw=0.8 )
        ax1.plot( [-0.14,-0.07],[-0.28,-0.185], ls='-', marker=' ', color='k', lw=0.8 )

        ax2.text( 0.2, 0.21,  '1.8', color='k',ha='left',va='bottom',fontweight='bold', fontsize='large')
        ax2.text( -0.18, -0.27, '2.0', color='k',ha='right',va='bottom',fontweight='bold', fontsize='large')
        ax2.text( 0.3, 0.01,  '2.2', color='k',ha='left',va='bottom',fontweight='bold', fontsize='large')




#------------------------------------------------------------

def plotMedBandFilters( z = 2, day=5 ):
    from hstsntools import stardust
    w1a, f1a = stardust.snsed.getsed( sedfile='/usr/local/SNDATA_ROOT/snsed/Hsiao07.dat', day=day )

    w1az = w1a * (1+z)
    f1az = f1a / f1a.max() / 2.
    pl.clf()

    ax1 = subplot(3,1,1)
    plot(w125, f125, 'b--', label='F125W')
    plot(w127, f127, 'b-', label='F127M')
    plot(w1az, f1az, 'r-', label='_nolegend_')
    ax1.legend( loc='upper right', frameon=False, numpoints=2, handlelen=0.2, labelspacing=0.1 )
    ax1.set_xlim( 9000, 20000 )
    ax1.text(9500,0.2, 'SNIa\nz=%.1f\nt=%i'%(z,day), color='r',ha='left',va='bottom')
    setp(ax1.get_xticklabels(), visible=False)
    setp(ax1.get_yticklabels(), visible=False)


    ax2 = subplot(3,1,2, sharex=ax1, sharey=ax1)
    plot(w140, f140, 'g--',label='F140W')
    plot(w139, f139, 'g-',label='F139M')
    plot(w1az, f1az, 'r-', label='_nolegend_')
    ax2.legend( loc='upper right', frameon=False, numpoints=2, handlelen=0.2, labelspacing=0.1 )
    ax2.set_xlim( 9000, 20000 )
    setp(ax2.get_xticklabels(), visible=False)
    setp(ax2.get_yticklabels(), visible=False)
    ax2.set_ylabel('Flux / Transmission (arbitrary units)')

    ax3= subplot(3,1,3, sharex=ax1, sharey=ax1)
    plot(w160, f160, 'm--',label='F160W')
    plot(w153, f153, 'm-',label='F153M')
    plot(w1az, f1az, 'r-',label='_nolegend_')
    ax3.legend( loc='upper right', frameon=False, numpoints=2, handlelen=0.2, labelspacing=0.1 )
    setp(ax3.get_yticklabels(), visible=False)

    ax1.set_xlim( 9000, 20000 )
    ax1.set_xlabel('observed wavelength (Angstroms)')

    fig = gcf()
    fig.subplots_adjust( wspace=0, hspace=0, left=0.05, bottom=0.15, right=0.95, top=0.95)



def plotbroadbandz(  zvals=[1,1.5,2.0], day=0 ):
    """ show how broad bands cover the SED at high z"""
    from hstsnpipe import tools
    from tools import snana
    #w1a, f1a = snana.snsed.getsed( sedfile='/usr/local/SNDATA_ROOT/snsed/Hsiao07.extrap.dat', day=day )
    print("SALT2")
    # w1a, f1a = snana.snsed.getsed( sedfile='/usr/local/SNDATA_ROOT/models/SALT2/SALT2.Guy10_UV2IR/salt2_template_0.dat', day=day )
    w1a, f1a = snana.snsed.getsed( sedfile='/usr/local/SNDATA_ROOT/models/SALT2/SALT2.Guy10_UV2IR/salt2_template_1.dat', day=day )
    #wII, fII = snana.snsed.getsed( sedfile='/usr/local/SNDATA_ROOT/snsed/non1a/SDSS-000018.DAT', day=0 )
    #wIb, fIb = snana.snsed.getsed( sedfile='/usr/local/SNDATA_ROOT/snsed/non1a/SDSS-000020.DAT', day=0 )

    clf()

    i = 0
    for z in zvals: 
        i+=1
        w1az = w1a * (1+z)
        f1az = f1a / f1a.max() / 2.

        #wII = wII * (1+z)
        #fII = fII / fII.max() / 2.

        #wIb = wIb * (1+z)
        #fIb = fIb / fIb.max() / 2.

        ax = subplot(3,1,i)
        plot(w350, f350, 'b--', label='F350LP(W)')
        plot(w125, f125, 'g--', label='F125W(J)')
        plot(w160, f160, 'r--', label='F160W(H)')

        plot(w1az, f1az, 'k-', label='_nolegend_')
        #ax.legend( loc='upper right', frameon=False, numpoints=2, handlelen=0.2, labelspacing=0.1 )
        ax.set_xlim( 3000, 20000 )
        ax.text(0.98,0.95, 'z=%.1f'%(z), color='k',ha='right',va='top',transform=ax.transAxes)
        setp(ax.get_yticklabels(), visible=False)

        if i==1 : 
            top = ax.get_ylim()[1]
            ax.text(16000,top, 'F160W(H)', color='r',ha='center',va='bottom')
            ax.text(12500,top, 'F125W(J)', color='g',ha='center',va='bottom')
            ax.text(3500,top, 'F350LP(W)', color='b',ha='left',va='bottom')
        if i<3 : 
            setp(ax.get_xticklabels(), visible=False)
        if i==2 : 
            ax.set_ylabel('Flux or Transmission (arbitrary units)')
        if i==3 : 
            ax.set_xlabel('observed wavelength (Angstroms)')

    fig = gcf()
    fig.subplots_adjust( wspace=0, hspace=0, left=0.05, bottom=0.12, right=0.95, top=0.95)
