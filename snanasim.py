# 2013.02.12
# S.Rodney
"""
Figure to demonstrate the use of medium band filters 
for efficient SN classification and redshift estimation
"""
import numpy as np
import os
import stardust

sndataroot = os.environ['SNDATA_ROOT']
topdir = os.path.abspath( '.' )

_medbanddat = 'ffsn_exampleIa_medband.dat'


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

    return readSimDataMC( simname )

def readSimDataMC(simname='simFF_medband'):
    simIa = stardust.SimTable( simname+'_Ia')
    simIbc = stardust.SimTable( simname+'_Ibc')
    simII = stardust.SimTable( simname+'_II')
    return [simIa, simIbc, simII]



def dosimGrid( sn, ngridz=50, bands='X7I8LYOJPNQH',
               x1range=[-2,2], crange=[-0.2,0.5],
               avrange=[0,0.7], trestrange = [-5,5],
               clobber=0, verbose=1 ) :
    """ set up and run a SNANA Grid simulation to show the
    medium band technique.
    """
    import stardust
    simname = 'sim_'+sn.nickname+'_medbandGrid'

    if 'ClassSim' in sn.__dict__ :
        # we allow a generous range of possible redshifts, up to dz=0.2 or 10-sigma
        zrange = [ sn.ClassSim.Ia.z_maxprob-min(0.2,10*sn.ClassSim.Ia.z_maxprob_err),
                   sn.ClassSim.Ia.z_maxprob+min(0.2,10*sn.ClassSim.Ia.z_maxprob_err) ]
    else :
        zrange = [ sn.z - sn.zerr, sn.z + sn.zerr ]

    stardust.simulate.mkgridsimlib( simname+'.simlib', survey='HST',
                                    field='default', bands=bands,
                                    clobber=clobber )

    stardust.simulate.mkinputGrid(
        simname+'_Ia', inputfile=simname+'_Ia.input', simlibfile=simname+'.simlib',
        simtype='Ia', clobber=clobber,
        GENFILTERS=bands,
        NGRID_TREST=10,    GENRANGE_TREST=trestrange,
        NGRID_LOGZ=ngridz,     GENRANGE_REDSHIFT=zrange,
        NGRID_LUMIPAR=10,  GENRANGE_SALT2X1=x1range,
        NGRID_COLORPAR=10, GENRANGE_SALT2C=crange,
        NGRID_COLORLAW=1,  GENRANGE_RV=[3.1,3.1] )

    # stardust.simulate.mkinput( simname+'_Ibc', inputfile=simname+'_Ibc.input',
    #                                simlibfile=simname+'.simlib',
    #                                simtype='Ibc', ratemodel='flat',
    #                                clobber=clobber, GENSOURCE='GRID',
    #                                GENFILTERS=bands,
    #                                GENPERFECT=2, GENRANGE_REDSHIFT=zrange,
    #                                 smear=False,
    #                                GENRANGE_AV=avrange, )

    # stardust.simulate.mkinput( simname+'_II', inputfile=simname+'_II.input',
    #                            simlibfile=simname+'.simlib',
    #                            simtype='II', rate='flat',
    #                            clobber=clobber,
    #                            GENFILTERS=bands, GENSOURCE='GRID',
    #                            GENPERFECT=2, GENRANGE_REDSHIFT=zrange,
    #                            smear=False,
    #                            GENRANGE_AV=avrange, )

    stardust.simulate.dosim('%s_Ia.input'%simname, verbose=verbose )
    # stardust.simulate.dosim('%s_Ibc.input'%simname, perfect=True, verbose=verbose )
    # stardust.simulate.dosim('%s_II.input'%simname, perfect=True, verbose=verbose )

    simdat = stardust.SimTable( '%s_Ia'%simname)
    return( simdat )





