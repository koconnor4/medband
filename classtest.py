
# Probability that a SN of a given type (Ia,Ibc,II) belongs to any
# given sub-class (Ia,Ib,Ic,IIP,IIL,IIn).   from Li et al 2011a
SubClassProbs = {
           'ia':{'Ia':1.0 },
           # For simulations containing only Type II sub-types
           'ii':{'IIP':0.7,'IIn':0.3 },
           # For simulations containing only Type Ib/c sub-types
           'ibc':{'Ib':0.46,'Ic':0.54},
           # For simulations containing all CC sub-types together
           'cc':{'Ib':0.19/0.76*0.46,'Ic':0.19/0.76*0.54,
                 'IIP':0.57/0.76*0.7,'IIn':0.57/0.76*0.3 },
           }

class SncosmoSim( object ):
    """ A class for holding the results of sncosmo simulations,
    including parameters, cosmology, and light curves.
    """
    def __init__(self, sntype, observations=None, z_range=[1.8,2.2],
                 t0_range=[0,0], nsim=100, perfect=True,
                 Om=0.3, H0=70, templateset='NOTPSNID' ):
        """ Run a monte carlo sim using sncosmo to simulate <nsim> SNe
        of the given <sntype> over the given <z_range>.

        Simulates Type Ia SNe with the SALT2 model, and CC SNe with
        the SNANA CC templates.

        Set perfect=True for noiseless "observations" of the simulated SNe.
        :return:
        """
        from astropy import cosmology
        import sncosmo
        from numpy.random import normal, uniform, choice
        import numpy as np
        # get dictionaries of sncosmo CCSN model names and their corresponding SN sub-type
        from .classify import SubClassDict_NOTPSNID, SubClassDict_PSNID, SubClassDict_SNANA

        self.sntype = sntype
        self.z_range = z_range
        self.nsim = nsim
        self.perfect = perfect

        if observations is None :
            observations = mk_obstable( )
        self.observations = observations

        if templateset.lower()=='psnid':
            SubClassDict = SubClassDict_PSNID[sntype.lower()]
        elif templateset.lower()=='snana':
            SubClassDict = SubClassDict_SNANA[sntype.lower()]
        elif templateset.lower()=='notpsnid':
            SubClassDict = SubClassDict_NOTPSNID[sntype.lower()]

        # Make a list of all the unique sncosmo source models available,
        # and assign a relative probability that any given simulated SN of this
        # primary type (Ia or Ibc or II) belongs to a given model.
        #
        # e.g. From Li+ 2011 we have that 46% of all Type Ib/c SN are in the Ib
        # subclass.  So if there are 10 type Ib models in our set, then each
        # simulated Ibc supernova has a 4.6% chance of being assigned to each
        # of the 10 Ib models.
        subClassProbs = SubClassProbs[sntype.lower()]
        self.SourcenameSet = np.array( SubClassDict.keys() )
        self.SubclassSet = np.array([ SubClassDict[source] for source in self.SourcenameSet ])
        self.SubclassCount = np.array([ len(np.where(self.SubclassSet==subclass)[0])
                                        for subclass in self.SubclassSet ], dtype=float)
        self.SourceprobSet = np.array([ subClassProbs[subclass]
                                        for subclass in self.SubclassSet ]) / self.SubclassCount
        self.SourceprobSet /= self.SourceprobSet.sum()

        # load the O'Donnell 1994 dust model
        self.dust = sncosmo.OD94Dust()

        # Define an sncosmo SN model for each available source
        modelset = np.array([ sncosmo.Model(source=source, effects=[self.dust],
                                    effect_names=['host'], effect_frames=['rest'])
                                 for source in self.SourcenameSet ])
        # Define a cosmology
        self.cosmo = cosmology.FlatLambdaCDM(Om0=Om, H0=H0)

        # For each simulated SN, draw random Av from distributions
        # as defined in Rodney et al 2014a :
        #   For SN Ia :  P(Av) = exp(-Av/0.33)
        #   For CC SN :  P(Av) = 4 * gauss(0.6) + exp(-Rv/1.7)
        if sntype=='Ia':
            tau,sigma,R0 = 0.33, 0, 0
        else :
            tau,sigma,R0 = 1.7, 0.6, 4
        self.Av = mcsample( pAv, nsim, tau=tau, sigma=sigma, R0=R0 )

        # For each simulated SN, draw a random Rv from a normal
        # distribution centered on 3.1 with width 0.5
        self.Rv = normal( 3.1, 0.5, nsim )
        self.Rv = np.where( self.Rv>0, self.Rv, 0.01 )

        # Convert Av and Rv to E(B-V) :
        # Rv = Av/EBV ==> EBV = Av/Rv
        self.EBV = self.Av / self.Rv

        # TODO : draw the redshifts with a metropolis-hastings sampler to match
        #  a distribution defined based on the expected SN rate

        # Disabled : uniform redshift spacing
        # zlist = np.linspace( z_range[0], z_range[1], nsim )

        # Draw a random redshift from a uniform distribution
        self.z = uniform( low=z_range[0], high=z_range[1], size=nsim )

        lightcurvelist = []
        peakabsmagRlist = []
        modelparamlist = []
        subclasslist = []
        modelindexlist = []
        sourcenamelist = []
        t0list = []
        if sntype=='Ia':
            x0list = []
            x1list = []
            clist = []
        else :
            amplitudelist = []
        for isim in range(self.nsim):
            # Randomly draw an sncosmo model from the available list, according to
            # the predefined probability list, setting the SN sub-class for this
            # simulated SN
            imodel = choice( np.arange(len(modelset)), replace=True, p=self.SourceprobSet )
            model =  modelset[imodel]
            subclass = self.SubclassSet[imodel]

            z = self.z[isim]
            EBV = self.EBV[isim]
            Rv = self.Rv[isim]

            # Set the peak absolute magnitude according to the observed
            # luminosity functions, as defined in Table 3 of Graur:2014a;
            # and set the host extinction according to the 'mid' dust model
            # of Rodney:2014a.
            if subclass == 'Ia' :
                MR = normal( -19.37, 0.47 )
            elif subclass == 'Ib' :
                MR = normal( -17.90, 0.90 )
            elif subclass == 'Ic' :
                MR = normal( -18.30, 0.60 )
            elif subclass == 'IIP' :
                MR = normal( -16.56, 0.80 )
            elif subclass == 'IIL' :
                MR = normal( -17.66, 0.42 )
            elif subclass == 'IIn' :
                MR = normal( -18.25, 1.00 )
            model.set(z=z)
            model.set_source_peakabsmag( MR, 'bessellr', 'vega', cosmo=self.cosmo)

            modelindexlist.append( imodel )
            subclasslist.append( subclass )
            peakabsmagRlist.append( MR )
            sourcenamelist.append( self.SourcenameSet[imodel] )
            if subclass =='Ia' :
                x0 = model.get('x0')
                # TODO : use bifurcated gaussians for more realistic x1,c dist'ns
                x1 = normal(0., 1.)
                c = normal(0., 0.1)
                t0 = uniform( t0_range[0], t0_range[1] )
                modelparams = {'z':z, 't0':t0, 'x0':x0, 'x1':x1, 'c':c, 'hostebv':EBV, 'hostr_v':Rv}
                t0list.append( t0 )
                x0list.append( x0 )
                x1list.append( x1 )
                clist.append( c )
                t0list.append( t0 )
            else :
                amplitude = model.get('amplitude')
                t0 = uniform( t0_range[0], t0_range[1] )
                modelparams = {'z':z, 't0':t0, 'amplitude':amplitude, 'hostebv':EBV, 'hostr_v':Rv }
                amplitudelist.append( amplitude )
                t0list.append( t0 )
            modelparamlist.append( modelparams )

            # Generate one simulated SN:
            snlc = sncosmo.realize_lcs(self.observations, model, [ modelparams ],
                                       thresh=None, perfect=perfect )
            lightcurvelist.append( snlc[0] )

        self.lightcurves = lightcurvelist
        self.t0 = np.array( t0list )
        self.modelindex = np.array( modelindexlist )
        self.sourcename = np.array( sourcenamelist )
        self.subclass = np.array( subclasslist )
        self.modelparam = np.array( modelparamlist )
        self.peakabsmagR = np.array( peakabsmagRlist )

        if sntype=='Ia':
            self.x0 = np.array( x0list )
            self.x1 = np.array( x1list )
            self.c  = np.array( clist )
        else :
            self.amplitude = np.array( amplitudelist )

        return




from astropy.io import ascii

testsnIadat = """
# time band     flux         fluxerr       zp  zpsys
 0.0 f127m 0.491947902265  0.017418231547 24.6412    ab
 0.0 f139m 0.513425670819 0.0168000764011 24.4793    ab
 0.0 f153m 0.486808758939 0.0167684488219 24.4635    ab
 0.0 f125w  2.14010106322 0.0649063974142   26.25    ab
 0.0 f140w  2.78151131439 0.0722039093523   26.46    ab
 0.0 f160w   1.6716457987 0.0594698101517   25.96    ab
"""

testsnCCdat = """
#time  band      flux         fluxerr       zp  zpsys
 0.0 f127m 0.9359 0.9674 26.47    ab
 0.0 f139m 0.8960 0.9466 26.49    ab
 0.0 f153m  1.004  1.002  26.7    ab
 0.0 f125w  3.937  1.984 28.02    ab
 0.0 f140w  5.606  2.367 28.48    ab
 0.0 f160w  3.978  1.994 28.19    ab
"""

testsnIa = ascii.read( testsnIadat )
testsnCC = ascii.read( testsnCCdat )


def do_classification_test( sntypelist=['Ia','Ibc','II'],
        outfileroot='classtest',
        # Parameters for making the obs table, giving survey parameters
        nepochs=10, cadence=15, medbands=True,
        bandlist=['f125w','f160w','f814w'], exptimelist=[2500.,2500.,2500.],
        # Parameters for the simulated SN sample
        nsim=1000, z_range=[1.5,2.5], t0_range=[10,50], fixSNR=False,
        # Parameters for the nested sampling bayesian classifer
        nobj=100, maxiter=10000,
        clobber=False, verbose=True ) :
    """  Simulate a bunch of SNe and then classify them, recording results
    into text files.
    """
    import cPickle
    import os
    import numpy as np
    from copy import deepcopy
    from numpy.random import normal, uniform
    from time import asctime,time
    from . import classify

    t0searchrange = [t0_range[0]*0.8,t0_range[1]*1.1]
    zsearchrange = [z_range[0]*0.8,z_range[1]*1.1]

    # Construct an observation table with the user-specified survey parameters
    obstable = mk_obstable( nepochs=nepochs, cadence=cadence, medbands=medbands,
                            bandlist=bandlist, exptimelist=exptimelist )

    simpicklelist = []
    outfilelist = []
    if isinstance( sntypelist, str ) : sntypelist = [sntypelist]
    for sntype in sntypelist :
        # check for existing output data files
        outfile = outfileroot + '_%s.%i.dat'%(sntype,nsim)
        if os.path.isfile( outfile ) and not clobber :
            print("%s exists. Use clobber to overwrite."%outfile)
            print( "Carrying on to the next sn type")
            continue

        # Load simulated SN from pickles, if possible. Otherwise, run the simulations
        classSimpkl = '%s.classSim%s.%i.pkl'%(outfileroot,sntype,nsim)
        if os.path.isfile( classSimpkl ) and not clobber>1 :
            if verbose: print("Loading %s simulation from pickle : %s"%(sntype,classSimpkl) )
            fin = open( classSimpkl, 'rb' )
            classSim = cPickle.load( fin )
            fin.close()
        else :
            if verbose: print("Running a new %s simulation, then saving to pickle : %s"%(sntype,classSimpkl) )
            classSim = SncosmoSim( sntype , observations=obstable, z_range=z_range,
                                   t0_range=t0_range, nsim=nsim, perfect=fixSNR, Om=0.3, H0=70 )
            fout = open( classSimpkl, 'wb' )
            cPickle.dump( classSim, fout, protocol=-1 )
            fout.close()
        simpicklelist.append( classSim )

        if verbose: print("Writing classification output to %s "%(outfile))

        # erase any existing output files, and start a new version with a header
        fout = open( outfile, 'w')
        print >> fout, '# %s'% asctime()
        print >> fout, '# Classifier Validation Test output generated with %s'% __file__
        print >> fout, '# Classifying %i simulated Type %s SN'%(nsim,sntype)
        print >> fout, '# survey: nepochs=%i ; cadence=%i; medbands=%s; '%( nepochs, cadence, medbands)
        print >> fout, '# survey: bands=%s; exptimes=%s'%(bandlist, exptimelist )
        print >> fout, '# nestlc: nobj=%i; maxiter=%i'%( nobj, maxiter )
        if sntype=='Ia' :
            print >> fout, '#index pIa  pIbc   pII     z     t0         x0     x1     c   zphot zphoterr snrmax zpost zposterr t0post t0posterr'
        else :
            print >> fout, '#index pIa  pIbc   pII     z     t0  amplitude       model     Av     Rv  zphot zphoterr snrmax zpost zposterr t0post t0posterr'
        fout.close()

        # Run the classifier
        for isn in range(nsim) :
            if verbose: print("Classifying Type %s SN %i  of  %i"%(sntype, isn+1, nsim))
            tstart = time()
            snlc = classSim.lightcurves[isn]

            # assign a photometric redshift
            zphoterr = uniform( 0.05, 0.5 )
            zphotmin  = zsearchrange[0]*1.1
            zphotmax = zsearchrange[1]*0.9
            zphot = classSim.z[isn] + normal( 0, zphoterr )
            zphot = max( zphotmin, min( zphot, zphotmax) )

            if fixSNR :
                # Fix the S/N ratio and add appropriate noise to the perfect flux
                snlc = deepcopy( classSim.lightcurves[isn])
                snlc['fluxerr'] = np.abs(snlc['flux']) / fixSNR
                snlc['flux'] = normal( snlc['flux'], snlc['fluxerr'] )

            # Skip if no points with S/N > 5
            pIa=pIbc=pII=0
            zout=dzout=t0out=dt0out=x0out=dx0out=x1out=dx1out=cout=dcout=0
            ampout=dampout=ebvout=debvout=rvout=drvout=0
            snrmax = np.abs( snlc['flux']/snlc['fluxerr'] ).max()
            if snrmax > 5 :
                # Run the classifier and append the result to the running output file.
                classout = classify.classify( snlc, zhost=zphot, zhosterr=zphoterr,
                                              t0_range=t0searchrange, zminmax=zsearchrange,
                                              nobj=nobj, maxiter=maxiter,
                                              templateset='PSNID', verbose=(verbose>1) )
                pIa  = classout['pIa']
                pIbc = classout['pIbc']
                pII  = classout['pII']

                bestmodelres = classout[ classout['bestmodel']]['res']
                bestpdf = classify.get_marginal_pdfs( bestmodelres, nbins=0, verbose=verbose>2 )
                zout,dzout = bestpdf['z'][2],bestpdf['z'][3]
                t0out,dt0out = bestpdf['t0'][2],bestpdf['t0'][3]
                # if sntype=='Ia':
                #     x0out,dx0out = bestpdf['x0'][2],bestpdf['x0'][3]
                #     x1out,dx1out = bestpdf['x1'][2],bestpdf['x1'][3]
                #     cout, dcout = bestpdf['c'][2],bestpdf['c'][3]
                # else :
                #     ampout,dampout = bestpdf['amplitude'][2],bestpdf['amplitude'][3]
                #     ebvout,debvout = bestpdf['hostebv'][2],bestpdf['hostebv'][3]
                #     rvout, drvout = bestpdf['hostr_v'][2],bestpdf['hostr_v'][3]

            if sntype == 'Ia' :
                outline = '%4i  %4.2f  %4.2f  %4.2f  %4.2f  %4.2f  %9.3e  %5.2f    %5.2f  %5.2f  %5.2f  %5.2f  %5.2f %5.2f   %7.1f %4.1f '% (
                        isn, pIa, pIbc, pII, classSim.z[isn],  classSim.t0[isn],
                        classSim.x0[isn], classSim.x1[isn], classSim.c[isn], zphot, zphoterr, snrmax,
                        zout, dzout, t0out, dt0out )#, x0out,dx0out, x1out,dx1out, cout, dcout )
            else :
                outline = '%4i  %4.2f  %4.2f  %4.2f  %4.2f  %4.2f  %9.3e  %10s  %5.2f  %5.2f  %5.2f  %5.2f    %5.2f %5.2f %5.2f   %7.1f %4.1f '% (
                    isn, pIa, pIbc, pII, classSim.z[isn],  classSim.t0[isn],
                    classSim.amplitude[isn], classSim.sourcename[isn],
                    classSim.Av[isn], classSim.Rv[isn],
                    zphot, zphoterr, snrmax,
                    zout, dzout, t0out, dt0out)#, ampout,dampout, ebvout, debvout, rvout, drvout )

            fout = open( outfile, 'a')
            print >> fout, outline
            fout.close()

            tend = time()
            if verbose>1 :
                print( "%i : %.1f sec"%(isn,tend-tstart))
                print("-------------------------")
        outfilelist.append( outfile )
        if verbose: print("Done classifying %i Type %s SN."%(nsim,sntype))

    return( simpicklelist, outfilelist )

def plotSNRtest01( datfileIa='colorColorClassify_Ia.dat',
                   datfileCC='colorColorClassify_CC.dat' ) :
    """  Plot the results of a classification validation test.
    :return:
    """
    from matplotlib import pyplot as pl
    from pytools import plotsetup
    import numpy as np
    from astropy.io import ascii

    datIa = ascii.read( datfileIa, format='commented_header', header_start=-1, data_start=0)
    datCC = ascii.read( datfileCC, format='commented_header', header_start=-1, data_start=0)

    pIaIa = datIa['P(Ia|D)']
    zIa = datIa['z']
    SNRIa = datIa['S/N']

    pIaCC = datCC['P(Ia|D)']
    zCC = datCC['z']
    SNRCC = datCC['S/N']

    fig1 = plotsetup.fullpaperfig( 1, figsize=[8,3])
    SNRlist = np.unique( SNRIa )
    ncol = len(SNRlist)
    icol=0
    axlist = []
    for SNR in SNRlist :
        icol+=1
        if icol == 1 :
            ax1 = pl.subplot(1,ncol,icol)
            ax = ax1
        else :
            ax = pl.subplot(1,ncol,icol,sharex=ax1, sharey=ax1)
        axlist.append( ax )
        out1 = ax.hist( pIaIa[np.where( SNRIa==SNR)[0]],
                        bins=np.arange(0,1.01,0.05),
                        color='darkmagenta', alpha=0.3 )
        out2 = ax.hist( pIaCC[np.where( SNRCC==SNR)[0]],
                        bins=np.arange(0,1.01,0.05),
                        color='teal', alpha=0.3 )

    ax1.set_xlim( 0.0, 1.05 )
    ax1.set_xticks( [0.2,0.4,0.6,0.8,1.0])

    ax1 = axlist[0]
    imid = int(round(len(axlist)/2))
    ax2 = axlist[imid]
    ax3 = axlist[-1]
    ax2.set_xlabel('P(Ia$|$D)')
    ax1.set_ylabel('Number of Test SN')

    pl.setp( ax2.get_yticklabels(),visible=False )
    pl.setp( ax2.get_ylabel(),visible=False )
    ax3.yaxis.set_ticks_position('right')
    ax3.yaxis.set_ticks_position('both')
    ax3.yaxis.set_label_position('right')
    ax1.set_title('S/N=%i'%SNRlist[0],color='0.3')
    ax2.set_title('S/N=%i'%SNRlist[imid],color='0.3')
    ax3.set_title('S/N=%i'%SNRlist[-1],color='0.3')
    ax3.text( 0.1,0.95,'CCSN\nTest\nSet',color='teal',ha='left',va='top',transform=ax3.transAxes)
    ax3.text( 0.88,0.95,'SN Ia\nTest\nSet',color='darkmagenta',ha='right',va='top',transform=ax3.transAxes)
    pl.subplots_adjust( wspace=0, hspace=0, left=0.09, bottom=0.18, right=0.92, top=0.9)
    pl.draw()


def plotSNRtest02( datfileIa='colorColorClassify_Ia.dat',
                   datfileCC='colorColorClassify_CC.dat' ) :
    """  Plot the results of a classification validation test.
    :return:
    """
    from matplotlib import pyplot as pl
    from pytools import plotsetup
    import numpy as np
    from astropy.io import ascii

    datIa = ascii.read( datfileIa, format='commented_header', header_start=-1, data_start=0)
    datCC = ascii.read( datfileCC, format='commented_header', header_start=-1, data_start=0)

    pIaIa = datIa['P(Ia|D)']
    zIa = datIa['z']
    SNRIa = datIa['S/N']

    pIaCC = datCC['P(Ia|D)']
    zCC = datCC['z']
    SNRCC = datCC['S/N']

    fig2 = plotsetup.fullpaperfig( 2, figsize=[8,3])

    def purity( nTrue, nFalse, W=3 ):
        return( float(nTrue) / ( float(nTrue) + float(W*nFalse) ) )

    def efficiency( nTrue, nTot ):
        return( float(nTrue) / float(nTot) )

    SNRlist = np.unique( SNRIa )
    ncol = len(SNRlist)
    icol=0
    axlist = []
    for SNR in SNRlist :
        icol+=1
        irow=-1
        iSNRIa = np.where( SNRIa==SNR )[0]
        iSNRCC = np.where( SNRCC==SNR )[0]
        zIaSNR = zIa[iSNRIa]
        zCCSNR = zCC[iSNRCC]
        for threshold in [0.5,0.75] :
            irow+=1
            if icol == 1 :
                ax1 = pl.subplot(2,ncol,icol+irow*ncol)
                ax = ax1
            else :
                ax = pl.subplot(2,ncol,icol+irow*ncol,sharex=ax1, sharey=ax1)
            axlist.append( ax )

            nIaTrue = np.count_nonzero( pIaIa[iSNRIa] > threshold )
            nIaFalse = np.count_nonzero( pIaCC[iSNRCC] > threshold )
            nIaTot = len( iSNRIa )
            ppIa = purity( nIaTrue, nIaFalse, W=3 )
            pIa = purity( nIaTrue, nIaFalse, W=1 )
            eIa = efficiency( nIaTrue, nIaTot )
            print( 'SNR=%i, Thresh=%.2f, Purity = %.2f, Pseudo-Purity = %.2f'%(
                SNR, threshold, pIa, ppIa ) )
            zbinmin = np.arange(1.8,2.2,0.05)
            ppIaz, pIaz, eIaz = [], [], []
            zbinmid = []
            for zmin in zbinmin :
                zmax = zmin+0.05
                izthisbinIa = np.where( (zIaSNR>=zmin) & (zIaSNR<zmax) )[0]
                izthisbinCC = np.where( (zCCSNR>=zmin) & (zCCSNR<zmax) )[0]
                if not len(izthisbinCC) : continue
                if not len(izthisbinIa) : continue
                nIaTrue = np.count_nonzero( pIaIa[iSNRIa][izthisbinIa]  > threshold )
                nIaFalse = np.count_nonzero( pIaCC[iSNRCC][izthisbinCC] > threshold )
                nIaTot = len(izthisbinIa)
                if nIaTrue == 0 and nIaFalse==0 : continue
                ppIaz.append( purity( nIaTrue, nIaFalse, W=3 ) )
                pIaz.append( purity( nIaTrue, nIaFalse, W=1 ) )
                eIaz.append( efficiency( nIaTrue, nIaTot ) )
                zbinmid.append( zmin+0.025 )
            ax.plot( zbinmid, ppIaz, color='darkorange', marker='d', ls='--', lw=1 )
            ax.plot( zbinmid, pIaz, color='darkcyan', marker='o', ls='-', lw=1 )
            ax.plot( zbinmid, eIaz, color='0.5', marker='^', ls=':', lw=1 )
            ax.set_xlim( 1.75, 2.25 )
            ax.set_ylim( -0.1, 1.1 )
            if icol==2 :
                pl.setp( ax.get_yticklabels(), visible=False )
            if icol==3 :
                ax.yaxis.set_ticks_position('right')
                ax.yaxis.set_ticks_position('both')
                ax.yaxis.set_label_position('right')
            if irow==0 :
                pl.setp( ax.get_xticklabels(), visible=False )
                ax.set_title( 'S/N=%i'%SNR, fontsize=12, color='0.5' )
            #if irow==1 and icol==2 :
            #    ax.set_xlabel('Redshift')
            #if icol==1 :
            #    ax.set_ylabel('Ia Sample Purity')
            if icol==2 :
                ax.text( 0.5, 0.05, 'thresh=%0.2f'%(threshold), ha='center',va='bottom',
                         transform=ax.transAxes, fontsize=12, color='0.5')

            if irow==0 and icol==ncol :
                ax.text( 0.85, 0.05, 'true\npurity', ha='center',va='bottom',color='darkcyan', fontsize=10, transform=ax.transAxes)
                ax.text( 0.5, 0.05, 'pseudo\npurity', ha='center',va='bottom',color='darkorange', fontsize=10, transform=ax.transAxes)
                ax.text( 0.05, 0.05, 'efficiency', ha='left',va='bottom',color='0.3', fontsize=10, transform=ax.transAxes)

    # make a background axis frame for x and y axis labels
    bgAxes = pl.axes([0.05, 0.1, 0.9, 0.85],frameon=False )
    bgAxes.set_xticks([])
    bgAxes.set_yticks([])
    bgAxes.set_xlabel('Redshift')
    bgAxes.set_ylabel('Ia Sample Purity')
    pl.subplots_adjust( wspace=0, hspace=0, left=0.08, bottom=0.18, right=0.92, top=0.9)

    pl.draw()


def testClassify( nsim=100 ):
    from . import classify
    simIa = SncosmoSim( 'Ia', nsim=nsim )
    simIbc = SncosmoSim( 'Ibc', nsim=nsim )
    simII = SncosmoSim( 'II', nsim=nsim )

    pIaArrayLog = []
    pIaArray = []
    for simset in [simIa, simIbc, simII ] :
        for simsn in simset :
            classout = classify.classify( simsn )

    return( )


def mk_obstable( nepochs=10, cadence=15, medbands=True,
                 bandlist=['f125w','f160w','f814w'],
                 exptimelist=[2500.,2500.,2500.]  ):
    """ construct a table of observations for simulating supernovae with sncosmo
    :param exptime:
    :return:
    """

    from simparam import gethstzp, gethstbgnoise
    import numpy as np
    from astropy.table import Table

    epochlist = range(0,nepochs+1)

    camlist = [ 'WFC3IR' if b.lower().startswith('f1') else 'ACSWFC' for b in bandlist ]

    band = np.ravel( [ bandlist for epoch in epochlist ])
    exptime = np.ravel( [ exptimelist for epoch in epochlist ] )
    time = np.ravel( [ [cadence*epoch for b in bandlist] for epoch in epochlist ] )
    zp = np.ravel( [ [gethstzp(b,'ab',c) for b,c in zip(bandlist,camlist)] for epoch in epochlist ] )
    zpsys = np.ravel( [ ['ab' for b in bandlist] for epoch in epochlist ] )

    if medbands :
        # add a single epoch of med band data, 3 filters, with 3 orbits per filter
        medbandlist = ['F127M','F139M', 'F153M']
        band = np.append( band, medbandlist  )
        exptime = np.append( exptime, np.ones(3) * (3*2500) )  # 3 orbits for each med band to get S/N~10 at z=2-2.5
        time = np.append( time, np.ones(3) * int(cadence*nepochs/3) )
        zp = np.append( zp,  [ gethstzp(b) for b in medbandlist ] )
        zpsys = np.append( zpsys, [ 'ab' for b in medbandlist ] )

    gain = exptime
    bgnoise = np.array( [ gethstbgnoise( b.lower(), et, 'WFC3IR' if b.lower().startswith('f1') else 'ACSWFC' ) for b,et in zip(band,exptime) ] )

    observations = Table({'band': band, 'time': time, 'zp': zp,
                          'zpsys': zpsys, 'gain': gain, 'exptime':exptime,
                          'skynoise': bgnoise,
                          })
    return( observations )



def mcsample( p, Ndraws, x0=None, mcsigma=0.05,
              Nburnin=100,  debug=False, *args, **kwargs ) :
    """ Crude metropolis-hastings monte carlo sampling funcion.

    The first argument is a callable function that defines
    the posterior probability at position x:  p(x).

    Positional arguments and optional keyword arguments for the function p
    may be provided at the end.  The function p will be called as
     p(x, *args, **kwargs).

    We construct a Markov Chain with  Ndraws  steps using the
    Metropolis-Hastings algorithm with a gaussian proposal distribution
    of stddev sigma.
    """
    from numpy import random
    if debug: import pdb; pdb.set_trace()

    # if user doesn't provide a starting point,
    # then draw an initial random position between 0 and 1
    if not x0 : x0 = random.uniform()
    xsamples = []
    istep = 0
    p0 = p(x0, *args, **kwargs)
    while len(xsamples) < Ndraws :
        # draw a new position from a Gaussian proposal dist'n
        x1 = random.normal( x0, mcsigma )
        p1 = p( x1, *args, **kwargs )
        # compare new against old position
        if p1>=p0 :
            # new position has higher probability, so
            # accept it unconditionally
            if istep>Nburnin : xsamples.append( x1 )
            p0=p1
            x0=x1
        else :
            # new position has lower probability, so
            # pick new or old based on relative probs.
            y = random.uniform( )
            if y<p1/p0 :
                if istep>Nburnin : xsamples.append( x1 )
                p0=p1
                x0=x1
            else :
                if istep>Nburnin : xsamples.append( x0 )
        istep +=1
    return( xsamples )



# TODO : This is not yet working with emcee.  I've fallen back to using my
#  own metropolis hastings sampler instead.
def mkAvSampler(type='Ia'):
    """ Define a monte carlo sampler to draw Av values from a host galaxy
    Extinction distribution as defined in Rodney:2014a.
    :param type:
    :return:
    """
    import numpy as np
    import emcee

    if type=='Ia' :
        # For SN Ia :
        # P(Av) = exp(-Av/0.33)
        tau = 0.33
        sigma=0
        R0=0
    else :
        # For CC SN :
        # P(Av) = 4 * gauss(0.6) + exp(-Rv/1.7)
        tau = 1.7
        sigma=0.6
        R0 = 4
    def lnprobfn(Av):
        return( np.log(pAv(Av,sigma,tau,R0) ))
    return( emcee.Sampler( 1, lnprobfn, args=[], kwargs={}) )


def pAv( Av, sigma=0, tau=0, R0=0, noNegativeAv=True ):
    """  Dust models:   P(Av)
    :param Av:
    :param sigma:
    :param tau:
    :param R0:
    :param noNegativeAv:
    :return:
    """
    import numpy as np
    if not np.iterable( Av ) : Av = np.array( [Av] )

    # gaussian core
    core = lambda sigma,av : np.exp( -av**2 / (2*sigma**2) )
    # Exponential tail
    tail = lambda tau,av : np.exp( -av/tau )

    if tau!=0 and noNegativeAv:
        tailOut = np.where( Av>=0, tail(tau,Av), 0 )
    elif tau!=0 :
        tailOut = tail(tau,Av)
    else :
        tailOut = np.zeros( len( Av ) )

    if sigma!=0 and noNegativeAv:
        coreOut = np.where( Av>=0, core(sigma,Av), 0 )
    elif sigma!=0 :
        coreOut = core(sigma,Av)
    else :
        coreOut = np.zeros( len( Av ) )

    if len(Av) == 1 :
        coreOut = coreOut[0]
        tailOut = tailOut[0]
    if sigma==0 : return( tailOut )
    elif tau==0 : return( coreOut )
    else : return( R0 * coreOut + tailOut )




def plotColorz( snsimset ):
    """ plot the med band pseudo-colors vs z
    :param snsimset:
    :return:
    """
    import numpy as np
    from matplotlib import pyplot as pl

    z = np.array([ sn.meta['z'] for sn in snsimset ])
    mag = np.array([ -2.5*np.log10( sn['flux'] ) + sn['zp'] for sn in snsimset ])
    # magerr = np.array( [2.5 * np.log10(np.e) *  sn['fluxerr']/ sn['flux'] for sn in sne])

    flt = np.array( [sn['band'] for sn in snsimset] )

    i127 = np.where( flt=='wfc3f127m' )
    i139 = np.where( flt=='wfc3f139m' )
    i153 = np.where( flt=='wfc3f153m' )
    i125 = np.where( flt=='wfc3f125w' )
    i140 = np.where( flt=='wfc3f140w' )
    i160 = np.where( flt=='wfc3f160w' )

    pl.clf()
    pl.plot( z, mag[i127]-mag[i125], marker='o', alpha=0.3, color='darkmagenta', ls=' ' )
    pl.plot( z, mag[i139]-mag[i140], marker='o', alpha=1, color='teal', ls=' ' )
    pl.plot( z, mag[i153]-mag[i160], marker='o', alpha=0.3, color='darkorange', ls=' ' )
    pl.draw()


def plotMagz( snsimset, **plotargs ):
    """ plot the broad-band mag vs z (hubble diagram)
    :param snsimset:
    :return:
    """
    import numpy as np
    from matplotlib import pyplot as pl


    z = np.array([ sn.meta['z'] for sn in snsimset ])

    # the magnitude if we had a zeropoint of 0
    mag = np.array([ -2.5*np.log10( sn['flux'] ) + sn['zp'] for sn in snsimset ])
    #magerr = np.array( [2.5 * np.log10(np.e) *  sn['fluxerr']/ sn['flux'] for sn in snsimset])

    flt = np.array( [sn['band'] for sn in snsimset] )

    i127 = np.where( flt=='wfc3f127m' )
    i139 = np.where( flt=='wfc3f139m' )
    i153 = np.where( flt=='wfc3f153m' )
    i125 = np.where( flt=='wfc3f125w' )
    i140 = np.where( flt=='wfc3f140w' )
    i160 = np.where( flt=='wfc3f160w' )


    plotdefaults = {'marker':'o', 'alpha':0.3, 'color':'darkred', 'ls':' '}
    plotdefaults.update( **plotargs )
    # pl.clf()
    #pl.plot( z, mag[i125], marker='o', alpha=0.3, color='darkmagenta', ls=' ' )
    #pl.plot( z, mag[i140], marker='o', alpha=1, color='teal', ls=' ' )
    #pl.errorbar( z, mag[i160], magerr[i160], marker='o', alpha=0.3, color='darkred', ls=' ' )
    #mag160 = mag0[i160] + (zptTrue['wfc3f160w'])
    pl.plot( z, mag[i160], **plotdefaults )
    pl.draw()


def plotColorColorSN( sn ) :
    """
    Plot a single SN onto color-color diagrams.
    :param sn:
    :return:
    """
    import numpy as np
    from matplotlib import pyplot as pl

    # TODO : pick out the med-wide pairs observed at the same time.

    mag = -2.5*np.log10( sn['flux'] ) + sn['zp']
    dmag = 2.5*np.log10(np.e) * sn['fluxerr'] / sn['flux']

    i127 = np.where( sn['band'] == 'f127m' )[0]
    i125 = np.where( sn['band'] == 'f125w' )[0]
    i139 = np.where( sn['band'] == 'f139m' )[0]
    i140 = np.where( sn['band'] == 'f140w' )[0]
    i153 = np.where( sn['band'] == 'f153m' )[0]
    i160 = np.where( sn['band'] == 'f160w' )[0]

    c127 = mag[i127]-mag[i125]
    c139 = mag[i139]-mag[i140]
    c153 = mag[i153]-mag[i160]

    dc127 = np.sqrt( dmag[i127]**2 + dmag[i125]**2)
    dc139 = np.sqrt( dmag[i139]**2 + dmag[i140]**2)
    dc153 = np.sqrt( dmag[i153]**2 + dmag[i160]**2)

    fig = pl.gcf()
    ax1 = fig.add_subplot( 1,2,1 )
    ax2 = fig.add_subplot( 1,2,2, sharex=ax1 )

    ax1.errorbar( c139, c127, dc127, dc139, marker='o', color='k', capsize=0, lw=2)
    ax2.errorbar( c139, c153, dc153, dc139, marker='o', color='k', capsize=0, lw=2)









def scumsum( a ):
    """
    Sorted Cumulative Sum function :
    Construct an array "sumabove" such that the cell at index i in sumabove
    is equal to the sum of all cells from the input array "a" that have a
    cell value higher than a[i]
    """
    # Collapse the array into 1 dimension
    sumabove = a.ravel()

    # Sort the raveled array by descending cell value
    iravelsorted = sumabove.argsort( axis=0 )[::-1]

    # Reassign each cell to be the cumulative sum of all
    # input array cells with a higher value :
    sumabove[iravelsorted] = sumabove[iravelsorted].cumsum()

    # Now unravel back into shape of original array and return
    return( sumabove.reshape( a.shape ) )

