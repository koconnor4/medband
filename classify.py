
# Dictionary of sncosmo CCSN model names and their corresponding SN sub-type
ccSubClassDict = {
    'ib.01':'Ib','ib.02':'Ib','ib.03':'Ib','ib.04':'Ib','ib.05':'Ib','ib.06':'Ib','ib.07':'Ib',
    'ic.01':'Ic','ic.02':'Ic','ic.03':'Ic','ic.04':'Ic','ic.05':'Ic','ic.06':'Ic','ic.07':'Ic','ic.08':'Ic','ic.09':'Ic',
    'iip.01':'IIP','iip.02':'IIP','iip.03':'IIP','iip.04':'IIP','iip.05':'IIP','iip.06':'IIP','iip.07':'IIP','iip.08':'IIP','iip.09':'IIP','iip.10':'IIP',
    'iip.11':'IIP','iip.12':'IIP','iip.13':'IIP','iip.14':'IIP','iip.15':'IIP','iip.16':'IIP','iip.17':'IIP','iip.18':'IIP','iip.19':'IIP','iip.20':'IIP',
    'iip.21':'IIP','iip.22':'IIP','iip.23':'IIP','iip.24':'IIP',
    'iin.01':'IIn','iin.02':'IIn',
    }
iaSubClassDict = {
    # 'salt2':'Ia',
    'salt2-extended':'Ia',
    }

# Probability that a CC SN belongs to any given CC sub-class (Ib,Ic,IIP,IIL,IIn)
# from Li et al 2011a
ccSubClassProbs = {'Ib':0.19/0.76*0.46,'Ic':0.19/0.76*0.54,
                   'IIP':0.57/0.76*0.7,'IIn':0.57/0.76*0.3 }

class SncosmoSim( object ):
    """ A class for holding the results of sncosmo simulations,
    including parameters, cosmology, and light curves.
    """
    def __init__(self, sntype, observations=None, z_range=[1.8,2.2],
                 t0_range=[0,0], nsim=100, perfect=True,
                 Om=0.3, H0=70 ):
        """ Run a monte carlo sim using sncosmo to simulate <nsim> SNe
        of the given <sntype> over the given <z_range>.

        Simulates Type Ia SNe with the SALT2 model, and CC SNe with
        the SNANA CC templates.

        Set perfect=True for noiseless "observations" of the simulated SNe.
        :return:
        """
        from astropy import cosmology
        import sncosmo
        from sncosmohst import hstbandpasses, ccsnmodels
        from numpy.random import normal, uniform, choice
        import numpy as np

        self.sntype = sntype
        self.z_range = z_range
        self.nsim = nsim
        self.perfect = perfect

        if observations is None :
            observations = mkobservationsTable( )
        self.observations = observations

        # Make a list of all the unique sncosmo source models available,
        # and assign a relative probability that any given simulated SN of this
        # type (CC or Ia) belongs to that subclass
        if sntype.lower()=='cc' :
            self.SourcenameSet = np.array( ccSubClassDict.keys() )
            self.SubclassSet = np.array([ ccSubClassDict[source] for source in self.SourcenameSet ])
            self.SubclassCount = np.array([ len(np.where(self.SubclassSet==subclass)[0])
                                         for subclass in self.SubclassSet ], dtype=float)
            self.SourceprobSet = np.array([ ccSubClassProbs[subclass]
                                          for subclass in self.SubclassSet ]) / self.SubclassCount
            self.SourceprobSet /= self.SourceprobSet.sum()
        else :
            # No sub-class divisions for SNIa
            self.SourcenameSet = np.array(['salt2-extended'])
            self.SubclassSet = np.array( ['Ia'] )
            self.SourceprobSet = np.array( [1] )
            self.SubclassCount = np.array( [1] )

        # load the O'Donnell 1994 dust model
        self.dust = sncosmo.OD94Dust()

        # Define an sncosmo SN model for each available source
        modelset = np.array([ sncosmo.Model(source=source, effects=[self.dust],
                                    effect_names=['host'], effect_frames=['rest'])
                                 for source in self.SourcenameSet ])
        # Define a cosmology
        # self.Om = Om
        # self.H0 = H0
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


def testClassifySNR( nsim=1000, nclass=None, outfileroot='colorColorClassify',
                     clobber=False, verbose=True ):
    import cPickle
    import os
    import numpy as np
    from copy import deepcopy
    from numpy.random import normal
    from time import asctime

    if nclass is None :
        nclass = nsim

    outfileIa = outfileroot + '_Ia.dat'
    outfileCC = outfileroot + '_CC.dat'
    if os.path.isfile( outfileIa ) and not clobber :
        print("%s exists. Use clobber to overwrite."%outfileIa)
    if os.path.isfile( outfileCC ) and not clobber :
        print("%s exists. Use clobber to overwrite."%outfileCC)

    classSimIapkl = 'colorColorClassify.classSimIa.%i.pkl'%nsim
    classSimCCpkl = 'colorColorClassify.classSimCC.%i.pkl'%nsim

    if os.path.isfile( classSimIapkl ) and not clobber>1 :
        if verbose: print("Loading Ia simulation from pickle : %s"%classSimIapkl)
        fin = open( classSimIapkl, 'rb' )
        classSimIa = cPickle.load( fin )
        fin.close()
    else :
        if verbose: print("Running a new Ia simulation, then saving to pickle : %s"%classSimIapkl)
        classSimIa = SncosmoSim( 'Ia' , nsim=nsim )
        fout = open( classSimIapkl, 'wb' )
        cPickle.dump( classSimIa, fout, protocol=-1 )
        fout.close()

    if os.path.isfile( classSimCCpkl ) and not clobber>1 :
        if verbose: print("Loading CC simulation from pickle : %s"%classSimCCpkl)
        fin = open( classSimCCpkl, 'rb' )
        classSimCC = cPickle.load(fin)
        fin.close()
    else :
        if verbose: print("Running a new CC simulation, then saving to pickle : %s"%classSimCCpkl)
        classSimCC = SncosmoSim( 'CC' , nsim=nsim )
        fout = open( classSimCCpkl, 'wb' )
        cPickle.dump( classSimCC, fout, protocol=-1 )
        fout.close()

    if verbose: print("Writing classification output to %s and %s "%(outfileIa, outfileCC))

    # erase any existing output files, and start a new version with a header line
    foutIa = open( outfileIa, 'w')
    print >> foutIa, '# %s'% asctime()
    print >> foutIa, '# Med Band Classifier Validation Test output generated with %s'% __file__
    print >> foutIa, '# Classification samples : %s  , %s '%(classSimIapkl,classSimCCpkl)
    print >> foutIa, '# Classifying %i simulated Type Ia SN'%nclass
    print >> foutIa, '# P(Ia|D)  S/N   z    t0   x0   x1   c  Av   Rv'
    foutIa.close()

    foutCC = open( outfileCC, 'w')
    print >> foutCC, '# %s'% asctime()
    print >> foutCC, '# Med Band Classifier Validation Test output generated with %s'% __file__
    print >> foutCC, '# Classification samples : %s  , %s '%(classSimIapkl,classSimCCpkl)
    print >> foutCC, '# Classifying %i simulated CC SN'%nclass
    print >> foutCC, '# P(Ia|D)   S/N   model  z  amplitude  t0  Av   Rv'
    foutCC.close()

    for isn in range(nclass) :
        if verbose: print("Classifying CCSN %i  and SNIa %i of  %i"%(isn, isn, nclass))
        for SNR in [10,20,30]:

            # -----------------------------------------------------------------
            # Classify one SNIa test particle :
            x1 = classSimIa.x1[isn]
            c = classSimIa.c[isn]
            dx1c = np.sqrt( (classSimIa.x1-x1)**2 + (classSimIa.c-c)**2 )
            iNotThisx1c = np.where( dx1c>0.05 )[0]

            # Fix the S/N ratio and add appropriate noise to the perfect flux
            testlcIa = deepcopy( classSimIa.lightcurves[isn] )
            testlcIa['fluxerr'] = testlcIa['flux'] / SNR
            testlcIa['flux'] = normal( testlcIa['flux'], testlcIa['fluxerr'] )

            # exclude from the classification set any simulated SN with (almost)
            # the same x1 and c as our test particle
            classSimIaLCsubset = [ classSimIa.lightcurves[inot] for inot in iNotThisx1c ]

            # Run the classifier and append the result to the running output file.
            testsnIapIa = colorColorClassify( testlcIa, classSimIaLCsubset, classSimCC.lightcurves,  )
            foutIa = open( outfileIa, 'a')
            print >> foutIa, '%4.2f  %i  %5.3f  %7.1f   %9.3e  %5.2f  %5.2f  %5.2f  %5.2f'%(
                testsnIapIa, SNR, classSimIa.z[isn],  classSimIa.t0[isn],
                classSimIa.x0[isn], classSimIa.x1[isn], classSimIa.c[isn],
                classSimIa.Av[isn], classSimIa.Rv[isn] )
            foutIa.close()

            # -----------------------------------------------------------------
            # Classify one CCSN, excluding from the classification set any
            #  simulated SN using the same CC model
            testlcCC = deepcopy( classSimCC.lightcurves[isn] )
            testlcCC['fluxerr'] = testlcCC['flux'] / SNR
            testlcCC['flux'] = normal( testlcCC['flux'], testlcCC['fluxerr'] )
            sourcename = classSimCC.sourcename[isn]
            iNotThisSource = np.where( classSimCC.sourcename != sourcename )[0]
            classSimCCLCsubset = [ classSimCC.lightcurves[inot] for inot in iNotThisSource ]
            testsnCCpIa = colorColorClassify( testlcCC, classSimIa.lightcurves, classSimCCLCsubset )

            # Append the classification result to the running output file.
            foutCC = open( outfileCC, 'a')
            print >> foutCC, '%4.2f  %i  %9s  %5.3f  %9.3e  %7.1f  %5.2f  %5.2f'%(
                testsnCCpIa, SNR, classSimCC.sourcename[isn], classSimCC.z[isn],
                classSimCC.amplitude[isn], classSimCC.t0[isn],
                classSimCC.Av[isn], classSimCC.Rv[isn] )
            foutCC.close()
    if verbose: print("Done classifying %i SN of each class."%(nclass))

    return( outfileIa, outfileCC )

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
    ax2.set_xlabel('P(Ia|D)')
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


def testClassify100():
    simIalist = [ mcSim( 'Ia', nsim=100 ) for i in range(10) ]
    simCClist = [ mcSim( 'Ia', nsim=100 ) for i in range(10) ]
    pIaArrayLog = []
    pIaArray = []
    for i in range(10) :
        simIa = simIalist[i]
        for j in range(10) :
            simCC = simCClist[j]
            pIaArrayLog.append( colorColorClassifyLogLike( testsnIa, simIa, simCC ) )
            pIaArray.append( colorColorClassify( testsnIa, simIa, simCC ) )
    return( pIaArrayLog, pIaArray )

def testClassify1000():
    simIalist = [ mcSim( 'Ia', nsim=1000 ) for i in range(10) ]
    simCClist = [ mcSim( 'CC', nsim=1000 ) for i in range(10) ]
    pIaArray = []
    for i in range(10) :
        simIa = simIalist[i]
        for j in range(10) :
            simCC = simCClist[j]
            pIaArray.append( colorColorClassify( testsnIa, simIa, simCC ) )
    return( simIalist, simCClist, pIaArray )

def colorColorClassifyLogLike( sn=testsnIa, simIa=None, simCC=None,
                               z_range=[1.8,2.2], nsim=5 ):
    """  Classify a SN with med-wide band colors by comparing
    measured flux ratios against those from sncosmo sims.
    :return:
    """
    import numpy as np

    # run monte carlo simulations if needed
    if simIa is None :
        simIa = SncosmoSim( 'Ia', z_range=z_range, nsim=nsim, perfect=True )
    if simCC is None :
        simCC = SncosmoSim( 'CC', z_range=z_range, nsim=nsim, perfect=True )

    # extract the flux ratios from the observed and simulated SN data
    fRsn, dfRsn = [], []
    fRsimIa, fRsimCC = [], []
    for bandpair in [['f127m','f125w'],['f139m','f140w'],['f153m','f160w']]:
        imed = np.where( sn['band']==bandpair[0])
        iwide = np.where( sn['band']==bandpair[1])
        fRatio = sn['flux'][imed] / sn['flux'][iwide]
        fRsn.append( fRatio )
        fRatioErr = np.sqrt(
            (sn['fluxerr'][imed]/sn['flux'][imed])**2 + \
            (sn['fluxerr'][iwide]/sn['flux'][iwide])**2 ) * np.abs( fRatio )
        dfRsn.append( fRatioErr )

        # Assuming all simulated SN have the same observations, as is the case
        # when the sim is generated with the mcSim() function.
        imedsim = np.where( simIa[0]['band']==bandpair[0] )
        iwidesim = np.where( simIa[0]['band']==bandpair[1] )
        fRatioSimIa = np.array( [ simIa[isim]['flux'][imedsim] / simIa[isim]['flux'][iwidesim]
                            for isim in range(len(simIa))] )
        fRsimIa.append( fRatioSimIa.T[0] )

        fRatioSimCC = np.array( [ simCC[isim]['flux'][imedsim] / simCC[isim]['flux'][iwidesim]
                            for isim in range(len(simCC))] )
        fRsimCC.append( fRatioSimCC.T[0] )

    fRsn = np.array( fRsn, dtype=float )
    dfRsn = np.array( dfRsn, dtype=float )
    fRsimIa = np.array( fRsimIa, dtype=float )
    fRsimCC = np.array( fRsimCC, dtype=float )

    def loglike( fRsn, fRsim, dfRsn ):
        """  Returns the natural log likelihood matrix, derived
        by comparing the observed flux ratios of the real SN (fRsn, dfRsn),
        against the flux ratios of a simulated SN (fRsim).

        :param fRsn: an N-dimensional array giving the ratios of medium to
               wide band fluxes from an observed SN, for N med+wide
               bandpass pairs.
        :param fRsim: an NxM-dimensional array, giving the same N flux ratios
               for each of M simulated SN.
        :param dfRsn: an N-dimensional array giving observed flux ratio uncertainties
        :return: M-dimensional array giving the log likelihood from comparing
               the observed SN to each of the M simulated SN.
        """
        # likelihood = 1/sqrt(2*pi*sigma^2)  *  exp( -(fRsn-fRsim)^2/dfRsn^2 )
        halfchi2 = 0.5 * (fRsn-fRsim)**2 / dfRsn**2  # 3x1k
        lognormfactor = np.log( np.sqrt(2*np.pi) * dfRsn ) # 3x1
        loglike =  -lognormfactor.sum(axis=0) - halfchi2.sum(axis=0) # 1k
        return(  loglike )

    loglikeIa = loglike( fRsn, fRsimIa, dfRsn )
    loglikeCC = loglike( fRsn, fRsimCC, dfRsn )
    likeIaSum = np.exp( loglikeIa ).sum()
    likeCCSum = np.exp( loglikeCC ).sum()
    pIaFromLog = likeIaSum / ( likeIaSum + likeCCSum )

    return( pIaFromLog )


def colorColorClassify( sn=testsnIa, simIa=None, simCC=None,
                        z_range=[1.8,2.2], nsim=5 ):
    """  Classify a SN with med-wide band colors by comparing
    measured flux ratios against those from sncosmo sims.
    :return:
    """
    import numpy as np

    # run monte carlo simulations if needed
    if simIa is None :
        simIa = SncosmoSim( 'Ia', z_range=z_range, nsim=nsim, perfect=True )
    if simCC is None :
        simCC = SncosmoSim( 'CC', z_range=z_range, nsim=nsim, perfect=True )

    # extract the flux ratios from the observed and simulated SN data
    fRsn, dfRsn = [], []
    fRsimIa, fRsimCC = [], []
    for bandpair in [['f127m','f125w'],['f139m','f140w'],['f153m','f160w']]:
        imed = np.where( sn['band']==bandpair[0])
        iwide = np.where( sn['band']==bandpair[1])
        fRatio = sn['flux'][imed] / sn['flux'][iwide]
        fRsn.append( fRatio )
        fRatioErr = np.sqrt(
            (sn['fluxerr'][imed]/sn['flux'][imed])**2 + \
            (sn['fluxerr'][iwide]/sn['flux'][iwide])**2 ) * np.abs( fRatio )
        dfRsn.append( fRatioErr )

        # Assuming all simulated SN have the same observations, as is the case
        # when the sim is generated with the mcSim() function.
        imedsim = np.where( simIa[0]['band']==bandpair[0] )
        iwidesim = np.where( simIa[0]['band']==bandpair[1] )
        fRatioSimIa = np.array( [ simIa[isim]['flux'][imedsim] / simIa[isim]['flux'][iwidesim]
                            for isim in range(len(simIa))] )
        fRsimIa.append( fRatioSimIa.T[0] )

        fRatioSimCC = np.array( [ simCC[isim]['flux'][imedsim] / simCC[isim]['flux'][iwidesim]
                            for isim in range(len(simCC))] )
        fRsimCC.append( fRatioSimCC.T[0] )

    fRsn = np.array( fRsn, dtype=float )
    dfRsn = np.array( dfRsn, dtype=float )
    fRsimIa = np.array( fRsimIa, dtype=float )
    fRsimCC = np.array( fRsimCC, dtype=float )

    def like( fRsn, fRsim, dfRsn ):
        """  Returns the likelihood matrix, derived
        by comparing the observed flux ratios of the real SN (fRsn, dfRsn),
        against the flux ratios of a simulated SN (fRsim).

        :param fRsn: an N-dimensional array giving the ratios of medium to
               wide band fluxes from an observed SN, for N med+wide
               bandpass pairs.
        :param fRsim: an NxM-dimensional array, giving the same N flux ratios
               for each of M simulated SN.
        :param dfRsn: an N-dimensional array giving observed flux ratio uncertainties
        :return: M-dimensional array giving the log likelihood from comparing
               the observed SN to each of the M simulated SN.
        """
        # likelihood = 1/sqrt(2*pi*sigma^2)  *  exp( -(fRsn-fRsim)^2/dfRsn^2 )
        halfchi2 = 0.5 * (fRsn-fRsim)**2 / dfRsn**2  # 3x1k
        normfactor = 1./(np.sqrt(2*np.pi)*dfRsn) # 3x1
        like =  np.prod( normfactor*np.exp(-halfchi2), axis=0) # 1k
        return(  like.sum() )

    likeIa = like( fRsn, fRsimIa, dfRsn )
    likeCC = like( fRsn, fRsimCC, dfRsn )
    pIa = likeIa / ( likeIa + likeCC )

    return( pIa )

def mkobservationsTable( orbits=10.  ):
    from sncosmohst.simparam import gethstzp, gethstbgnoise
    import numpy as np
    from astropy.table import Table

    # medium band at peak observation set :
    bandlist = ['f127m','f139m','f153m','f125w','f140w','f160w']
    if not np.iterable( orbits ) : orbits = np.ones(len(bandlist)) * orbits
    exptimelist = orbits * 2500

    timelist = np.zeros( len( bandlist ) )
    zplist = [ gethstzp(band) for band in bandlist ]
    zpsyslist = [ 'ab' for i in range(len(bandlist)) ]
    # gainlist = [ gethstgain( et ) for et in exptimelist ]
    gainlist = [ et for et in exptimelist ]
    bgnoiselist = [ gethstbgnoise( band, et ) for band,et in zip(bandlist,exptimelist) ]
    # bgnoiselist = np.zeros( len(bandlist) )
    observations = Table({'band': bandlist, 'time': timelist, 'zp': zplist,
                          'zpsys': zpsyslist, 'gain': gainlist, 'exptime':exptimelist,
                          'skynoise': bgnoiselist,
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



def plotColorColor( snsim, plotstyle='points', nbins=None, **plotargs ):
    """ plot the med band pseudo-colors in a color-color diagram
    :param snsim:
    :return:
    """
    import numpy as np
    from matplotlib import pyplot as pl
    from matplotlib import ticker

    z = np.array([ sn.meta['z'] for sn in snsim ])
    mag = np.array([ -2.5*np.log10( sn['flux'] ) + sn['zp'] for sn in snsim ])
    # magerr = np.array( [2.5 * np.log10(np.e) *  sn['fluxerr']/ sn['flux'] for sn in sne])

    flt = np.array( [sn['band'] for sn in snsim] )

    i127 = np.where( flt=='f127m' )
    i139 = np.where( flt=='f139m' )
    i153 = np.where( flt=='f153m' )
    i125 = np.where( flt=='f125w' )
    i140 = np.where( flt=='f140w' )
    i160 = np.where( flt=='f160w' )

    fig = pl.gcf()
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2, sharex=ax1)

    plotdefaults = {'marker':'o', 'alpha':0.3, 'color':'darkorange', 'ls':' '}
    plotdefaults.update( **plotargs )
    plotargs = plotdefaults

    if plotstyle=='points':
        ax1.plot( mag[i139]-mag[i140], mag[i127]-mag[i125], **plotdefaults )
        ax2.plot( mag[i139]-mag[i140], mag[i153]-mag[i160], **plotdefaults )

    elif plotstyle.startswith('contour'):
        xarray = mag[i139]-mag[i140]
        nsim = len(xarray)
        if nbins is None : nbins = int( np.sqrt( nsim  ) )

        plotdefaults = {'levels':[0.68,0.0],'colors':['r','g'],'ls':'-',
                        'alpha':0.5, 'extend':'neither'}
        plotdefaults.update( **plotargs )

        # Plot filled contours, showing  the full extent of the population,
        # and contour lines containing 68% of the population.
        # First, bin the points into a 2-d histogram:
        # (Note that we reverse the x-y order here to get the binned arrays
        #  plotted in the correct direction )
        for yarray,ax in zip( [mag[i127]-mag[i125], mag[i153]-mag[i160]], [ax1,ax2] ):

            count,y,x = np.histogram2d( yarray, xarray, bins=nbins,
                                        range=[[-0.5,0.5],[-0.5,0.5]] )

            # Renormalize relative to the sum of all SNe in this class :
            count /= count.sum()

            # Now set up an array 'cabove' such that  the cell value in cabove[i,j]
            # is equal to the sum of all cells that have a value higher than c[i,j]
            cabove = scumsum( count )

            # solid lines give probability contours at specified levels
            # (defaults to 0.68 for "1-sigma contours")
            ax.contour( x[:-1], y[:-1], cabove, **plotdefaults )
            if plotstyle=='contourf' :
                ax.contourf( x[:-1], y[:-1], cabove, **plotdefaults )

    ax1.set_xlabel( 'F139M-F140W' )
    ax1.set_ylabel( 'F127M-F125W' )
    ax2.set_xlabel( 'F139M-F140W' )
    ax2.set_ylabel( 'F153M-F160W' , rotation=-90)
    ax2.yaxis.set_label_position('right')
    ax2.yaxis.set_ticks_position('right')
    ax2.yaxis.set_ticks_position('both')

    ax1.xaxis.set_major_locator( ticker.MultipleLocator( 0.1 ) )
    ax1.xaxis.set_minor_locator( ticker.MultipleLocator( 0.05 ) )
    ax2.xaxis.set_major_locator( ticker.MultipleLocator( 0.1 ) )
    ax2.xaxis.set_minor_locator( ticker.MultipleLocator( 0.05 ) )
    ax1.yaxis.set_major_locator( ticker.MultipleLocator( 0.1 ) )
    ax1.yaxis.set_minor_locator( ticker.MultipleLocator( 0.05 ) )
    ax2.yaxis.set_major_locator( ticker.MultipleLocator( 0.1 ) )
    ax2.yaxis.set_minor_locator( ticker.MultipleLocator( 0.05 ) )

    ax1.set_xlim( -0.25, 0.35 )
    ax1.set_ylim( -0.49, 0.25 )
    ax2.set_ylim( -0.19, 0.35 )


    fig.subplots_adjust( left=0.12, bottom=0.12, right=0.88, top=0.95, wspace=0.1 )

    pl.draw()


def mkCircleFig( simIa=None, simCC=None, z_range=[1.8,2.2], nsim=100,
                 nbins=None, **plotargs ):
    """ Make a circle diagram, i.e. a med-wide band pseudo-color-color plot,
    showing both Type Ia and CC simulations over the given redshift range.
    :param snsim:
    :return:
    """
    from matplotlib import pyplot as pl
    from pytools import plotsetup
    import numpy as np

    if simIa is None :
        simIa = mcSim( 'Ia', z_range=z_range, nsim=nsim )
    if simCC is None :
        simCC = mcSim( 'CC', z_range=z_range, nsim=nsim )

    fig = plotsetup.fullpaperfig( figsize=[8,4] )

    plotColorColor(simIa,nbins=nbins, plotstyle='contourf',levels=[0.68,0.], colors=['darkmagenta'] )
    plotColorColor(simCC,nbins=nbins, plotstyle='contourf',levels=[0.68,0.], colors=['teal'] )

    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2, sharex=ax1)

    ax1.set_xlim( -0.25, 0.3 )
    ax1.set_ylim( -0.45, 0.2 )
    ax2.set_ylim( -0.2, 0.3 )

    ax1.text( 0.1, 0.1, 'SN Ia', color='darkmagenta',
              ha='left',va='bottom')#,transform=ax1.transAxes,fontsize=14)
    ax1.text( 0., -0.05, 'CC SN', color='teal',
              ha='center',va='center')#,transform=ax1.transAxes,fontsize=14)
    pl.draw()


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


