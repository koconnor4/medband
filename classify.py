
# Dictionary of sncosmo CCSN model names and their corresponding SN sub-type
SubClassDict_SNANA = {
           'ii':{
               'iip.01':'IIP','iip.02':'IIP','iip.03':'IIP','iip.04':'IIP','iip.05':'IIP','iip.06':'IIP','iip.07':'IIP','iip.08':'IIP','iip.09':'IIP','iip.10':'IIP',
               'iip.11':'IIP','iip.12':'IIP','iip.13':'IIP','iip.14':'IIP','iip.15':'IIP','iip.16':'IIP','iip.17':'IIP','iip.18':'IIP','iip.19':'IIP','iip.20':'IIP',
               'iip.21':'IIP','iip.22':'IIP','iip.23':'IIP','iip.24':'IIP',
               'iin.01':'IIn','iin.02':'IIn',
               },
           'ibc':{
               'ib.01':'Ib','ib.02':'Ib','ib.03':'Ib','ib.04':'Ib','ib.05':'Ib','ib.06':'Ib','ib.07':'Ib',
               'ic.01':'Ic','ic.02':'Ic','ic.03':'Ic','ic.04':'Ic','ic.05':'Ic','ic.06':'Ic','ic.07':'Ic','ic.08':'Ic','ic.09':'Ic',
               },
           'ia': {'salt2':'Ia'},
}


SubClassDict_PSNID = {
           'ii':{ 's11-2004hx':'II','s11-2005lc':'IIP','s11-2005gi':'IIP','s11-2006jl':'IIP' },
           'ibc':{ 's11-2005hl':'Ib','s11-2005hm':'Ib','s11-2006fo':'Ic', 's11-2006jo':'Ib'},
           'ia': {'salt2':'Ia'},
}

SubClassDict_NOTPSNID = {
           'ii':{  # 'iip.01':'IIP','iip.02':'IIP','iip.05':'IIP',
               'iip.03':'IIP','iip.04':'IIP','iip.06':'IIP','iip.07':'IIP','iip.08':'IIP','iip.09':'IIP','iip.10':'IIP',
               'iip.11':'IIP','iip.12':'IIP','iip.13':'IIP','iip.14':'IIP','iip.15':'IIP','iip.16':'IIP','iip.17':'IIP','iip.18':'IIP','iip.19':'IIP','iip.20':'IIP',
               'iip.21':'IIP','iip.22':'IIP','iip.23':'IIP','iip.24':'IIP',
               'iin.01':'IIn','iin.02':'IIn',
               },
           'ibc':{
               'ib.01':'Ib','ib.02':'Ib','ib.03':'Ib','ib.04':'Ib','ib.07':'Ib', # 'ib.05':'Ib' 'ib.06':'Ib',
               'ic.01':'Ic','ic.02':'Ic','ic.03':'Ic', 'ic.05':'Ic','ic.06':'Ic','ic.07':'Ic','ic.08':'Ic','ic.09':'Ic', # 'ic.04':'Ic'
               },
           'ia': {'salt2':'Ia'},
}




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

def gauss( x, mu, sigma, range=None):
    """ Return values from a (bifurcated) gaussian.
    If sigma is a scalar, then this function returns a  symmetric
    normal distribution.

    If sigma is a 2-element iterable, then we define a bifurcated
    gaussian (i.e. two gaussians with different widths that meet with a
    common y value at x=mu)
    In this case, sigma must contain a positive value
    giving sigma for the right half gaussian, and a negative value giving
    sigma for the left half gaussian.

    If range is specified, then we include a normalization factor to
    ensure that the function integrates to unity over the given interval.
    """
    import numpy as np
    from scipy.special import erf

    if np.iterable( sigma ) :
        assert np.sign(sigma[0])!=np.sign(sigma[1]), \
            "sigma must be [+sigmaR,-sigmaL] or [-sigmaL,+sigmaR] :  " \
            "i.e. components must have opposite signs"
        sigmaL = - np.min( sigma )
        sigmaR = np.max( sigma )
    else :
        sigmaL = sigmaR = sigma

    if range is not None :
        normfactor = 2. / ( np.abs( erf( (range[0]-mu)/(np.sqrt(2)*sigmaL) ) ) + \
                            np.abs( erf( (range[1]-mu)/(np.sqrt(2)*sigmaR) ) ) )
    else :
        normfactor = 1.

    if np.iterable(x) and type(x) != np.ndarray :
        x = np.asarray( x )

    normaldist = lambda x,mu,sig : np.exp(-(x-mu)**2/(2*sig**2))/(np.sqrt(2*np.pi)*sig)
    gaussL = normfactor * 2*sigmaL/(sigmaL+sigmaR) * normaldist( x, mu, sigmaL )
    gaussR = normfactor * 2*sigmaR/(sigmaL+sigmaR) * normaldist( x, mu, sigmaR )

    if not np.iterable( x ) :
        if x <= mu :
            return( gaussL )
        else :
            return( gaussR )
    else :
        return( np.where( x<=mu, gaussL, gaussR ) )



def get_evidence( sn=testsnIa, modelsource='salt2',
                 zhost=None, zhosterr=None, t0_range=None,
                 zminmax=[0.1,2.8],
                 nobj=100, maxiter=1000, verbose=True ):
    """  compute the Bayesian evidence (and likelihood distributions)
    for the given SN class using the sncosmo nested sampling algorithm.
    :return:
    """
    import os
    import numpy as np
    from scipy import interpolate, integrate
    import sncosmo
    import time
    tstart = time.time()

    # standardize the data column names and normalize to zpt=25 AB
    sn = sncosmo.fitting.standardize_data( sn )
    sn = sncosmo.fitting.normalize_data( sn )

    # Define parameter bounds and priors for z, x1, c, Rv, etc
    if zhost is None :
        zhost = None
    elif isinstance(zhost,str) :
        # read in the z prior from a file giving z and p(z)
        assert os.path.isfile( zhost ), "If zprior is a string, it must be a filename"
        z,pdf = np.loadtxt( zhost, unpack=True )
        # normalize so that it integrates to unity over the allowed z range
        izgood = np.where( (zminmax[0]<z) & (z<zminmax[1]) )[0]
        pdfint = integrate.simps( pdf[izgood], z[izgood] )
        pdf = pdf / pdfint
        zprior = interpolate.interp1d( z, pdf, bounds_error=False, fill_value=0)
        import pdb; pdb.set_trace()
    else :
        if zhosterr is None :
            zhosterr = 0.1
        if np.iterable( zhosterr ) :
            assert np.sign(zhosterr[0])!=np.sign(zhosterr[1]), \
                "zphoterr must be [+err,-err] or [-err,+err] :  " \
                "i.e. components must have opposite signs"
            zhostminus = - np.min( zhosterr )
            zhostplus = np.max( zhosterr )
        else :
            zhostminus = zhostplus = zhosterr
        zmin, zmax = zminmax
        zminmax = [ max( zmin, zhost-zhostminus*5), min(zmax,zhost+zhostplus*5) ]
        def zprior( z ) :
            return( gauss( z, zhost, [-zhostminus,zhostplus], range=zminmax ) )

    if t0_range is None :
        t0_range = [sn['time'].min()-20,sn['time'].max()+20]

    if zhosterr>0.01 :
        bounds={'z':(zminmax[0],zminmax[1]),'t0':(t0_range[0],t0_range[1]) }
    else :
        bounds={'t0':(t0_range[0],t0_range[1]) }

    if modelsource.lower().startswith('salt2') :
        # define the Ia SALT2 model parameter bounds and priors
        model = sncosmo.Model( source=modelsource)
        if zhosterr>0.01 :
            param_names = ['z','t0','x0','x1','c']
        else :
            param_names = ['t0','x0','x1','c']
        bounds['x1'] = (-5.,5.)
        bounds['c'] = (-0.5,3.0)
        def x1prior( x1 ) :
            return( gauss( x1, 0, [-1.5,0.9], range=bounds['x1'] ) )
        def cprior( c ) :
            return( gauss( c, 0, [-0.08,0.14], range=bounds['c'] ) )
        if zhost :
            priorfn = {'z':zprior, 'x1':x1prior, 'c':cprior}
        else :
            priorfn = { 'x1':x1prior, 'c':cprior }

    else :
        # define a host-galaxy dust model
        dust = sncosmo.CCM89Dust( )

        # Define the CC model, parameter bounds and priors
        model = sncosmo.Model( source=modelsource, effects=[dust],
                               effect_names=['host'], effect_frames=['rest'])

        if zhosterr>0.01 :
            param_names = ['z','t0','amplitude','hostebv','hostr_v']
        else :
            param_names = ['t0','amplitude','hostebv','hostr_v']
        bounds['hostebv'] = (0.0,1.0)
        bounds['hostr_v'] = (2.0,4.0)
        def rvprior( rv ) :
            return( gauss( rv, 3.1, 0.3, range=bounds['host_rv'] ) )
        # TODO : include a proper Av or E(B-V) prior for CC models
        if zhost and zhosterr>0.01:
            priorfn = {'z':zprior, 'rv':rvprior }
        else :
            priorfn = { 'rv':rvprior }

    model.set(z=np.mean(zminmax))

    res, fit = sncosmo.fitting.nest_lc( sn, model, param_names, bounds,
                                        guess_amplitude_bound=True,
                                        priors=priorfn,
                                        nobj=nobj, maxiter=maxiter,
                                        verbose=verbose )
    tend = time.time()
    if verbose : print("  Total Fitting time = %.1f sec"%(tend-tstart))
    return( sn, res, fit, priorfn )

def get_marginal_pdfs( res, nbins=51, verbose=True ):
    """ Given the results <res> from a nested sampling chain, determine the
    marginalized posterior probability density functions for each of the
    parameters in the model.

    :param res:  the results of a nestlc run
    :param nbins: number of bins (steps along the x axis) for sampling
       each parameter's marginalized posterior probability
    :return: a dict with an entry for each parameter, giving a 2-tuple containing
       NDarrays of length nbins.  The first array in each pair gives the parameter
       value that defines the left edge of each bin along the parameter axis.
       The second array gives the posterior probability density integrated
       across that bin.
    """
    import numpy as np
    import sncosmo
    param_names = res.param_names
    weights = res.weights
    samples = res.samples
    pdfdict = {}

    for param in param_names :
        ipar = param_names.index( param )
        paramvals = samples[:,ipar]

        if nbins>1:
            if param in res.bounds :
                parvalmin, parvalmax = res.bounds[param]
            else :
                parvalmin, parvalmax = 0.99*paramvals.min(), 1.01*paramvals.max()
            parambins = np.linspace( parvalmin, parvalmax, nbins, endpoint=True )
            binindices = np.digitize( paramvals, parambins )

            # we estimate the marginalized pdf by summing the weights of all points in the bin,
            # where the weight of each point is the prior volume at that point times the
            # likelihood, divided by the total evidence
            pdf = np.array( [ weights[np.where( binindices==ibin )].sum() for ibin in range(len(parambins)) ] )
        else :
            parambins = None
            pdf = None

        mean = (weights  * samples[:,ipar]).sum()
        std = np.sqrt( (weights * (samples[:,ipar]-mean)**2 ).sum() )

        pdfdict[param] = (parambins,pdf,mean,std)

        if verbose :
            if np.abs(std)>=0.1:
                print( '  <%s> =  %.2f +- %.2f'%( param, np.round(mean,2), np.round(std,2))  )
            elif np.abs(std)>=0.01:
                print( '  <%s> =  %.3f +- %.3f'%( param, np.round(mean,3), np.round(std,3)) )
            elif np.abs(std)>=0.001:
                print( '  <%s> =  %.4f +- %.4f'%( param, np.round(mean,4), np.round(std,4)) )
            else :
                print( '  <%s> = %.3e +- %.3e'%( param, mean, std) )


        if param == 'x0' :
            salt2 = sncosmo.Model( source='salt2')
            salt2.source.set_peakmag( 0., 'bessellb', 'ab' )
            x0_AB0 = salt2.get('x0')
            mBmean = -2.5*np.log10(  mean / x0_AB0 )
            mBstd = 2.5*np.log10( np.e ) *  std / mean
            mBbins = -2.5*np.log10(  parambins / x0_AB0 )

            pdfdict['mB'] = ( mBbins, pdf, mBmean, mBstd )
            print( '  <%s> =  %.3f +- %.3f'%( 'mB', np.round(mBmean,3), np.round(mBstd,3)) )

    return( pdfdict )




def plot_marginal_pdfs( res, nbins=101, **kwargs):
    """ plot the results of a classification run
    :return:
    """
    from matplotlib import pyplot as pl
    import numpy as np

    nparam = len(res.param_names)
    # nrow = np.sqrt( nparam )
    # ncol = nparam / nrow + 1
    nrow, ncol = 1, nparam

    pdfdict = get_marginal_pdfs( res, nbins )

    fig = pl.gcf()
    for parname in res.param_names :
        iax = res.param_names.index( parname )+1
        ax = fig.add_subplot( nrow, ncol, iax )

        parval, pdf, mean, std = pdfdict[parname]
        ax.plot(  parval, pdf, **kwargs )
        if np.abs(std)>=0.1:
            ax.text( 0.95, 0.95, '%s  %.1f +- %.1f'%( parname, np.round(mean,1), np.round(std,1)),
                     ha='right',va='top',transform=ax.transAxes )
        elif np.abs(std)>=0.01:
            ax.text( 0.95, 0.95, '%s  %.2f +- %.2f'%( parname, np.round(mean,2), np.round(std,2)),
                     ha='right',va='top',transform=ax.transAxes )
        elif np.abs(std)>=0.001:
            ax.text( 0.95, 0.95, '%s  %.3f +- %.3f'%( parname, np.round(mean,3), np.round(std,3)),
                     ha='right',va='top',transform=ax.transAxes )
        else :
            ax.text( 0.95, 0.95, '%s  %.3e +- %.3e'%( parname, mean, std),
                     ha='right',va='top',transform=ax.transAxes )

    pl.draw()


def classify( sn, zhost=1., zhosterr=1., t0_range=None, zminmax=[0.1,2.8],
              nobj=100, maxiter=1000, templateset='PSNID',
              nsteps_pdf=0, verbose=True ):
    """  Collect the bayesian evidence for all SN sub-classes.
    :param sn:
    :param zhost:
    :param zhosterr:
    :param t0_range:
    :param zminmax:
    :param nobj:
    :param maxiter:
    :param verbose:
    :return:
    """
    import numpy as np
    import time

    tstart = time.time()

    if templateset.lower()=='psnid':
        SubClassDict = SubClassDict_PSNID
    elif templateset.lower()=='snana':
        SubClassDict = SubClassDict_SNANA
    elif templateset.lower()=='notpsnid':
        SubClassDict = SubClassDict_NOTPSNID

    iimodelnames = SubClassDict['ii'].keys()
    ibcmodelnames = SubClassDict['ibc'].keys()
    iamodelnames = SubClassDict['ia'].keys()

    outdict = {}
    allmodelnames = np.append(np.append(iimodelnames, ibcmodelnames), iamodelnames)

    logpriordict = {
        'ia':np.log(0.24/len(iamodelnames)),
        'ibc':np.log(0.19/len(ibcmodelnames)),
        'ii':np.log(0.57/len(iimodelnames)),
        }

    logz = { 'Ia':[], 'II':[], 'Ibc':[] }
    bestlogz = -np.inf
    for modelsource in allmodelnames :
        if verbose >1 :
            dt = time.time()-tstart
            print('------------------------------')
            print( "model: %s  dt=%i sec"%( modelsource, dt))
        sn, res, fit, priorfn = get_evidence(
            sn, modelsource=modelsource, zhost=zhost, zhosterr=zhosterr,
            t0_range=t0_range, zminmax=zminmax,
            nobj=nobj, maxiter=maxiter, verbose=max(0,verbose-1) )

        if nsteps_pdf :
            pdf =  get_marginal_pdfs( res, nbins=nsteps_pdf, verbose=max(0,verbose-1) )
        else :
            pdf = None
        outdict[modelsource] = {'sn':sn, 'res':res, 'fit':fit, 'pdf':pdf, 'priorfn':priorfn }

        if res.logz>bestlogz :
            outdict['bestmodel'] = modelsource
            bestlogz = res.logz

        # multiply the model evidence by the sub-type prior
        if modelsource in iimodelnames :
            logprior = logpriordict['ii']
            logz['II'].append( logprior + res.logz )
        elif modelsource in ibcmodelnames :
            logprior = logpriordict['ibc']
            logz['Ibc'].append( logprior + res.logz )
        elif modelsource in iamodelnames :
            logprior = logpriordict['ia']
            logz['Ia'].append( logprior + res.logz )

    # sum up the evidence from all models for each sn type
    logztype = {}
    for modelsource in ['II','Ibc','Ia']:
        logztype[modelsource] = logz[modelsource][0]
        for i in range( 1,len(logz[modelsource]) ):
            logztype[modelsource] = np.logaddexp( logztype[modelsource], logz[modelsource][i] )

    # define the total evidence (final denominator in Bayes theorem) and then
    # the classification probabilities
    logzall = np.logaddexp( np.logaddexp( logztype['Ia'], logztype['Ibc']), logztype['II'] )
    pIa = np.exp( logztype['Ia'] - logzall )
    pIbc = np.exp( logztype['Ibc'] - logzall )
    pII = np.exp( logztype['II'] - logzall )

    outdict['pIa'] = pIa
    outdict['pIbc'] = pIbc
    outdict['pII'] = pII

    outdict['logztype'] = logztype
    outdict['logzall'] = logzall

    return( outdict )



def medbandClassTestStone( sfxlist=['J','H','JNH','WIDE','MB','MBJ','MBH','MBIR','ALL'],
                           templateset='SNANA', usehostprior=False):
    """  Run the classifier on a series of versions of the SN light curve
    holding a subset of the data.
    :param sn: name of the SN to evaluate
    :return:
    for sfx in  :

    """

    outdictlist = []
    for sfx in sfxlist :
        datfile = 'HST_CANDELS4_stone%s.sncosmo.dat'%(sfx)
        sn = ascii.read( datfile, format='commented_header', header_start=-1, data_start=0 )
        if usehostprior :
            outdict = classify( sn,  zhost=1.8, zhosterr=0.02, t0_range=[56460,56500],zminmax=[1.7,1.9], verbose=False, templateset=templateset )
            fout = open('stone_classtest_withHostPrior.dat','a')
        else :
            outdict = classify( sn,  zhost=2.0, zhosterr=0.8, zminmax=[1.2,2.8], t0_range=[56460,56500], verbose=False, templateset=templateset )
            fout = open('stone_classtest_noHostPrior.dat','a')
        print >> fout, 'stone %s : PIa=%.3e'%(sfx, outdict['pIa'])
        fout.close()

        outdictlist.append( outdict )
    return( outdictlist )




def medbandClassTestColfax( sfxlist=['J','H','JNH','WIDE','MB','MBJ','MBH','MBIR','ALL'],
                           templateset='SNANA', usehostprior=True):
    """  Run the classifier on a series of versions of the SN light curve
    holding a subset of the data.
    :param sn: name of the SN to evaluate
    :return:
    for sfx in  :

    """

    outdictlist = []
    for sfx in sfxlist :
        datfile = 'HST_CANDELS4_colfax%s.sncosmo.dat'%(sfx)
        sn = ascii.read( datfile, format='commented_header', header_start=-1, data_start=0 )
        if usehostprior :
            outdict = classify( sn,  zhost=2.26, zhosterr=0.09, zminmax=[1.8,2.6], t0_range=[56055,56095], verbose=False, templateset=templateset )
            fout = open('colfax_classtest_withHostPrior.dat','a')
        else :
            outdict = classify( sn,  zhost=2.0, zhosterr=0.8, zminmax=[1.2,2.8], t0_range=[56055,56095], verbose=False, templateset=templateset )
            fout = open('colfax_classtest_noHostPrior.dat','a')
        print >> fout, 'colfax %s : PIa=%.3e'%(sfx, outdict['pIa'])
        fout.close()

        outdictlist.append( outdict )
    return( outdictlist )



def medbandClassTestBush( sfxlist=['J','H','IR','WIDE','MB','MBJ','MBH','MBIR','ALL'],
                           templateset='SNANA', usehostprior=False):
    """  Run the classifier on a series of versions of the SN light curve
    holding a subset of the data.
    :param sn: name of the SN to evaluate
    :return:
    for sfx in  :

    """

    outdictlist = []
    for sfx in sfxlist :
        datfile = 'HST_CANDELS4_bush%s.sncosmo.dat'%(sfx)
        sn = ascii.read( datfile, format='commented_header', header_start=-1, data_start=0 )
        if usehostprior :
            outdict = classify( sn,  zhost=1.15, zhosterr=0.4, t0_range=[55760,55840],zminmax=[0.7,1.6], verbose=False, templateset=templateset )
            fout = open('bush_classtest_withHostPrior.dat','a')
        else :
            outdict = classify( sn,  zhost=1.5, zhosterr=1.0, zminmax=[0.5,2.0], t0_range=[55760,55840], verbose=False, templateset=templateset )
            fout = open('bush_classtest_noHostPrior.dat','a')
        print >> fout, 'bush %s : PIa=%.3e'%(sfx, outdict['pIa'])
        fout.close()

        outdictlist.append( outdict )
    return( outdictlist )





