
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


def mcSim( sntype, z_range=[1.8,2.2], nsim=100, perfect=True ):
    """ Run a monte carlo sim using sncosmo to simulate <nsim> SNe
    of the given <sntype> over the given <z_range>.

    Simulates Type Ia SNe with the SALT2 model, and CC SNe with
    the SNANA CC templates.

    Set perfect=True for noiseless "observations" of the simulated SNe.
    :return:
    """
    from astropy.table import Table
    from astropy import cosmology
    import sncosmo
    from numpy.random import normal, uniform, choice
    import numpy as np

    # List of all available sncosmo source models for this SN type :

    # Make a list of all the unique sncosmo source models available,
    # and assign a relative probability that any given simulated SN of this
    # type (CC or Ia) belongs to that subclass
    if sntype.lower()=='cc' :
        sourcelist = ccSubClassDict.keys()
        subclasslist = np.array([ ccSubClassDict[source] for source in sourcelist ])
        ninthissubclass = np.array([ len(np.where(subclasslist==subclass)[0])
                                     for subclass in subclasslist ], dtype=float)
        psourcelist = np.array([ ccSubClassProbs[subclass]
                                 for subclass in subclasslist ]) / ninthissubclass
        psourcelist /= psourcelist.sum()
    else :
        # No sub-class divisions for SNIa
        sourcelist = np.array(['salt2-extended'])
        subclasslist = np.array( ['Ia'] )
        psourcelist = np.array( [1] )

    # load the O'Donnell 1994 dust model
    dust = sncosmo.OD94Dust()

    # Define an sncosmo SN model for each available source
    modellist = [ sncosmo.Model(source=source, effects=[dust],
                                effect_names=['host'], effect_frames=['rest'])
                  for source in sourcelist ]

    # medium band at peak observation set :
    observations = Table({'band': ['wfc3f127m','wfc3f139m','wfc3f153m',
                                   'wfc3f125w','wfc3f140w','wfc3f160w',
                                   ],
                          'time': [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ],
                          'zp': [26.47,26.49,26.70, 28.02,28.48,28.19,],
                          'zpsys': ['ab','ab','ab','ab','ab','ab'],
                          'gain': [1., 1., 1.,1.,1.,1.],
                          'skynoise': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                          })

    # Define a cosmology
    cosmo = cosmology.FlatLambdaCDM(Om0=0.3, H0=70.)

    # For each simulated SN, draw random Av from distributions
    # as defined in Rodney et al 2014a :
    #   For SN Ia :  P(Av) = exp(-Av/0.33)
    #   For CC SN :  P(Av) = 4 * gauss(0.6) + exp(-Rv/1.7)
    if sntype=='Ia':
        tau,sigma,R0 = 0.33, 0, 0
    else :
        tau,sigma,R0 = 1.7, 0.6, 4
    Avlist = mcsample( pAv, nsim, tau=tau, sigma=sigma, R0=R0 )

    # For each simulated SN, draw a random Rv from a normal
    # distribution centered on 3.1 with width 0.5
    Rvlist = normal( 3.1, 0.5, nsim )
    Rvlist = np.where( Rvlist>0, Rvlist, 0.01 )

    # Convert Av and Rv to E(B-V) :
    # Rv = Av/EBV ==> EBV = Av/Rv
    EBVlist = Avlist / Rvlist

    # TODO : draw the redshifts with a metropolis-hastings sampler to match
    #  a distribution defined based on the expected SN rate

    # Disabled : uniform redshift spacing
    # zlist = np.linspace( z_range[0], z_range[1], nsim )

    # Draw a random redshift from a uniform distribution
    zlist = uniform( low=z_range[0], high=z_range[1], size=nsim )

    simset = []
    for isim in range(nsim):
        # Randomly draw an sncosmo model from the available list, according to
        # the predefined probability list, setting the SN sub-class for this
        # simulated SN
        imodel = choice( np.arange(len(modellist)), replace=True, p=psourcelist )
        model = modellist[imodel]
        subclass = subclasslist[imodel]

        z=zlist[isim]
        EBV = EBVlist[isim]
        Rv = Rvlist[isim]

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
        model.set_source_peakabsmag( MR, 'bessellr', 'vega', cosmo=cosmo)

        if subclass =='Ia' :
            x0 = model.get('x0')
            # TODO : use bifurcated gaussians for more realistic x1,c dist'ns
            x1 = normal(0., 1.)
            c = normal(0., 0.1)
            # TODO : draw a random t0
            modelparams = {'z':z, 't0':0, 'x0':x0, 'x1':x1, 'c':c, 'hostebv':EBV, 'hostr_v':Rv}
        else :
            amplitude = model.get('amplitude')
            modelparams = {'z':z, 't0':0, 'amplitude':amplitude }

        # Generate one simulated SN:
        sn = sncosmo.realize_lcs(observations, model, [ modelparams ],
                                 thresh=None, perfect=perfect )[0]
        simset.append( sn )
    return( simset )

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

    i127 = np.where( flt=='wfc3f127m' )
    i139 = np.where( flt=='wfc3f139m' )
    i153 = np.where( flt=='wfc3f153m' )
    i125 = np.where( flt=='wfc3f125w' )
    i140 = np.where( flt=='wfc3f140w' )
    i160 = np.where( flt=='wfc3f160w' )

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


