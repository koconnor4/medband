from pylab import *
from numpy import *

def gethubbledat( datfile='union2.dat') :
    from astropy.io import ascii
    from astropy.table import Table
    from pytools import cosmo
    import exceptions
    import os

    datfile = os.path.abspath(datfile)
    if not os.path.isfile( datfile ):
        thisfile = sys.argv[0]
        if 'ipython' in thisfile :
            thisfile = __file__
        thisdir = os.path.dirname( thisfile )
        datfile = os.path.join( thisdir, datfile )
    if not os.path.isfile( datfile ):
        raise exceptions.RuntimeError("No such file %s"%datfile)

    # read in redshifts and  distances
    hubbledat = ascii.read( datfile, format='commented_header', header_start=-1, data_start=0)
    if os.path.basename(datfile).startswith('union'):
        H0 = 70.0
        Om = 0.295
        MB = -19.175  # absolute SNIa mag assuming H0=70
        muLCDM = cosmo.mu( hubbledat['z'], H0=H0, Om=Om, Ode=(1-Om) )
        mures = hubbledat['mu'] - muLCDM
        hubbledat['muLCDM'] = muLCDM
        hubbledat['mures'] = mures
        return( hubbledat, H0, Om, MB  )

    elif os.path.basename(datfile).startswith('ps1') :
        #H0 = 70.0
        #Om = 0.283671
        #MB = -19.363

        H0 = 70.88
        Om = 0.272
        MB = -19.34
        alpha=0.147
        beta=3.13

        z = hubbledat['z']
        x1 = hubbledat['x1']
        c = hubbledat['c']
        mB = hubbledat['mB']
        fitprob = hubbledat['FITPROB']
        mupull = hubbledat['MUPULL']
        igood = np.where( (np.abs(x1)<3) & (np.abs(c)<0.3) & (fitprob > 0.001) & (np.abs(mupull)<4) )

        # igood = np.where( (np.abs(x1)<3) & (np.abs(c)<0.3) & (np.abs(mupull)<4) )

        code2survey = {1:'sdss',4:'snls',15:'ps1',
                       5:'low-z',50:'low-z',53:'low-z',54:'low-z',
                       101:'hstacs', 111:'hstacs', 116:'hstacs',
                       106:'hstacs',107:'hstwfc3'}

        hubbledatgood = Table()
        hubbledatgood['name'] = hubbledat['CID'][igood]
        # hubbledatgood['mu'] = hubbledat['MU'][igood]
        # hubbledatgood['mures'] = hubbledat['MURES'][igood]
        hubbledatgood['mberr'] = hubbledat['mBERR'][igood]
        hubbledatgood['x1err'] = hubbledat['x1ERR'][igood]
        hubbledatgood['cerr'] = hubbledat['cERR'][igood]
        hubbledatgood['x1'] = hubbledat['x1'][igood]
        hubbledatgood['c'] = hubbledat['c'][igood]
        muerr = np.sqrt( hubbledat['mBERR'][igood]**2 + (alpha*hubbledatgood['x1err'])**2 + \
                       (beta*hubbledat['cERR'][igood])**2 + 0.12**2 + (0.093*z[igood])**2 )
        hubbledatgood['muerr'] = muerr
        hubbledatgood['muerr2'] = muerr

        hubbledatgood['mu'] = mB[igood] + alpha*x1[igood] - beta*c[igood] - MB

        muLCDM = cosmo.mu( z[igood], H0=H0, Om=Om, Ode=(1-Om) )
        hubbledatgood['mures'] = hubbledatgood['mu'] - muLCDM

        # hubbledatgood['muerrDS'] = hubbledat['MUERR'][igood]
        # hubbledatgood['muerr'] = hubbledat['MUERR'][igood]

        hubbledatgood['z'] = hubbledat['z'][igood]
        hubbledatgood['zerr'] = hubbledat['zERR'][igood]
        hubbledatgood['survey'] = np.array([code2survey[idsur] for idsur in hubbledat['IDSURVEY'][igood]])
        hubbledatgood['spec'] = np.ones( len(hubbledatgood) )
        return( hubbledatgood, H0, Om, MB )

    elif os.path.basename(datfile).startswith('jla') :
        H0 = 70.0
        Om = 0.295
        muLCDM = cosmo.mu( hubbledat['z'], H0=H0, Om=Om, Ode=(1-Om) )

        alpha = 0.141
        beta = 3.101
        MB = -19.05
        muJLA = hubbledat['mb'] - MB + alpha * hubbledat['x1'] - beta * hubbledat['c']

        mures = muJLA - muLCDM
        hubbledat['muLCDM'] = muLCDM
        hubbledat['mures'] = mures
        hubbledat['muerr'] = np.sqrt( hubbledat['mberr']**2 + (alpha*hubbledat['x1err'])**2 +
                                      (beta*hubbledat['cerr'])**2 )
        hubbledat['spec'] = np.ones( len( mures ))
        hubbledat['zerr'] = np.zeros( len( mures ))
        code2survey = {2:'sdss',1:'snls',3:'low-z',
                       4:'hstacs',5:'hstwfc3'}
        hubbledat['survey'] = np.array([code2survey[idsur] for idsur in hubbledat['set']])
        return( hubbledat, H0, Om, MB )
    elif os.path.basename(datfile).startswith('rodney'):
        H0 = 70.0
        Om = 0.277
        MB = -19.223  # absolute SNIa mag assuming H0=70
        alpha = 0.141
        beta = 3.101
        muSR = hubbledat['mb'] - MB + alpha * hubbledat['x1'] - beta * hubbledat['c']
        muLCDM = cosmo.mu( hubbledat['z'], H0=H0, Om=Om, Ode=(1-Om) )
        hubbledat['mu'] = muSR
        mures = muSR - muLCDM
        hubbledat['muLCDM'] = muLCDM
        hubbledat['mures'] = mures
        hubbledat['muerr'] = np.sqrt( hubbledat['mberr']**2 + (alpha*hubbledat['x1err'])**2 +
                                      (beta*hubbledat['cerr'])**2 )
        return( hubbledat, H0, Om, MB  )


def get_visibility_time_grid( deteff_threshold=0.8, x1=-3, c=0.3 ):
    """ Define the line in a hubble residuals diagram where our control
    time is < 50%
    :param m50pct:
    :return:
    """
    import sncosmo
    from pytools import cosmo

    if deteff_threshold == 0.8 :
        mlim125 = 25.84
        mlim160=25.85
    elif deteff_threshold == 0.5 :
        mlim125 = 25.5 + 0.9
        mlim160 = 25.2 + 1.25

    H0 = 72.0
    Om = 0.277
    MB = -19.175  # absolute SNIa mag assuming H0=70

    model = sncosmo.Model(source='salt2')
    model.source.set_peakmag( 0, 'bessellb', 'ab' )
    x0_MB0 = model.get('x0')

    zlist = np.arange( 1.0, 2.7, 0.1 )

    dmlist = np.arange( -1.5,1.5,0.1 )
    tdetgrid = []
    for z in zlist :

        mB = MB + cosmo.mu(z, H0=H0, Om=Om)
        x0 = x0_MB0 * 10**(-0.4*mB)

        model.set( z=z, x1=x1, c=c, x0=x0 )

        tobs =  np.arange( -14*(1+z), 60*(1+z), 1.0 )
        magH = model.bandmag('f160w', 'ab', tobs )
        magJ = model.bandmag('f125w', 'ab', tobs )

        tdetlist = [ len(np.where( (magH+deltamag<=mlim160) & (magJ+deltamag<=mlim125) )[0]) for deltamag in dmlist ]
        tdetgrid.append( tdetlist )
    tdetgrid = np.array(tdetgrid).reshape( [len(zlist), len(dmlist)] ).swapaxes(0,1)

    return( zlist, dmlist, tdetgrid )

def get_visibility_time_line( tvislim=50, deteff_threshold=0.8, x1=-3, c=0.3,
                              zrange=[1.0,2.7], tstep=1.0 ):
    """ Define the line in a hubble residuals diagram where our control
    time is < 50%
    :param m50pct:
    :return:
    """
    import sncosmo
    from pytools import cosmo

    if deteff_threshold == 0.8 :
        mlim125 = 25.84
        mlim160=25.85
    elif deteff_threshold == 0.5 :
        mlim125 = 26.4 # 25.5 + 0.9
        mlim160 = 26.45 # 25.2 + 1.25

    H0 = 72.0
    Om = 0.277
    MB = -19.175  # absolute SNIa mag assuming H0=70

    model = sncosmo.Model(source='salt2')
    model.source.set_peakmag( 0, 'bessellb', 'ab' )
    x0_MB0 = model.get('x0')

    zlist = np.arange(zrange[0], zrange[1], 0.1 )

    dmlist = np.arange( -1.5,1.5,0.1 )
    dmmaxlist = []
    for z in zlist :

        mB = MB + cosmo.mu(z, H0=H0, Om=Om)
        x0 = x0_MB0 * 10**(-0.4*mB)

        model.set( z=z, x1=x1, c=c, x0=x0 )

        tobs =  np.arange( -14*(1+z), 50*(1+z), tstep )
        magH = model.bandmag('f160w', 'ab', tobs )
        magJ = model.bandmag('f125w', 'ab', tobs )
        for deltamag in dmlist :
            tvis = tstep * len(np.where((magH+deltamag<=mlim160) & (magJ+deltamag<=mlim125))[0])
            if tvis < tvislim : break
        dmmaxlist.append( deltamag )
    return( zlist, dmmaxlist )



def mkresidplot( datfile='ps1union.dat', ageaxis=True,
                 highz='Malmquist', highzlabels=True,
                 lensfactor=0.093, presfig=False, bigtext=False ) :
    """ construct the hubble residual plot
    """
    from pytools import plotsetup, cosmo, colorpalette as cp
    from matplotlib import patches
    from astropy.io import ascii
    import numpy as np
    from matplotlib import pyplot as pl
    from scipy.interpolate import interp1d
    if presfig :
        plotsetup.presfig( wide=True )
    else :
        plotsetup.fullpaperfig( [8.5,4.5], bigtext=bigtext)
    clf()

    highz= highz.lower()
    if highz=='malmquist':
        zlim=[0,2.78]
    elif highz=='rodney':
        zlim=[0,3.49]
        if presfig : zlim=[0,3.78]
    else :
        zlim=[0,2.45]

    fig = gcf()
    if presfig :
        bottomAxeslim = [ 0.13, 0.18, 0.85, 0.65 ]
    else :
        bottomAxeslim = [ 0.09, 0.12, 0.88, 0.75 ]


    ax1 = fig.add_axes( bottomAxeslim, frameon=True, alpha=0.0, label='axes1')
    ax1top = ax1.twiny()

    # read in redshifts and  distances
    hubbledat, H0, Om, MB = gethubbledat( datfile )

    if lensfactor :
        # muerr = 0.093z : Holz & Linder, 2005
        # muerr = 0.055z : Jonsson+ 2012
        hubbledat['muerr'] = np.sqrt( hubbledat['muerr']**2 + (lensfactor*hubbledat['z'])**2 )

    # plot ground data
    colordict = {'low-z':cp.purple,
                 'sdss':cp.lightblue,
                 'ps1':cp.darkblue,
                 'snls':cp.darkgreen,
                 'hstacs':cp.darkorange,
                 'hstwfc3':cp.medred,
                 'hstwfc3a':cp.medred,
                 'hstwfc3b':cp.darkred,
                 'jwst':cp.darkgrey,
    }

    elw=0.75
    cs=0
    for survey in ['low-z','sdss','snls','ps1','hstacs']:
        isurvey = where( (hubbledat['survey']==survey) & (hubbledat['spec']>0) )

        ax1.errorbar( hubbledat['z'][isurvey], hubbledat['mures'][isurvey],
                      hubbledat['muerr'][isurvey], hubbledat['zerr'][isurvey],
                      marker='o', ls=' ',mew=0, markersize=5,
                      color=colordict[survey], mfc=colordict[survey],
                      ecolor=colordict[survey],
                      elinewidth=elw, capsize=cs, alpha=0.7 )

    if highz :
        # get the high-z SN fits from Dan
        labeldict = {#'primo':'Rodney+ 2012',
                     #'wilson':'Jones, Rodney+ 2013',
                     'stone':'GND13Sto',
                     'colfax':'GND12Col',
                     }
        if highz=='rodney':
            labeldict = {'primo':'Rodney+ 2012',
                         'wilson':'Jones, Rodney+ 2013',
                         'stone':'Rodney+ in prep',
                         'colfax':'Rodney+ in prep',
                         }

        for survey, namelist in zip( ['hstwfc3a','hstwfc3b'],
                                     [['primo','wilson'],['stone','colfax']] ):
            color=colordict[survey]
            for name in namelist :
                if name not in hubbledat['name'] : continue
                iname = np.where( hubbledat['name'] == name )[0][0]
                z,zerr = hubbledat['z'][iname], hubbledat['zerr'][iname]
                mu,muerr = hubbledat['mures'][iname],hubbledat['muerr'][iname]
                if name=='colfax' :
                    # Not applying a lensing correction for mu=1.04 +-0.03  :  0.0426 +- 0.031 mag
                    icol = np.where( [n.startswith('col') for n in hubbledat['name']] )[0]
                    mucol = interp1d( hubbledat['z'][icol], hubbledat['mures'][icol] )
                    ax1.plot( [2.16,2.19,2.21,2.24,2.26,2.28],
                              mucol([2.16,2.19,2.21,2.24,2.26,2.28]),
                              color=color, lw=0.75, ls='-' )
                    ax1.errorbar( 2.26, mucol(2.26), muerr,
                                  marker='D', ls=' ',mew=0, markersize=8,
                                  color=color, mfc=color, ecolor=color,
                                  elinewidth=elw, capsize=cs, alpha=1 )
                    z = 2.26
                    mu = mucol(2.26)
                else :
                    # Not applying a lensing correction for Stone
                    #if name=='stone':
                    #    mu = mu + 0.0215
                    #    muerr = np.sqrt( muerr**2 + 0.02**2 )
                    ax1.errorbar( z, mu, muerr, zerr,
                                  marker='D', ls=' ',mew=0, markersize=8,
                                  color=color, mfc=color, ecolor=color,
                                  elinewidth=elw, capsize=cs, alpha=1 )

                if highzlabels and name in labeldict :
                    # ax1.text( z-0.035, mu+0.07, labeldict[name],
                    if name=='primo':
                        ax1.text( z-0.05, mu+0.035, labeldict[name],
                                  color=color, # backgroundcolor='w',
                                  va='bottom', ha='center', rotation=90 )
                    else :
                        ax1.text( z, mu+muerr+0.07, labeldict[name],
                                  color=color, # backgroundcolor='w',
                                  va='bottom', ha='center', rotation=90 )



    # plot the lcdm cosmology
    ax1.axhline( 0, color='0.4', ls='-' )

    # plot extreme w'
    zLCDM = np.arange( 0.01, 5, 0.05)
    muLCDM = cosmo.mu( zLCDM, H0=H0, Om=Om, Ode=(1-Om), w0=-1., wa=0)
    #muextremeDE = cosmo.mu( z, H0=_H0, Om=_Om, Ode=(1-_Om), w0=-0.7, wa=-2)
    mu1 = cosmo.mu( zLCDM, H0=H0, Om=Om-0.015, Ode=(1-Om+0.015), w0=-1.2, wa=0.7)
    mu2 = cosmo.mu( zLCDM, H0=H0, Om=Om+0.055, Ode=(1-Om-0.055), w0=-0.65, wa=-2.2)
    ax1.fill_between( zLCDM, mu1-muLCDM, mu2-muLCDM, color='k', alpha=0.2, lw=0, zorder=10 )

    if not presfig :
        ax1.text( 0.08, -0.8,  'low-z',    fontsize=(presfig and 'small') or 10.5, ha='center', color=colordict['low-z']     )
        ax1.text( 0.25, -0.85, 'SDSS',     fontsize=(presfig and 'small') or 10.5, ha='center', color=colordict['sdss']      )
        ax1.text( 0.5 , -0.9,  'PS1',      fontsize=(presfig and 'small') or 10.5, ha='center', color=colordict['ps1']       )
        ax1.text( 0.75, -0.95, 'SNLS',     fontsize=(presfig and 'small') or 10.5, ha='center', color=colordict['snls']      )
        ax1.text( 1.25, -0.8, 'HST+ACS',   fontsize=(presfig and 'small') or 10.5, ha='center', color=colordict['hstacs']    )
        ax1.text( 1.95, -0.95, 'HST+WFC3', fontsize=(presfig and 'small') or 10.5, ha='center', color=colordict['hstwfc3b']  )

        lowzRect = patches.Rectangle( [0.0,-0.85], 0.08, 0.03, fill=True, lw=0, alpha=0.3, color=colordict['low-z']     )
        sdssRect = patches.Rectangle( [0.06,-0.9], 0.34, 0.03, fill=True, lw=0, alpha=0.3, color=colordict['sdss']      )
        ps1Rect = patches.Rectangle(  [0.1,-0.95], 0.5,  0.03, fill=True, lw=0, alpha=0.3, color=colordict['ps1']       )
        snlsRect = patches.Rectangle( [0.15,-1.0], 0.9,  0.03, fill=True, lw=0, alpha=0.3, color=colordict['snls']      )
        acsRect = patches.Rectangle(  [0.60,-0.85], 1.15, 0.03, fill=True, lw=0, alpha=0.3, color=colordict['hstacs']    )
        wfc3Rect = patches.Rectangle( [1.20,-1.0], 1.1,  0.03, fill=True, lw=0, alpha=0.3, color=colordict['hstwfc3b']  )

        ax1.add_patch( lowzRect )
        ax1.add_patch( sdssRect )
        ax1.add_patch( ps1Rect )
        ax1.add_patch( snlsRect )
        ax1.add_patch( acsRect )
        ax1.add_patch( wfc3Rect )

    if highz=='rodney':
        if not presfig :
            #ffsnRect = patches.Rectangle( [2.30,-1.0], 0.7,  0.03, fill=False,lw=1, alpha=1.0, color=colordict['hstwfc3b']  )
            jwstRect = patches.Rectangle( [2.00,-0.85],2.0,  0.03, fill=False,lw=1, alpha=1.0, color=colordict['jwst']  )
            #ax1.add_patch( ffsnRect )
            ax1.add_patch( jwstRect )
            #ax1.text( 2.5,  -0.95, '+ lensing',fontsize=10.5, ha='center', color=colordict['hstwfc3b']  )
            ax1.text( 3.0,  -0.78, 'JWST', fontsize=(presfig and 'medium') or 10.5, ha='center', color=colordict['jwst']  )

        color = cp.bluegray
        ax1.arrow( 2.5, 0.15, 0.3, 0.0, fc=color, ec=color, head_width=0.05, head_length=0.0, overhang=0.0 )
        ax1.plot( 2.5, 0.15, marker='o', mfc='w', mec=color )
        ax1.text( 2.49, 0.19, '0.3 Gyr', ha='left', va='bottom', color=color)
        ax1.text( 2.82, 0.15, 'z=2.8', ha='left', va='center', color=color)

        color = cp.teal
        ax1.arrow( 2.5, 0.45, 0.55, 0.0, fc=color, ec=color, head_width=0.05, head_length=0.08, overhang=0.25, length_includes_head=True )
        ax1.plot( 2.5, 0.45, marker='o', mfc='w', mec=color )
        ax1.text( 2.6, 0.48, '1 Gyr', ha='left', va='bottom', color=color)
        ax1.text( 3.08, 0.45, 'z=3.8', ha='left', va='center', color=color)

        color = cp.coral
        ax1.arrow( 2.5, 0.7, 0.75, 0.0, fc=color, ec=color, head_width=0.05, head_length=0.08, overhang=0.25, length_includes_head=True  )
        ax1.plot( 2.5, 0.7, marker='o', mfc='w', mec=color )
        ax1.text( 2.7, 0.72, '2 Gyr', ha='left', va='bottom', color=color)
        ax1.text( 3.28, 0.7, 'z=8.3', ha='left', va='center', color=color)

        if not presfig :
            ax1.text( 2.5, 0.7, 'z=2.5', ha='center', va='bottom', color='0.4' )
            ax1.text( 2.5, 0.8, 'explosion', ha='center', va='bottom', color='0.3' )#, fontsize=10.5 )
            ax1.text( 2.85, 0.8, 'delay',     ha='center', va='bottom', color='0.3' )#, fontsize=10.5 )
            ax1.text( 3.15, 0.8, 'formation', ha='left',   va='bottom', color='0.3' )#, fontsize=10.5 )

    if presfig :
        ax1.text( 0.82,  0.46, 'flat $\Lambda$CDM',
                  color='0.1', ha='right', va='top',
                  transform=ax1.transAxes )# , fontsize=10.5, backgroundcolor='w' )
    else :
        ax1.text( 0.99,  0.56, 'flat $\Lambda$CDM',
                  color='0.1', ha='right', va='top',
                  transform=ax1.transAxes )# , fontsize=10.5, backgroundcolor='w' )

    ax1.text( 0.99,  0.46, '$\pm$95\% conf. \non $w_0 , w_a$',
              color='0.4', ha='right', va='top',
              transform=ax1.transAxes )# , fontsize=10.5,backgroundcolor='w' )

    # ax1.text( 0.95, 0.95, 'Open Symbols:\n no spec. confirmation', color='0.5', transform=ax1.transAxes, ha='right',va='top' )

    ax1.set_xlabel('Redshift',fontsize='large')
    # ax1.set_ylabel('$\mu_{observed}$ - $\mu_{\Lambda CDM}$')
    if highz=='rodney' :
        ax1.set_ylabel(r'$m_B^{\prime} - m_{\Lambda CDM}$', fontsize='large')
    else :
        ax1.set_ylabel(r'$m_B + \alpha x_1 - \beta c - m_{\Lambda CDM}$', fontsize='large')
    if presfig :
        ax1.yaxis.set_label_coords( -0.08, 0.65 )

    # mulim = array( [33,46.5] )
    muresidlim = array( [-1.05,1.05] )
    # zlim  = array( [0.0,3.49] )
    ax1.xaxis.set_ticks( [0.1,] + np.arange(0.5,zlim[-1],0.5).tolist() )
    ax1.set_xlim( zlim )
    ax1.set_ylim( muresidlim )


    from pytools import colorpalette as cp

    if highz.lower()=='malmquist' :
        # z, dm, tvis = get_visibility_time_grid( 0.5, x1=0., c=0. )
        # contourlevels = [50]
        # ax1.contour( z, dm, tvis, levels=contourlevels, colors=cp.darkbluegray, linestyles='dashed' )
        z5050, dm5050 = get_visibility_time_line( tvislim=50, deteff_threshold=0.5, x1=0, c=0, zrange=[1.0,3.0] )
        ax1.fill_between( z5050, dm5050, np.ones(len(z5050))*muresidlim[1]+0.05, color=cp.bluegray, zorder=-200, alpha=0.3)
        # ax1.plot( z5050, dm5050, color='0.8')

        z050, dm050 = get_visibility_time_line( tvislim=0.01, deteff_threshold=0.5, x1=0, c=0, zrange=[1.7,3.0] )
        ax1.fill_between( z050, dm050, np.ones(len(z050))*muresidlim[1]+0.05, color=cp.darkbluegray, zorder=-100, alpha=0.3)

        ax1.text( 1.53, 0.92, 't$_{\\rm vis}<$50 days', ha='left', va='top',  color=cp.darkgrey )
        ax1.text( 2.75, 0.92, 'undetectable', ha='right', va='top', color=cp.darkgrey )

    if highz=='evolution' :
        # plot SNIa evolution
        minmassdict = progenitorMinMassLines( zform=11 )

        z = minmassdict['z']
        for key,color in zip( ['mm25','mm15','mm04'], ['r','g','b'] ):
            minmass = minmassdict[key]
            mavez = Mave( _MAXMASS, minmass )
            ax1.plot( z, (mavez-mavez[0]) * _DMUDMASS, ls='-', color=color )



        ax1.text( 2.5,  -0.95, '+ lensing',fontsize=10.5, ha='center', color=colordict['hstwfc3b']  )
        ax1.text( 3.0,  -0.78, 'JWST',     fontsize=10.5, ha='center', color=colordict['jwst']  )

        color = cp.bluegray
        ax1.arrow( 2.5, 0.25, 0.3, 0.0, fc=color, ec=color, head_width=0.05, head_length=0.0, overhang=0.0 )
        ax1.plot( 2.5, 0.25, marker='o', mfc='w', mec=color )
        ax1.text( 2.55, 0.27, '0.3 Gyr', ha='left', va='bottom', color=color)
        ax1.text( 2.82, 0.25, 'z=2.8', ha='left', va='center', color=color)

        color = cp.teal
        ax1.arrow( 2.5, 0.45, 0.55, 0.0, fc=color, ec=color, head_width=0.05, head_length=0.08, overhang=0.25, length_includes_head=True )
        ax1.plot( 2.5, 0.45, marker='o', mfc='w', mec=color )
        ax1.text( 2.65, 0.47, '1 Gyr', ha='left', va='bottom', color=color)
        ax1.text( 3.08, 0.45, 'z=3.8', ha='left', va='center', color=color)

        color = cp.coral
        ax1.arrow( 2.5, 0.65, 0.75, 0.0, fc=color, ec=color, head_width=0.05, head_length=0.08, overhang=0.25, length_includes_head=True  )
        ax1.plot( 2.5, 0.65, marker='o', mfc='w', mec=color )
        ax1.text( 2.75, 0.67, '2 Gyr', ha='left', va='bottom', color=color)
        ax1.text( 3.28, 0.65, 'z=8.3', ha='left', va='center', color=color)

        ax1.text( 2.5, 0.7, 'z=2.5', ha='center', va='bottom', color='0.4' )
        ax1.text( 2.5, 0.8, 'explosion', ha='center', va='bottom', color='0.3', fontsize=10.5 )
        ax1.text( 2.85, 0.8, 'delay', ha='center', va='bottom', color='0.3',    fontsize=10.5  )
        ax1.text( 3.15, 0.8, 'formation', ha='left', va='bottom', color='0.3',  fontsize=10.5 )


    if ageaxis :
        if highz=='rodney':
            # ax1top.set_xlabel( 'age of Universe (Gyr)',fontsize='large',x=0.15,ha='left')
            ageticks = np.array( [13,  9,  6, 4, 3, 2 ] )
        else :
            ageticks = np.array( [13,  9,  6, 4, 3 ] )
        ax1top.set_xlabel( 'Age of Universe (Gyr)',fontsize='large')
        ax1top.set_xlim( zlim )
        ztickloctop = cosmo.zfromt( ageticks )
        ax1top.xaxis.set_ticks( ztickloctop )
        ax1top.xaxis.set_ticklabels( ageticks )

    pl.draw()

def mkhubblefig( datfile='ps1union.dat', ageaxis=True,
                 highz='rodney', bigtext=False ) :
    """ construct the hubble residual plot
    """
    from pytools import plotsetup, cosmo, colorpalette as cp
    from matplotlib import patches, ticker
    from mpl_toolkits.axes_grid1 import host_subplot
    from astropy.io import ascii
    import numpy as np
    from matplotlib import pyplot as pl
    from scipy.interpolate import interp1d
    plotsetup.presfig(wide=True)
    pl.clf()
    fig = pl.gcf()

    if highz=='Malmquist':
        zlim=[0,2.58]
    elif highz=='rodney':
        zlim=[0,3.49]
        if bigtext : zlim=[0,3.53]
    else :
        zlim=[0,2.45]

    # bottomAxeslim = [ 0.14, 0.14, 0.72, 0.72 ]
    # ax1 = fig.add_axes( bottomAxeslim, frameon=True, alpha=0.0, label='axes1')
    ax1 = host_subplot(111)
    fig.subplots_adjust( left=0.11, bottom=0.16, right=0.89, top=0.84 )

    # read in redshifts and  distances
    hubbledat, H0, Om, MB = gethubbledat( datfile )

    lensfactor=0.093
    hubbledat['muerr'] = np.sqrt( hubbledat['muerr']**2 + (lensfactor*hubbledat['z'])**2 )

    # plot ground data
    colordict = {'low-z':cp.purple,
                 'sdss':cp.lightblue,
                 'ps1':cp.darkblue,
                 'snls':cp.darkgreen,
                 'hstacs':cp.darkorange,
                 'hstwfc3a':cp.medred,
                 'hstwfc3b':cp.darkred,
                 'jwst':cp.darkgrey,
    }

    elw=0.75
    cs=0
    for survey in ['low-z','sdss','snls','ps1','hstacs','hstwfc3']:
        isurvey = where( (hubbledat['survey']==survey) & (hubbledat['spec']>0) )

        ax1.errorbar( hubbledat['z'][isurvey], hubbledat['mu'][isurvey],
                      hubbledat['muerr'][isurvey], hubbledat['zerr'][isurvey],
                      marker='o', ls=' ',mew=0, markersize=5,
                      color=colordict[survey], mfc=colordict[survey],
                      ecolor=colordict[survey],
                      elinewidth=elw, capsize=cs, alpha=0.7 )

    if highz :
        # get the high-z SN fits from Dan
        labeldict = {#'primo':'Rodney+ 2012',
                     #'wilson':'Jones, Rodney+ 2013',
                     'stone':'GND13Sto',
                     'colfax':'GND12Col',
                     }
        if highz=='rodney':
            labeldict = {'primo':'Rodney+ 2012',
                         'wilson':'Jones, Rodney+ 2013',
                         'stone':'Rodney+ in prep',
                         'colfax':'Rodney+ in prep',
                         }

        for survey, namelist in zip( ['hstwfc3a','hstwfc3b'],
                                     [['primo','wilson'],['stone','colfax']] ):
            color=colordict[survey]
            for name in namelist :
                if name not in hubbledat['name'] : continue
                iname = np.where( hubbledat['name'] == name )[0][0]
                z,zerr = hubbledat['z'][iname], hubbledat['zerr'][iname]
                mu,muerr = hubbledat['mu'][iname],hubbledat['muerr'][iname]
                if name=='colfax' :
                    icol = np.where( [n.startswith('col') for n in hubbledat['name']] )[0]
                    # lensing correction for mu=1.04 +-0.03  :  0.0426 +- 0.031 mag
                    mucol = interp1d( hubbledat['z'][icol], hubbledat['mu'][icol]+0.0426)
                    muerr = np.sqrt( muerr**2 + 0.031**2 )
                    ax1.plot( [2.16,2.19,2.25,2.28], mucol([2.16,2.19,2.25,2.28]),
                              color=color, lw=0.75, ls='-' )
                    ax1.errorbar( 2.22, mucol(2.22), muerr,
                                  marker='D', ls=' ',mew=0, markersize=8,
                                  color=color, mfc=color, ecolor=color,
                                  elinewidth=elw, capsize=cs, alpha=1 )
                    z = 2.1
                    mu = mucol(2.1)
                else :
                    if name=='stone':
                        mu = mu + 0.0215
                        muerr = np.sqrt( muerr**2 + 0.02**2 )
                    ax1.errorbar( z, mu, muerr, zerr,
                                  marker='D', ls=' ',mew=0, markersize=8,
                                  color=color, mfc=color, ecolor=color,
                                  elinewidth=elw, capsize=cs, alpha=1 )

                if highzlabels and name in labeldict :
                    ax1.text( z+0.02, mu-0.6, labeldict[name],
                              color=color, # backgroundcolor='w',
                              fontsize='small', va='top', ha='center', rotation=90 )



    # plot the lcdm cosmology
    ax1.axhline( 0, color='0.4', ls='-' )

    # plot extreme w'
    zLCDM = np.arange( 0.01, 5, 0.05)
    muLCDM = cosmo.mu( zLCDM, H0=H0, Om=Om, Ode=(1-Om), w0=-1., wa=0)
    #muextremeDE = cosmo.mu( z, H0=_H0, Om=_Om, Ode=(1-_Om), w0=-0.7, wa=-2)
    mu1 = cosmo.mu( zLCDM, H0=H0, Om=Om-0.015, Ode=(1-Om+0.015), w0=-1.2, wa=0.7)
    mu2 = cosmo.mu( zLCDM, H0=H0, Om=Om+0.055, Ode=(1-Om-0.055), w0=-0.65, wa=-2.2)
    ax1.fill_between( zLCDM, mu1, mu2, color='k', alpha=0.2, lw=0, zorder=10 )


    ax1.text( 0.3, 35,  'low-z',   fontsize='large', ha='center', color=colordict['low-z']     )
    ax1.text( 0.4, 37.5, 'SDSS',     fontsize='large', ha='center', color=colordict['sdss']      )
    ax1.text( 0.7, 40,  'PS1',     fontsize='large', ha='center', color=colordict['ps1']       )
    ax1.text( 1.0, 42, 'SNLS',     fontsize='large', ha='center', color=colordict['snls']      )
    ax1.text( 1.4, 46, 'HST+ACS',  fontsize='large', ha='right', color=colordict['hstacs']    )
    ax1.text( 1.83, 47, 'HST+WFC3', fontsize='large', ha='center', color=colordict['hstwfc3b']  )
    ax1.text( 3.0, 46, '( JWST )',  fontsize='large', ha='center', color=colordict['jwst']  )

    ax1.set_xlabel('Redshift',fontsize='large')
    ax1.set_ylabel(r'Distance Modulus', fontsize='large')

    mulim = array( [33, 49.5] )
    muresidlim = array( [-1.05,1.05] )
    ax1.xaxis.set_ticks( [0.1,] + np.arange(0.5,zlim[-1],0.5).tolist() )
    ax1.set_xlim( zlim )
    ax1.set_ylim( mulim )

    ax1.yaxis.set_major_locator( ticker.MultipleLocator( 4 ) )
    ax1.yaxis.set_minor_locator( ticker.MultipleLocator( 1 ) )
    if ageaxis :
        ax2 = ax1.twin()
        if highz=='rodney':
            # ax1top.set_xlabel( 'age of Universe (Gyr)',fontsize='large',x=0.15,ha='left')
            ageticks = np.array( [13,  9,  6, 4, 3, 2 ] )
        else :
            ageticks = np.array( [13,  9,  6, 4, 3 ] )
        ax2.set_xlim( zlim )
        ztickloctop = cosmo.zfromt( ageticks )
        ax2.xaxis.set_ticks( ztickloctop )
        ax2.xaxis.set_ticklabels( ageticks )
        ax2.yaxis.set_major_locator( ticker.MultipleLocator( 4 ) )
        ax2.yaxis.set_minor_locator( ticker.MultipleLocator( 1 ) )
        ax2.yaxis.set_ticklabels( [ '%i'%round(float(tick)+MB) for tick in ax1.yaxis.get_majorticklocs()] )
        ax2.set_xlabel( 'Age of Universe (Gyr)',fontsize='large')
        ax2.set_ylabel( 'Peak Magnitude',fontsize='large',rotation=-90, labelpad=27)

    pl.draw()





_Nzsteps = 100 # number of redshift steps for progenitor mass effect curves
_MINMASS = 0.6 # absolute minimum MS mass of a SNIa progenitor
_MAXMASS = 8  # maximum MS mass of a SNIa progenitor (above this it becomes a CC SN)
_DMUDMASS = 0.03 # change in apparent distance modulus as MS mass increases

# Average SNIa progenitor mass as a function of redshift
# assumes a Kroupa IMF :  N(M)dM ~ M^{-alpha}   with alpha = 2.7  for M > 1 Msun
# progenitors are drawn from stars with a MS mass between Mmin and Mmax [Msun]
Mave = lambda Mmax,Mmin : 2.43 * (Mmax**(-0.7) - Mmin**(-0.7)) / (Mmax**(-1.7) - Mmin**(-1.7))


def mkEvolPlot( minmassdict={} ) :
    from matplotlib import pyplot as pl

    if not minmassdict :
        minmassdict = progenitorMinMassLines( zform=11 )

    ax1 = pl.subplot( 311 )
    ax2 = pl.subplot( 312, sharex=ax1 )
    ax3 = pl.subplot( 313, sharex=ax1 )

    z = minmassdict['z']
    for key,color in zip( ['mm3','mm2','mm1','mm05'], ['r','g','b','k'] ):
        minmass = minmassdict[key]
        mavez = Mave( _MAXMASS, minmass )
        ax1.plot( z, minmass, ls='-', color=color )
        ax2.plot( z, mavez, ls='-', color=color )
        ax3.plot( z, (mavez-mavez[0]) * _DMUDMASS, ls='-', color=color )
    pl.setp( ax1.get_xticklabels(), visible=False )
    pl.setp( ax2.get_xticklabels(), visible=False )
    ax1.set_ylabel('Min Progenitor Mass')
    ax2.set_ylabel('Mean Progenitor Mass')
    ax3.set_xlabel('Redshift')
    ax3.set_ylabel('Systematic bias in $\mu$')
    #ax1.grid()
    #ax2.grid()
    #ax3.grid()
    ax1.text( 2.0, 5.5, r'$\tau$=2.5 Gyr', color='r', ha='right',va='bottom')
    ax1.text( 2.8, 4, r'$\tau$=1.5 Gyr', color='g', ha='right',va='bottom')
    ax1.text( 3.3, 2.5, r'$\tau$=0.4 Gyr', color='b', ha='right',va='bottom')
    ax2.text( 0.5, 5.0, 'M$_{max}$ = %.1f M$_{\odot}$\n M$_{min}$ = %.1f M$_{\odot}$'%(_MAXMASS,_MINMASS), color='k', ha='left',va='bottom', fontsize='large')
    ax3.text( 0.4, 0.06, r'$d\mu/dM_{MS}$=0.03 mag M$_{\odot}^{-1}$', color='k', ha='left',va='bottom', fontsize='large')
    fig = pl.gcf()
    fig.subplots_adjust( hspace=0, right=0.98, top=0.98)

    ax3.plot( [ 2.24,3.33], [ 0.14,0.176], 'k--' )
    ax1.set_xlim( 0, 3.6 )
    ax3.set_ylim( -0.01, 0.20 )
    ax3.text( 2.4,0.16,'no SNIa',fontsize='small',ha='left',va='bottom',rotation=5)
    ax2.set_ylim( 1.01, 7.9 )
    ax1.set_ylim( 1.01, 7.9 )

    pl.draw()

    return( minmassdict )


def progenitorMinMassLines( zform=11, nz=_Nzsteps ):
    """ returns a dictionary object
    containing a redshift vector and
    a set of vectors giving the mass of a
    SNIa progenitor for each of 4 total
    delay time assumptions :
      z  = redshift vector
      mm3 = min prog. mass for 3.0 Gyr delay from star formation to explosion
      mm2.5 = 2.5 Gyr
      mm2 = 2.0 Gyr
      mm1 = 1.0 Gyr
      mm04 = 0.4 Gyr
    zform gives the redshift at which the SN progenitor star formed.
    Use zform=11 (the redshift of reionization) to get the absolute
    minimum possible progenitor mass allowed by the Universe at
    each redshift.
    """
    from numpy import linspace,where
    from pytools import cosmo

    # our redshift range of interest, with tighter sampling at the top
    z = linspace( 0.2, 4.5, nz )

    # initial the dictionary of minimum mass curves that will be returned
    mmdict = {'z':z }
    _H0 = 72.0
    _Om = 0.277
    _Ode = 1 - _Om
    _MB = -19.175  # absolute SNIa mag assuming H0=70

    # age of Universe at time of progenitor star formation
    ageAtForm = cosmo.agez( zform, H0=_H0, Om=_Om, Ode=_Ode )

    # age of Univ. at this z :
    ageUniverse = cosmo.agez( z, H0=_H0, Om=_Om, Ode=_Ode )

    # time since formation (i.e. time since first epoch of star formation
    tSinceForm = ageUniverse - ageAtForm

    #for tdelay, key in zip( [ 2.5, 1.5, 1.0, 0.4], ['mm25','mm20','mm15','mm04'] ):
    for tdelay, key in zip( [3, 2, 1, 0.5], ['mm3','mm2','mm1','mm05'] ):

        # The min prog. mass is the mass of a star that is born at zform,
        # evolves off the main sequence and then instantly reaches the point of
        # explosion at this redshift.
        # So this star has the maximum possible main sequence lifetime
        # for a Ia progenitor at this redshift :
        tMSmax = tSinceForm - tdelay

        itMSpos = where( tMSmax>0 )

        # main sequ. lifetime as a function of Mass, in Gyrs:
        # tMS = 8 * Mass**(-2.5)
        # invert this to get the mass as a function of lifetime,
        # restricting to positive lifetimes
        # and restricting to plausible SNIa progenitor masses (3-8 Msun)
        minmass = where(  tMSmax>0,  (tMSmax/8)**(-0.4), _MAXMASS )
        minmass = where( minmass<_MAXMASS, minmass, _MAXMASS )
        minmass = where( minmass>_MINMASS, minmass, _MINMASS )
        mmdict[ key ] = minmass

    return( mmdict )
