__author__ = 'rodney'

from astropy import io
import numpy as np

def readascii( filename ):
    from astropy.io import ascii
    return( ascii.read( filename, format='commented_header', header_start=-1, data_start=0 ))

def do_colfax_fit( nbins=51 ):
    """  run three rounds of the nested sampling to get redshift posteriors with
      * just the wide band data, no host prior
      * just the med band epoch, no host prior
      * all data, including host prior
    :return:
    """
    from medband import classify

    #TODO chdir to .dat dir
    colfaxMB = readascii('HST_CANDELS4_colfaxMB.sncosmo.dat')
    colfaxMBout = classify.get_evidence( colfaxMB, zhost=2.26, zhosterr=5, zminmax=[1.3, 2.8], nobj=100, maxiter=10000 )
    colfaxMBpdfs = classify.get_marginal_pdfs( colfaxMBout[1], nbins=nbins, verbose=True )
    io.misc.pickle_helpers.fnpickle( colfaxMBpdfs, '../colfax_posterior_pdfs_MB.pkl')

    colfaxWB = readascii('HST_CANDELS4_colfaxWB.sncosmo.dat')
    colfaxWBout = classify.get_evidence( colfaxWB, zhost=2.26, zhosterr=5, zminmax=[1.3, 2.8], nobj=100, maxiter=10000 )
    colfaxWBpdfs = classify.get_marginal_pdfs( colfaxWBout[1], nbins=nbins, verbose=True )
    io.misc.pickle_helpers.fnpickle( colfaxWBpdfs, '../colfax_posterior_pdfs_WB.pkl')

    colfaxALL = readascii('HST_CANDELS4_colfaxALL.sncosmo.dat')
    colfaxALLout = classify.get_evidence( colfaxALL, zhost=2.26, zhosterr=0.09, zminmax=[1.3, 2.8], nobj=100, maxiter=10000 )
    colfaxALLpdfs = classify.get_marginal_pdfs( colfaxALLout[1], nbins=101, verbose=True )
    io.misc.pickle_helpers.fnpickle( colfaxALLpdfs, '../colfax_posterior_pdfs_ALL.pkl')
    return( colfaxMBout, colfaxWBout, colfaxALLout )


def do_stone_fit():
    """  run three rounds of the nested sampling to get redshift posteriors with
      * just the wide band data, no host prior
      * just the med band epoch, no host prior
      * all data, including host prior
    :return:
    """
    from medband import classify

    #TODO chdir to .dat dir
    stoneMB = readascii('HST_CANDELS4_stoneMB.sncosmo.dat')
    stoneMBout = classify.get_evidence( stoneMB, zhost=1.8, zhosterr=5.0, zminmax=[1.3, 2.8], nobj=100, maxiter=10000 )
    stoneMBpdfs = classify.get_marginal_pdfs( stoneMBout[1], nbins=35, verbose=True )
    io.misc.pickle_helpers.fnpickle( stoneMBpdfs, '../stone_posterior_pdfs_MB.pkl')

    stoneWB = readascii('HST_CANDELS4_stoneWB.sncosmo.dat')
    stoneWBout = classify.get_evidence( stoneWB, zhost=1.8, zhosterr=5.0, zminmax=[1.3, 2.8], nobj=100, maxiter=10000 )
    stoneWBpdfs = classify.get_marginal_pdfs( stoneWBout[1], nbins=35, verbose=True )
    io.misc.pickle_helpers.fnpickle( stoneWBpdfs, '../stone_posterior_pdfs_WB.pkl')

    stoneALL = readascii('HST_CANDELS4_stoneALL.sncosmo.dat')
    stoneALLout = classify.get_evidence( stoneALL, zhost=1.80, zhosterr=0.02, zminmax=[1.3, 2.8], nobj=100, maxiter=10000 )
    stoneALLpdfs = classify.get_marginal_pdfs( stoneALLout[1], nbins=35, verbose=True )
    io.misc.pickle_helpers.fnpickle( stoneALLpdfs, '../stone_posterior_pdfs.pkl')
    return( stoneMBout, stoneWBout, stoneALLout )


def do_bush_fit( verbose=True, nobj=100, maxiter=10000 ):
    """  run three rounds of the nested sampling to get redshift posteriors with
      * just the wide band data, no host prior
      * just the med band epoch, no host prior
      * all data, including host prior
    :return:
    """
    from medband import classify

    #TODO chdir to .dat dir
    bushClassOutDict = {}
    for datversion in ['MB','WB','ALL'] :
        if datversion == 'ALL' :
            zhost=1.15
            zhosterr=[-0.38,+0.19]
        else :
            zhost=None
            zhosterr=None
        bush = readascii('HST_CANDELS4_bush%s.sncosmo.dat'%datversion )
        bushClassOut = classify.classify( bush, zhost=1.15, zhosterr=[-0.06,0.13], t0_range=[55803.1-20,55803.1+20], zminmax=[0.6,1.6], nobj=100, maxiter=10000, nsteps_pdf=101, verbose=3 )


        for source in ['s11-2004hx','s11-2005lc','s11-2005gi','s11-2006jl','s11-2005hl','s11-2005hm','s11-2006fo', 's11-2006jo','salt2','salt2-extended']:


            classify.get_evidence( bush, modelsource='salt2',zhost=zhost, zhosterr=zhosterr, zminmax=[0.26,2.21], nobj=nobj, maxiter=maxiter, verbose=verbose )
            io.misc.pickle_helpers.fnpickle( bushClassOut['pdf'], '../bush_posterior_pdfs_%s_MB.pkl'%source)

            bushClassOutDict[ 'bush%s+%s'%(datversion,source) ] = bushClassOut
    return( bushClassOutDict )

def mk_composite_pdf( classifyOutput, snid='bush', suffix='MB',
                      templateset='PSNID', nsteps_pdf=101, showplots=False ):
    from scipy.interpolate import interp1d
    from .classify import SubClassDict_NOTPSNID, SubClassDict_PSNID, SubClassDict_SNANA
    from . import classify
    logztot = classifyOutput['logzall']

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

    z = np.arange( 0, 3, 0.01 )
    pdfz = np.zeros( len(z) )
    for modelsource in ['s11-2004hx','s11-2005lc','s11-2005gi','s11-2006jl','s11-2005hl','s11-2005hm','s11-2006fo', 's11-2006jo','salt2','salt2-extended']:
        if modelsource not in classifyOutput.keys() : continue
        logz = classifyOutput[modelsource]['res']['logz']

        # multiply the model evidence by the sub-type prior
        if modelsource in iimodelnames :
            logprior = logpriordict['ii']
        elif modelsource in ibcmodelnames :
            logprior = logpriordict['ibc']
        elif modelsource in iamodelnames :
            logprior = logpriordict['ia']
        logzpost = logprior + logz

        scale = np.exp( logzpost - logztot )

        if classifyOutput[modelsource]['pdf'] is None :
            res = classifyOutput[modelsource]['res']
            pdf =  classify.get_marginal_pdfs( res, nbins=nsteps_pdf )
            classifyOutput[modelsource]['pdf'] = pdf

        pdfzsamp = classifyOutput[modelsource]['pdf']['z']
        pdfzinterp = interp1d( pdfzsamp[0], pdfzsamp[1], bounds_error=False, fill_value=0)
        pdfz += pdfzinterp( z ) * scale
        if showplots :
            from matplotlib import pyplot as pl
            print( "%s : %.3e"%(modelsource, scale))
            pl.plot( z, pdfzinterp(z)*scale, label=modelsource )
    if showplots :
        pl.legend()

    zmax = z[ np.argmax( pdfz ) ]
    ihalfmax = np.where( np.abs( pdfz - pdfz.max()/2.) < 0.1*pdfz.max() )[0]
    zerr = np.min( zmax - z[ihalfmax[0]], z[ihalfmax[-1]] - zmax )

    pdf = {'z':[ z, pdfz, zmax, zerr ]}
    io.misc.pickle_helpers.fnpickle( pdf, '../%s_posterior_pdfs_%s.pkl'%(snid,suffix) )


def mk_fig( sn, legend=False, nsmooth=3 ):
    """ construct the redshift pdf(z) figure(s) for the paper
    :param sn:
    :return:
    """
    from matplotlib import pyplot as pl
    from medband import classify
    from pytools import plotsetup, colorpalette as cp
    from scipy.interpolate import interp1d

    fillcolors = [cp.bluegrey, cp.red, cp.brown, cp.black]
    textcolors = [cp.darkbluegrey, cp.red, cp.brown, cp.grey]

    plotsetup.halfpaperfig( )
    pl.clf()
    fig = pl.gcf()
    ax = pl.gca()

    pdfWB = io.misc.pickle_helpers.fnunpickle('../%s_posterior_pdfs_WB.pkl'%sn)
    pdfMB = io.misc.pickle_helpers.fnunpickle('../%s_posterior_pdfs_MB.pkl'%sn)
    pdfALL = io.misc.pickle_helpers.fnunpickle('../%s_posterior_pdfs_ALL.pkl'%sn)

    if sn == 'colfax' :
        zhost = 2.26
        zhosterr = [-0.09,+0.09]
        zsamp, pdfsamp = np.loadtxt("../HOSTS/SEDFITS/colfaxA_photoz_pdf.dat", unpack=True )
        zprior = interp1d( zsamp, pdfsamp, bounds_error=False, fill_value=0 )
        #z = np.linspace( 0, 3., 301 )
        z = np.linspace( 1.3,2.8,101 )
        def zprior( z ) :
            return( classify.gauss( z, zhost, zhosterr ) )
    elif sn == 'stone' :
        zhost = 1.80
        zhosterr = [-0.02,+0.02]
        z = np.linspace( 1.3,2.8,101 )
        def zprior( z ) :
            return( classify.gauss( z, zhost, zhosterr ) )
    elif sn == 'bush' :
        zhost = 1.15
        zhosterr = [-0.38,+0.19]
        zsamp, pdfsamp = np.loadtxt("../HOSTS/SEDFITS/bushA_photoz_pdf.dat", unpack=True )
        z = np.linspace( 0, 3., 301 )
        zprior = interp1d( zsamp, pdfsamp, bounds_error=False, fill_value=0 )


    zstrlist = []
    for pdf in [ pdfWB, pdfMB, pdfALL ] :
        zpk = round( pdf['z'][0][ pdf['z'][1].argmax() ], 2 )
        zerrplus = round( pdf['z'][2] + pdf['z'][3] - zpk, 2 )
        zerrminus = round( pdf['z'][2] - pdf['z'][3] - zpk, 2 )
        zstrlist.append( '$z = %.2f ^{%+.2f}_{%+.2f}$'%( zpk, zerrplus, zerrminus ) )

    #zWB_str = '$z = %.2f \pm %.2f$'%(pdfWB['z'][2],pdfWB['z'][3])
    #zMB_str = '$z = %.2f \pm %.2f$'%(pdfMB['z'][2],pdfMB['z'][3])
    #zAll_str= '$z = %.2f \pm %.2f$'%(pdfALL['z'][2],pdfALL['z'][3])
    #zhost_str='$z = %.2f ^{%+.2f}_{%+.2f}$'%( zhost, zhosterr[1], zhosterr[0] )
    if sn in [ 'stone', 'colfax' ]:
        zhost_str='$z = %.2f \pm %.2f$'%( zhost, zhosterr[1] )

    medsmooth = lambda f,N : np.array( [ np.mean( f[max(0,i-N):min(len(f),max(0,i-N)+2*N)]) for i in range(len(f)) ] )
    for x,y,c,fill in zip( [pdfWB['z'][0],pdfMB['z'][0],z,pdfALL['z'][0]],
                         [pdfWB['z'][1],pdfMB['z'][1],zprior(z),pdfALL['z'][1]],
                         fillcolors,
                         [1,1,1,0]) :
        y = y /y.max()
        if fill :
            pl.fill_between( x, np.zeros(len(y)), y, color=c, alpha=0.3 )
        else :
            if nsmooth>1 :
                y = medsmooth(y,nsmooth) / medsmooth(y,nsmooth).max()
            pl.plot( x, y, color=c, lw=1.7,  ls='-' )

    pl.xlabel('Redshift')
    pl.ylabel('p($z|$D) (normalized)')

    arrowlw = 0.5
    if legend :
        pl.legend(loc='upper left', ncol=2, frameon=False, fontsize=10, handlelength=0.4, numpoints=3 )
    elif sn=='stone' :
        ax.text( 2.2, 0.74, 'wide band\n light curve\n '+zstrlist[0], ha='center', va='bottom', color=textcolors[0], fontsize=10 )
        ax.plot( [2.18,2.0], [0.72,0.57], ls='-', color=textcolors[0], lw=arrowlw )

        ax.text( 2.5, 0.56, 'med. band colors\n '+zstrlist[1], ha='center', va='bottom', color=textcolors[1], fontsize=10 )
        ax.plot( [2.48,2.28], [0.54,0.41], ls='-', color=textcolors[1], lw=arrowlw )

        ax.text( 1.50, 0.82, 'host prior\n '+zhost_str, ha='center', va='bottom', color=textcolors[2], fontsize=10 )
        ax.plot( [1.52,1.8], [0.8,0.33], ls='-', color=textcolors[2],lw=arrowlw )

        ax.text( 1.84, 1.05, 'combined constraint\n'+zstrlist[2], ha='center', va='bottom', color='k', fontsize=10 )
        ax.plot( [1.82,1.8], [1.04,0.98], ls='-', color='0.2',lw=arrowlw )

        ax.text( 0.95, 0.95,  'GND13Sto', fontsize='large', ha='right', va='top', transform=ax.transAxes )
    elif sn=='colfax' :
        ax.text( 2.60, 0.84, 'wide band\n light curve\n '+zstrlist[0], ha='center', va='bottom', color=textcolors[0], fontsize=10 )
        ax.plot( [2.56,2.4], [0.82,0.42], ls='-', color=textcolors[0], lw=arrowlw )

        ax.text( 1.57, 0.61, 'med. band colors\n '+zstrlist[1], ha='center', va='bottom', color=textcolors[1], fontsize=10 )
        ax.plot( [1.62,2.02], [0.59,0.36], ls='-', color=textcolors[1], lw=arrowlw )

        ax.text( 1.82, 0.86, 'host prior\n '+zhost_str, ha='center', va='bottom', color=textcolors[2], fontsize=10 )
        ax.plot( [1.86,2.18], [0.84,0.61], ls='-', color=textcolors[2],lw=arrowlw )

        ax.text( 1.96, 1.02, 'combined constraint\n'+zstrlist[2], ha='center', va='bottom', color='k', fontsize=10 )
        ax.plot( [2.0,2.25], [1.00,0.81], ls='-', color='0.2',lw=arrowlw )

        #ax.text( 1.84, 1.05, 'host prior\n '+zhost_str, ha='center', va='bottom', color=cp.darkblue, fontsize=10 )
        #ax.plot( [1.82,2.02], [1.04,0.73], ls='-', color=cp.darkblue,lw=arrowlw )

        ax.text( 0.95, 0.95,  'GND12Col', fontsize='large', ha='right', va='top', transform=ax.transAxes )
        ax.set_ylim( -0.003, 1.19 )
        ax.set_xlim( 1.31,2.79 )


    elif sn=='bush' :
        ax.text( 2.05, 0.84, 'wide band\nlight curve\n '+zWB_str, ha='center', va='bottom', color=cp.green, fontsize=10 )
        # ax.plot( [2.56,2.45], [0.82,0.48], ls='-', color=cp.green, lw=arrowlw )

        ax.text( 2.02, 0.45, 'med. band\ncolors\n '+zMB_str, ha='center', va='bottom', color=cp.red, fontsize=10 )
        # ax.plot( [1.62,2.0], [0.3,0.12], ls='-', color=cp.red, lw=arrowlw )

        ax.text( 1.35, 1.05, 'host prior\n '+zhost_str, ha='center', va='bottom', color=cp.darkblue, fontsize=10 )
        # ax.plot( [1.82,2.02], [1.04,0.73], ls='-', color=cp.darkblue,lw=arrowlw )

        ax.text( 0.83, 1.01, 'combined\nconstraint\n'+zAll_str, ha='center', va='bottom', color='k', fontsize=10 )
        # ax.plot( [1.72,2.1], [0.68,0.33], ls='-', color='0.2',lw=arrowlw )

        ax.text( 0.95, 0.95,  'GSD11Bus', fontsize='large', ha='right', va='top', transform=ax.transAxes )
        ax.set_ylim( -0.003, 1.19 )
        ax.set_xlim( 0.2,2.4 )

    fig.subplots_adjust( left=0.12,bottom=0.12,right=0.97, top=0.97 )
    pl.draw()


def get_pdf_err( x, y ):
    """ measure the 68% confidence region of a probability dist'n function y(x)
    :param pdf:
    :return:
    """
    from scipy.integrate import simps
    import numpy as np

    ytot = simps( y, x )
    ipeak = [ np.argmax(y) ]

    nsteps = 0
    #import pdb; pdb.set_trace()
    while nsteps<len(y):
        if simps( y[ipeak], x[ipeak] )/ytot >= 0.68 : break
        iright = ipeak[-1]
        ileft = ipeak[0]
        if iright+1>=len(y):
            ipeak = [ileft-1] + ipeak
        elif ileft-1<0:
            ipeak = ipeak + [iright+1]
        elif y[iright+1]<y[ileft-1]:
            ipeak = [ileft-1] + ipeak
        elif y[iright+1]>y[ileft-1] :
            ipeak = ipeak + [iright+1]
        else :
            ipeak = [ileft-1] + ipeak + [iright+1]
        print( "%.2f : %.3e %.3e: %i  %i" %( simps( y[ipeak], x[ipeak] )/ytot, y[ileft], y[iright], ileft, iright ) )
        print('continue')
    ipk = np.argmax(y)
    ipkL = ipeak[0]
    ipkR = ipeak[-1]
    print( "%.2f  %+.2f  %+.2f"%( x[ipk], x[ipk]-x[ipkL], x[ipk]-x[ipkR]) )



