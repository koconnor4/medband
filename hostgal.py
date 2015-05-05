__author__ = 'rodney'


# SN positions, measured on the goods-n/s 60mas mosaic images by hand in ds9
xColfax, yColfax = 12349.177, 14478.117
xStone, yStone = 7793.1145, 12619.196
xBush, yBush = 6255.20,   11186.55

def host_separation( ):
    """  Colfax has two possible host galaxies.  This function
    is just a record of an ipython session, in which I extracted the
    parameters of sextractor Kron ellipses from the CANDELS hsex catalog
    and converted them into the R parameter, following Sullivan et al 2006.
    I get R_A = 1.6  and R_B = 4.4, indicating that the SN is most likely
    associated with host candidate A.
    :return:
    """
    import numpy as np
    from astropy.io import ascii
    cat1 = ascii.read('/store/goods-n/catalogs/gn_all_candels_sx_h_130829_hphotom_comb_merge_psfmatch2h.cat')

    # identify the catalog indices for hosts A and B
    ra = cat1['ALPHA_J2000']
    dec = cat1['DELTA_J2000']
    raA,decA = 189.15635, 62.30910  # measured from the images in ds9
    raB,decB = 189.15681, 62.30839
    distA = np.sqrt( (ra-raA)**2 + (dec-decA)**2 )
    distB = np.sqrt( (ra-raB)**2 + (dec-decB)**2 )

    # check that the distances are small
    assert distA.min() < 1e-4
    assert distB.min() < 1e-4
    iA = distA.argmin()
    iB = distB.argmin()

    # extract the ellipse parameters :
    CxxA = cat1['CXX_IMAGE'][iA]
    CyyA = cat1['CYY_IMAGE'][iA]
    CxyA = cat1['CXY_IMAGE'][iA]
    CxxB = cat1['CXX_IMAGE'][iB]
    CyyB = cat1['CYY_IMAGE'][iB]
    CxyB = cat1['CXY_IMAGE'][iB]
    xA = cat1['X_IMAGE'][iA]
    xB = cat1['X_IMAGE'][iB]
    yA = cat1['Y_IMAGE'][iA]
    yB = cat1['Y_IMAGE'][iB]

    # SN position, measured on the goods-n mosaic image by hand in ds9
    xSN, ySN = 12349.177, 14478.117

    xrA = xSN - xA
    yrA = ySN - yA
    xrB = xSN - xB
    yrB = ySN - yB
    RA = np.sqrt(  CxxA*xrA**2 + CyyA*yrA**2 + CxyA*xrA*yrA )
    RB = np.sqrt(  CxxB*xrB**2 + CyyB*yrB**2 + CxyB*xrB*yrB )
    thetaB = cat1['THETA_IMAGE'][iB]
    thetaA = cat1['THETA_IMAGE'][iA]
    rA_A = cat1['A_IMAGE'][iA]
    rB_A = cat1['B_IMAGE'][iA]
    rA_B = cat1['A_IMAGE'][iB]
    rB_B = cat1['B_IMAGE'][iB]
    print( "host A : ", RA, rA_A, rB_A, thetaA, xA, yA )
    print( "host B : ", RB, rA_B, rB_B, thetaB, xB, yB )


def get_host_ellipses( sn='colfax', cat1=None, rclose=10 ):
    """
    get ellipse parameters for all sources within <rclose> arcsec
    :return:
    """
    import numpy as np
    from astropy.io import ascii

    if sn == 'colfax' : 
        xSN, ySN = xColfax, yColfax
    elif sn == 'stone' : 
        xSN, ySN = xStone, yStone
    elif sn == 'bush' : 
        xSN, ySN = xBush, yBush
        raSN,decSN = 53.178232, -27.801973 

    if cat1 is None and sn in ['colfax','stone']:
        cat1 = ascii.read('/store/goods-n/catalogs/gn_all_candels_sx_h_130829_hphotom_comb_merge_psfmatch2h.cat')

    elif cat1 is None and sn in ['bush']:
        cat1 = ascii.read('/store/goods-s/catalogs/gs_all_sx_jh_120605_hphotom_comb_merge_psfmatch2h.cat')

    ra = cat1['ALPHA_J2000']
    dec = cat1['DELTA_J2000']

    x = cat1['X_IMAGE']
    y = cat1['Y_IMAGE']

    dist = np.sqrt( ((ra-raSN)*np.cos(np.deg2rad(dec))) **2 + (dec-decSN)**2 ) * 3600.
    iclose = np.where( dist< rclose )[0]

    if len(iclose)<1: 
        import pdb; pdb.set_trace()

    print("#idx           RA        Dec          d          x          y    theta        A       B       R")

    for i in iclose :
        # extract the ellipse parameters :

        Cxx = cat1['CXX_IMAGE'][i]
        Cyy = cat1['CYY_IMAGE'][i]
        Cxy = cat1['CXY_IMAGE'][i]
        xgal = x[i]
        ygal = y[i]

        ragal = ra[i]
        decgal = dec[i]

        xr = xSN - xgal
        yr = ySN - ygal
        R = np.sqrt(  Cxx*xr**2 + Cyy*yr**2 + Cxy*xr*yr )
        theta = cat1['THETA_IMAGE'][i]
        A = cat1['A_IMAGE'][i]
        B = cat1['B_IMAGE'][i]
        darcsec = dist[i]

        print( "%5i %11.6f %10.6f %10.3f %10.2f %10.2f %8.3f %8.3f %8.3f %6.2f"%( i,ragal,decgal,darcsec,xgal,ygal,theta,A,B,R) )


band2wave = {'u':3600, 'k':21900, 'ks':21500,
             'ch1':36000, 'ch2':45000, 'ch3':58000,'ch4':80000 }

def mk_host_figure( sn='colfax', showpdf=True ):
    """ figure showing galaxy SED fit and pdf(z) for colfaxA
    :return:
    """
    from astropy.io import ascii
    import os
    import numpy as np
    from matplotlib import pyplot as pl, ticker
    from pytools import plotsetup
    import sncosmo

    fig = plotsetup.halfpaperfig( 1, figsize=[5,3])
    ax = fig.add_subplot(1,1,1)

    datadir = "/Users/rodney/Dropbox/MEDBAND/HOSTS/SEDFITS/"

    # read in the observed  SED from the catalog
    # hostphot = ascii.read( os.path.join(datadir,'bushA_ABmag.dat'), format='fixed_width_two_line')
    hostphot = ascii.read( os.path.join(datadir,sn+'A_ABmag.dat'), format='fixed_width_two_line')
    wavedict = {}
    ibest = np.where( [ hostname.endswith('best') for hostname in hostphot['SNHOST_NAME'] ])[0][0]

    for colname in hostphot.colnames :
        if not colname.endswith( 'MAG' ) : continue
        bandname = colname.split('_')[1].lower()
        if bandname.lower() in band2wave.keys() :
            wave = band2wave[bandname]
        else :
            if not bandname.startswith('f') :
                bandname = 'bessell' + bandname
            try :
                band = sncosmo.get_bandpass( bandname )
                wave = band.wave_eff
            except :
                print("no wave_eff for %s"%bandname)
                continue
        wavedict[ colname ] = wave

        mag = hostphot[colname][ibest]
        magerr = hostphot[colname+'ERR'][ibest]
        ax.errorbar( wave, mag, magerr, marker='o', color='k', ls=' ')

    # read in the best-fit SED
    # sedfit = ascii.read( os.path.join(datadir,sn+'A_bestfit_sed.dat'), format='commented_header')
    # ax.plot( sedfit['wavelength'], sedfit['ABmag'] )

    sedfit = ascii.read( os.path.join(datadir,sn+'A_bestfit_sed.dat'), format='commented_header', header_start=-1, data_start=0)
    ax.plot( sedfit['wavelength'], sedfit['ABmag'], color='0.7', zorder=-100 )
    ax.set_xscale('log')

    xlim = [2500,90000]
    if sn=='colfax':
        ylim = [27.9,21.1]
        ax.text( 0.05, 0.92, 'GND12Col Host SED', ha='left', va='top',transform=ax.transAxes )

    elif sn=='bush':
        ylim = [30.05,26.4]
        ax.text( 0.05, 0.92, 'GSD11Bus Host SED', ha='left', va='top',transform=ax.transAxes )

    ax.set_xlim( xlim )
    ax.set_ylim( ylim )

    if showpdf :
        # read in the pdf(z)
        pdfdat = ascii.read( os.path.join(datadir,sn+'A_photoz_pdf.dat'), format='commented_header', header_start=-1, data_start=0)

        ax2 = pl.axes( [0.55,0.07,0.4,0.28] )
        ax2.plot( pdfdat['z'], pdfdat['pdf']/pdfdat['pdf'].max(), color='0.5' )
        ax2.text( 0.1, 0.85, 'P(z)', ha='left', va='top',transform=ax2.transAxes, color='0.5' )

        ax2.set_xlim( 0.01, 2.9 )
        ax2.set_ylim( 0.01, 1.1 )
        ax2.xaxis.set_tick_params( which='both', pad=0)
        ax2.xaxis.set_ticks_position('top')
        ax2.xaxis.set_ticks_position('both')
        ax2.xaxis.set_label_position('top')
        ax2.yaxis.set_ticklabels( [] )
        ax2.xaxis.set_major_locator( ticker.MultipleLocator(1.0) )
        ax2.xaxis.set_minor_locator( ticker.MultipleLocator(0.2) )
        ax2.yaxis.set_major_locator( ticker.MultipleLocator(0.5) )
        ax2.yaxis.set_minor_locator( ticker.MultipleLocator(0.1) )
        ax2.set_xlabel( 'Redshift')#, labelpad=12)

    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_ticks_position('both')
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticklabels( [] )
    ax.set_xlabel( 'Wavelength', labelpad=7)
    ax.set_ylabel( 'AB mag')
    ax.text( 0.0, 1.02, '300 nm', ha='left', va='bottom', fontsize='small', transform=ax.transAxes)
    ax.text( 0.4, 1.02, '1 $\mu$m', ha='center', va='bottom', fontsize='small', transform=ax.transAxes)
    ax.text( 1.0, 1.02, '8 $\mu$m', ha='right', va='bottom', fontsize='small', transform=ax.transAxes)

    fig.subplots_adjust( left=0.12,bottom=0.07,right=0.95, top=0.85)
    pl.draw()


def colfax_photoz_combine(savedat=False,showplot=False):
    """  combine the photoz pdfs from Bahram/Shooby and Alberto
    :return:
    """
    from astropy.io import ascii
    from scipy.interpolate import interp1d
    from scipy.integrate import simps
    import numpy as np

    datadir = "/Users/rodney/Dropbox/MEDBAND/HOSTS/SEDFITS/"

    ucr = ascii.read(datadir+'colfaxA_photoz_pdf_UCR.dat',format='commented_header',header_start=-1,data_start=0)
    bpz = ascii.read(datadir+'colfaxA_photoz_pdf_BPZ.dat',format='commented_header',header_start=-1,data_start=0)

    ucrint = interp1d( ucr['z'], ucr['p(z)'], bounds_error=False, fill_value=0)
    bpzint = interp1d( bpz['z'], bpz['p(z)'], bounds_error=False, fill_value=0)


    znew = np.arange( 0.01, 3.0, 0.01 )
    ucrnew = ucrint( znew )
    bpznew = bpzint( znew )

    mu,sigma = 2.26,0.09
    medgauss = 1/np.sqrt(2*np.pi*sigma**2 ) * np.exp( -(mu-znew)**2/(2*sigma**2))


    ucrnewsum = simps( ucrnew, znew )
    bpznewsum = simps( bpznew, znew )

    combopdf = (ucrnew/ucrnewsum + bpznew/bpznewsum) / 2.


    fout = open(datadir+"colfaxA_photoz_pdf.dat",'w')
    print >> fout, '#  z   pdf  '
    for zz,pp in zip( znew, combopdf ) :
        print >> fout, '%.3f  %.5e'%(zz,pp)
    fout.close()

    if showplot :
        from matplotlib import pyplot as pl
        pl.plot( znew, ucrnew/ucrnewsum, 'r--', label='UCR')
        pl.plot( znew, bpznew/bpznewsum, 'g--', label='BPZ')
        pl.plot( znew, combopdf, 'b-', label='UCR+BPZ')
        pl.plot( znew, medgauss, 'k-', label='2.26 +- 0.09')
        pl.legend( loc='upper right', numpoints=2 )
    return( znew, combopdf )



