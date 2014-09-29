__author__ = 'rodney'

def synphotSNR( filtname, mag, etime=1200, nexp=4, aperture=0.4, pixscale=0.128):
    import pysynphot as syn
    from numpy import pi, sqrt
    bandpass = syn.ObsBandpass('wfc3,ir,%s,aper#%.1f'%(filtname,aperture))
    sn1a = syn.FileSpectrum('sn1a_flux_max.dat').redshift( 1.5 ).renorm( mag, 'abmag', bandpass )
    sourcecps = syn.Observation( sn1a, bandpass ).countrate()
    srccounts = sourcecps * etime


    skybandpass = syn.ObsBandpass('wfc3,ir,%s'%filtname)
    earthshine = syn.FileSpectrum( 'earthshine_average.dat')
    earthcps_arcsec2 = syn.Observation( earthshine, skybandpass, force='extrap' ).countrate()
    V = syn.ObsBandpass('landolt,V')
    zodiacal = syn.FileSpectrum( 'zodiacal_high.dat').renorm( 22.7, 'vegamag',V)
    zodcps_arcsec2 = syn.Observation( zodiacal, skybandpass, force='extrap').countrate()

    aparea_arcsec2 = pi * aperture**2
    aparea_pix = aparea_arcsec2 / pixscale**2
    npix = (aparea_arcsec2 / pixscale**2)
    sharpness = 3.3 / npix

    skycounts = (earthcps_arcsec2+zodcps_arcsec2) * aparea_arcsec2 * etime
    skycounts_pix = skycounts / npix

    # TODO : get thermal noise, dark current and readnoise from SYNPHOT
    darkcps_arcsec2 = 3.05
    darkcounts = darkcps_arcsec2 * aparea_arcsec2 * etime
    darkcounts_pix = darkcounts / npix

    thermcps_arcsec2 = 8.15
    thermcounts = thermcps_arcsec2 * aparea_arcsec2 * etime
    thermcounts_pix = thermcounts / npix

    rdnoise_pix = 12.5 # read noise per pixel (per exposure)
    rdnoiseTot = rdnoise_pix * aparea_pix  # 161.61

    bgcounts_pix = skycounts_pix + darkcounts_pix + thermcounts_pix


    noise = sqrt( srccounts + npix*bgcounts_pix + npix*nexp*rdnoise_pix**2 )

    SNR = srccounts / noise

    onoise = sqrt(srccounts  +  (bgcounts_pix + rdnoise_pix**2)/sharpness )
    oSNR = srccounts / onoise

    return SNR, oSNR


def update_psfdat():
    import os
    from astropy.io import ascii
    from astropy.table import Column
    import pyfits
    # read in the psf mags computed by Dan
    psfdat = ascii.read( 'psfphot.dat')

    # get exposure time and number of exposures from the header
    etimelist, nexplist, etcSNRlist, etcOSNRlist, psfSNRlist = [], [], [], [], []
    for row in psfdat :
        imdir = '/store/snreproc/%s.090mas/'%(row['sn'])
        suffix = '-e00_sub_masked' if row['image']=='sub' else '_reg_drz_sci'
        imfile = os.path.join( imdir,
            '%s_%s_%s%s.fits'%(
                row['sn'], row['band'], row['epoch'], suffix ) )
        etime = pyfits.getval( imfile, 'EXPTIME' )
        nexp = pyfits.getval( imfile, 'NCOMBINE' )
        etimelist.append( etime )
        nexplist.append( nexp )
        etcSNR, etcSNRopt = synphotSNR( row['band'], row['stack'], etime=etime, nexp=nexp)
        etcSNRlist.append( etcSNR )
        etcOSNRlist.append( etcSNRopt )
        psfSNR = 1.08574 / row['err']
        psfSNRlist.append( psfSNR )

    nexpCol = Column( nexplist, 'nexp', dtype=int, format='%i')
    etimeCol = Column( etimelist, 'etime', dtype=float, format='%8.2f' )
    etcSNRCol = Column( etcSNRlist, 'SNRetc', dtype=float, format='%6.2f' )
    etcOSNRCol = Column( etcOSNRlist, 'optSNRetc', dtype=float, format='%6.2f' )
    psfSNRCol = Column( psfSNRlist, 'SNRpsf', dtype=float, format='%6.2f' )

    psfdat.add_columns( [ nexpCol, etimeCol, etcSNRCol, etcOSNRCol, psfSNRCol], indexes=[4,4,9,9,9])
    psfdat.write( 'psfphot2.dat', format='ascii.commented_header' )
    psfdat.write( 'psfphot3.dat', format='ascii.fixed_width' )

def plotSNRcompare() :
    from astropy.io import ascii
    import numpy as np
    from matplotlib import pyplot as pl
    from pytools import plotsetup

    fig = plotsetup.presfig()
    psfdat = ascii.read( 'psfphot3.dat', format='fixed_width')

    isub = np.where(psfdat['image']=='sub')[0]

    oSNRetc = psfdat['optSNRetc'][isub]
    SNRetc = psfdat['SNRetc'][isub]
    SNRpsf = psfdat['SNRpsf'][isub]
    psfmag = psfdat['stack'][isub]
    band = psfdat['band'][isub]

    iwide = np.array([ i for i in range(len(band)) if band[i].endswith('w')] )
    imed = np.array([ i for i in range(len(band)) if band[i].endswith('m')] )

    pl.clf()
    ax = pl.subplot( 2, 2,1 )
    ax.plot( oSNRetc[iwide], SNRpsf[iwide], color='darkorange', marker='d', ls=' ', alpha=0.5, mew=0.5, ms=8, label='wide band subs' )
    ax.plot( oSNRetc[imed], SNRpsf[imed], 'mo', ls=' ', alpha=0.7, mew=0.5, ms=8, label='med band subs', )
    ax.plot( [0,200], [0,200], 'k-' )
    ax.set_xlim( 0, 45 )
    ax.set_ylim( 0, 45 )
    ax.set_xlabel( 'Optimal S/N from ETC' )
    ax.set_ylabel( 'Observed S/N from PSF drops' )
    ax.legend( numpoints=1, loc='lower right' )

    ax = pl.subplot( 2,2,2 )
    ax.plot( oSNRetc[iwide], SNRpsf[iwide], color='darkorange', marker='d', ls=' ', alpha=0.5, mew=0.5, ms=8 )
    ax.plot( oSNRetc[imed], SNRpsf[imed], 'mo', ls=' ', alpha=0.7, mew=0.5, label='med band subs', ms=8)
    ax.plot( [0,100], [0,100], 'k-' )
    ax.set_xlim( 0, 15 )
    ax.set_ylim( 0, 15 )
    ax.set_xlabel( 'Optimal S/N from ETC' )
    ax.set_ylabel( 'Observed S/N from PSF drops' )

    ax = pl.subplot( 2, 2, 3 )
    ax.plot( SNRetc[iwide], SNRpsf[iwide], color='darkred', marker='d', ls=' ', alpha=0.5, mew=0.5, ms=8, label='wide band subs' )
    ax.plot( SNRetc[imed], SNRpsf[imed], 'co', ls=' ', alpha=0.7, mew=0.5, ms=8, label='med band subs' )
    ax.plot( [0,200], [0,200], 'k-' )
    ax.set_xlim( 0, 45 )
    ax.set_ylim( 0, 45 )
    ax.set_xlabel( 'Aperture S/N from ETC' )
    ax.set_ylabel( 'Observed S/N from PSF drops' )
    ax.legend( numpoints=1, loc='lower right' )

    ax = pl.subplot( 2,2,4 )
    ax.plot( SNRetc[iwide], SNRpsf[iwide], color='darkred', marker='d', ls=' ', alpha=0.5, mew=0.5,  ms=8)
    ax.plot( SNRetc[imed], SNRpsf[imed], 'co', ls=' ', alpha=0.7, mew=0.5, ms=8)
    ax.plot( [0,100], [0,100], 'k-' )
    ax.set_xlim( 0, 15 )
    ax.set_ylim( 0, 15 )
    ax.set_xlabel( 'Aperture S/N from ETC' )
    ax.set_ylabel( 'Observed S/N from PSF drops' )

    pl.draw()




def plotSNR_vs_mag() :
    """ make a plot comparing S/N vs mag for wide and med bands
    :return:
    """
    from astropy.io import ascii
    import numpy as np
    from matplotlib import pyplot as pl
    from pytools import plotsetup

    fig = plotsetup.presfig()
    psfdat = ascii.read( 'psfphot3.dat', format='fixed_width')

    isub = np.where(psfdat['image']=='sub')[0]

    oSNRetc = psfdat['optSNRetc'][isub]
    SNRetc = psfdat['SNRetc'][isub]
    SNRpsf = psfdat['SNRpsf'][isub]
    psfmag = psfdat['stack'][isub]
    band = psfdat['band'][isub]

    iwide = np.array([ i for i in range(len(band)) if band[i].endswith('w')] )
    imed = np.array([ i for i in range(len(band)) if band[i].endswith('m')] )

    pl.clf()
    ax = pl.subplot( 1,1,1 )
    ax.plot( psfmag[iwide], SNRpsf[iwide], color='darkorange', marker='d', ls=' ', alpha=0.8, mew=0.5, ms=8, label='wide band subs' )
    ax.plot( psfmag[imed], SNRpsf[imed], 'mo', ls=' ', alpha=0.8, mew=0.5, ms=8, label='med band subs' )
    ax.set_xlim( 25, 30)
    ax.set_ylim( 0, 30)
    ax.set_xlabel( 'Observed PSF-fitting mag (AB)' )
    ax.set_ylabel( 'Observed S/N from PSF drops' )
    ax.legend( numpoints=1, loc='upper right' )

    pl.draw()



def plothostless_compare():
    """
    make a plot comparing S/N and observed mag for reg and sub in hostless
    """
    from astropy.io import ascii
    import numpy as np
    from matplotlib import pyplot as pl
    from pytools import plotsetup

    psfdat = ascii.read( 'psfphot3.dat', format='fixed_width')
    snlist = psfdat['sn']
    bandlist = psfdat['band']
    imagelist = psfdat['image']
    epochlist = psfdat['epoch']

    fig = plotsetup.presfig()
    pl.clf()
    ax1 = pl.subplot( 1,2,1 )
    ax2 = pl.subplot( 1,2,2 )

    for thissn,color in zip( ['stone'], ['red'] ) :

        SNRsub, magsub, SNRreg, magreg, bandmatch = [],[],[],[], []
        for row in psfdat :
            sn = row['sn']
            if thissn=='hosted' and sn in ['stone','dawes'] : continue
            if thissn!='hosted' and sn != thissn : continue
            band = row['band']
            epoch = row['epoch']
            image = row['image']
            if image=='reg' : continue
            imatch = np.where( (bandlist==band) & (snlist==sn) & (epochlist==epoch) & (imagelist=='reg'))[0]
            if len(imatch) == 1 :
                SNRsub.append( row['SNRpsf'])
                # magsub.append( row['stack'])
                magsub.append( row['aperture'])
                SNRreg.append( psfdat['SNRpsf'][imatch[0]])
                # magreg.append( psfdat['stack'][imatch[0]])
                magreg.append( psfdat['aperture'][imatch[0]])
                bandmatch.append( band )
        SNRsub, magsub, SNRreg, magreg, bandmatch =  np.array(SNRsub), np.array(magsub), np.array(SNRreg), np.array(magreg), np.array(bandmatch)

        # imed = np.array([ i for i in range(len(bandmatch)) if bandmatch[i].endswith('m')])
        # iwide = np.array([ i for i in range(len(bandmatch)) if bandmatch[i].endswith('w')])

        imed = np.array([ i for i in range(len(bandmatch)) if bandmatch[i]=='f125w'])
        iwide = np.array([ i for i in range(len(bandmatch)) if bandmatch[i]!='f125w'])


        ax1.plot( magsub[iwide], magreg[iwide]-magsub[iwide], color=color, marker='d', ls=' ', alpha=0.8, mew=0.5, ms=8, label='%s other'%thissn )
        if len(imed) :
            ax1.plot( magsub[imed], magreg[imed]-magsub[imed], color=color, marker='o', ls=' ', alpha=0.8, mew=0.5, ms=12, label='%s f125w'%thissn )
        ax1.set_xlim( 24.8, 28.1)
        ax1.set_ylim( 1.1, -1.1)
        ax1.set_ylabel( 'observed difference :  reg mag - sub mag' )
        ax1.set_xlabel( 'sub mag (AB)' )

        ax2.plot( SNRsub[iwide], SNRreg[iwide]-SNRsub[iwide], color=color, marker='d', ls=' ', alpha=0.8, mew=0.5, ms=8, label='other' )
        if len(imed) :
            ax2.plot( SNRsub[imed], SNRreg[imed]-SNRsub[imed], color=color, marker='o', ls=' ', alpha=0.8, mew=0.5, ms=12, label='f125w' )
        ax2.set_xlim( 0, 25)
        ax2.set_ylim( -5, 11)
        ax2.set_ylabel( 'observed difference :  reg S/N - sub S/N' )
        ax2.set_xlabel( 'sub S/N' )

    ax1.text( 28, -0.02, 'SN appears brighter in reg than in sub', color='0.3', ha='right', va='bottom',fontsize=12)
    ax1.text( 28, 0.02, 'SN appears fainter in reg than in sub', color='0.3', ha='right', va='top',fontsize=12)
    ax1.axhline( 0, color='0.5', ls='--', )

    ax2.axhline( 0, color='0.5', ls='--', )
    ax2.text( 25, 0.5, 'S/N is better in reg than in sub', color='0.3', ha='right', va='bottom',fontsize=12)
    ax2.text( 25, -0.5, 'S/N is worse in reg than in sub', color='0.3', ha='right', va='top',fontsize=12)
    ax1.legend( numpoints=1, loc='lower left')
    pl.draw()


def plot_stone_newap():
    """
    make a plot comparing S/N and observed mag for reg and sub in hostless
    """
    from astropy.io import ascii
    import numpy as np
    from matplotlib import pyplot as pl
    from pytools import plotsetup

    psfdat = ascii.read( 'psfphot_stonemod.dat', format='fixed_width')
    snlist = psfdat['sn']
    bandlist = psfdat['band']
    imagelist = psfdat['image']
    epochlist = psfdat['epoch']

    aplist = psfdat['aperture']
    aperrlist = psfdat['aperr']

    fig = plotsetup.presfig()
    isub = np.where( imagelist=='sub' )[0]
    ireg = np.where( imagelist=='reg' )[0]

    ax1 = pl.subplot( 1,1,1 )
    ax1.errorbar(  aplist[isub], aplist[ireg], aperrlist[ireg], aperrlist[isub], ls=' ', marker='o', color='cyan', ms=10, capsize=0, alpha=0.5 )

    ax1.set_xlabel('sub image')
    ax1.set_ylabel('reg image')
    ax1.plot( [25,28],[25,28], 'k-')



def plot_stone_ap_vs_psf():
    """
    make a plot comparing S/N and observed mag for reg and sub photometry
    for SN stone, by DS and SR

    """
    from astropy.io import ascii
    import numpy as np
    from matplotlib import pyplot as pl
    from pytools import plotsetup
    plotsetup.presfig( figsize=[10,10])

    dat = ascii.read( 'photcompare_stone.dat', format='commented_header', header_start=-1, data_start=0)
    iapSRbigreg = np.where( (dat['PHOTTYPE']=='apSR') &
                            (dat['SKYANN']=='big') &
                            (dat['IMTYPE']=='reg')  )
    iapSRsmallreg = np.where( (dat['PHOTTYPE']=='apSR') &
                            (dat['SKYANN']=='small') &
                            (dat['IMTYPE']=='reg')  )

    iapSRbigsub = np.where( (dat['PHOTTYPE']=='apSR') &
                            (dat['SKYANN']=='big') &
                            (dat['IMTYPE']=='sub')  )
    iapSRsmallsub = np.where( (dat['PHOTTYPE']=='apSR') &
                            (dat['SKYANN']=='small') &
                            (dat['IMTYPE']=='sub')  )

    iapDSreg = np.where( (dat['PHOTTYPE']=='apDS') & (dat['IMTYPE']=='reg') )
    iapDSsub = np.where( (dat['PHOTTYPE']=='apDS') & (dat['IMTYPE']=='sub') )

    ipsfDSreg = np.where( (dat['PHOTTYPE']=='psfDS') & (dat['IMTYPE']=='reg') )
    ipsfDSsub = np.where( (dat['PHOTTYPE']=='psfDS') & (dat['IMTYPE']=='sub') )

    mag = dat['MAG']
    magerr = dat['MAGERR']
    flux = dat['FLUX']
    fluxerr = dat['FLUXERR']
    SNR = np.abs( flux / fluxerr )

    plotargs = {'alpha':0.8, 'mew':0 }
    pl.plot( mag[iapSRbigsub], mag[iapSRbigreg] - mag[iapSRbigsub], 'ro', label='ap big sky ann. (SR)', **plotargs )
    pl.plot( mag[iapSRsmallsub], mag[iapSRsmallreg] - mag[iapSRsmallsub], 'go', label='ap small sky ann. (SR)', **plotargs  )

    pl.plot( mag[iapDSsub], mag[iapDSreg] - mag[iapDSsub], 'mD', label='ap (DS)', **plotargs  )
    pl.plot( mag[ipsfDSsub], mag[ipsfDSreg] - mag[ipsfDSsub], 'cD', label='psf (DS)', **plotargs  )

    ax1 = pl.gca()
    ax1.set_xlabel('sub mag [AB]')
    ax1.set_ylabel('observed difference : reg mag - sub mag')
    ax1.axhline( 0.0, ls='--', color='0.5')
    ax1.legend( loc='upper left')
    ax1.invert_yaxis()
    ax1.text( 26.8, -0.02, 'SN appears brighter in reg than in sub', color='0.3', ha='right', va='bottom',fontsize=12)
    ax1.text( 26.8, 0.02, 'SN appears fainter in reg than in sub', color='0.3', ha='right', va='top',fontsize=12)
    ax1.set_ylim( 0.4, -0.4)
    ax1.set_xlim( 25, 27.7)


def plot_stone_template_tests():
    """
    make a plot comparing S/N and observed mag for reg and sub photometry
    for SN stone, by DS and SR

    """
    from astropy.io import ascii
    import numpy as np
    from matplotlib import pyplot as pl
    from pytools import plotsetup
    plotsetup.presfig( figsize=[10,10])

    dat = ascii.read( 'stone_template_phot.dat', format='commented_header', header_start=-1, data_start=0)
    iapSRbigreg = np.where( (dat['PHOTTYPE']=='apSR') &
                            (dat['SKYANN']=='big') &
                            (dat['IMTYPE']=='reg')  )
    iapSRsmallreg = np.where( (dat['PHOTTYPE']=='apSR') &
                            (dat['SKYANN']=='small') &
                            (dat['IMTYPE']=='reg')  )

    iapSRbigsub = np.where( (dat['PHOTTYPE']=='apSR') &
                            (dat['SKYANN']=='big') &
                            (dat['IMTYPE']=='sub')  )
    iapSRsmallsub = np.where( (dat['PHOTTYPE']=='apSR') &
                            (dat['SKYANN']=='small') &
                            (dat['IMTYPE']=='sub')  )


    iapSRbigreg10 = np.where( (dat['EPOCH']<20) & (dat['SKYANN']=='big') & (dat['IMTYPE']=='reg')  )
    iapSRsmallreg10 = np.where( (dat['EPOCH']<20) & (dat['SKYANN']=='small') & (dat['IMTYPE']=='reg')  )
    iapSRbigreg20 = np.where( (dat['EPOCH']>=20) & (dat['EPOCH']<30) & (dat['SKYANN']=='big') & (dat['IMTYPE']=='reg')  )
    iapSRsmallreg20 = np.where( (dat['EPOCH']>=20) &  (dat['EPOCH']<30) & (dat['SKYANN']=='small') & (dat['IMTYPE']=='reg')  )
    iapSRbigreg30 = np.where(   (dat['EPOCH']<40) & (dat['EPOCH']>=30) & (dat['SKYANN']=='big') & (dat['IMTYPE']=='reg')  )
    iapSRsmallreg30 = np.where( (dat['EPOCH']<40) & (dat['EPOCH']>=30) & (dat['SKYANN']=='small') & (dat['IMTYPE']=='reg')  )
    iapSRbigreg40 = np.where( (dat['EPOCH']>=40) & (dat['SKYANN']=='big') & (dat['IMTYPE']=='reg')  )
    iapSRsmallreg40 = np.where( (dat['EPOCH']>=40) & (dat['SKYANN']=='small') & (dat['IMTYPE']=='reg')  )

    iapSRbigsub10 = np.where( (dat['EPOCH']<20) & (dat['SKYANN']=='big') & (dat['IMTYPE']=='sub')  )
    iapSRsmallsub10 = np.where( (dat['EPOCH']<20) & (dat['SKYANN']=='small') & (dat['IMTYPE']=='sub')  )
    iapSRbigsub20 = np.where( (dat['EPOCH']>=20) & (dat['EPOCH']<30) & (dat['SKYANN']=='big') & (dat['IMTYPE']=='sub')  )
    iapSRsmallsub20 = np.where( (dat['EPOCH']>=20) &  (dat['EPOCH']<30) & (dat['SKYANN']=='small') & (dat['IMTYPE']=='sub')  )
    iapSRbigsub30 = np.where(   (dat['EPOCH']<40) & (dat['EPOCH']>=30) & (dat['SKYANN']=='big') & (dat['IMTYPE']=='sub')  )
    iapSRsmallsub30 = np.where( (dat['EPOCH']<40) & (dat['EPOCH']>=30) & (dat['SKYANN']=='small') & (dat['IMTYPE']=='sub')  )
    #iapSRbigsub40 = np.where( (dat['EPOCH']>=40) & (dat['SKYANN']=='big') & (dat['IMTYPE']=='sub')  )
    #iapSRsmallsub40 = np.where( (dat['EPOCH']>=40) & (dat['SKYANN']=='small') & (dat['IMTYPE']=='sub')  )


    sky = dat['SKY']
    skyerr = dat['SKYERR']

    flux = dat['FLUX']
    fluxerr = dat['FLUXERR']
    SNR = np.abs( flux / fluxerr )

    plotargs = {'alpha':0.8, 'mew':0,  'ls':' '}

    #fig1 = pl.figure(1)
    #pl.clf()

    # ax1 = pl.subplot(1,1,1)
    # pl.errorbar( sky[iapSRsmallreg10], sky[iapSRbigreg10] , skyerr[iapSRbigreg10], skyerr[iapSRsmallreg10],'bo', label='1 visit', **plotargs )
    # pl.errorbar( sky[iapSRsmallreg20], sky[iapSRbigreg20] , skyerr[iapSRbigreg20], skyerr[iapSRsmallreg20],'go', label='2 visits', **plotargs )
    # pl.errorbar( sky[iapSRsmallreg30], sky[iapSRbigreg30] , skyerr[iapSRbigreg30], skyerr[iapSRsmallreg30],'ro', label='3 visits', **plotargs )
    # pl.errorbar( sky[iapSRsmallreg40], sky[iapSRbigreg40] , skyerr[iapSRbigreg40], skyerr[iapSRsmallreg40],'kD', label='4 visits', **plotargs )
    # pl.plot( [-10,10],[-10,10],'k--')
    # ax1.axhline(0.0,color='k',ls=':')
    # ax1.set_xlim(0.45,0.7)
    # ax1.set_ylim(0.45,0.7)
    # ax1.set_xlabel('sky in alt template reg image (big annulus)')
    # ax1.set_ylabel('sky in alt template reg image (small annulus)')
    # ax1.legend( loc='upper left')
    # ax1.set_title('alternate drizzle combinations for stone template')
    # 
    # fig2 = pl.figure(2)
    # pl.clf()
    # pl.errorbar( sky[iapSRsmallreg10] - sky[iapSRbigreg10] , sky[iapSRsmallsub10]-sky[iapSRbigsub10], skyerr[iapSRsmallreg10], skyerr[iapSRbigreg10], 'bo', label='1 visit', **plotargs )
    # pl.errorbar( sky[iapSRsmallreg20] - sky[iapSRbigreg20] , sky[iapSRsmallsub20]-sky[iapSRbigsub20], skyerr[iapSRsmallreg20], skyerr[iapSRbigreg20], 'go', label='2 visits', **plotargs )
    # pl.errorbar( sky[iapSRsmallreg30] - sky[iapSRbigreg30] , sky[iapSRsmallsub30]-sky[iapSRbigsub30], skyerr[iapSRsmallreg30], skyerr[iapSRbigreg30], 'ro', label='3 visits', **plotargs )
    # # pl.errorbar( sky[iapSRsmallreg40] - sky[iapSRbigreg40] , sky[iapSRsmallsub40]-sky[iapSRbigsub40], skyerr[iapSRsmallreg40], skyerr[iapSRbigreg40], 'kD', label='4 visits', **plotargs )
    # pl.ylabel( 'Fsky[small] - Fsky[big] from alt template diff images' )
    # pl.xlabel( 'Fsky[small] - Fsky[big] from alt template reg images' )

    import numpy as np
    A  = np.pi * (0.3/0.09)**2
    trueflux = flux + sky * A 
    fig3 = pl.figure(3)
    pl.clf()
    # pl.errorbar( sky[iapSRsmallreg10] - sky[iapSRbigreg10] , trueflux[iapSRsmallreg10]-trueflux[iapSRbigreg10] ,  fluxerr[iapSRsmallreg10], skyerr[iapSRsmallreg10], color='b',marker='o', label='1 visit', **plotargs )
    # pl.errorbar( sky[iapSRsmallreg20] - sky[iapSRbigreg20] , trueflux[iapSRsmallreg20]-trueflux[iapSRbigreg20] ,  fluxerr[iapSRsmallreg20], skyerr[iapSRsmallreg20], color='g',marker='o', label='2 visits', **plotargs )
    # pl.errorbar( sky[iapSRsmallreg30] - sky[iapSRbigreg30] , trueflux[iapSRsmallreg30]-trueflux[iapSRbigreg30] ,  fluxerr[iapSRsmallreg30], skyerr[iapSRsmallreg30], color='r',marker='o', label='3 visits', **plotargs )
    # pl.errorbar( sky[iapSRsmallreg40] - sky[iapSRbigreg40] , trueflux[iapSRsmallreg40]-trueflux[iapSRbigreg40] ,  fluxerr[iapSRsmallreg40], skyerr[iapSRsmallreg40], color='k',marker='D', label='4 visits', **plotargs )

    pl.plot(  -(sky[iapSRsmallreg10] - sky[iapSRbigreg10])*A, (flux[iapSRsmallreg10]-flux[iapSRbigreg10])/1.22074468 ,  color='b',marker='o', label='1 visit', **plotargs )
    pl.plot(  -(sky[iapSRsmallreg20] - sky[iapSRbigreg20])*A, (flux[iapSRsmallreg20]-flux[iapSRbigreg20])/1.22074468 ,  color='g',marker='o', label='2 visits', **plotargs )
    pl.plot(  -(sky[iapSRsmallreg30] - sky[iapSRbigreg30])*A, (flux[iapSRsmallreg30]-flux[iapSRbigreg30])/1.22074468 ,  color='r',marker='o', label='3 visits', **plotargs )
    pl.plot(  -(sky[iapSRsmallreg40] - sky[iapSRbigreg40])*A, (flux[iapSRsmallreg40]-flux[iapSRbigreg40])/1.22074468 ,  color='k',marker='D', label='4 visits', **plotargs )

    pl.plot([-0.35,-0.05],[-0.35,-0.05])
    # pl.xlabel( 'Total flux in 0.3" aperture ' )
    # pl.ylabel( 'Ftot[smallsky] - Ftot[bigsky] \n after adding back in the sky flux' )




def plot_stone_deltaflux():
    """
    make a plot showing sky value variations for sn stone
    """
    from astropy.io import ascii
    import numpy as np
    from matplotlib import pyplot as pl
    from pytools import plotsetup
    plotsetup.presfig( figsize=[10,10])

    dat = ascii.read( 'stone_skyMMM_phot.dat', format='commented_header', header_start=-1, data_start=0)


    fig1 = pl.figure(1)
    pl.clf()
    # ax1 = pl.subplot(1,1,1)

    for filt,color,apcor ,irow in zip( ['F125W','F160W'],['b','r'],[0.185,0.194],[0,1])  :
        fapcor = 10**(-0.4*apcor)
        for marker,skyann in zip( ['^','s','o','D'], [ 0.5,6.0 ] ): # ,0.6,1.0, 6.0 ]) : # 0.5,0.6,1.0,6.0]) :
            #ibig = np.where( (dat['SKYANN']==6.0)  & (dat['FILTER']==filt) & (dat['IMTYPE']=='reg'))
            #ismall = np.where( (dat['SKYANN']==skyann) & (dat['FILTER']==filt)  & (dat['IMTYPE']=='reg') )

            if skyann==0.5 :
                fig1 = pl.figure(1)
                fig = fig1
            else :
                fig2 = pl.figure(2)
                fig = fig2

            ireg = np.where( (dat['SKYANN']==skyann) & (dat['FILTER']==filt) & (dat['IMTYPE']=='reg') & (dat['MJD']>56400) )
            isub = np.where( (dat['SKYANN']==skyann) & (dat['FILTER']==filt) & (dat['IMTYPE']=='sub') & (dat['MJD']>56400) )
            itmp = np.where( (dat['SKYANN']==skyann) & (dat['FILTER']==filt) & (dat['IMTYPE']=='reg') & (dat['MJD']<56400) )[0]

            mjd = dat['MJD'][ireg]
            skyreg = dat['SKY'][ireg]
            skysub = dat['SKY'][isub]
            skytmp = dat['SKY'][itmp]

            skyerrreg = dat['SKYERR'][ireg]
            skyerrsub = dat['SKYERR'][isub]
            skyerrtmp = dat['SKYERR'][itmp]

            fluxreg = dat['FLUX'][ireg]
            fluxsub = dat['FLUX'][isub]
            fluxtmp = dat['FLUX'][itmp]

            fluxerrreg = dat['FLUXERR'][ireg]
            fluxerrsub = dat['FLUXERR'][isub]
            fluxerrtmp = dat['FLUXERR'][itmp]

            # area = np.pi * ( 0.4 / 0.09 )**2
            area = 62.0562


            ftottmp = fluxtmp + skytmp * area
            ftotreg = fluxreg + skyreg * area

            ftotsub = fluxsub + skysub * area

            deltaftot = ftotsub + ftottmp - ftotreg

            #ftoterr = fluxerrtmp + skyerrtmp * area
            # import pdb; pdb.set_trace()

            ftotsub = fluxsub + (skysub*area) + fluxtmp +  skytmp * area
            ftotreg = fluxreg + (skyreg*area)

            deltaflux = fluxsub - fluxreg
            deltafluxerr = np.sqrt( fluxerrsub**2  + fluxerrreg**2 )

            plotargs = {'alpha':0.5, 'mew':0 }

            print( ' "source" Flux in template = %.3f (filter=%s, skyann=%.1f)'%(fluxtmp,filt,skyann))
            print( ' "sky" Flux in template = %.3f (filter=%s, skyann=%.1f)'%(skytmp,filt,skyann))

            # pl.errorbar( mjd+skyann*10, deltaflux, deltafluxerr,

            # Difference between reg and (sub+template) source flux:
            #pl.errorbar( mjd+skyann, fluxreg-fluxsub-fluxtmp, fluxerrreg, ls='-',
            #             marker=marker,color=color,
            #             label=filt+' %.1f'%skyann, **plotargs )
            mjdbins = np.arange(len(mjd))*3 # np.append( mjd, mjd[-1]+5 )

            ax1 = fig.add_subplot(2,2,irow*2 + 1 )

            # total flux in the reg images
            ax1.bar( mjdbins-1, fluxreg, width=1, color='darkorange', label='reg', alpha=0.5 )
            # total pre-sub flux in the diff images
            ax1.bar( mjdbins, fluxsub, width=1, label='sub', color='b', alpha=0.5 )
            ax1.bar( mjdbins, np.zeros(len(mjd))+fluxtmp, bottom=fluxsub, width=1, color='r', label='tmp', alpha=0.5 )
            ax1.set_title( 'measured "source" flux in %s'%filt )

            ax2 = fig.add_subplot(2,2,irow*2 + 2 )
            ax2.bar( mjdbins-1, skyreg, width=1, color='darkorange', label='reg', alpha=0.5 )
            ax2.bar( mjdbins, skysub, width=1, color='b', label='sub', alpha=0.5 )
            ax2.bar( mjdbins, np.zeros(len(mjd))+skytmp, bottom=skysub, width=1, color='r', label='tmp', alpha=0.5 )
            ax2.set_title( 'measured sky in %s'%filt )
    for fig in [fig1,fig2]:
        ax = fig.add_subplot( 2, 2, 1 )
        ax.legend( loc='upper left', frameon=False, fontsize='small')

    fig1.suptitle('Sky annulus = [0.5",1.0"]')
    fig1.subplots_adjust(left=0.12, bottom=0.12, top=0.88, right=0.95, wspace=0.25, hspace=0.25)
    fig2.suptitle('Sky annulus = [6.0",12.0"]')
    fig2.subplots_adjust(left=0.12, bottom=0.12, top=0.88, right=0.95, wspace=0.25, hspace=0.25)

    # pl.xlabel( 'MJD' )
    # pl.ylabel( '$\Delta$Flux [sub-reg] / TemplateFlux' )
    pl.draw()
