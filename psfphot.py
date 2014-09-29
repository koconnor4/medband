# 2014.08.29
# S.Rodney
# read in psf photometry and write it out to .dat files for SNANA or sncosmo


bushheader = """
SURVEY: CANDELS
SNID: 1 # bush
PHOTOMETRY_VERSION: highzIaLCs%i
FILTERS: YJNHLOPQWVIZ
MAGTYPE: LOG10
MAGREF: AB
RA: 53.178232  deg
DECL: -27.8019730  deg
MWEBV: 0.0080    MW E(B-V)
REDSHIFT_FINAL: 1.76 +- 0.53
SEARCH_PEAKMJD: 55803.1
SEARCH_PEAKMJDERR: 15.0

NOBS: %i
NVAR: 7
VARLIST:  MJD  FLT FIELD   FLUXCAL   FLUXCALERR    MAG     MAGERR   ZPT
"""

colfaxheader = """
SURVEY: CANDELS
SNID: 2  # colfax
PHOTOMETRY_VERSION: highzIaLCs%i
FILTERS: HJNOPQWVIZ
MAGTYPE: LOG10
MAGREF: AB
RA: 189.15659    deg
DECL: 62.309172  deg
MWEBV: 0.012    MW E(B-V)
REDSHIFT_FINAL: 2.13 +- 0.2
HOST_GALAXY_PHOTO-Z: 2.13 +- 0.2
REDSHIFT_STATUS: OK
SEARCH_PEAKMJD: 56074.3
SEARCH_PEAKMJDERR: 15.0

NOBS: %i
NVAR: 7
VARLIST:  MJD  FLT FIELD   FLUXCAL   FLUXCALERR    MAG     MAGERR   ZPT
"""

stoneheader = """
SURVEY: CANDELS
SNID: 3  # stone
PHOTOMETRY_VERSION: highzIaLCs%i
FILTERS: HJNPQWIZ
RA: 189.31989    deg
DECL: 62.278179  deg
MWEBV: 0.0104
REDSHIFT_FINAL: 1.90 +- 0.3
SEARCH_PEAKMJD: 56467.5
SEARCH_PEAKMJDERR: 15.0

NOBS: %i
NVAR: 7
VARLIST:  MJD  FLT FIELD   FLUXCAL   FLUXCALERR    MAG     MAGERR   ZPT
"""

zpt = {'F105W':26.0974,'F110W':26.6424,'F125W':26.0449,'F140W':26.2608,
       'F160W':25.7551,'F098M':25.5041,'F127M':24.4545,'F139M':24.2880,
       'F153M':24.2725,
       'F435W': 25.65777,'F475W': 26.05923,'F555W': 25.7184,'F606W': 26.49113,
       'F625W': 25.90667,'F775W': 25.66504,'F814W': 25.94333,'F850L':24.84245,
       'F218W':22.7776,'F225W':23.8629,'F275W':23.9740,'F336W':24.5377,
       'F350L':26.8413,'F390W':25.2389,'F438W':24.7097,'F600L':25.7681,
       'F763M':24.2070,'F845M':23.7811,
       }

FilterAlpha = { 'unknown':'?',
                'F225W':'S','F275W':'T','F336W':'U','F390W':'C',
                'F350L':'W',
                'F435W':'B','F475W':'G','F606W':'V','F625W':'R',
                'F775W':'X','F814W':'I','F850L':'Z',
                'F125W':'J','F160W':'H','F125W+F160W':'A',
                'F105W':'Y','F110W':'M','F140W':'N',
                'F098M':'L','F127M':'O','F139M':'P','F153M':'Q',
                'G141':'4','G102':'2','blank':'0',
                'F763M':'7','F845M':'9',
                }

def mk_all_datfiles( ):
    masterphotfile=mk_master_phot_file()
    mk_datfiles('bush', masterphotfile=masterphotfile)
    mk_datfiles('colfax', masterphotfile=masterphotfile)
    mk_datfiles('stone', masterphotfile=masterphotfile)


def mk_datfiles( sn='colfax', masterphotfile='psfphot_medband.dat') :
    """  Write out 4 SNANA-style .dat files for the given SN, using the four
    flavors of photometry in the master .dat file.
    :param sn:
    :return:
    """
    from astropy.io import ascii
    import numpy as np

    psfdat = ascii.read(masterphotfile, format='fixed_width')

    isn = np.where( (psfdat['sn']==sn) )[0]
    nobs = len(isn)

    if sn == 'bush' : header = bushheader
    elif sn == 'colfax' : header = colfaxheader
    elif sn == 'stone' : header = stoneheader

    for iphot,sfx in zip( [2,4], ['_ap','_psf']) :
        snana = open('HST_CANDELS%i_%s.dat'%(iphot,sn),'w')
        sncosmo = open('HST_CANDELS%i_%s.sncosmo.dat'%(iphot,sn),'w')
        hdr = header.lstrip()%(iphot,nobs)
        print >> snana, hdr
        for hdrline in hdr.split('\n') :
            if len(hdrline.strip()) == 0 : break
            print >>sncosmo,  "#" + hdrline
        print >> sncosmo, '#%8s  %6s   %8s  %8s   %7s  %7s  %6s  %6s'%(
            'mjd','filter','flux','fluxerr','mag','magerr','magsys','zpt')

        for i in isn :
            mjd = psfdat['mjd'][i]
            filter = psfdat['band'][i]
            band = FilterAlpha[filter.upper()]
            if filter.endswith('l') : filter += 'p'
            zpt = psfdat['zpt_AB'][i]
            mag = psfdat['mag%s'%sfx][i]
            magerr = psfdat['magerr%s'%sfx][i]
            flux = psfdat['flux%s'%sfx][i]
            fluxerr = psfdat['fluxerr%s'%sfx][i]
            fluxcal = flux * 10**(0.4*(27.5-zpt))
            fluxcalerr = fluxerr * 10**(0.4*(27.5-zpt))
            snanaline = 'OBS: %9.3f  %s  %s  %8.4f  %8.4f   %7.3f %7.3f  %.3f'%(
                mjd, band, sn, fluxcal, fluxcalerr, mag, magerr, zpt )
            print >> snana, snanaline
            sncosmoline = '%9.3f  %6s   %8.4f  %8.4f   %7.3f  %7.3f      AB  %6.3f'%(
                mjd, filter, flux, fluxerr, mag, magerr,  zpt )
            print >> sncosmo, sncosmoline
        print >> snana, 'END:'
        snana.close()
        sncosmo.close()



def mk_master_phot_file( infile='round2sn_force_all.txt', outfile='psfphot_medband.dat',
                        doapphot=True, acsfile='apphot_medband_acsuvis.dat' ):
    """  read in a .txt file from Dan and convert it to a human-readable
    formatted ascii file
    :return:
    """
    from astropy.io import ascii
    import pyfits
    import numpy as np
    from astropy.table import Column, vstack
    import os
    from hstphot import hstapphot

    irtab = ascii.read( infile, format='commented_header', header_start=-1, data_start=0 )
    fileroots = [ os.path.splitext(os.path.basename(fn))[0]  for fn in  irtab['file'] ]
    iremove = np.where( [ (fr.split('_')[0] not in ['bush','stone','colfax']) or (fr.split('_')[3] == 'reg') for fr in fileroots ] )[0]
    irtab.remove_rows( iremove )

    fileroots = [ os.path.splitext(os.path.basename(fn))[0]  for fn in  irtab['file'] ]
    sncol = Column( name='sn', data=[ fr.split('_')[0] for fr in fileroots ] )
    bandcol = Column( name='band', data=[ fr.split('_')[1] for fr in fileroots ] )
    epochcol = Column( name='epoch', data=[ fr.split('_')[2].split('-')[0] for fr in fileroots ] )
    imtypecol  = Column( name='imtype', data=[ fr.split('_')[3] for fr in fileroots ] )


    imfilelist =  [ os.path.join('../'+sn,fr+'.fits') for sn,fr in zip( sncol, fileroots ) ]
    nexpcol = Column( name='nexp', data=[ pyfits.getval( imfile, 'NCOMBINE') for imfile in imfilelist ], format='%i' )
    exptimecol = Column( name='exptime', data=[ pyfits.getval( imfile, 'EXPTIME') for imfile in imfilelist ], format='%i' )
    mjdcol = Column( name='mjd', data=[ pyfits.getval( imfile, 'EXPSTART')+(et/2./3600./24.) for imfile,et in zip( imfilelist, exptimecol ) ], format='%.2f' )
    irtab.add_columns( [ sncol, bandcol, epochcol, imtypecol, nexpcol, exptimecol, mjdcol ], indexes=[1,1,1,1,1,1,1] )

    irtab.rename_column( 'apphot', 'mag_apDS' )
    irtab.rename_column( 'apphoterr', 'magerr_apDS' )
    irtab.rename_column( 'psfphotstack_zp', 'zpt_psf' )
    irtab.rename_column( 'drop_error_mag', 'magerr_drop' )
    irtab.rename_column( 'sky', 'fsky_psf' )
    irtab.rename_column( 'skyerr', 'fskyerr_psf' )

    # psfphotstack_mag = psfphotstack_zp + 2.5*alog10(psfphotstack_psf )
    zptAB = { 'f098m':25.6672,'f105w':26.2688,'f125w':26.2302,'f127m': 24.641,
              'f139m':24.4793,'f140w':26.4523,'f153m':24.4634,'f160w':25.9462, }
    zptlist = np.array( [ zptAB[band] for band in irtab['band'] ] )

    # psfphotstack_mag = psfphotstack_zp - 2.5*alog10(psfphotstack_psf )

    fluxpsfcol = Column( name='flux_psf', data= irtab['psfphotstack_psf'] * 10**(-0.4*(irtab['zpt_psf'] - zptlist )), format='%.5f')
    fluxerrpsfcol = Column( name='fluxerr_psf', data= irtab['drop_error_flux'] * 10**(-0.4*(irtab['zpt_psf'] - zptlist )), format='%.5f')
    # fluxerrpsfcol = Column( name='fluxerr_psf', data= np.abs(irtab['fluxerr_psf']), format='%.5f')
    # fluxerrdropcol = Column( name='fluxerr_drop', data= np.abs(0.92103 * irtab['magerr_psf'] * 10**(-0.4*(irtab['drop_median']-zptlist))), format='%.5f')

    # flux_HST =  flux_psf  * 10**(-0.4 * ( zpt_psf - zpt_AB ) )
    # dflux = 0.4 * ln(10) * dmag * flux

    magpsfdat = np.where( fluxpsfcol>1e-10, -2.5*np.log10(fluxpsfcol.data) + zptlist, -2.5*np.log10( 3 * np.abs(fluxerrpsfcol.data) ) + zptlist )
    magpsfcol = Column( name='mag_psf', data=magpsfdat , format='%.4f')

    magerrpsfdat = np.where( fluxpsfcol>1e-10,(fluxerrpsfcol.data / fluxpsfcol.data ) / 0.92103403719761845, 1 / 0.92103403719761845 )
    magerrpsfcol = Column( name='magerr_psf', data=magerrpsfdat, format='%.4f')

    # import pdb; pdb.set_trace()

    irtab.add_column( magerrpsfcol, index=12 )
    irtab.add_column( magpsfcol, index=12 )
    irtab.add_column( fluxerrpsfcol, index=12 )
    irtab.add_column( fluxpsfcol, index=12 )

    if doapphot :
        apflux, apfluxerr, apmag, apmagerr, apzpt, apsky, apskyerr = [], [], [], [], [], [], []
        for imfile,ra,dec in zip( imfilelist, irtab['ra'], irtab['dec']) :
            x,y = hstapphot.radec2xy( imfile, ra,dec)
            apphotlist = hstapphot.dophot( imfile, x, y, aparcsec=0.3, system='AB', skyannarcsec=[6.0,12.0], skyalgorithm='sigmaclipping' )
            apphot = apphotlist[0].split()
            #apmjdobs = float(apphot[0])
            #apfiltname = apphot[1]
            #apaparcsec = float( apphotlist[2] )
            apflux.append( float(apphot[3]) )
            apfluxerr.append( float(apphot[4]) )
            apmag.append( float(apphot[5]) )
            apmagerr.append( float(apphot[6]) )
            apzpt.append( float(apphot[8]) )
            apsky.append( float(apphot[9]) )
            apskyerr.append( float(apphot[10]) )
        for apdat,apcolname in zip( [apflux, apfluxerr, apmag, apmagerr, apzpt, apsky, apskyerr],
                                    ['flux_ap', 'fluxerr_ap', 'mag_ap', 'magerr_ap', 'zpt_AB', 'fsky_ap', 'fskyerr_ap'] ) :
            newcol = Column( name=apcolname, data= apdat )
            irtab.add_column( newcol, index=8 )

    irtab.remove_columns( [ 'mag_apDS', 'magerr_apDS', 'file', 'filt', 'scale', 'aper', 'psfphotstack_psf', 'psf_model_mag', 'id', 'syntemp', 'x1','x2','y1','y2','ra','dec', 'psfphotidl3','psfphotidl_zp','psfphoterr3','psfphotidl_psf', 'drop_median',  'fsky_psf', 'fskyerr_psf',   'fsky_ap', 'fskyerr_ap', 'zpt_psf',]  )
    neworder = ('sn', 'band', 'epoch', 'mjd',   'mag_psf', 'magerr_psf', 'magerr_drop', 'mag_ap',  'magerr_ap', 'flux_psf', 'fluxerr_psf', 'flux_ap', 'fluxerr_ap',  'zpt_AB',  'nexp', 'exptime', 'imtype')
    newirtab = irtab[neworder]

    acstab = ascii.read( acsfile, format='fixed_width' )
    acsmagerrdropcol = Column( name='magerr_drop', data=acstab['magerr_ap'] )
    acstab.add_column( acsmagerrdropcol )
    outtab = vstack( [newirtab, acstab], join_type='outer' )
    outtab.sort( ['sn','band','epoch'] )

    fixtab = fixnans(outtab)
    fixtab.write( outfile, format='ascii.fixed_width')
    # outtab.write( outfile, format='ascii.fixed_width')
    return( outfile )


def fixnans( table ):
    import numpy as np
    for row in table :
        if np.isnan( row['mag_psf'] ) :
            row['mag_psf'] = row['mag_ap']
            row['magerr_psf'] = row['magerr_ap']
        if not np.isfinite( row['flux_psf'] ) :
            row['flux_psf'] = row['flux_ap']
            row['fluxerr_psf'] = row['fluxerr_ap']
        #if row['fluxerr_psf']/row['fluxerr_ap'] < 0.5 :
        #    row['fluxerr_psf'] = row['fluxerr_ap']

    return( table )



def add_flux_columns( ):
    """  read in a commentedheader .dat file hand-prepared by SR, and write
    out a pair of new ascii data files that includes new columns with the
    psf-fitting magnitudes converted to fluxes.
    :return:
    """
    from astropy.io import ascii
    import numpy as np
    from astropy.table import Column

    psfdat = ascii.read('psfphot_medband_commentedheader.dat', format='commented_header', header_start=-1, data_start=0)

    magerrpsf = psfdat['magerrpsf']
    magerrap = psfdat['magerrapSR']
    magap = psfdat['magapSR']
    fluxap = psfdat['fluxapSR']
    fluxerrap = psfdat['fluxerrapSR']
    zpt = psfdat['zpt']
    for col in ['magpsfstack', 'magpsfidl3','magpsfidl5'] :
        magpsf = psfdat[col]
        fluxpsf = np.where( magerrpsf>-1, 10**(-0.4*(magpsf.data-zpt.data)), fluxap )
        fluxpsfcol = Column( name=col.replace('mag','flux'), data=fluxpsf, format='%.4f' )
        fluxerrpsf = np.where( magerrpsf>-1, 0.92103 * magerrpsf * fluxpsf, fluxerrap )
        fluxerrpsfcol = Column( name=col.replace('mag','fluxerr'), data=fluxerrpsf, format='%.4f' )
        psfdat.add_column( fluxpsfcol )
        psfdat.add_column( fluxerrpsfcol )
    psfdat.write( 'psfphot_medband_fixedwidth.dat', format='ascii.fixed_width')
    psfdat.write( 'psfphot_medband_ascii.dat', format='ascii.commented_header')




def plotComparison( name='colfax' ):
    """  Make a comparison plot showing the light curves from each of the
    four photometry flavors.
    :param sn:
    :return:
    """
    import sncosmo
    # from sncosmost import hstbandpasses
    from astropy.io import ascii
    import numpy as np
    from matplotlib import pyplot as pl
    from pytools import plotsetup

    plotsetup.presfig( figsize = [20,12])
    fig = pl.gcf()
    pl.clf()

    for phottype, iphot, color, marker in zip( ['ap','psf'], [2,4], ['r','g'], ['D','o'] ):
        sndat = ascii.read( 'HST_CANDELS%i_%s.sncosmo.dat'%(iphot,name), format='commented_header', header_start=-1, data_start=0 )
        filterlist = np.unique(  [ filt for filt in sndat['filter'] if filt.lower().startswith('f1') ] )
        naxes = len(filterlist)
        nrow = int( np.sqrt( naxes ) )
        ncol = ( naxes / nrow ) + 1

        for filt, iax in zip( filterlist, range(1,naxes+1) ) :
            ax = fig.add_subplot( nrow, ncol, iax )
            ifilt = np.where( sndat['filter'] == filt )
            mjd = sndat['mjd'][ifilt] + iphot
            flux = sndat['flux'][ifilt]
            fluxerr = sndat['fluxerr'][ifilt]
            ax.errorbar( mjd, flux, fluxerr, marker=marker, color=color, ls='-', label=phottype )
            ax.set_xlabel('mjd')
            ax.set_ylabel('flux')
            ax.text( 0.02, 0.95, filt.upper(), fontsize=14, ha='left', va='top',transform=ax.transAxes )
            ax.axhline( 0, ls='--', color='0.5')

        fig.subplots_adjust( wspace=0.2, hspace=0.2, left=0.08, right=0.95, bottom=0.12, top=0.92 )
        ax.legend( loc='upper right', numpoints=1, bbox_to_anchor=[1.05,1.05] )
    fig.suptitle( name.upper() )
    pl.draw()




def plotLightCurve( name='colfax'):
    """ plot the ACS, UVIS and IR light curves
    :param name:
    :return:
    """
    import sncosmo
    # from sncosmost import hstbandpasses
    from astropy.io import ascii
    import numpy as np
    from matplotlib import pyplot as pl
    from pytools import plotsetup

    plotsetup.presfig( figsize = [20,12])
    fig = pl.gcf()
    pl.clf()


    sndat = ascii.read( 'HST_CANDELS4_%s.sncosmo.dat'%(name), format='commented_header', header_start=-1, data_start=0 )
    filterlist = np.unique(  [ filt for filt in sndat['FILTER'] ] )
    naxes = len(filterlist)
    nrow = int( np.sqrt( naxes ) )
    ncol = ( naxes / nrow ) + 1

    for filt, iax in zip( filterlist, range(1,naxes+1) ) :
        ax = fig.add_subplot( nrow, ncol, iax )
        ifilt = np.where( sndat['FILTER'] == filt )
        mjd = sndat['MJD'][ifilt]
        flux = sndat['FLUX'][ifilt]
        fluxerr = sndat['FLUXERR'][ifilt]
        mag = sndat['MAG'][ifilt]
        magerr = sndat['MAGERR'][ifilt]

        if filt.startswith('f1') :
            color='r'
            marker='D'
        elif filt.startswith('f3') :
            color='k'
            marker='o'
        else :
            color='g'
            marker='s'

        ax.errorbar( mjd, flux, fluxerr, marker=marker, color=color, ls='-' )
        ax.set_xlabel('MJD')
        ax.set_ylabel('FLUX')
        ax.text( 0.02, 0.95, filt.upper(), fontsize=20, ha='left', va='top',transform=ax.transAxes )
        # ax.invert_yaxis()
        ax.axhline( 0, ls='--', color='0.5')

    fig.subplots_adjust( wspace=0.2, hspace=0.2, left=0.08, right=0.95, bottom=0.12, top=0.92 )
    fig.suptitle( name.upper() )
    pl.draw()

