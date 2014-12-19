__author__ = 'rodney'


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


def get_host_ellipses( cat1=None, xSN=7793.1145, ySN=12619.196):
    """
    get ellipse parameters for all sources within 10 arcsec = 170 pixels
    :return:
    """
    import numpy as np
    from astropy.io import ascii

    if cat1 is None :
        cat1 = ascii.read('/store/goods-n/catalogs/gn_all_candels_sx_h_130829_hphotom_comb_merge_psfmatch2h.cat')


    # identify the catalog indices for hosts A and B
    x = cat1['X_IMAGE']
    y = cat1['Y_IMAGE']
    dist = np.sqrt( (x-xSN)**2 + (y-ySN)**2 )

    iclose = np.where( dist<170 )[0]


    for i in iclose :
        # extract the ellipse parameters :
        Cxx = cat1['CXX_IMAGE'][i]
        Cyy = cat1['CYY_IMAGE'][i]
        Cxy = cat1['CXY_IMAGE'][i]
        x = cat1['X_IMAGE'][i]
        y = cat1['Y_IMAGE'][i]

        xr = xSN - x
        yr = ySN - y
        R = np.sqrt(  Cxx*xr**2 + Cyy*yr**2 + Cxy*xr*yr )
        theta = cat1['THETA_IMAGE'][i]
        A = cat1['A_IMAGE'][i]
        B = cat1['B_IMAGE'][i]
        ra = cat1['ALPHA_J2000'][i]
        dec = cat1['DELTA_J2000'][i]

        darcsec = dist[i] * 0.06
        print( "index :  %i"%    cat1['NUMBER'][i] )
        print( "  RA    :  %.6f"%ra )
        print( "  Dec   :  %.6f"%dec )
        print( '  d (") :  %.3f'%darcsec )
        print( "  x_gal :  %.3f"%x )
        print( "  y_gal :  %.3f"%y )
        print( "  theta :  %.3f"%theta )
        print( "  A     :  %.3f"%A )
        print( "  B     :  %.3f"%B )
        print( "  R     :  %.3f"%R )
        print( "-------------------------")
        print("")

