#! /usr/bin/env python
# S.Rodney
# 2014.05.04


def mkTemplates( clobber=False, verbose=False):
    """ Make med-band templates from sndrizpipe epoch 00 broad-band images.
    """
    from hstsntools import imageops, filters

    f125to127 = filters.computeScaling( 'F125W','F127M')
    print("F127M ~ F125W * %.3f"%f125to127)
    f160to153 = filters.computeScaling( 'F160W','F153M')
    print("F153M ~ F160W * %.3f"%f160to153)

    f125and160to140 = filters.computeScaling2to1( 'F125W','F160W','F140W')
    print("F140W ~ (F125W+F160W) * %.3f"%f125and160to140)

    f140to139 = filters.computeScaling( 'F140W','F139M')
    print("F139M ~ F140W * %.3f"%f140to139)

    # make epoch 00 med-band templates
    imageops.imscaleflux( "colfax.e00/colfax_f125w_e00_reg_drz_sci.fits",
                          f125to127,
                          outfile='colfax.e00/colfax_f127m_e00_reg_drz_sci.fits',
                          clobber=clobber, verbose=verbose)
    imageops.imscaleflux( "colfax.e00/colfax_f160w_e00_reg_drz_sci.fits",
                          f160to153,
                          outfile='colfax.e00/colfax_f153m_e00_reg_drz_sci.fits',
                          clobber=clobber, verbose=verbose)

    imageops.imsum( "colfax.e00/colfax_f125w_e00_reg_drz_sci.fits",
                    "colfax.e00/colfax_f160w_e00_reg_drz_sci.fits",
                    "colfax.e00/colfax_f125w+f160w_e00_reg_drz_sci.fits",
                    clobber=clobber, verbose=verbose )

    imageops.imscaleflux( "colfax.e00/colfax_f125w+f160w_e00_reg_drz_sci.fits",
                          f125and160to140,
                          outfile='colfax.e00/colfax_f140w_e00_reg_drz_sci.fits',
                          clobber=clobber, verbose=verbose)
    imageops.imscaleflux( "colfax.e00/colfax_f140w_e00_reg_drz_sci.fits",
                          f140to139,
                          outfile='colfax.e00/colfax_f139m_e00_reg_drz_sci.fits',
                          clobber=clobber, verbose=verbose)

def mkSubs( clobber=False, verbose=False ):
    """ Make med-band diff images.
    """
    from hstsntools import imageops

    for band in ['f127m','f139m','f153m']:
        template='colfax.e00/colfax_%s_e00_reg_drz_sci.fits'%band
        snimage='colfax.e03/colfax_%s_e03_reg_drz_sci.fits'%band
        diffimage='colfax.e03/colfax_%s_e03-e00_sub_sci.fits'%band

        imageops.imsubtract( template, snimage, diffimage,
                             clobber=clobber, verbose=verbose )

    template='colfax.e00/colfax_f140w_e00_reg_drz_sci.fits'
    for epoch in [3,4,5]:
        snimage='colfax.e%02i/colfax_f140w_e%02i_reg_drz_sci.fits'%(epoch,epoch)
        diffimage='colfax.e%02i/colfax_f140w_e%02i-e00_sub_sci.fits'%(epoch,epoch)
        imageops.imsubtract( template, snimage, diffimage,
                             clobber=clobber, verbose=verbose )
