import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
import os

from Functions import *
os.chdir('/project/galaxies') 



def continuum_subtract(feature_file, cont1_file, cont2_file):
    
    # ---------------------------
    # Load feature image
    # ---------------------------
    with fits.open(feature_file) as hdul:
        feature_data = hdul['SCI'].data.astype(float)
        feature_header = hdul['SCI'].header
        feature_wcs = WCS(feature_header)

    shape_out = feature_data.shape

    # ---------------------------
    # Load continuum images
    # ---------------------------
    with fits.open(cont1_file) as hdul:
        cont1_data = hdul['SCI'].data.astype(float)
        cont1_wcs = WCS(hdul['SCI'].header)

    with fits.open(cont2_file) as hdul:
        cont2_data = hdul['SCI'].data.astype(float)
        cont2_wcs = WCS(hdul['SCI'].header)

    # ---------------------------
    # Reproject continuum images
    # ---------------------------
    cont1_regrid, _ = reproject_interp(
        (cont1_data, cont1_wcs),
        feature_wcs,
        shape_out=shape_out
    )

    cont2_regrid, _ = reproject_interp(
        (cont2_data, cont2_wcs),
        feature_wcs,
        shape_out=shape_out
    )

    # ---------------------------
    # Get pivot wavelengths
    # ---------------------------
    lam_feature = jwst_pivots[extract_filter_name(feature_file)].value
    lam_cont1 = jwst_pivots[extract_filter_name(cont1_file)].value
    lam_cont2 = jwst_pivots[extract_filter_name(cont2_file)].value

    # Ensure continuum wavelengths are ordered
    if lam_cont1 > lam_cont2:
        lam_cont1, lam_cont2 = lam_cont2, lam_cont1
        cont1_regrid, cont2_regrid = cont2_regrid, cont1_regrid

    # ---------------------------
    # Estimate continuum at feature wavelength
    # ---------------------------
    continuum = cont1_regrid + (
        (cont2_regrid - cont1_regrid)
        * (lam_feature - lam_cont1)
        / (lam_cont2 - lam_cont1)
    )

    # ---------------------------
    # Subtract continuum
    # ---------------------------
    line_image = feature_data - continuum

    # Preserve NaN mask
    line_image[np.isnan(feature_data)] = np.nan

    # ---------------------------
    # Save output FITS
    # ---------------------------
    base = feature_file.split('.fi')[-2]
    output_file = base + "_cont_subtracted.fits"

    fits.writeto(
        output_file,
        line_image,
        feature_header,
        overwrite=True
    )

    print(f"Saved continuum-subtracted image: {output_file}")

names = ['ngc1433', 'ngc1512', 'ngc1672']
for name in names:
    files = glob.glob(f'tjuchau/data_files/JWST/images/{name.upper()}/*')
    continuum_subtract(files[1], files[0], files[2])
m51_f150w = '/project/galaxies/tjuchau/data_files/JWST/images/v0p3p2/ngc5194/ngc5194_nircam_lv3_f150w_i2d_anchor.fits'
m51_f187n = '/project/galaxies/tjuchau/data_files/JWST/images/v0p3p2/ngc5194/ngc5194_nircam_lv3_f187n_i2d_anchor.fits'
m51_f300m = '/project/galaxies/tjuchau/data_files/JWST/images/v0p3p2/ngc5194/ngc5194_nircam_lv3_f300m_i2d_anchor.fits'
continuum_subtract(m51_f187n, m51_f150w, m51_f300m)