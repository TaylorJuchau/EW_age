import os
from Functions import *
os.chdir('/project/galaxies') #TJ change working directory to be the parent directory
from astropy.table import Table
import glob
from matplotlib.path import Path
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
#TJ To create table of galaxies, run tjuchau/projects/EW_vs_Age/Create_cluster_table.py
table = Table.read('/project/galaxies/tjuchau/data_files/Kiana_Cluster_Files/full_table.csv')

for row in table:
    loc = [row['ra'], row['dec']]
    radius = row['radius']
    galaxy_name = row['galaxy']
    if galaxy_name == "M51":
        hst_file = '/project/galaxies/tjuchau/data_files/HST/ngc5194/F555W_NGC5194_ACS_WFC_drc.fits'
        f150_file = '/project/galaxies/tjuchau/data_files/JWST/images/v0p3p2/ngc5194/ngc5194_nircam_lv3_f150w_i2d_anchor.fits'
        f187_file = '/project/galaxies/tjuchau/data_files/JWST/images/v0p3p2/ngc5194/ngc5194_nircam_lv3_f187n_i2d_anchor.fits'
        f300_file = '/project/galaxies/tjuchau/data_files/JWST/images/v0p3p2/ngc5194/ngc5194_nircam_lv3_f300m_i2d_anchor.fits'
        cont_sub_file = '/project/galaxies/tjuchau/data_files/JWST/images/v0p3p2/ngc5194/ngc5194_nircam_lv3_f187n_i2d_anchor_cont_subtracted.fits'
    elif galaxy_name == 'ngc1433':
        hst_file = '/project/galaxies/tjuchau/data_files/HST/ngc1433/hlsp_phangs-hst_hst_wfc3-uvis_ngc1433_f555w_v1_exp-drc-sci.fits'
        f150_file = '/project/galaxies/tjuchau/data_files/JWST/images/NGC1433/ngc1433_nircam_lv3_f150w_i2d_anchor.fits'
        f187_file = '/project/galaxies/tjuchau/data_files/JWST/images/NGC1433/ngc1433_nircam_lv3_f187n_i2d_anchor.fits'
        f300_file = '/project/galaxies/tjuchau/data_files/JWST/images/NGC1433/ngc1433_nircam_lv3_f300m_i2d_anchor.fits'
        cont_sub_file = '/project/galaxies/tjuchau/data_files/JWST/images/NGC1433/ngc1433_nircam_lv3_f187n_i2d_anchor_cont_subtracted.fits'
    elif galaxy_name == 'ngc1512':
        hst_file = '/project/galaxies/tjuchau/data_files/HST/ngc1512/hlsp_phangs-hst_hst_wfc3-uvis_ngc1512mosaic_f555w_v1_exp-drc-sci.fits'
        f150_file = '/project/galaxies/tjuchau/data_files/JWST/images/NGC1512/ngc1512_nircam_lv3_f150w_i2d_anchor.fits'
        f187_file = '/project/galaxies/tjuchau/data_files/JWST/images/NGC1512/ngc1512_nircam_lv3_f187n_i2d_anchor.fits'
        f300_file = '/project/galaxies/tjuchau/data_files/JWST/images/NGC1512/ngc1512_nircam_lv3_f300m_i2d_anchor.fits'
        cont_sub_file = '/project/galaxies/tjuchau/data_files/JWST/images/NGC1512/ngc1512_nircam_lv3_f187n_i2d_anchor_cont_subtracted.fits'
    elif galaxy_name == 'ngc1672':
        hst_file = '/project/galaxies/tjuchau/data_files/HST/ngc1672/hlsp_phangs-hst_hst_wfc3-uvis_ngc1672mosaic_f555w_v1_exp-drc-sci.fits'
        f150_file = '/project/galaxies/tjuchau/data_files/JWST/images/NGC1672/ngc1672_nircam_lv3_f150w_i2d_anchor.fits'
        f187_file = '/project/galaxies/tjuchau/data_files/JWST/images/NGC1672/ngc1672_nircam_lv3_f187n_i2d_anchor.fits'
        f300_file = '/project/galaxies/tjuchau/data_files/JWST/images/NGC1672/ngc1672_nircam_lv3_f300m_i2d_anchor.fits'
        cont_sub_file = '/project/galaxies/tjuchau/data_files/JWST/images/NGC1672/ngc1672_nircam_lv3_f187n_i2d_anchor_cont_subtracted.fits'
    ew = get_EW_using_filters(f187_file, [f150_file, f300_file], loc, radius)
    row['EW_187'] = ew
try:
    table.write(output_path, format='csv', overwrite=True)
    print(f"Table successfully saved to {output_path}")
except Exception as e:
    print(f"Error saving table: {e}")

    