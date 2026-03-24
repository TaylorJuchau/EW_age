import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from matplotlib.patches import Circle
from astropy.table import Table

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
import astropy.units as u
from matplotlib.patches import Circle
from astropy.table import Table

def interactive_aperture_editor(table, i, image_files, zoom_factor=5):
    
    row = table[i]
    ra, dec = row['ra'], row['dec']
    radius_arcsec = row['radius']*u.arcsec

    coord = SkyCoord(ra, dec, unit='deg')

    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    axes = axes.flatten()

    circles = []
    
    for ax, file in zip(axes, image_files):
        
        hdu = fits.open(file)['SCI']
        data = hdu.data
        wcs = WCS(hdu.header)

        # pixel position
        x, y = wcs.world_to_pixel(coord)

        # pixel scale (arcsec/pixel)
        pixscale = np.mean(np.abs(wcs.pixel_scale_matrix.diagonal())) * 3600

        # cutout size in pixels
        size_pix = int((zoom_factor * radius_arcsec.value) / pixscale)

        # make cutout
        cutout = Cutout2D(data, (x, y), (size_pix, size_pix), wcs=wcs)

        cut_data = cutout.data
        cut_wcs = cutout.wcs

        # center of cutout
        cy, cx = np.array(cut_data.shape) / 2

        ax.imshow(
            cut_data,
            origin='lower',
            cmap='gray',
            vmin=np.nanpercentile(cut_data, 5),
            vmax=np.nanpercentile(cut_data, 99)
        )

        ax.set_title(file.split('/')[-1], fontsize=8)

        # initial radius in pixels
        r_pix = radius_arcsec.value / pixscale

        circ = Circle((cx, cy), r_pix, edgecolor='red', facecolor='none', lw=2)
        ax.add_patch(circ)

        circles.append((circ, pixscale))

    dragging = {'active': False}

    def on_press(event):
        if event.inaxes is None:
            return
        dragging['active'] = True

    def on_release(event):
        dragging['active'] = False

    def on_motion(event):
        if not dragging['active'] or event.inaxes is None:
            return
        
        ax = event.inaxes
        idx = list(axes).index(ax)

        circ, pixscale = circles[idx]

        x0, y0 = circ.center
        dx = event.xdata - x0
        dy = event.ydata - y0

        new_r = np.sqrt(dx**2 + dy**2)

        # update ALL circles
        for c, _ in circles:
            c.set_radius(new_r)

        fig.canvas.draw_idle()

    def on_key(event):
        if event.key == 'enter':
            r_pix = circles[0][0].get_radius()
            pixscale = circles[0][1]
            new_radius = r_pix * pixscale

            table[i]['radius'] = new_radius
            print(f"Updated row {i} radius → {new_radius:.3f} arcsec")
            plt.close(fig)

    fig.canvas.mpl_connect('button_press_event', on_press)
    fig.canvas.mpl_connect('button_release_event', on_release)
    fig.canvas.mpl_connect('motion_notify_event', on_motion)
    fig.canvas.mpl_connect('key_press_event', on_key)

    plt.tight_layout()
    plt.show()

    
hst_file = '/project/galaxies/tjuchau/data_files/HST/ngc5194/F555W_NGC5194_ACS_WFC_drc.fits'
f187_file = '/project/galaxies/tjuchau/data_files/JWST/images/v0p3p2/ngc5194/ngc5194_nircam_lv3_f187n_i2d_anchor.fits'
f300_file = '/project/galaxies/tjuchau/data_files/JWST/images/v0p3p2/ngc5194/ngc5194_nircam_lv3_f300m_i2d_anchor.fits'
cont_sub_file = '/project/galaxies/tjuchau/data_files/JWST/images/v0p3p2/ngc5194/ngc5194_nircam_lv3_f187n_i2d_anchor_cont_subtracted.fits'
table = Table.read('/project/galaxies/tjuchau/data_files/Kiana_Cluster_Files/full_table.csv')
table = table[table['galaxy'] == 'M51']
for i in range(len(table)):
    interactive_aperture_editor(table, i, image_files=[
    hst_file,
    f187_file,
    f300_file,
    cont_sub_file
])
