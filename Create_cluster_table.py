output_path = '/project/galaxies/tjuchau/data_files/Kiana_Cluster_Files/full_table.csv'
import os
os.chdir('/project/galaxies') #TJ change working directory to be the parent directory
import glob
from astropy.table import Table
import numpy as np
import astropy.units as u
import re
from matplotlib.path import Path as draw_path
from astropy.coordinates import SkyCoord
from Functions import *
os.chdir('/project/galaxies') #TJ change working directory to be the parent directory
regions = glob.glob('/project/galaxies/tjuchau/data_files/misc_data/*.deg.reg')
m51_clusters = Table.read('/project/galaxies/tjuchau/data_files/misc_data/clusters.csv')

def get_EW_using_filters(feature_filter_file, continuum_filter_files, location, radius):
    #TJ load files

    #TJ get all 3 image fluxes
    Fnu_feature = get_image_flux(feature_filter_file, location, radius, replace_negatives=False)
    Fnu_cont = [get_image_flux(f, location, radius, replace_negatives=False) for f in continuum_filter_files]
    feature_filter = extract_filter_name(feature_filter_file)
    continuum_filters = [extract_filter_name(x) for x in continuum_filter_files]
    if Fnu_feature.unit != Fnu_cont[0].unit:
        print('units are not the same in the feature image and continuum image!')
    elif Fnu_feature.unit == u.W/(u.m**2*u.Hz):
            
        #TJ look up pivot wavelengths
        pivot_feat = jwst_pivots[feature_filter]
        pivot_cont = [jwst_pivots[extract_filter_name(f)] for f in continuum_filter_files]
        
        #TJ convert continuum levels into F_lambda using pivot wavelengths, still need to multiply by dlamda
        fλ_cont = [(Fnu * c / pivot**2).to(u.W / u.m**2 / u.m)
                   for Fnu, pivot in zip(Fnu_cont, pivot_cont)]
        
        #TJ get mean wavelengths
        cont_wls = [jwst_means[f] for f in continuum_filters]
        line_wl = jwst_means[feature_filter]
    
        #TJ interpolate continuum values to the feature wavelength
        feature_continuum = np.interp(
            line_wl.value,
            [w.value for w in cont_wls],
            [f.value for f in fλ_cont]
        ) * u.W / u.m**2 / u.m
        
        #TJ print continuum if needed
        #print("F_lamda of photo continuum : ", feature_continuum)
        #TJ get filter transmission curve info
        wl, T = get_filter_data(feature_filter)
    
        #TJ multiply feature F_lambda by dlambda to complete unit conversion
        norm = np.trapezoid(T, wl) / np.max(T)
        cont_in_filter = feature_continuum * norm
        
        #TJ convert feature filter's F_nu into F_lamda
        fλ_feature = ((Fnu_feature * c / pivot_feat**2).to(u.W / u.m**2 / u.m))*norm 
    
        #TJ feature area is only the area above the continuum
        feature_only = fλ_feature - cont_in_filter
    
        #TJEquivalent width is this area divided by the continuum level
        EW = (feature_only / feature_continuum).to(u.m)
    
        return EW

kiana_files = np.concatenate([glob.glob('/project/galaxies/tjuchau/data_files/Kiana_Cluster_Files/ngc1433*nircam*.csv'),
glob.glob('/project/galaxies/tjuchau/data_files/Kiana_Cluster_Files/ngc1512*nircam*.csv'),
glob.glob('/project/galaxies/tjuchau/data_files/Kiana_Cluster_Files/ngc1672*nircam*.csv')])

names = ['ngc1433', 'ngc1512', 'ngc1672']
table = Table.read(kiana_files[0])
table.add_column([names[0]]*len(table), index=0, name = 'galaxy')
for i, file in enumerate(kiana_files):
    if i ==0:
        continue
    tab = Table.read(file)
    tab.add_column([names[i]]*len(tab), index=0, name = 'galaxy')
    for row in tab:
        table.add_row(row)

paEWs = []
radii = []
for row in table:
    loc = [row['ra'], row['dec']]
    distance = (row['best.universe.luminosity_distance']*u.m).to(u.pc)

    gal_name = row['galaxy']
    if gal_name == 'ngc1433':
        physical_size = 14.7*u.pc
    elif gal_name == 'ngc1512':
        physical_size = 13.9*u.pc
    elif gal_name == 'ngc1672':
        physical_size = 15*u.pc

    radius = np.arcsin(physical_size/distance).to(u.arcsec)
    radii.append(radius)

    files = glob.glob(f'tjuchau/data_files/JWST/images/{gal_name.upper()}/*')
    f150 = [x for x in files if extract_filter_name(x)=='F150W'][0]
    f187 = [x for x in files if extract_filter_name(x)=='F187N'][0]
    f300 = [x for x in files if extract_filter_name(x)=='F300M'][0]
    pa_EW = get_EW_using_filters(f187, [f150,f300], loc, radius)
    paEWs.append(pa_EW)
table.add_column(paEWs, name = 'EW_187')
table.add_column(radii, name = 'radius')
table = table[np.isfinite(table['EW_187'])]

keep_cols = ['col0', 'galaxy', 'best.attenuation.A550', 'best.nebular.logU', 
'best.sfh.age', 'best.stellar.age_m_star', 'best.universe.luminosity_distance',
'best.universe.redshift', 'best.stellar.m_gas', 'best.stellar.m_star', 'ra',
'dec', 'EW_658', 'EW_187', 'radius']
my_table = table[keep_cols]
new_col = [None]*len(my_table)
my_table.add_column(new_col, name = 'spec_data')



def read_ds9_region(filename):
    shapes = []

    with open(filename) as f:
        for line in f:
            line = line.strip()

            if line.startswith("polygon"):
                coords = line[8:-1]
                vals = np.array(coords.split(","), dtype=float)
                ra = vals[0::2]
                dec = vals[1::2]
                shapes.append(("polygon", np.column_stack((ra, dec))))

            if line.startswith("box"):
                vals = line[4:-1].split(",")

                ra = float(vals[0])
                dec = float(vals[1])

                width = float(vals[2].replace('"','')) / 3600
                height = float(vals[3].replace('"','')) / 3600
                angle = float(vals[4])

                shapes.append(("box", (ra, dec, width, height, angle)))

    return shapes


def points_in_box(points, ra_c, dec_c, width, height, angle):
    """
    points: Nx2 array of RA,Dec
    width/height in degrees
    angle in degrees
    """

    # offsets
    dx = (points[:,0] - ra_c) * np.cos(np.radians(dec_c))
    dy = (points[:,1] - dec_c)

    theta = np.radians(angle)

    xr = dx*np.cos(theta) + dy*np.sin(theta)
    yr = -dx*np.sin(theta) + dy*np.cos(theta)

    return (np.abs(xr) <= width/2) & (np.abs(yr) <= height/2)


# Load regions
all_shapes = []
for f in regions:
    all_shapes.extend(read_ds9_region(f))


points = np.column_stack((m51_clusters["ra_gaia"], m51_clusters["dec_gaia"]))
inside_any = np.zeros(len(points), dtype=bool)


for shape in all_shapes:

    if shape[0] == "polygon":
        path = draw_path(shape[1])
        inside_any |= path.contains_points(points)

    if shape[0] == "box":
        ra, dec, w, h, ang = shape[1]
        inside_any |= points_in_box(points, ra, dec, w, h, ang)


region_id = np.full(len(points), -1)

for i, shape in enumerate(all_shapes):

    if shape[0] == "polygon":
        mask = draw_path(shape[1]).contains_points(points)

    if shape[0] == "box":
        ra, dec, w, h, ang = shape[1]
        mask = points_in_box(points, ra, dec, w, h, ang)

    region_id[mask] = i

m51_clusters["region_id"] = region_id
clusters_in_regions = m51_clusters[inside_any]
m51_f150w = '/project/galaxies/tjuchau/data_files/JWST/images/v0p3p2/ngc5194/ngc5194_nircam_lv3_f150w_i2d_anchor.fits'
m51_f187n = '/project/galaxies/tjuchau/data_files/JWST/images/v0p3p2/ngc5194/ngc5194_nircam_lv3_f187n_i2d_anchor.fits'
m51_f300m = '/project/galaxies/tjuchau/data_files/JWST/images/v0p3p2/ngc5194/ngc5194_nircam_lv3_f300m_i2d_anchor.fits'
for i, row in enumerate(clusters_in_regions):
    loc = [row['ra_gaia'], row['dec_gaia']]
    pa_EW = get_EW_using_filters(m51_f187n, [m51_f150w, m51_f300m], loc, 0.3*u.arcsec)
    new_row = [i, "M51", None, None, row['age_best_yr']/1e6, 
    None, None, None, row['mass_best_msun']/2, row['mass_best_msun']/2, 
    row['ra_gaia'], row['dec_gaia'], None, pa_EW, 0.3, regions[row['region_id']]]
    my_table.add_row(new_row)
try:
    my_table.write(output_path, format='csv', overwrite=True)
    print(f"Table successfully saved to {output_path}")
except Exception as e:
    print(f"Error saving table: {e}")