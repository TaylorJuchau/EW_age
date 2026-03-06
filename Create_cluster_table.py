output_path = '/project/galaxies/tjuchau/data_files/Kiana_Cluster_Files/full_table.csv'
import os
os.chdir('/project/galaxies') #TJ change working directory to be the parent directory
import glob
from astropy.table import Table
import numpy as np
import astropy.units as u
import re
from Functions import *
os.chdir('/project/galaxies') #TJ change working directory to be the parent directory

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
try:
    table.write(output_path, format='csv', overwrite=True)
    print(f"Table successfully saved to {output_path}")
except Exception as e:
    print(f"Error saving table: {e}")
