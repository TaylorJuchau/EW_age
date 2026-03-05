
import os
os.chdir('/project/galaxies') #TJ change working directory to be the parent directory
from Functions import *
M51images, filter_files = collect_M51_image_and_filter_files(filter_directory, image_directory)

import os
os.chdir('/cluster/medbow/project/galaxies/tjuchau/') #TJ change working directory to be the parent directory
import re
from astropy.table import Table
import numpy as np
import pickle
M51_pa_alpha_file = [x for x in M51images if extract_filter_name(x) == 'F187N'][0]
M51_cont_files = [x for x in M51images if extract_filter_name(x) in ['F150W', "F300M"]]

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
        
    elif Fnu_feature.unit == u.W/u.m**2:
        print('test')

def get_Pa_a_continuum(continuum_filter_files, location, radius):

    #TJ get the two image fluxes
    app_sum = [get_image_flux(f, location, radius, replace_negatives=False) for f in continuum_filter_files]
    continuum_filters = [extract_filter_name(x) for x in continuum_filter_files]
    if app_sum[0].unit != app_sum[1].unit:
        print('units are not the same in the two images!')
        
    elif app_sum[0].unit == u.W/(u.m**2*u.Hz):
        #TJ look up pivot wavelengths
        pivot_feat = jwst_pivots['F187N']
        pivot_cont = [jwst_pivots[extract_filter_name(f)] for f in continuum_filter_files]
        
        #TJ convert continuum levels into F_lambda using pivot wavelengths, still need to multiply by dlamda
        fλ_cont = [(Fnu * c / pivot**2).to(u.W / u.m**2 / u.m)
                   for Fnu, pivot in zip(app_sum, pivot_cont)]

    elif app_sum[0].unit == u.W/u.m**2:
        fλ_cont = []
    
        for Fband, filt in zip(app_sum, continuum_filters):

            wl, T = get_filter_data(filt)

            width = np.trapezoid(T, wl)  # effective width

            fλ = (Fband / width).to(u.W/u.m**2/u.m)

            fλ_cont.append(fλ)
    #TJ interpolate continuum values to the feature wavelength
    feature_continuum = np.interp(
        pivot_feat.value,
        [w.value for w in pivot_cont],
        [f.value for f in fλ_cont]
    ) * u.W / u.m**2 / u.m
    
    #TJ print continuum if needed
    #print("F_lamda of photo continuum : ", feature_continuum)
    #TJ get filter transmission curve info
    wl, T = get_filter_data("F187N")

    #TJ multiply feature F_lambda by dlambda to complete unit conversion
    norm = np.trapezoid(T, wl) / np.max(T)
    cont_in_filter = feature_continuum * norm
    return cont_in_filter

def get_filter_Weff(filter_name):
    wl, T = get_filter_data(filter_name)
    return np.trapezoid(T,wl) / np.max(T)
        
        
def get_galaxy_ID_from_file(file):
    try:
        gal = file.split('/')[-1].split('_')[0]
        return gal
    except:
        return None

from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import glob

# -------------------------------------------------
# Gather file info per galaxy
# -------------------------------------------------
kiana_files = glob.glob('data_files/Kiana_Cluster_Files/*nircam_all_clusters_results.csv')
radius = 0.3*u.arcsec
data = {}
print("Found Kiana's csv files...")
for file in kiana_files:
    name = get_galaxy_ID_from_file(file)
    data.setdefault(name, {})
    data[name]['files'] = {}
'''
# Danny files (images etc.)
danny_files = glob.glob('/d/crow2/phangs/cycle*/**/*', recursive=True)

for file in danny_files:
    name = get_galaxy_ID_from_file(file)

    if name not in data:
        continue

    if ("f187n_continuum.fits" in file) or ('paa_native_anchor_cont.fits' in file):
        data[name]['files']['F187N_continuum'] = file

    elif ('f187n_line.fits' in file) or ('paa_native_anchor_line.fits' in file):
        data[name]['files']['F187N_line'] = file

    elif extract_filter_name(file) == "F300M":
        data[name]['files']["F300M"] = file

    elif extract_filter_name(file) == "F150W":
        data[name]['files']["F150W"] = file
    elif "f187n_i2d_anchor" in file:
        data[name]['files']["F187N"] = file

'''
data['ngc1433']['files'] = glob.glob('data_files/JWST/images/NGC1433/*')
data['ngc1512']['files'] = glob.glob('data_files/JWST/images/NGC1512/*')
data['ngc1672']['files'] = glob.glob('data_files/JWST/images/NGC1672/*')
galaxy_tables = {}

for file in kiana_files:

    name = get_galaxy_ID_from_file(file)

    # Require needed Pa-alpha files
    if data[name]['files'] == {}:
        continue
    print(f"collecting data for {name} from Kiana's files...")
    f187n_file = [x for x in data[name]['files'] if extract_filter_name(x) == 'F187N'][0]
    f150w_file = [x for x in data[name]['files'] if extract_filter_name(x) == 'F150W'][0]
    f300m_file = [x for x in data[name]['files'] if extract_filter_name(x) == 'F300M'][0]

    t = Table.read(file, format="csv")

    locations = []
    ages = []
    stellar_masses = []
    gas_masses = []
    Avs = []
    EWs = []
    mass = []
    for row in t:
        
        loc = SkyCoord(ra=row['ra']*u.deg,
                       dec=row['dec']*u.deg)

        EW = get_EW_using_filters(f187n_file, [f150w_file, f300m_file], loc, radius)
        
        # --- Append ---
        locations.append(loc)
        ages.append(row['best.sfh.age'])
        stellar_masses.append(row['best.stellar.m_star'])
        gas_masses.append(row['best.stellar.m_gas'])
        Avs.append(row['best.attenuation.A550'])
        EWs.append(EW.to(u.AA))
        mass.append(row['best.stellar.m_star'] + row['best.stellar.m_gas'])
    # --- Create Table ---
    print(f"Data collected, adding table to dictionary")
    tab = Table()
    tab['location'] = locations
    tab['radius'] = [0.3]*len(locations)
    tab['age_Myr'] = ages
    tab['stellar_mass_Msun'] = stellar_masses
    tab['gas_mass_Msun'] = gas_masses
    tab['Av'] = Avs
    tab['PaA_EW'] = EWs
    tab['mass'] = mass
    galaxy_tables[name] = tab
print('Data saved as dictionary of astropy tables to /project/galaxies/tjuchau/data_files/misc_data/saved_data/Kiana_galaxies.pkl')
with open("/project/galaxies/tjuchau/data_files/misc_data/saved_data/Kiana_galaxies.pkl", "wb") as file:
    pickle.dump(galaxy_tables, file)