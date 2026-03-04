output_path = '/project/galaxies/tjuchau/data_files/Kiana_Cluster_Files/full_table.csv'
import os
os.chdir('/project/galaxies') #TJ change working directory to be the parent directory
from glob import glob
from astropy.table import Table
import numpy as np
kiana_files = np.concatenate([glob('/project/galaxies/tjuchau/data_files/Kiana_Cluster_Files/ngc1433*nircam*.csv'),
glob('/project/galaxies/tjuchau/data_files/Kiana_Cluster_Files/ngc1512*nircam*.csv'),
glob('/project/galaxies/tjuchau/data_files/Kiana_Cluster_Files/ngc1672*nircam*.csv')])

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
try:
    table.write(output_path, format='csv', overwrite=True)
    print(f"Table successfully saved to {output_path}")
except Exception as e:
    print(f"Error saving table: {e}")


