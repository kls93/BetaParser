
import os
import shutil
import pandas as pd
from Subroutines import beta_structure_coords
prompt = '> '

# Determines whether user wants to analyse beta-barrels or beta-sandwiches
print('Beta-sandwiches or beta-barrels?')
cathcode = input(prompt).lower()
if cathcode in ['beta-barrels', 'barrels']:
    run = '2.40'
elif cathcode in ['beta-sandwiches', 'sandwiches']:
    run = '2.60'
else:
    for char in cathcode:
        if char.isdigit():
            run = cathcode
            break

# Determines the raw file path of the domain descriptions list
print('Provide directory of domains description file')
directory = input(prompt)
os.chdir('{}/CATH_{}'.format(directory, cathcode))

# Specifies the resolution and Rfactor (working value) cutoffs used in previous
# steps
print('Select resolution cutoff')
resn = float(input(prompt))
print('Select Rfactor (working value) cutoff')
rfac = float(input(prompt))

# Loads the dataframe generated in previous steps
filtered_domain_dict = pd.read_pickle('CATH_{}_resn_{}_rfac_{}_pre_cd_hit.pkl'.format(run, resn, rfac))
filtered_domain_dict = filtered_domain_dict.drop('index', axis=1)
filtered_domain_dict = filtered_domain_dict.reset_index(drop=True)

# Obtains xyz coordinates for the sequences output from CD-HIT
beta_structure = beta_structure_coords(run=run, resn=resn, rfac=rfac)
cd_hit_domain_dict = beta_structure.gen_cd_hit_dict(filtered_domain_dict)

if os.path.isdir('PDB_files'):
    shutil.rmtree('PDB_files')
os.mkdir('PDB_files')

beta_structure.get_xyz_coords(cd_hit_domain_dict)
