
import os
import shutil
import pandas as pd
from Subroutines import extract_beta_structure_coords
prompt = '> '

# Determines whether user wants to analyse beta-barrels or beta-sandwiches
print('Beta-sandwiches or beta-barrels?')
structure_type = input(prompt).lower()
if structure_type in ['beta-barrels', 'barrels']:
    run = '2.40'
elif structure_type in ['beta-sandwiches', 'sandwiches']:
    run = '2.60'
else:
    for char in structure_type:
        if char.isdigit():
            run = structure_type
            break

# Determines the absolute file path of the domain directory
print('Specify absolute file path of working directory')
directory = input(prompt)
os.chdir('{}'.format(directory))

# Specifies the resolution and Rfactor (working value) cutoffs used in previous
# steps
print('Select resolution cutoff')
resn = float(input(prompt))
print('Select Rfactor (working value) cutoff')
rfac = float(input(prompt))

# Loads the dataframe generated in previous steps
filtered_domain_dict = pd.read_pickle(
    'CATH_{}_resn_{}_rfac_{}_pre_cd_hit.pkl'.format(run, resn, rfac)
    )

# Obtains xyz coordinates for the sequences output from CD-HIT
beta_structure = exttract_beta_structure_coords(run=run, resn=resn, rfac=rfac)
cd_hit_domain_dict = beta_structure.gen_cd_hit_dict(filtered_domain_dict)
if os.path.isdir('PDB_files'):
    shutil.rmtree('PDB_files')
os.mkdir('PDB_files')
beta_structure.get_xyz_coords(cd_hit_domain_dict)
