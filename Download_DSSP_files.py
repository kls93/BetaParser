
import os
import pandas as pd
from Subroutines import filter_dssp_database
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

# Specifies the resolution and R_factor (working value) cutoffs used in
# previous steps
print('Select resolution cutoff')
resn = float(input(prompt))
print('Select Rfactor (working value) cutoff')
rfac = float(input(prompt))

# Determines the absolte file path of the domain directory
print('Specify absolute file path of working directory')
directory = input(prompt)
os.chdir('{}'.format(directory))

# Copies required DSSP files from database (on hard drive) to local machine
cd_hit_domain_dict_xyz = pd.read_pickle(
    'CATH_{}_resn_{}_rfac_{}_filtered_xyz.pkl'.format(run, resn, rfac)
    )
pdb_list = cd_hit_domain_dict_xyz['PDB_CODE'].tolist()

filtered_files = filter_dssp_database(run=run, resn=resn, rfac=rfac)
filtered_files.copy_files_from_dssp_database(pdb_list)
