
# Code for running on Tombstone

import os
import shutil
import pandas as pd
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

cd_hit_domain_dict_xyz = pd.read_pickle('CATH_{}_resn_{}_rfac_{}_filtered_xyz.pkl'.format(run, resn, rfac))
pdb_list = cd_hit_domain_dict_xyz['PDB_CODE'].tolist()
unprocessed_list = []

os.makedirs('CATH_{}/DSSP_files'.format(run))
os.chdir('../shared/structural_bioinformatics/data/')
for pdb in pdb_list:
    middle_characters = pdb[1:len(pdb)-1]
    try:
        shutil.copy2(
            '{}/{}/dssp/{}_1.mmol.dssp'.format(middle_characters, pdb, pdb),
            '../../../ks17361/CATH_{}/DSSP_files/{}.dssp'.format(run, pdb))
    except FileNotFoundError:
        unprocessed_list.append(pdb)

os.chdir('../../../ks17361/CATH_{}/DSSP_files'.format(run))
with open('Unprocessed_pdb_files.txt', 'w') as unprocessed_file:
    for pdb in unprocessed_list:
        unprocessed_file.write('{}\n'.format(pdb))
