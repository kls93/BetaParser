
import os
import pandas as pd
from Subroutines import beta_structure_dssp_classification
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

# Specifies the resolution and Rfactor (working value) cutoffs used in previous
# steps
print('Select resolution cutoff')
resn = float(input(prompt))
print('Select Rfactor (working value) cutoff')
rfac = float(input(prompt))

os.mkdir('DSSP_PDB_files_{}'.format(run))
os.chdir('CD_HIT_PDB_files_{}'.format(run))

for pdb in os.listdir('.'):
    with open('{}'.format(pdb), 'r') as pdb_file:
        pdb_file_lines = [line.strip('\n') for line in pdb_file]
    pdb_code = pdb.split('_')[0]
    beta_structure = beta_structure_dssp_classification(run=run, resn=resn,
                                                        rfac=rfac, pdb_name=pdb,
                                                        pdb_code=pdb)
    beta_structure.extract_dssp_secondary_structure(pdb_file_lines)
