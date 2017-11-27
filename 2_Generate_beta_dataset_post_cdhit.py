
import os
import shutil
import pandas as pd
from Subroutines import extract_beta_structure_coords
prompt = '> '

# Determines whether user wants to analyse beta-barrels or beta-sandwiches
print('Specify CATHCODE')
run = input(prompt)

# Determines the absolute file path of the domain directory
print('Specify absolute file path of working directory:')
directory = input(prompt)
os.chdir('/{}/CATH_{}'.format(directory.strip('/'), run))

# Determines the absolute file path of the (locally saved) PDB database
print('Specify absolute file path of PDB database:')
pdb_database = '{}'.format(input(prompt).strip('/'))

# Specifies the resolution and Rfactor (working value) cutoffs used in previous
# steps
print('Select resolution cutoff:')
resn = float(input(prompt))
print('Select Rfactor (working value) cutoff:')
rfac = float(input(prompt))

# Loads the dataframe generated in previous steps
filtered_domain_df = pd.read_pickle(
    'CATH_{}_resn_{}_rfac_{}_pre_cd_hit.pkl'.format(run, resn, rfac)
    )

# Obtains xyz coordinates for the sequences output from CD-HIT
beta_structure = extract_beta_structure_coords(run=run, resn=resn, rfac=rfac,
                                               pdb_database=pdb_database)
cd_hit_domain_df = beta_structure.gen_cd_hit_dict(filtered_domain_df)
if os.path.isdir('CD_HIT_DSEQS_PDB_files'):
    shutil.rmtree('CD_HIT_DSEQS_PDB_files')
os.mkdir('CD_HIT_DSEQS_PDB_files')
beta_structure.get_xyz_coords(cd_hit_domain_df)
