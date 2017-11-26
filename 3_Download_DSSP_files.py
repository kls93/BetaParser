
import os
import shutil
import pandas as pd
from Subroutines import filter_dssp_database, beta_structure_dssp_classification
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

# Determines the absolte file path of the domain directory
print('Specify absolute file path of working directory')
directory = input(prompt)
os.chdir('{}'.format(directory))

# Determines the absolute file path of the (locally saved) DSSP database
print('Specify absolute file path of DSSP database')
dssp_database = input(prompt)

# Specifies the resolution and R_factor (working value) cutoffs used in
# previous steps
print('Select resolution cutoff')
resn = float(input(prompt))
print('Select Rfactor (working value) cutoff')
rfac = float(input(prompt))

# Makes a list of PDB accession codes from the CD-HIT filtered dataframe
cd_hit_domain_df_xyz = pd.read_pickle(
    'CATH_{}_resn_{}_rfac_{}_filtered_xyz.pkl'.format(run, resn, rfac)
    )
pdb_list = cd_hit_domain_df_xyz.PDB_CODE.tolist()

# Copies required DSSP files from database (on hard drive) to local machine
filtered_files = filter_dssp_database(run=run, resn=resn, rfac=rfac,
                                      dssp_database=dssp_database)
dssp_domain_df = filtered_files.copy_files_from_dssp_database(pdb_list)

# Extracts beta-strands (as classified by DSSP) from the beta-structure domains
beta_structure = beta_structure_dssp_classification(
    run=run, resn=resn, rfac=rfac
    )
dssp_domain_df, dssp_residues_dict = beta_structure.extract_dssp_file_lines(
    dssp_domain_df=dssp_domain_df
    )
if os.path.isdir('DSSP_filtered_DSEQS'):
    shutil.rmtree('DSSP_filtered_DSEQS')
os.mkdir('DSSP_filtered_DSEQS')
beta_structure.get_dssp_sec_struct(dssp_residues_dict=dssp_residues_dict)
beta_structure.write_dssp_sec_struct_pdb(dssp_residues_dict=dssp_residues_dict)
