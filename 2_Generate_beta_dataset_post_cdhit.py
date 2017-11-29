
import os
import shutil
import pandas as pd
from Subroutines import (extract_beta_structure_coords, filter_dssp_database,
    beta_structure_dssp_classification, manipulate_beta_structure,
    manipulate_beta_network)
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

# Determines the absolute file path of the (locally saved) DSSP database
print('Specify absolute file path of DSSP database:')
dssp_database = '{}'.format(input(prompt).strip('/'))

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
cd_hit_domain_df = beta_structure.gen_cd_hit_dict(
    filtered_domain_df=filtered_domain_df
    )

if os.path.isdir('CD_HIT_DSEQS'):
    shutil.rmtree('CD_HIT_DSEQS')
os.mkdir('CD_HIT_DSEQS')

cd_hit_domain_df, pdb_dfs_dict = beta_structure.get_xyz_coords(
    cd_hit_domain_df=cd_hit_domain_df
    )

# Copies required DSSP files from database (on hard drive) to local machine
if os.path.isdir('DSSP_files'):
    shutil.rmtree('DSSP_files')
os.mkdir('DSSP_files')

filtered_files = filter_dssp_database(run=run, resn=resn, rfac=rfac,
                                      dssp_database=dssp_database)
dssp_domain_df = filtered_files.copy_files_from_dssp_database(
    cd_hit_domain_df=cd_hit_domain_df
    )

# Extracts beta-strands (as classified by DSSP) from the beta-structure domains
beta_structure = beta_structure_dssp_classification(
    run=run, resn=resn, rfac=rfac
    )
dssp_residues_dict = beta_structure.extract_dssp_file_lines(
    dssp_domain_df=dssp_domain_df
    )

shutil.rmtree('DSSP_files')

if os.path.isdir('DSSP_filtered_DSEQS'):
    shutil.rmtree('DSSP_filtered_DSEQS')
os.mkdir('DSSP_filtered_DSEQS')

dssp_dfs_dict = beta_structure.get_dssp_sec_struct_df(
    dssp_residues_dict=dssp_residues_dict, pdb_dfs_dict=pdb_dfs_dict
    )
beta_structure.write_dssp_sec_struct_pdb(
    dssp_residues_dict=dssp_residues_dict, dssp_dfs_dict=dssp_dfs_dict
    )

# Combines the beta-strands into sheets and translates the identified
# beta-strand interactions into a network
beta_structure = manipulate_beta_structure()
domain_networks_dict = beta_structure.identify_strand_interactions(
    dssp_residues_dict=dssp_residues_dict, dssp_dfs_dict=dssp_dfs_dict
    )

# Idetifies strand interactions from the networks generated in previous steps
# beta_structure = manipulate_beta_network(domain_networks_dict)
