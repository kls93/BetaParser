
import os
import shutil
from Subroutines import select_beta_structures, filter_beta_structure
prompt = '> '

# Determines whether user wants to analyse beta-barrels or beta-sandwiches
print('Specify CATHCODE:')
run = input(prompt)

# Determines the absolute file path of the domain descriptions file
print('Specify absolute file path of domains description file:')
directory = input(prompt)
directory = 'Users/ks17361/Lab_work_DW/Beta_structure/Bioinformatics/CATH_database'
os.chdir('{}'.format(directory))

# Determines the absolute file path of the (locally saved) PDB database
print('Specify absolute file path of PDB database:')
pdb_database = '{}'.format(input(prompt).strip('/'))
pdb_database = 'Volumes/Seagate_Backup_Plus_Drive/pdb'

# Generates a list of the domain descriptions provided in
# CATH_domain_description_v_4_2_0.txt. Then filters the domain descriptions
# list for beta-structures (the type dependent upon the earlier user input),
# picking out PDB accession codes and sequences (whose values are stored in the
# 'domain_df' dataframe).
beta_structure = select_beta_structures(run=run)
domains_description = beta_structure.domain_desc_list()
domain_df = beta_structure.domain_desc_filter(domains_description)

# Filters the domain_df for X-ray structures with resolution < 1.6 Angstroms
# (to allow distinction of hydrogen bonds) and R_factor (working value) < 0.20.
# Writes a file listing all PDB ids that meet these criteria suitable for
# uploading to the cd_hit web server.
if os.path.isdir('CATH_{}'.format(run)):
    shutil.rmtree('CATH_{}'.format(run))
os.mkdir('CATH_{}'.format(run))
os.chdir('CATH_{}'.format(run))

print('Select resolution cutoff:')
resn = float(input(prompt))
print('Select Rfactor (working value) cutoff:')
rfac = float(input(prompt))

beta_structure = filter_beta_structure(run=run, resn=resn, rfac=rfac,
                                       domain_df=domain_df,
                                       pdb_database=pdb_database)
filtered_domain_df = beta_structure.resn_rfac_filter()
beta_structure.gen_cd_hit_list(filtered_domain_df)
