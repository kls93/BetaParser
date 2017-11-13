
import os
from Subroutines import gen_beta_structure_df, beta_structure_df
prompt = '> '

# Determines whether user wants to analyse beta-barrels or beta-sandwiches
print('Beta-sandwiches or beta-barrels?')
owChoice = input(prompt).lower()
if owChoice in ['beta-barrels', 'barrels']:
    run = '2.40'
elif owChoice in ['beta-sandwiches', 'sandwiches']:
    run = '2.60'
else:
    for char in owChoice:
        if char.isdigit():
            run = owChoice
            break

# Determines the raw file path of the domain descriptions list
print('Provide directory of domains description file')
owChoice = input(prompt)
os.chdir('{}'.format(owChoice))

# Generates a list of the domain descriptions provided in
# CATH_domain_description_v_4_2_0.txt. Then filters the domain descriptions
# list for beta-structures (the type dependent upon the earlier user input),
# picking out PDB accession codes and sequences (whose values are stored in the
# 'domain_dict' dataframe).
structure = gen_beta_structure_df(run=run)
domains_description = structure.domain_desc_list()
domain_dict = structure.domain_desc_filter(domains_description)

# Filters the domain_dict for X-ray structures with resolution < 1.6 Angstroms
# (to allow distinction of hydrogen bonds) and R_factor (working value) < 0.25.
# Writes a file listing all PDB ids that meet these criteria suitable for
# uploading to the PISCES web server.
print('Select resolution cutoff')
resn = float(input(prompt))
print('Select Rfactor (working value) cutoff')
rfac = float(input(prompt))
structure = beta_structure_df(run=run, resn=resn, rfac=rfac, domain_dict=domain_dict)
filtered_domain_dict = structure.resn_rfac_filter()
structure.gen_PISCES_list(filtered_domain_dict)
