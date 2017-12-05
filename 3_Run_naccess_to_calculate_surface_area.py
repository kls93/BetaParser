
import pickle
from Subroutines import calculate_solvent_accessible_surface_area
prompt = '> '

# Determines whether user wants to analyse beta-barrels or beta-sandwiches
print('Specify CATHCODE')
run = input(prompt)

# Specifies the resolution and Rfactor (working value) cutoffs used in previous
# steps
print('Select resolution cutoff:')
# resn = float(input(prompt))
resn = 1.6
print('Select Rfactor (working value) cutoff:')
#rfac = float(input(prompt))
rfac = 0.20

# Loads the domain networks dictionaries generated in previous steps
with open(
    'CATH_{}_resn_{}_rfac_{}_domain_networks_dict.pkl'.format(
        run, resn, rfac
        ), 'rb'
    ) as pickle_file:
    (dssp_dfs_dict, domain_networks_dict, domain_sheets_dict) = pickle.load(pickle_file)

beta_structure = calculate_solvent_accessible_surface_area(
    run=run, resn=resn, rfac=rfac, dssp_dfs_dict=dssp_dfs_dict,
    domain_networks_dict=domain_networks_dict,
    domain_sheets_dict=domain_sheets_dict
    )
beta_structure.run_naccess()
