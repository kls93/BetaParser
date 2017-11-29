
import os
import shutil
from Subroutines import select_beta_structures, filter_beta_structure
prompt = '> '

# Determines whether user wants to analyse beta-barrels or beta-sandwiches
print('Specify SCOPe code:')
run = input(prompt)

# Determines the absolute file path of the domain descriptions file
print('Specify absolute file path of domains description file:')
directory = input(prompt)
os.chdir('{}'.format(directory))

# Determines the absolute file path of the (locally saved) PDB database
print('Specify absolute file path of PDB database:')
pdb_database = '{}'.format(input(prompt).strip('/'))
