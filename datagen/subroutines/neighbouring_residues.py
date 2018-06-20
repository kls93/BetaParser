
import isambard
import pandas as pd
import numpy as np
from collections import OrderedDict
if __name__ == 'subroutines.neighbouring_residues':
    from subroutines.run_stages import run_stages
else:
    from datagen.subroutines.run_stages import run_stages


class nearest_neighbours(run_stages):

    def __init__(self, run_parameters, radius):
        run_stages.__init__(self, run_parameters)
        self.radius = radius

    def calculate_nearest_neighbours(self, sec_struct_dfs_dict):
        # For each C_alpha atom in the domain, makes a list of all residues
        # within a user-specified radius in the context of the biological
        # assembly (TODO: optimise this value)
        for domain_id in list(sec_struct_dfs_dict.keys()):
            print('Calculating neighbouring residues for each atom in '
                  '{}'.format(domain_id))

            dssp_df = sec_struct_dfs_dict[domain_id]
            neighbours_dict = OrderedDict()

            # Creates AMPAL object
            domain = isambard.ampal.convert_pdb_to_ampal(
                'Parent_assemblies/{}.pdb'.format(domain_id[0:4])
            )

            # Creates list of res_ids ('neighbouring residues') within a
            # user-specified radius of the centre of mass of every residue in
            # the input domain
            res_ids = dssp_df[dssp_df['ATMNAME'] == 'CA']['RES_ID'].tolist()
            for res_id in res_ids:
                res_df = dssp_df[(dssp_df['RES_ID'] == res_id) &
                                 (~dssp_df['ATMNAME'].isin(['N', 'C', 'O']))]
                res_df = res_df.reset_index(drop=True)

                try:
                    atom_df = res_df[res_df['ATMNAME'] == 'CA']
                    atom_df = atom_df.reset_index(drop=True)
                except:
                    pass

                if atom_df.shape[0] == 1:
                    coordinates = np.zeros((res_df.shape[0], 3))
                    for row in range(res_df.shape[0]):
                        coordinates[row][0] = res_df['XPOS'][row]
                        coordinates[row][1] = res_df['YPOS'][row]
                        coordinates[row][2] = res_df['ZPOS'][row]
                    com = np.average(coordinates, axis=0)
                    dummy_atom = isambard.ampal.Atom([com[0], com[1], com[2]], 'C')

                    neighbouring_res = []
                    for neighbour in domain.is_within(self.radius, dummy_atom, ligands=False):
                        parent_res = (neighbour.ampal_parent.ampal_parent.id
                                      + neighbour.ampal_parent.id
                                      + neighbour.ampal_parent.insertion_code)
                        neighbouring_res.append(parent_res)

                    # Creates ordered set of neighbouring residues
                    neighbours_set = []
                    for res in neighbouring_res:
                        if not res in neighbours_set:
                            neighbours_set.append(res)

                    # Removes central residue from the list of its neighbours,
                    # then stores neighbours list in dictionary
                    neighbours_set.remove(res_id)
                    neighbours_dict[res_id] = neighbours_set

            # Updates dataframe with neighbouring residue information
            neighbours_list = ['']*dssp_df.shape[0]
            for row in range(dssp_df.shape[0]):
                if (
                    dssp_df['RES_ID'][row] in list(neighbours_dict.keys())
                    and dssp_df['ATMNAME'][row] == 'CA'
                ):
                    neighbours_list[row] = neighbours_dict[dssp_df['RES_ID'][row]]
            neighbours_df = pd.DataFrame({'NEIGHBOURS': neighbours_list})
            dssp_df = pd.concat([dssp_df, neighbours_df], axis=1)
            sec_struct_dfs_dict[domain_id] = dssp_df

        return sec_struct_dfs_dict
