
import isambard
import pandas as pd
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
        # within a user-specified radius (TODO: optimise this value)
        for domain_id in list(sec_struct_dfs_dict.keys()):
            print('Calculating neighbouring residues for each atom in '
                  '{}'.format(domain_id))

            dssp_df = sec_struct_dfs_dict[domain_id]

            neighbours_dict = OrderedDict()
            domain = isambard.ampal.convert_pdb_to_ampal(
                'Beta_strands/{}.pdb'.format(domain_id)
            )

            atom_list = list(domain.get_atoms())
            atom_ca_list = [atom for atom in atom_list if atom.res_label == 'CA']
            for atom_ca in atom_ca_list:
                neighbouring_res = []
                for atom in domain.is_within(self.radius, atom_ca, ligands=False):
                    parent_res = (atom.ampal_parent.ampal_parent.id
                                  + atom.ampal_parent.id
                                  + atom.ampal_parent.insertion_code)
                    neighbouring_res.append(parent_res)

                neighbours_set = []
                for res in neighbouring_res:
                    if not res in neighbours_set:
                        neighbours_set.append(res)

                parent_res = (atom_ca.ampal_parent.ampal_parent.id
                              + atom_ca.ampal_parent.id
                              + atom_ca.ampal_parent.insertion_code)
                neighbours_set.remove(parent_res)
                neighbours_dict[parent_res] = neighbours_set

            neighbours_list = ['']*dssp_df.shape[0]
            for row in range(dssp_df.shape[0]):
                if (
                    dssp_df['ATMNAME'][row] == 'CA'
                        and dssp_df['RES_ID'][row] in list(neighbours_dict.keys())
                ):
                    neighbours_list[row] = neighbours_dict[dssp_df['RES_ID'][row]]
            neighbours_df = pd.DataFrame({'NEIGHBOURS': neighbours_list})
            dssp_df = pd.concat([dssp_df, neighbours_df], axis=1)
            sec_struct_dfs_dict[domain_id] = dssp_df

        return sec_struct_dfs_dict
