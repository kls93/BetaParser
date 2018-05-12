

import os
import pandas as pd
if __name__ == 'subroutines.RING':
    from subroutines.run_stages import run_stages
else:
    from datagen.subroutines.run_stages import run_stages


class calculate_residue_interaction_network(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    def parse_RING_output(self, sec_struct_dfs_dict):
        # Appends residue interaction network information to dssp_df.
        # **TODO** Currently treats main chain / main chain, main chain / side
        # chain and side chain / side chain interactions identically, might
        # want to update this in the future.
        for domain_id in list(sec_struct_dfs_dict.keys()):
            dssp_df = sec_struct_dfs_dict[domain_id]

            interactions = {'VDW': {},
                            'HBOND': {},
                            'IONIC': {},
                            'SSBOND': {},
                            'PIPISTACK': {},
                            'PICATION': {}}

            for res_id in dssp_df['RES_ID'].tolist():
                for interaction_type in list(interactions.keys()):
                    interactions[interaction_type][res_id] = []

            with open('{}{}.ring'.format(self.ring_database, domain_id[0:4]), 'r') as ring_file:
                ring_output = ring_file.readlines()

            ring_output = ''.join(ring_output)
            ring_output = ring_output.split('NodeID')[1]
            ring_output = ring_output.split('\n')[1:]
            for line in ring_output:
                aa_1 = line.split()[0].split(':')
                aa_1 = (aa_1[0].replace('_', '') + aa_1[1].replace('_', '')
                        + aa_1[2].replace('_', ''))
                aa_2 = line.split()[2].split(':')
                aa_2 = (aa_2[0].replace('_', '') + aa_2[1].replace('_', '')
                        + aa_2[2].replace('_', ''))

                interaction_type = line.split()[1].split(':')[0]
                if aa_2 not in interactions[interaction_type][aa_1]:
                    interactions[interaction_type][aa_1].append(aa_2)
                if aa_1 not in interactions[interaction_type][aa_2]:
                    interactions[interaction_type][aa_2].append(aa_1)

            ring_df_dict = {'VDW': ['']*len(dssp_df.shape[0]),
                            'HBOND': ['']*len(dssp_df.shape[0]),
                            'IONIC': ['']*len(dssp_df.shape[0]),
                            'SSBOND': ['']*len(dssp_df.shape[0]),
                            'PIPISTACK': ['']*len(dssp_df.shape[0]),
                            'PICATION': ['']*len(dssp_df.shape[0])}
            for row in range(dssp_df.shape[0]):
                res_id = dssp_df['RES_ID'][row]
                for interaction_type in list(ring_df_dict.keys()):
                    ring_df_dict[interaction_type][row] = interactions[interaction_type][res_id]
            ring_df = pd.DataFrame(ring_df_dict)
            dssp_df = pd.concat([dssp_df, ring_df], axis=1)
            sec_struct_dfs_dict[domain_id] = dssp_df

        return sec_struct_dfs_dict
