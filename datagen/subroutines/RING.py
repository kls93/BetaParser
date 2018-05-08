

import os
import pandas as pd
if __name__ == 'subroutines.RING':
    from subroutines.run_stages import run_stages
else:
    from datagen.subroutines.run_stages import run_stages


class RING(run_stages):

    def __init__(self, run_parameters, ring_exe_path):
        run_stages.__init__(self, run_parameters)
        self.ring_exe = ring_exe_path

    def run_RING(self, domain_id):
        # Calculates interaction network of every residue in the domain
        path_to_bio_assembly = '{}/{}/{}.pdb'.format(
            self.pdb_database, domain_id[1:3], domain_id
        )
        if self.resn > 2.0:
            res_approx = 'lollipop'
        elif self.resn <= 2.0:
            res_approx = 'closest'
        ring_run_command = './{} -i {} --all -n {} -g 1'.format(
            self.ring_exe, path_to_bio_assembly, res_approx
        )
        os.system('{} > {}_RING_output.txt'.format(ring_run_command, domain_id))

    def parse_RING_output(self, domain_id, dssp_df):
        # Appends residue interaction network information to dssp_df
        van_der_waals_dict = {}
        hydrogen_bonds_dict = {}
        ionic_bonds_dict = {}
        pi_pi_stacking_dict = {}
        cation_pi_interactions_dict = {}
        for res_id in dssp_df['RES_ID'].tolist():
            van_der_waals_dict[res_id] = []
            hydrogen_bonds_dict[res_id] = []
            ionic_bonds_dict[res_id] = []
            pi_pi_stacking_dict[res_id] = []
            cation_pi_interactions_dict[res_id] = []

        with open('{}_RING_output.txt'.format(domain_id), 'r') as ring_file:
            ring_output = ring_file.readlines()
        os.remove('{}_RING_output.txt'.format(domain_id))

        ring_output = ring_output[1:]
        for line in ring_output:
            aa_1 = line.split()[0].split(':')
            aa_1 = (aa_1[0].replace('_', '') + aa_1[1].replace('_', '')
                    + aa_1[2].replace('_', ''))
            aa_2 = line.split()[2].split(':')
            aa_2 = (aa_2[0].replace('_', '') + aa_2[1].replace('_', '')
                    + aa_2[2].replace('_', ''))
            interaction = line.split()[1]

            van_der_waals_dict[aa_1].append(aa_2)
            hydrogen_bonds_dict[aa_1].append(aa_2)
            ionic_bonds_dict[aa_1].append(aa_2)
            pi_pi_stacking_dict[aa_1].append(aa_2)
            cation_pi_interactions_dict[aa_1].append(aa_2)

            van_der_waals_dict[aa_2].append(aa_1)
            hydrogen_bonds_dict[aa_2].append(aa_1)
            ionic_bonds_dict[aa_2].append(aa_1)
            pi_pi_stacking_dict[aa_2].append(aa_1)
            cation_pi_interactions_dict[aa_2].append(aa_1)

        van_der_waals = ['']*dssp_df.shape[0]
        hydrogen_bonds = ['']*dssp_df.shape[0]
        ionic_bonds = ['']*dssp_df.shape[0]
        pi_pi_stacking = ['']*dssp_df.shape[0]
        cation_pi_interactions = ['']*dssp_df.shape[0]
        for row in range(dssp_df.shape[0]):
            res_id = dssp_df['RES_ID'][row]
            van_der_waals[row] = van_der_waals_dict[res_id]
            hydrogen_bonds[row] = hydrogen_bonds_dict[res_id]
            ionic_bonds[row] = ionic_bonds_dict[res_id]
            pi_pi_stacking[row] = pi_pi_stacking_dict[res_id]
            cation_pi_interactions[row] = cation_pi_interactions_dict[res_id]
        ring_df = pd.DataFrame({'VAN_DER_WAALS': van_der_waals,
                                'HYDROGEN_BONDS': hydrogen_bonds,
                                'IONIC_BONDS': ionic_bonds,
                                'PI_PI': pi_pi_stacking,
                                'CATION_PI': cation_pi_interactions})
        dssp_df = pd.concat([dssp_df, ring_df], axis=1)
        return dssp_df


class calc_residue_interaction_network(run_stages):

    def __init__(self, run_parameters, ring_exe_path):
        run_stages.__init__(self, run_parameters)
        self.ring_exe = ring_exe_path

    def calc_rin(self):
        for domain_id in list(sec_struct_dfs_dict.keys()):
            dssp_df = sec_struct_dfs_dict[domain_id]

            res_int_net = RING()
            res_int_net.run_RING(domain_id)
            dssp_df = res_int_net.parse_RING_output(domain_id, dssp_df)

            sec_struct_dfs_dict[domain_id] = dssp_df
