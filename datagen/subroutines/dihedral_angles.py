
import pandas as pd
import isambard
from collections import OrderedDict
if __name__ == 'subroutines.dihedral_angles':
    from subroutines.run_stages import run_stages
else:
    from datagen.subroutines.run_stages import run_stages


class calc_torsion_angles():

    def __init__(self, pdb):
        self.pdb = pdb

    def calc_omega_phi_psi_angles(self, domain_id, dssp_df):
        # Calculates omega, phi and psi dihedral angles for each residue in the
        # input beta_structure
        print('Calculating backbone torsion angles in {}'.format(domain_id))

        # Initialises lists of backbone torsion angle values
        omega = ['']*dssp_df.shape[0]
        phi = ['']*dssp_df.shape[0]
        psi = ['']*dssp_df.shape[0]

        # Generates ordered set of res_ids
        res_ids = dssp_df['RES_ID'].tolist()
        res_ids_filtered = []
        for res_id in res_ids:
            if res_id not in res_ids_filtered:
                res_ids_filtered.append(res_id)
        res_ids = res_ids_filtered

        # Uses ISAMBARD to calculate the omega, phi and psi angles of every
        # residue
        residues = list(self.pdb.get_monomers())
        angles = isambard.ampal.analyse_protein.measure_torsion_angles(residues)
        angles = OrderedDict(zip(res_ids, angles))

        # Updates domain dataframe (dssp_df)
        for row in range(dssp_df.shape[0]):
            if (
                dssp_df['RES_ID'][row] in list(angles.keys())
                and dssp_df['ATMNAME'][row] == 'CA'
            ):
                omega_val = angles[dssp_df['RES_ID'][row]][0]
                if type(omega_val) == float:
                    omega_val = round(omega_val, 1)
                else:
                    omega_val = ''

                phi_val = angles[dssp_df['RES_ID'][row]][1]
                if type(phi_val) == float:
                    phi_val = round(phi_val, 1)
                else:
                    phi_val = ''

                psi_val = angles[dssp_df['RES_ID'][row]][2]
                if type(psi_val) == float:
                    psi_val = round(psi_val, 1)
                else:
                    psi_val = ''

                omega[row] = omega_val
                phi[row] = phi_val
                psi[row] = psi_val

        angle_df = pd.DataFrame(OrderedDict({'PHI': phi,
                                             'PSI': psi,
                                             'OMEGA': omega}))

        dssp_df = pd.concat([dssp_df, angle_df], axis=1)

        return dssp_df

    def calc_chi_angles(self, domain_id, dssp_df):
        # Calculates side chain torsion angles for each residue in the input
        # beta_structure
        print('Calculating side chain torsion angles in {}'.format(domain_id))

        # Initialises list of side chain torsion angle values
        chi = ['']*dssp_df.shape[0]

        # Uses ISAMBARD to calculate the chi angles of every residue
        chi_angles = OrderedDict()
        residues = list(self.pdb.get_monomers())
        chain_res_num = [str(res.ampal_parent.id) + str(res.id) +
                         str(res.insertion_code) for res in residues]
        for index, res in enumerate(residues):
            angles = isambard.ampal.analyse_protein.measure_sidechain_torsion_angles(
                res
            )
            angles_rounded = []
            for angle in angles:
                if type(angle) == float:
                    angles_rounded.append(round(angle, 1))
                else:
                    angles_rounded.append('')
            chi_angles[chain_res_num[index]] = angles_rounded

        # Updates domain dataframe (dssp_df)
        for row in range(dssp_df.shape[0]):
            if (
                dssp_df['RES_ID'][row] in list(chi_angles.keys())
                and dssp_df['ATMNAME'][row] == 'CA'
            ):
                chi[row] = chi_angles[dssp_df['RES_ID'][row]]

        angle_df = pd.DataFrame({'CHI': chi})
        dssp_df = pd.concat([dssp_df, angle_df], axis=1)

        return dssp_df


class backbone_geometry():
    # TODO Maybe calculate rise per residue, residues per turn and radius of
    # helix (for beta-strands)

    def __init__(self, pdb):
        self.pdb = pdb


class dihedral_angles(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    def calc_dihedral_angles(self, sec_struct_dfs_dict):
        # Pipeline function to calculate the omega, phi, psi and chi angles in
        # every residue in the input beta-barrel/sandwich domain

        for domain_id in list(sec_struct_dfs_dict.keys()):
            dssp_df = sec_struct_dfs_dict[domain_id]

            # Creates AMPAL object
            pdb = isambard.ampal.convert_pdb_to_ampal(
                'Beta_strands/{}.pdb'.format(domain_id)
            )

            # Measures backbone torsion angles
            beta_structure = calc_torsion_angles(pdb)
            dssp_df = beta_structure.calc_omega_phi_psi_angles(domain_id, dssp_df)
            dssp_df = beta_structure.calc_chi_angles(domain_id, dssp_df)

            sec_struct_dfs_dict[domain_id] = dssp_df

        return sec_struct_dfs_dict
