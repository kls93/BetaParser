
import os
import networkx as nx
import pandas as pd
from collections import OrderedDict
if __name__ == 'subroutines.output_dataframe':
    from subroutines.run_stages import run_stages
    from subroutines.variables import gen_amino_acids_dict
else:
    from datagen.subroutines.run_stages import run_stages
    from datagen.subroutines.variables import gen_amino_acids_dict

class gen_output(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    def identify_edge_central(self, domain_sheets_dict, sec_struct_dfs_dict):
        # Uses domain_networks_dict to identify whether strands are edge or
        # central
        for domain_id in list(sec_struct_dfs_dict.keys()):
            networks = [network for key, network in domain_sheets_dict.items()
                        if domain_id in key]
            G = networks[0]
            for num in range(1, len(networks)):
                H = networks[1]
                G = nx.compose(G, H)
            sec_struct_df = sec_struct_dfs_dict[domain_id]

            edge_or_central = OrderedDict()
            strands = set(sec_struct_df['STRAND_NUM'].tolist())
            for strand_num in strands:
                if strand_num in G:
                    num_of_edges = len(G.neighbors(strand_num))
                    if num_of_edges == 1:
                        edge_id = 'edge'
                    elif num_of_edges > 1:
                        edge_id = 'central'
                    edge_or_central[strand_num] = edge_id

            df_strands = sec_struct_df['STRAND_NUM'].tolist()
            df_edge_or_central = ['']*len(df_strands)
            for index, value in enumerate(df_strands):
                if value in edge_or_central:
                    df_edge_or_central[index] = edge_or_central[value]
            df_edge_or_central = pd.DataFrame({'EDGE_OR_CNTRL': df_edge_or_central})
            sec_struct_df = pd.concat([sec_struct_df, df_edge_or_central], axis=1)
            sec_struct_dfs_dict[domain_id] = sec_struct_df

        return sec_struct_dfs_dict

    def write_beta_strand_dataframe(self, sec_struct_dfs_dict, opm_database,
                                    tilt_angles, strand_numbers, shear_numbers):
        # Generates dataframe of beta-strands
        amino_acids_dict = gen_amino_acids_dict()

        domain_strand_ids = []
        tilt_angle = []
        strand_number = []
        shear_number = []
        res_ids = []
        edge_or_central = []
        residues = []
        int_ext = []
        tm_ext = []
        bckbn_phi_psi = []
        solv_ascblty = []

        for domain_id, sec_struct_df in sec_struct_dfs_dict.items():
            strands_list = set(sec_struct_df['STRAND_NUM'].tolist())
            strands_list = [strand_num for strand_num in strands_list if
                            strand_num != '']

            for strand_num in strands_list:
                strand_df = sec_struct_df[sec_struct_df['STRAND_NUM']==strand_num]
                strand_df = strand_df.reset_index(drop=True)

                # Determines strand orientation
                pdb_code = domain_id[0:4]
                res_ids_list = strand_df['RES_ID'].tolist()
                strand_coordinates = []
                upper_bound = 0
                lower_bound = 0
                reverse = False
                if os.path.isfile('{}/{}.pdb'.format(opm_database, pdb_code)):
                    with open('{}/{}.pdb'.format(opm_database, pdb_code), 'r') as opm_file:
                        for line in opm_file:
                            if line[17:20].strip() != 'DUM':
                                if line[0:6].strip() in ['ATOM', 'HETATM']:
                                    chain_res_num = (line[21:22].strip()
                                                     + line[22:26].strip()
                                                     + line[26:27].strip())
                                    if (line[12:16].strip() == 'CA'
                                        and chain_res_num in res_ids_list
                                        ):
                                        strand_coordinates.append(float(line[46:54]))
                            else:
                                if float(line[46:54]) > upper_bound:
                                    upper_bound = float(line[46:54])
                                elif float(line[46:54]) < lower_bound:
                                    lower_bound = float(line[46:54])

                    tm_ext_list = ['']* len(res_ids_list)
                    if len(strand_coordinates) != len(res_ids_list):
                        print('ERROR: Failed to locate all strand coordinates')
                    else:
                        for index, z_coord in enumerate(strand_coordinates):
                            if z_coord < lower_bound or z_coord > upper_bound:
                                tm_ext_list[index] = 'exterior'
                            else:
                                tm_ext_list[index] = 'transmembrane'

                        if strand_coordinates[0] < strand_coordinates[len(strand_coordinates)-1]:
                            reverse = True

                # Gives strand a unique ID
                domain_strand_ids.append('{}_strand_{}'.format(domain_id, strand_num))

                # Lists tilt angle of strand
                if self.code[0:4] in ['2.40']:
                    tilt_angle.append(tilt_angles[domain_id])

                # Lists number of strands in barrel
                if self.code[0:4] in ['2.40']:
                    if pdb_code in strand_numbers:
                        strand_number.append(strand_numbers[pdb_code])
                    else:
                        strand_number.append('')

                # Lists barrel shear number
                if self.code[0:4] in ['2.40']:
                    shear_number.append('')

                # Lists residue IDs in strand
                if reverse is True:
                    res_ids_list.reverse()
                res_ids.append(res_ids_list)

                # Determines whether the strand is an edge or central strand
                if self.code in ['2.60']:
                    edge_or_central_list = strand_df['EDGE_OR_CNTRL'].tolist()[0]
                    if reverse is True:
                        edge_or_central_list.reverse()
                    edge_or_central.append(edge_or_central_list)

                # Generates FASTA sequence of strand
                res_list = strand_df['RESNAME'].tolist()
                sequence = ''
                for res in res_list:
                    if res in list(amino_acids_dict.keys()):
                        sequence += amino_acids_dict[res]
                    else:
                        sequence += 'X'
                if reverse is True:
                    sequence = sequence[::-1]
                residues.append(sequence)

                # Generates list of interior and exterior facing residues in
                # the strand
                int_ext_list = strand_df['INT_EXT'].tolist()
                if reverse is True:
                    int_ext_list.reverse()
                int_ext.append(int_ext_list)

                # Generates list of transmembrane and exterior residues in the
                # strand
                if self.code[0:4] in ['2.40']:
                    if reverse is True:
                        tm_ext_list.reverse()
                    tm_ext.append(tm_ext_list)

                # Generates list of phi and psi angles along the strand
                phi = strand_df['PHI'].tolist()
                if reverse is True:
                    phi.reverse()
                psi = strand_df['PSI'].tolist()
                if reverse is True:
                    psi.reverse()
                bckbn_geom = [[phi[index], psi[index]] for index, value in
                              enumerate(phi)]
                bckbn_phi_psi.append(bckbn_geom)

                # Generates list of residue solvent accessibility along the
                # strand
                solv_acsblty_list = strand_df['SOLV_ACSBLTY'].tolist()
                if reverse is True:
                    solv_acsblty_list.reverse()
                solv_ascblty.append(solv_acsblty_list)

        if self.code[0:4] in ['2.60']:
            beta_strands_df = pd.DataFrame({'STRAND_ID': domain_strand_ids,
                                            'TILT_ANGLE': tilt_angle,
                                            'RES_ID': res_ids,
                                            'EDGE_OR_CNTRL': edge_or_central,
                                            'RESIDUES': residues,
                                            'INT_EXT': int_ext,
                                            'BCKBN_GEOM': bckbn_phi_psi,
                                            'SOLV_ACSBLTY': solv_ascblty})
            beta_strands_df.to_pickle('Beta_strands_dataframe.pkl')
            beta_strands_df.to_csv('Beta_strands_dataframe.csv')

        elif self.code[0:4] in ['2.40']:
            beta_strands_df = pd.DataFrame({'STRAND_ID': domain_strand_ids,
                                            'TILT_ANGLE': tilt_angle,
                                            'STRAND_NUMBER': strand_number,
                                            'SHEAR_NUMBER': shear_number,
                                            'RES_ID': res_ids,
                                            'TM_OR_EXT': tm_ext,
                                            'RESIDUES': residues,
                                            'INT_EXT': int_ext,
                                            'BCKBN_GEOM': bckbn_phi_psi,
                                            'SOLV_ACSBLTY': solv_ascblty})
            beta_strands_df.to_pickle('Beta_strands_dataframe.pkl')
            beta_strands_df.to_csv('Beta_strands_dataframe.csv')

        else:
            return
