
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


class output_calcs():

    def determine_strand_orientation(domain_id, strand_id, strand_df,
                                     opm_database, unprocessed_list):
        # Determines orientation of input strand in barrel. If the N-terminus
        # is in the periplasm (Z-coordinate = negative) and C-terminus is
        # extracellular (Z-coordinate = positive), the strand orientation is
        # reversed in the output dataframe.
        print('Determining strand orientation in {} strand '
              '{}'.format(domain_id, strand_id))

        pdb_code = domain_id[0:4]
        res_ids_list = strand_df['RES_ID'].tolist()
        strand_coordinates = OrderedDict()
        lower_bound = 0
        upper_bound = 0
        reverse = False

        if os.path.isfile('{}/{}.pdb'.format(opm_database, pdb_code)):
            with open('{}/{}.pdb'.format(opm_database, pdb_code), 'r') as opm_file:
                for line in opm_file:
                    if line[17:20].strip() != 'DUM':  # Don't put as 'and'
                        # statement with line below (see corresponding 'else' statement)
                        if line[0:6].strip() in ['ATOM', 'HETATM']:
                            chain_res_num = line[21:27].replace(' ', '')
                            if (line[12:16].strip() == 'CA'
                                    and chain_res_num in res_ids_list
                                    ):
                                strand_coordinates[chain_res_num] = float(line[46:54])
                    else:
                        if float(line[46:54]) > upper_bound:
                            upper_bound = float(line[46:54])
                        elif float(line[46:54]) < lower_bound:
                            lower_bound = float(line[46:54])

            z_coordinates = list(strand_coordinates.values())
            if len(z_coordinates) != len(res_ids_list):
                unprocessed_list.append('{}_strand_{}'.format(domain_id, strand_id))
            else:
                if z_coordinates[0] < z_coordinates[-1]:
                    reverse = True

        return (reverse, res_ids_list, strand_coordinates, lower_bound,
                upper_bound, unprocessed_list)

    def determine_tm_or_ext(domain_id, strand_id, res_ids_list,
                            strand_coordinates, lower_bound, upper_bound,
                            unprocessed_list):
        # If the parent structure is in the OPM database, labels residues in an
        # input strand as either 'transmembrane' or 'external'
        print('Calculating positions of strands in membrane in '
              '{}'.format(domain_id))

        tm_ext_list = [''] * len(res_ids_list)

        if '{}_strand_{}'.format(domain_id, strand_id) not in unprocessed_list:
            for index, res_id in enumerate(res_ids_list):
                z_coord = strand_coordinates[res_id]
                if z_coord < lower_bound or z_coord > upper_bound:
                    tm_ext_list[index] = 'external'
                else:
                    tm_ext_list[index] = 'transmembrane'

        return tm_ext_list


class gen_output(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    def identify_edge_central(self, domain_sheets_dict, sec_struct_dfs_dict):
        # Uses domain_networks_dict to identify whether strands are edge or
        # central
        for domain_id in list(sec_struct_dfs_dict.keys()):
            print('Identifying edge and central strands in {}'.format(domain_id))

            networks = [network for key, network in domain_sheets_dict.items()
                        if domain_id in key]

            G = networks[0]
            for num in range(1, len(networks)):
                H = networks[num]
                G = nx.compose(G, H)

            dssp_df = sec_struct_dfs_dict[domain_id]

            edge_or_central_dict = OrderedDict()
            df_strands = [strand for strand in set(dssp_df['STRAND_NUM'].tolist())
                          if strand != '']

            for strand_num in df_strands:
                num_of_edges = len(G.neighbors(strand_num))
                if num_of_edges == 1:
                    edge_id = 'edge'
                elif num_of_edges > 1:
                    edge_id = 'central'
                edge_or_central_dict[strand_num] = edge_id

            edge_or_central_list = ['']*dssp_df.shape[0]
            for row in range(dssp_df.shape[0]):
                if (dssp_df['STRAND_NUM'][row] in list(edge_or_central_dict.keys())
                        and dssp_df['ATMNAME'][row] == 'CA'
                        ):
                    edge_or_central_list[row] = edge_or_central_dict[dssp_df['STRAND_NUM'][row]]
            df_edge_or_central = pd.DataFrame({'EDGE_OR_CNTRL': edge_or_central_list})
            dssp_df = pd.concat([dssp_df, df_edge_or_central], axis=1)
            sec_struct_dfs_dict[domain_id] = dssp_df

        return sec_struct_dfs_dict

    def write_beta_strand_dataframe(self, strand_or_res, sec_struct_dfs_dict,
                                    opm_database, tilt_angles, strand_numbers,
                                    shear_numbers):
        # Generates dataframes of residues and strands in the retained domains.
        if __name__ == 'subroutines.output_dataframe':
            from subroutines.output_dataframe import output_calcs
        else:
            from datagen.subroutines.output_dataframe import output_calcs

        # Generates dictionary of amino acid 1 and 3 letter codes
        amino_acids_dict = gen_amino_acids_dict()

        # Initialises lists of properties to be displayed in the dataframe
        # columns
        domain_strand_ids = []
        tilt_angle = []
        strand_number = []
        shear_number = []
        res_ids = []
        edge_or_central = []
        fasta_seq = []
        int_ext = []
        core_ext = []
        tm_ext = []
        bckbn_phi_psi = []
        solv_acsblty = []

        unprocessed_list = []

        # Extracts property values for each domain_id
        for domain_id, dssp_df in sec_struct_dfs_dict.items():
            strands_list = [strand_num for strand_num in set(dssp_df['STRAND_NUM'].tolist())
                            if strand_num != '']

            for strand_num in strands_list:
                strand_df = dssp_df[dssp_df['STRAND_NUM'] == strand_num]
                strand_df = strand_df.reset_index(drop=True)

                # If the strand's parent structure is in the OPM database,
                # determines the orientation of the strand in the OM
                (reverse, res_ids_list, strand_coordinates, lower_bound, upper_bound,
                 unprocessed_list) = output_calcs.determine_strand_orientation(
                    domain_id, strand_num, strand_df, opm_database, unprocessed_list
                )

                # Gives strand a unique ID
                if strand_or_res == 'strand':
                    domain_strand_ids.append('{}_strand_{}'.format(domain_id, strand_num))
                elif strand_or_res == 'res':
                    domain_strand_ids += ['{}_strand_{}'.format(
                        domain_id, strand_num)]*len(res_ids_list)

                # Lists tilt angle of strand
                if self.code[0:4] in ['2.40']:
                    if domain_id in list(tilt_angles.keys()):
                        value = tilt_angles[domain_id]
                    else:
                        value = ''

                    if strand_or_res == 'strand':
                        tilt_angle.append(value)
                    elif strand_or_res == 'res':
                        tilt_angle += [value]*len(res_ids_list)

                # Lists number of strands in barrel
                if self.code[0:4] in ['2.40']:
                    if strand_or_res == 'strand':
                        strand_number.append(strand_numbers[domain_id])
                    elif strand_or_res == 'res':
                        strand_number += [strand_numbers[domain_id]]*len(res_ids_list)

                # Lists barrel shear number
                if self.code[0:4] in ['2.40']:
                    # shear_number.append(shear_numbers[domain_id])
                    if strand_or_res == 'strand':
                        shear_number.append('')
                    elif strand_or_res == 'res':
                        shear_number += ['']*len(res_ids_list)

                # Determines whether the strand is an edge or central strand
                if self.code[0:4] in ['2.60']:
                    edge_or_central_list = [label for label in
                                            set(strand_df['EDGE_OR_CNTRL'].tolist())]

                    if strand_or_res == 'strand':
                        edge_or_central.append(edge_or_central_list[0])
                    elif strand_or_res == 'res':
                        edge_or_central += [edge_or_central_list[0]]*len(res_ids_list)

                # Lists residue IDs in strand
                if reverse is True:
                    res_ids_list.reverse()

                if strand_or_res == 'strand':
                    res_ids.append(res_ids_list)
                elif strand_or_res == 'res':
                    res_ids += res_ids_list

                # Generates FASTA sequence of strand
                res_list = strand_df['RESNAME'].tolist()
                sequence = ''
                for res in res_list:
                    if res in list(amino_acids_dict.keys()):
                        sequence += amino_acids_dict[res]
                    else:
                        sequence += 'X'
                if reverse is True:
                    sequence = sequence[::-1]  # Reverses string

                if strand_or_res == 'strand':
                    fasta_seq.append(sequence)
                elif strand_or_res == 'res':
                    fasta_seq += sequence

                # Generates list of interior and exterior facing residues in
                # the strand
                if self.code[0:4] in ['2.40']:
                    int_ext_list = strand_df['INT_EXT'].tolist()
                    if reverse is True:
                        int_ext_list.reverse()

                    if strand_or_res == 'strand':
                        int_ext.append(int_ext_list)
                    elif strand_or_res == 'res':
                        int_ext += int_ext_list

                # Generates list of residues that form the beta-sandwich core
                if self.code[0:4] in ['2.60']:
                    core_ext_list = strand_df['CORE_OR_EXT'].tolist()
                    if reverse is True:
                        core_ext_list.reverse()

                    if strand_or_res == 'strand':
                        core_ext.append(core_ext_list)
                    elif strand_or_res == 'res':
                        core_ext += core_ext_list

                # Generates list of transmembrane and external residues in the
                # strand - MUST COME AFTER REVERSAL (OR NOT) OF RES_IDS_LIST
                if self.code[0:4] in ['2.40']:
                    tm_ext_list = output_calcs.determine_tm_or_ext(
                        domain_id, strand_num, res_ids_list, strand_coordinates,
                        lower_bound, upper_bound, unprocessed_list
                    )
                    # Don't need to reverse list as will match res_ids_list
                    # order (which has already been reversed)
                    if strand_or_res == 'strand':
                        tm_ext.append(tm_ext_list)
                    elif strand_or_res == 'res':
                        tm_ext += tm_ext_list

                # Generates list of phi and psi angles along the strand
                phi = strand_df['PHI'].tolist()
                psi = strand_df['PSI'].tolist()
                if reverse is True:
                    phi.reverse()
                    psi.reverse()
                bckbn_geom = [[phi[index], psi[index]] for index, value in
                              enumerate(phi)]

                if strand_or_res == 'strand':
                    bckbn_phi_psi.append(bckbn_geom)
                elif strand_or_res == 'res':
                    bckbn_phi_psi += bckbn_geom

                # Generates list of residue solvent accessibility along the
                # strand
                solv_acsblty_list = strand_df['SOLV_ACSBLTY'].tolist()
                if reverse is True:
                    solv_acsblty_list.reverse()

                if strand_or_res == 'strand':
                    solv_acsblty.append(solv_acsblty_list)
                elif strand_or_res == 'res':
                    solv_acsblty += solv_acsblty_list

        if self.code[0:4] in ['2.40']:
            beta_strands_df = pd.DataFrame({'STRAND_ID': domain_strand_ids,
                                            'TILT_ANGLE(DEGREES)': tilt_angle,
                                            'TOTAL_STRAND_NUMBER': strand_number,
                                            'SHEAR_NUMBER': shear_number,
                                            'RES_ID': res_ids,
                                            'FASTA': fasta_seq,
                                            'INT_EXT': int_ext,
                                            'TM_OR_EXT': tm_ext,
                                            'BCKBN_GEOM': bckbn_phi_psi,
                                            'SOLV_ACSBLTY': solv_acsblty})
            cols = beta_strands_df.columns.tolist()
            cols = ([cols[6]] + [cols[7]] + [cols[9]] + [cols[4]] + [cols[3]] +
                    [cols[1]] + [cols[2]] + [cols[8]] + [cols[0]] + [cols[5]])
            beta_strands_df = beta_strands_df[cols]
            beta_strands_df.to_pickle('Beta_{}_dataframe.pkl'.format(strand_or_res))
            beta_strands_df.to_csv('Beta_{}_dataframe.csv'.format(strand_or_res))

        elif self.code[0:4] in ['2.60']:
            beta_strands_df = pd.DataFrame({'STRAND_ID': domain_strand_ids,
                                            'EDGE_OR_CNTRL': edge_or_central,
                                            'RES_ID': res_ids,
                                            'FASTA': fasta_seq,
                                            'CORE_OR_EXT': core_ext,
                                            'BCKBN_GEOM': bckbn_phi_psi,
                                            'SOLV_ACSBLTY': solv_acsblty})
            cols = beta_strands_df.columns.tolist()
            cols = ([cols[6]] + [cols[2]] + [cols[4]] + [cols[3]] + [cols[1]]
                    + [cols[0]] + [cols[5]])
            beta_strands_df = beta_strands_df[cols]
            beta_strands_df.to_pickle('Beta_{}_dataframe.pkl'.format(strand_or_res))
            beta_strands_df.to_csv('Beta_{}_dataframe.csv'.format(strand_or_res))

        with open('Unprocessed_domains.txt', 'a') as unprocessed_file:
            unprocessed_file.write('\n\nFailed to determine strand orientation:\n')
            for domain_id in unprocessed_list:
                unprocessed_file.write('{}\n'.format(domain_id))
