
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
        # Determines orientation of input strand in domain if that domain is in
        # the OPM database. If the N-terminus is in the periplasm
        # (Z-coordinate = negative) and C-terminus is extracellular
        # (Z-coordinate = positive), the strand orientation is reversed in the
        # output dataframe.

        pdb_code = domain_id[0:4]
        res_ids_list = strand_df['RES_ID'].tolist()
        strand_coordinates = OrderedDict()
        in_database = False
        lower_bound = 0
        upper_bound = 0
        reverse = False

        if os.path.isfile('{}/{}.pdb'.format(opm_database, pdb_code)):
            print('Determining strand orientation in {} strand '
                  '{}'.format(domain_id, strand_id))
            in_database = True

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

        return (in_database, reverse, res_ids_list, strand_coordinates,
                lower_bound, upper_bound, unprocessed_list)

    def determine_tm_or_ext(domain_id, strand_id, res_ids_list, in_database,
                            strand_coordinates, lower_bound, upper_bound,
                            unprocessed_list):
        # If the parent structure is in the OPM database, labels residues in an
        # input strand as either 'transmembrane' or 'external'
        print('Calculating positions of strands in membrane in '
              '{}'.format(domain_id))

        tm_ext_list = [''] * len(res_ids_list)

        if (
            in_database is True
            and '{}_strand_{}'.format(domain_id, strand_id) not in unprocessed_list
        ):
            for index, res_id in enumerate(res_ids_list):
                z_coord = strand_coordinates[res_id]
                if z_coord < lower_bound or z_coord > upper_bound:
                    tm_ext_list[index] = 'external'
                else:
                    tm_ext_list[index] = 'transmembrane'

        return tm_ext_list

    def determine_strand_position_sandwich(strand_df):
        # Calculates displacement (both absolute and percentage), measured in
        # number of residues, of each residue from the centre of the strand
        strand_abs_displacement = ['']*strand_df.shape[0]
        strand_percentage_displacement = ['']*strand_df.shape[0]
        strand_centre = (strand_df.shape[0] - 1) / 2

        for row in range(strand_df.shape[0]):
            diff = abs(row - strand_centre)
            strand_abs_displacement[row] = diff

            percentage = (abs(row - strand_centre) / strand_centre) * 100
            strand_percentage_displacement[row] = round(percentage, 1)

        return strand_abs_displacement, strand_percentage_displacement

    def determine_strand_position_barrel(strand_df, tm_ext_list):
        # Calculates position of residues within the barrel TM region (position
        # is measured as a percentage of the total number of TM residues in the
        # strand, with the extracellular environment and the periplasm
        # corresponding to 0% and 100%, respectively)
        membrane_loc = ['']*strand_df.shape[0]

        if 'transmembrane' in tm_ext_list:
            max_depth = tm_ext_list.count('transmembrane')

            count = 0
            for row in range(strand_df.shape[0]):
                if tm_ext_list[row] == 'transmembrane':
                    count += 1
                    percentage = (count / max_depth) * 100
                    membrane_loc[row] = round(percentage, 1)

        return membrane_loc


def reverse_strand_lists(sub_list, reverse):
    if reverse is True:
        sub_list = sub_list[::-1]
    return sublist()


def append_to_output_lists(sub_list, long_list, res_ids_list, strand_or_res):
    if strand_or_res == 'strand':
        long_list.append(sub_list)
    elif strand_or_res == 'res':
        if isinstance(sub_list, list):
            long_list += sub_list
        else:
            long_list += [sub_list]*len(res_ids_list)

    return long_list


def reverse_and_append(residue_property, strand_df, long_list, reverse,
                       res_ids_list, strand_or_res):
    sub_list = strand_df[residue_property].tolist()
    sub_list = reverse_strand_lists(sub_list, reverse)
    long_list = append_to_output_lists(sub_list, long_list, res_ids_list, strand_or_res)

    return long_list


class gen_output(run_stages):

    def __init__(self, run_parameters, radius):
        run_stages.__init__(self, run_parameters)
        self.radius = radius

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
                num_of_edges = len(list(G.neighbors(strand_num)))
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
        sheet_number = []
        tilt_angle = []
        strand_number = []
        shear_number = []
        res_ids = []
        edge_or_central = []
        fasta_seq = []
        z_coords = []
        strand_abs_pos = []
        strand_percentage_pos = []
        tm_pos = []
        int_ext = []
        core_surf = []
        buried_surface_area = []
        tm_ext = []
        omega = []
        phi = []
        psi = []
        chi = []
        solv_acsblty = []
        neighbours = []
        brige_pairs = []
        van_der_waals = []
        h_bonds = []
        ionic = []
        ss_bonds = []
        pi_pi_stacking = []
        cation_pi = []
        minus_1_list = []
        plus_1_list = []
        minus_2_list = []
        plus_2_list = []

        unprocessed_list = []

        # Extracts property values for each domain_id
        for domain_id in list(sec_struct_dfs_dict.keys()):
            dssp_df = sec_struct_dfs_dict[domain_id]

            strands_list = [strand_num for strand_num in set(dssp_df['STRAND_NUM'].tolist())
                            if strand_num != '']

            for strand_num in strands_list:
                strand_df = dssp_df[dssp_df['STRAND_NUM'] == strand_num]
                strand_df = strand_df.reset_index(drop=True)

                # If the strand's parent structure is in the OPM database,
                # determines the orientation of the strand in the OM
                (in_database, reverse, res_ids_list, strand_coordinates,
                 lower_bound, upper_bound, unprocessed_list
                 ) = output_calcs.determine_strand_orientation(
                    domain_id, strand_num, strand_df, opm_database, unprocessed_list
                )

                # Gives strand a unique ID
                if strand_or_res == 'strand':
                    domain_strand_ids.append('{}_strand_{}'.format(domain_id, strand_num))
                elif strand_or_res == 'res':
                    domain_strand_ids += ['{}_strand_{}'.format(
                        domain_id, strand_num)]*len(res_ids_list)

                # Lists parent sheet of each strand
                sheet = strand_df['SHEET_NUM'].tolist()[0]
                if strand_or_res == 'strand':
                    sheet_number.append(sheet)
                elif strand_or_res == 'res':
                    sheet_number += ['{}'.format(sheet)]*len(res_ids_list)

                # Lists tilt angle of strand
                if self.code[0:4] in ['2.40']:
                    if domain_id in list(tilt_angles.keys()):
                        value = tilt_angles[domain_id]
                    else:
                        value = ''

                    tilt_angle = append_to_output_lists(
                        value, tilt_angle, res_ids_list, strand_or_res
                    )

                # Lists number of strands in barrel
                if self.code[0:4] in ['2.40']:
                    strand_number = append_to_output_lists(
                        strand_numbers[domain_id], strand_number, res_ids_list,
                        strand_or_res
                    )

                # Lists barrel shear number
                if self.code[0:4] in ['2.40']:
                    # shear_number.append(shear_numbers[domain_id])
                    shear_number = append_to_output_lists(
                        '', shear_number, res_ids_list, strand_or_res
                    )

                # Determines whether the strand is an edge or central strand
                if self.code[0:4] in ['2.60']:
                    edge_or_central_list = [label for label in
                                            set(strand_df['EDGE_OR_CNTRL'].tolist())]

                    edge_or_central = append_to_output_lists(
                        edge_or_central_list[0], edge_or_central, res_ids_list,
                        strand_or_res
                    )

                # Lists residue IDs in strand
                res_ids_list = reverse_strand_lists(res_ids_list, reverse)
                res_ids = append_to_output_lists(
                    res_ids_list, res_ids, res_ids_list, strand_or_res
                )

                # Generates FASTA sequence of strand
                res_list = strand_df['RESNAME'].tolist()
                sequence = ''
                for res in res_list:
                    if res in list(amino_acids_dict.keys()):
                        sequence += amino_acids_dict[res]
                    else:
                        sequence += 'X'

                sequence = reverse_strand_lists(sequence, reverse)
                fasta_seq = append_to_output_lists(
                    sequence, fasta_seq, res_ids_list, strand_or_res
                )

                # Generates list of strand z-coordinates
                if self.code[0:4] in ['2.40']:
                    if strand_coordinates:
                        strand_z_coords = list(strand_coordinates.values())
                    else:
                        strand_z_coords = ['']*len(res_ids_list)

                    strand_z_coords = reverse_strand_lists(strand_z_coords, reverse)
                    z_coords = append_to_output_lists(
                        strand_z_coords, z_coords, res_ids_list, strand_or_res
                    )

                # Generates list of interior and exterior facing residues in
                # the strand
                int_ext = reverse_and_append(
                    'INT_EXT', strand_df, int_ext, reverse, res_ids_list,
                    strand_or_res
                )

                # Generates list of residues that form the beta-sandwich core
                if self.code[0:4] in ['2.60']:
                    core_surf = reverse_and_append(
                        'CORE_OR_SURFACE', strand_df, core_surf, reverse,
                        res_ids_list, strand_or_res
                    )

                # Generates list of percentage values for the surface area of
                # each residue that is buried in the sandwich as compared with
                # the individual parent sheet
                if self.code[0:4] in ['2.60']:
                    buried_surface_area = reverse_and_append(
                        'BURIED_SURFACE_AREA(%)', strand_df, buried_surface_area,
                        reverse, res_ids_list, strand_or_res
                    )

                # Generates list of transmembrane and external residues in the
                # strand - MUST COME AFTER REVERSAL (OR NOT) OF RES_IDS_LIST
                if self.code[0:4] in ['2.40']:
                    tm_ext_list = output_calcs.determine_tm_or_ext(
                        domain_id, strand_num, res_ids_list, in_database,
                        strand_coordinates, lower_bound, upper_bound,
                        unprocessed_list
                    )
                    # Don't need to reverse list as will match res_ids_list
                    # order (which has already been reversed)
                    tm_ext = append_to_output_lists(
                        tm_ext_list, tm_ext, res_ids_list, strand_or_res
                    )

                # Determines positions of residues relative to centre of strand
                # in sandwiches
                if self.code[0:4] in ['2.60']:
                    (strand_abs_displacement, strand_percentage_displacement
                     ) = output_calcs.determine_strand_position_sandwich(strand_df)

                    strand_abs_displacement = reverse_strand_lists(
                        strand_abs_displacement, reverse
                    )
                    strand_abs_pos = append_to_output_lists(
                        strand_abs_displacement, strand_abs_pos, res_ids_list,
                        reverse
                    )

                    strand_percentage_displacement = reverse_strand_lists(
                        strand_percentage_displacement, reverse
                    )
                    strand_percentage_pos = append_to_output_lists(
                        strand_percentage_displacement, strand_percentage_pos,
                        res_ids_list, strand_or_res
                    )

                # Determines positions of residues in TM region of strand in
                # beta-barrel - MUST COME AFTER TM / EXTERNAL CALCULATION
                if self.code[0:4] in ['2.40']:
                    membrane_loc = output_calcs.determine_strand_position_barrel(
                        strand_df, tm_ext_list
                    )

                    # Don't need to reverse list as will match tm_ext_list
                    # order (which has already been reversed)
                    tm_pos = append_to_output_lists(
                        membrane_loc, tm_pos, res_ids_list, strand_or_res
                    )

                # Generates list of torsion angles along the strand
                omega = reverse_and_append(
                    'OMEGA', strand_df, omega, reverse, res_ids_list, strand_or_res
                )
                phi = reverse_and_append(
                    'PHI', strand_df, phi, reverse, res_ids_list, strand_or_res
                )
                psi = reverse_and_append(
                    'PSI', strand_df, psi, reverse, res_ids_list, strand_or_res
                )
                chi = reverse_and_append(
                    'CHI', strand_df, chi, reverse, res_ids_list, strand_or_res
                )

                # Generates list of residue solvent accessibility along the
                # strand
                solv_acsblty = reverse_and_append(
                    'SOLV_ACSBLTY', strand_df, solv_acsblty, reverse,
                    res_ids_list, strand_or_res
                )

                # Generates list of neighbouring residues
                neighbours = reverse_and_append(
                    'NEIGHBOURS', strand_df, neighbours, reverse, res_ids_list,
                    strand_or_res
                )

                # Generates list of residues forming backbone hydrogen-bonds
                # with the residue in questio
                bridge_pairs_list = []
                dssp_to_pdb = dict(zip(dssp_df['DSSP_NUM'].tolist(),
                                       dssp_df['RES_ID'].tolist()))
                for row in range(strand_df.shape[0]):
                    bridges = strand_df['BRIDGE_PAIRS'][row]
                    bridges_renum = []
                    for dssp_num in bridges:
                        if dssp_num in list(dssp_to_pdb.keys()):
                            bridges_renum.append(dssp_to_pdb[dssp_num])
                    bridge_pairs_list.append(bridges_renum)

                bridge_pairs_list = reverse_strand_lists(bridge_pairs_list, reverse)
                bridge_pairs = append_to_output_lists(
                    bridge_pairs_list, bridge_pairs, res_ids_list, strand_or_res
                )

                # Generates list of residues in van der Waals contact
                van_der_waals = reverse_and_append(
                    'VDW', strand_df, van_der_waals, reverse, res_ids_list,
                    strand_or_res
                )

                # Generates list of residues that form hydrogen bonds with the
                # residue in question
                h_bonds = reverse_and_append(
                    'HBOND', strand_df, h_bonds, reverse, res_ids_list,
                    strand_or_res
                )

                # Generates list of residues that form ionic bonds with the
                # residue in question
                ionic = reverse_and_append(
                    'IONIC', strand_df, ionic, reverse, res_ids_list, strand_or_res
                )

                # Generates list of disulfide bonds
                ss_bonds = reverse_and_append(
                    'SSBOND', strand_df, ss_bonds, reverse, res_ids_list,
                    strand_or_res
                )

                # Generates list of residues that form pi-pi stacking
                # interactions with the residue in question
                pi_pi_stacking = reverse_and_append(
                    'PIPISTACK', strand_df, pi_pi_stacking, reverse,
                    res_ids_list, strand_or_res
                )

                # Generates list of residues that form cation_pi interactions
                # with the residue in question
                cation_pi = reverse_and_append(
                    'PICATION', strand_df, cation_pi, reverse, res_ids_list,
                    strand_or_res
                )

                # Generates list of residues in +/-1 and +/2 positions to the
                # residue in question
                minus_1 = []
                plus_1 = []
                minus_2 = []
                plus_2 = []

                with (open('{}{}/{}.pdb1'.format(self.ring_database, domain_id[1:3], domain_id), 'r')
                      as pdb_file
                      ):
                    pdb_file_lines = pdb_file.readlines()
                pdb_file_lines = [line for line in pdb_file_lines
                                  if line[0:6].strip() in ['ATOM', 'HETATM']
                                  and line[12:16].strip() == 'CA']
                consec_res_id_list = [line[21:27].replace(' ', '') for line in
                                      pdb_file_lines]
                fasta_list = [amino_acids_dict[line[17:20]] for line in
                              pdb_file_lines]

                for row in range(strand_df.shape[0]):
                    res_id = strand_df['RES_ID'][row]
                    res_id_index = consec_res_id_list.index(res_id)
                    try:
                        m_1 = fasta_list[res_id_index-1]
                        p_1 = fasta_list[res_id_index+1]
                        m_2 = fasta_list['FASTA'].tolist()[res_id_index-2]
                        p_2 = fasta_list['FASTA'].tolist()[res_id_index+2]
                        print(m_1)
                        print(p_1)
                        print(m_2)
                        print(p_2)
                    except KeyError:
                        m_1 = ''
                        p_1 = ''
                        m_2 = ''
                        p_2 = ''
                    minus_1.append(m_1)
                    plus_1.append(p_1)
                    minus_2.append(m_2)
                    plus_2.append(p_2)

                if reverse is True:
                    minus_1.reverse()
                    plus_1.reverse()
                    minus_2.reverse()
                    plus_2.reverse()

                if strand_or_res == 'strand':
                    minus_1_list.append(minus_1)
                    plus_1_list.append(plus_1)
                    minus_2_list.append(minus_2)
                    plus_2_list.append(plus_2)
                elif strand_or_res == 'res':
                    minus_1_list += minus_1
                    plus_1_list += plus_1
                    minus_2_list += minus_2
                    plus_2_list += plus_2

        # Generates csv file of beta-barrel dataset
        if self.code[0:4] in ['2.40']:
            beta_strands_df_dict = OrderedDict({'STRAND_ID': domain_strand_ids,
                                                'TILT_ANGLE(DEGREES)': tilt_angle,
                                                'TOTAL_STRAND_NUMBER': strand_number,
                                                'SHEAR_NUMBER': shear_number,
                                                'RES_ID': res_ids,
                                                'FASTA': fasta_seq,
                                                'Z_COORDS': z_coords,
                                                'STRAND_POS(%)': tm_pos,
                                                'INT_EXT': int_ext,
                                                'TM_OR_EXT': tm_ext,
                                                'OMEGA': omega,
                                                'PHI': phi,
                                                'PSI': psi,
                                                'CHI_ANGLES': chi,
                                                'SOLV_ACSBLTY': solv_acsblty,
                                                'NEIGHBOURING_RESIDUES(<{}A)'.format(
                                                    self.radius
                                                ): neighbours,
                                                'BRIDGE_PAIRS': bridge_pairs_list,
                                                'VDW': van_der_waals,
                                                'HBOND': h_bonds,
                                                'IONIC': ionic,
                                                'SSBOND': ss_bonds,
                                                'PIPISTACK': pi_pi_stacking,
                                                'PICATION': cation_pi,
                                                'MINUS_1_POS': minus_1_list,
                                                'PLUS_1_POS': plus_1_list,
                                                'MINUS_2_POS': minus_2_list,
                                                'PLUS_2_POS': plus_2_list})
            beta_strands_df = pd.DataFrame(beta_strands_df_dict)
            beta_strands_df.to_pickle('Beta_{}_dataframe.pkl'.format(strand_or_res))
            beta_strands_df.to_csv('Beta_{}_dataframe.csv'.format(strand_or_res))

        # Generates csv file of beta-sandwich dataset
        elif self.code[0:4] in ['2.60']:
            beta_strands_df_dict = OrderedDict({'STRAND_ID': domain_strand_ids,
                                                'SHEET_ID': sheet_number,
                                                'EDGE_OR_CNTRL': edge_or_central,
                                                'RES_ID': res_ids,
                                                'FASTA': fasta_seq,
                                                'STRAND_POS(ABS)': strand_abs_pos,
                                                'STRAND_POS(%)': strand_percentage_pos,
                                                'INT_EXT': int_ext,
                                                'CORE_OR_SURFACE': core_surf,
                                                'BURIED_SURFACE_AREA(%)': buried_surface_area,
                                                'OMEGA': omega,
                                                'PHI': phi,
                                                'PSI': psi,
                                                'CHI_ANGLES': chi,
                                                'SOLV_ACSBLTY': solv_acsblty,
                                                'NEIGHBOURING_RESIDUES(<{}A)'.format(
                                                    self.radius
                                                ): neighbours,
                                                'BRIDGE_PAIRS': bridge_pairs_list,
                                                'VDW': van_der_waals,
                                                'HBOND': h_bonds,
                                                'IONIC': ionic,
                                                'SSBOND': ss_bonds,
                                                'PIPISTACK': pi_pi_stacking,
                                                'PICATION': cation_pi,
                                                'MINUS_1_POS': minus_1_list,
                                                'PLUS_1_POS': plus_1_list,
                                                'MINUS_2_POS': minus_2_list,
                                                'PLUS_2_POS': plus_2_list})
            beta_strands_df = pd.DataFrame(beta_strands_df_dict)
            beta_strands_df.to_pickle('Beta_{}_dataframe.pkl'.format(strand_or_res))
            beta_strands_df.to_csv('Beta_{}_dataframe.csv'.format(strand_or_res))

        with open('Unprocessed_domains.txt', 'a') as unprocessed_file:
            unprocessed_file.write('\n\nFailed to determine strand orientation:\n')
            for domain_id in unprocessed_list:
                unprocessed_file.write('{}\n'.format(domain_id))
