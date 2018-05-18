
import os
import networkx as nx
import pandas as pd
import numpy as np
from collections import OrderedDict
if __name__ == 'subroutines.output_dataframe':
    from subroutines.run_stages import run_stages
    from subroutines.variables import gen_amino_acids_dict
else:
    from datagen.subroutines.run_stages import run_stages
    from datagen.subroutines.variables import gen_amino_acids_dict


class output_calcs():

    def make_res_id_to_fasta_dict(domain_id, amino_acids_dict):
        # Creates dictionary of residue ids and their corresponding 1 letter
        # amino acid codes from the input dataframe (dssp_df)

        # Extracts list of consecutive res_ids from biological assembly PDB file
        with open('Biological_assemblies/{}.pdb'.format(domain_id[0:4]), 'r') as pdb_file:
            pdb_file_lines = pdb_file.readlines()
        pdb_file_lines = [line for line in pdb_file_lines if line[0:6].strip()
                          in ['ATOM', 'HETATM']]
        consec_res_id_list = []
        consec_res_name_list = []
        for line in pdb_file_lines:
            chain_res_num = line[21:27].replace(' ', '')
            res_name = line[17:20].strip()
            if chain_res_num not in consec_res_id_list:
                consec_res_id_list.append(chain_res_num)
                consec_res_name_list.append(res_name)

        # Generates res_id to FASTA code dictionary
        fasta_list = ['']*len(consec_res_id_list)
        for index, res_id in enumerate(consec_res_id_list):
            res_name = consec_res_name_list[index]
            if res_name in list(amino_acids_dict.keys()):
                fasta = amino_acids_dict[res_name]
            else:
                fasta = 'X'
            fasta_list[index] = fasta

        res_id_to_fasta_dict = dict(zip(consec_res_id_list, fasta_list))

        return res_id_to_fasta_dict, consec_res_id_list

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
                            if (
                                line[12:16].strip() == 'CA'
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

        tm_ext_list = ['']*len(res_ids_list)

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

    def convert_dssp_num_to_res_id(strand_df, dssp_df, dssp_to_pdb_dict):
        # Converts bridge pairs extracted from DSSP file to res_ids. NOTE: will
        # only list bridge partners that are also in the retained
        # beta-sandwich/barrel strands.
        bridge_pairs_list = []
        for row in range(strand_df.shape[0]):
            bridges = strand_df['BRIDGE_PAIRS'][row]
            bridges_renum = []
            for dssp_num in bridges:
                if dssp_num in list(dssp_to_pdb_dict.keys()):
                    bridges_renum.append(dssp_to_pdb_dict[dssp_num])
            bridge_pairs_list.append(bridges_renum)

        return bridge_pairs_list

    def find_minus_plus_residues(domain_id, strand_df, displacement,
                                 consec_res_id_list, res_id_to_fasta_dict):
        # Generates lists of residues +/- x (where x is a user-specified number
        # of residues) from each input residue
        minus_list = []
        plus_list = []
        for row in range(strand_df.shape[0]):
            res_id = strand_df['RES_ID'][row]
            res_id_index = consec_res_id_list.index(res_id)
            try:
                minus_res = consec_res_id_list[res_id_index - displacement]
                plus_res = consec_res_id_list[res_id_index + displacement]
                if res_id[0:1] != minus_res[0:1]:
                    minus_res = ''
                else:
                    minus_res = res_id_to_fasta_dict[minus_res]
                if res_id[0:1] != plus_res[0:1]:
                    plus_res = ''
                else:
                    plus_res = res_id_to_fasta_dict[plus_res]
            except IndexError:
                minus_res = ''
                plus_res = ''

            minus_list.append(minus_res)
            plus_list.append(plus_res)

        return minus_list, plus_list

    def convert_res_id_to_fasta(sub_list, res_id_to_fasta_dict, nested):
        # Converts input list of residue ids (= chain + residue number +
        # insertion code) into a list of 1 letter amino acid codes
        fasta_list = []
        for item in sub_list:
            if nested is False:
                try:
                    fasta = res_id_to_fasta_dict[item]
                except KeyError:
                    fasta = ''

            elif nested is True:
                fasta = []
                for res in item:
                    try:
                        fasta_indv = res_id_to_fasta_dict[res]
                    except KeyError:
                        fasta_indv = ''
                    fasta.append(fasta_indv)

            fasta_list.append(fasta)

        return fasta_list


def reverse_strand_lists(sub_list, reverse):
    if reverse is True:
        sub_list = sub_list[::-1]
    return sub_list


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

    return long_list, sub_list


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
                if (
                    dssp_df['STRAND_NUM'][row] in list(edge_or_central_dict.keys())
                    and dssp_df['ATMNAME'][row] == 'CA'
                ):
                    edge_or_central_list[row] = edge_or_central_dict[dssp_df['STRAND_NUM'][row]]
            df_edge_or_central = pd.DataFrame({'EDGE_OR_CNTRL': edge_or_central_list})
            dssp_df = pd.concat([dssp_df, df_edge_or_central], axis=1)
            sec_struct_dfs_dict[domain_id] = dssp_df

        return sec_struct_dfs_dict

    def write_beta_strand_dataframe(self, strand_or_res, sec_struct_dfs_dict,
                                    opm_database, dssp_to_pdb_dict, tilt_angles,
                                    strand_numbers, shear_numbers):
        # Generates dataframes of residues and strands in the retained domains.
        if __name__ == 'subroutines.output_dataframe':
            from subroutines.output_dataframe import output_calcs
        else:
            from datagen.subroutines.output_dataframe import output_calcs

        # Generates dictionary of amino acid 3 and 1 letter codes
        amino_acids_dict = gen_amino_acids_dict()

        # Initialises lists of properties to be displayed in the dataframe
        # columns
        domain_strand_ids = []
        sheet_number = []
        tilt_angle = []
        total_strand_number = []
        indv_strand_number = []
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
        bridge_pairs = []
        van_der_waals = []
        h_bonds = []
        ionic = []
        ss_bonds = []
        pi_pi_stacking = []
        cation_pi = []
        minus_1 = []
        plus_1 = []
        minus_2 = []
        plus_2 = []
        neighbours_fasta = []
        bridge_pairs_fasta = []
        van_der_waals_fasta = []
        h_bonds_fasta = []
        ionic_fasta = []
        ss_bonds_fasta = []
        pi_pi_stacking_fasta = []
        cation_pi_fasta = []

        unprocessed_list = []

        # Extracts property values for each domain_id
        for domain_id in list(sec_struct_dfs_dict.keys()):
            dssp_df = sec_struct_dfs_dict[domain_id]

            # Generates dictionary of residue ids (= chain + residue number +
            # insertion code) and names
            res_id_to_fasta_dict, consec_res_id_list = output_calcs.make_res_id_to_fasta_dict(
                domain_id, amino_acids_dict
            )

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
                sheet_number = append_to_output_lists(
                    sheet, sheet_number, res_ids_list, strand_or_res
                )

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
                    total_strand_number = append_to_output_lists(
                        strand_numbers[domain_id], total_strand_number, res_ids_list,
                        strand_or_res
                    )

                # Lists parent strand number
                indv_strand_number = append_to_output_lists(
                    strand_df['STRAND_NUM'].tolist()[0], indv_strand_number,
                    res_ids_list, strand_or_res
                )

                # Lists barrel shear number
                if self.code[0:4] in ['2.40']:
                    shear_number = append_to_output_lists(
                        '', shear_number, res_ids_list, strand_or_res
                    )

                # Determines whether the strand is an edge or central strand
                if self.code[0:4] in ['2.60']:
                    edge_or_central_list = [label for label in
                                            set(strand_df['EDGE_OR_CNTRL'].tolist())]
                    if len(edge_or_central_list) != 1:
                        unprocessed_list.append(domain_id)

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
                sequence = list(sequence)

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
                int_ext, int_ext_sub_list = reverse_and_append(
                    'INT_EXT', strand_df, int_ext, reverse, res_ids_list,
                    strand_or_res
                )

                # Generates list of residues that form the beta-sandwich core
                if self.code[0:4] in ['2.60']:
                    core_surf, core_surf_sub_list = reverse_and_append(
                        'CORE_OR_SURFACE', strand_df, core_surf, reverse,
                        res_ids_list, strand_or_res
                    )

                # Generates list of percentage values for the surface area of
                # each residue that is buried in the sandwich as compared with
                # the individual parent sheet
                if self.code[0:4] in ['2.60']:
                    buried_surface_area, buried_surface_area_sub_list = reverse_and_append(
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
                        strand_or_res
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
                omega, omega_sub_list = reverse_and_append(
                    'OMEGA', strand_df, omega, reverse, res_ids_list, strand_or_res
                )
                phi, phi_sub_list = reverse_and_append(
                    'PHI', strand_df, phi, reverse, res_ids_list, strand_or_res
                )
                psi, psi_sub_list = reverse_and_append(
                    'PSI', strand_df, psi, reverse, res_ids_list, strand_or_res
                )
                chi, chi_sub_list = reverse_and_append(
                    'CHI', strand_df, chi, reverse, res_ids_list, strand_or_res
                )

                # Generates list of residue solvent accessibility along the
                # strand
                solv_acsblty, solv_acsblty_sub_list = reverse_and_append(
                    'SOLV_ACSBLTY', strand_df, solv_acsblty, reverse,
                    res_ids_list, strand_or_res
                )

                # Generates list of neighbouring residues
                neighbours, neighbours_sub_list = reverse_and_append(
                    'NEIGHBOURS', strand_df, neighbours, reverse, res_ids_list,
                    strand_or_res
                )

                neighbours_fasta_sub_list = output_calcs.convert_res_id_to_fasta(
                    neighbours_sub_list, res_id_to_fasta_dict, True
                )
                neighbours_fasta = append_to_output_lists(
                    neighbours_fasta_sub_list, neighbours_fasta, res_ids_list,
                    strand_or_res
                )

                # Generates list of residues forming backbone hydrogen-bonds
                # with the residue in question
                bridge_pairs_list = output_calcs.convert_dssp_num_to_res_id(
                    strand_df, dssp_df, dssp_to_pdb_dict
                )
                bridge_pairs_list = reverse_strand_lists(bridge_pairs_list, reverse)
                bridge_pairs = append_to_output_lists(
                    bridge_pairs_list, bridge_pairs, res_ids_list, strand_or_res
                )

                bridge_pairs_fasta_sub_list = output_calcs.convert_res_id_to_fasta(
                    bridge_pairs_list, res_id_to_fasta_dict, True
                )
                bridge_pairs_fasta = append_to_output_lists(
                    bridge_pairs_fasta_sub_list, bridge_pairs_fasta,
                    res_ids_list, strand_or_res
                )

                # Generates list of residues in van der Waals contact
                van_der_waals, van_der_waals_sub_list = reverse_and_append(
                    'VDW', strand_df, van_der_waals, reverse, res_ids_list,
                    strand_or_res
                )

                van_der_waals_fasta_sub_list = output_calcs.convert_res_id_to_fasta(
                    van_der_waals_sub_list, res_id_to_fasta_dict, True
                )
                van_der_waals_fasta = append_to_output_lists(
                    van_der_waals_fasta_sub_list, van_der_waals_fasta,
                    res_ids_list, strand_or_res
                )

                # Generates list of residues that form hydrogen bonds with the
                # residue in question
                h_bonds, h_bonds_sub_list = reverse_and_append(
                    'HBOND', strand_df, h_bonds, reverse, res_ids_list,
                    strand_or_res
                )

                h_bonds_fasta_sub_list = output_calcs.convert_res_id_to_fasta(
                    h_bonds_sub_list, res_id_to_fasta_dict, True
                )
                h_bonds_fasta = append_to_output_lists(
                    h_bonds_fasta_sub_list, h_bonds_fasta, res_ids_list,
                    strand_or_res
                )

                # Generates list of residues that form ionic bonds with the
                # residue in question
                ionic, ionic_sub_list = reverse_and_append(
                    'IONIC', strand_df, ionic, reverse, res_ids_list, strand_or_res
                )

                ionic_fasta_sub_list = output_calcs.convert_res_id_to_fasta(
                    ionic_sub_list, res_id_to_fasta_dict, True
                )
                ionic_fasta = append_to_output_lists(
                    ionic_fasta_sub_list, ionic_fasta, res_ids_list,
                    strand_or_res
                )

                # Generates list of disulfide bonds
                ss_bonds, ss_bonds_sub_list = reverse_and_append(
                    'SSBOND', strand_df, ss_bonds, reverse, res_ids_list,
                    strand_or_res
                )

                ss_bonds_fasta_sub_list = output_calcs.convert_res_id_to_fasta(
                    ss_bonds_sub_list, res_id_to_fasta_dict, True
                )
                ss_bonds_fasta = append_to_output_lists(
                    ss_bonds_fasta_sub_list, ss_bonds_fasta, res_ids_list,
                    strand_or_res
                )

                # Generates list of residues that form pi-pi stacking
                # interactions with the residue in question
                pi_pi_stacking, pi_pi_stacking_sub_list = reverse_and_append(
                    'PIPISTACK', strand_df, pi_pi_stacking, reverse,
                    res_ids_list, strand_or_res
                )

                pi_pi_stacking_fasta_sub_list = output_calcs.convert_res_id_to_fasta(
                    pi_pi_stacking_sub_list, res_id_to_fasta_dict, True
                )
                pi_pi_stacking_fasta = append_to_output_lists(
                    pi_pi_stacking_fasta_sub_list, pi_pi_stacking_fasta,
                    res_ids_list, strand_or_res
                )

                # Generates list of residues that form cation_pi interactions
                # with the residue in question
                cation_pi, cation_pi_sub_list = reverse_and_append(
                    'PICATION', strand_df, cation_pi, reverse, res_ids_list,
                    strand_or_res
                )

                cation_pi_fasta_sub_list = output_calcs.convert_res_id_to_fasta(
                    cation_pi_sub_list, res_id_to_fasta_dict, True
                )
                cation_pi_fasta = append_to_output_lists(
                    cation_pi_fasta_sub_list, cation_pi_fasta, res_ids_list,
                    strand_or_res
                )

                # Generates list of residues in +/-1 and +/2 positions to the
                # residue in question
                minus_1_list, plus_1_list = output_calcs.find_minus_plus_residues(
                    domain_id, strand_df, 1, consec_res_id_list, res_id_to_fasta_dict
                )
                minus_2_list, plus_2_list = output_calcs.find_minus_plus_residues(
                    domain_id, strand_df, 2, consec_res_id_list, res_id_to_fasta_dict
                )

                minus_1_list = reverse_strand_lists(minus_1_list, reverse)
                minus_1 = append_to_output_lists(
                    minus_1_list, minus_1, res_ids_list, strand_or_res
                )
                plus_1_list = reverse_strand_lists(plus_1_list, reverse)
                plus_1 = append_to_output_lists(
                    plus_1_list, plus_1, res_ids_list, strand_or_res
                )
                minus_2_list = reverse_strand_lists(minus_2_list, reverse)
                minus_2 = append_to_output_lists(
                    minus_2_list, minus_2, res_ids_list, strand_or_res
                )
                plus_2_list = reverse_strand_lists(plus_2_list, reverse)
                plus_2 = append_to_output_lists(
                    plus_2_list, plus_2, res_ids_list, strand_or_res
                )

        for index, strand_id in enumerate(domain_strand_ids):
            if strand_id.split('_')[0] in unprocessed_list:
                domain_strand_ids[index] = np.nan
                sheet_number[index] = np.nan
                tilt_angle[index] = np.nan
                total_strand_number[index] = np.nan
                indv_strand_number[index] = np.nan
                shear_number[index] = np.nan
                res_ids[index] = np.nan
                edge_or_central[index] = np.nan
                fasta_seq[index] = np.nan
                z_coords[index] = np.nan
                strand_abs_pos[index] = np.nan
                strand_percentage_pos[index] = np.nan
                tm_pos[index] = np.nan
                int_ext[index] = np.nan
                core_surf[index] = np.nan
                buried_surface_area[index] = np.nan
                tm_ext[index] = np.nan
                omega[index] = np.nan
                phi[index] = np.nan
                psi[index] = np.nan
                chi[index] = np.nan
                solv_acsblty[index] = np.nan
                neighbours[index] = np.nan
                bridge_pairs[index] = np.nan
                van_der_waals[index] = np.nan
                h_bonds[index] = np.nan
                ionic[index] = np.nan
                ss_bonds[index] = np.nan
                pi_pi_stacking[index] = np.nan
                cation_pi[index] = np.nan
                minus_1[index] = np.nan
                plus_1[index] = np.nan
                minus_2[index] = np.nan
                plus_2[index] = np.nan
                neighbours_fasta[index] = np.nan
                bridge_pairs_fasta[index] = np.nan
                van_der_waals_fasta[index] = np.nan
                h_bonds_fasta[index] = np.nan
                ionic_fasta[index] = np.nan
                ss_bonds_fasta[index] = np.nan
                pi_pi_stacking_fasta[index] = np.nan
                cation_pi_fasta[index] = np.nan

                # Generates csv file of beta-barrel dataset
        if self.code[0:4] in ['2.40']:
            beta_strands_df_dict = OrderedDict({'STRAND_ID': domain_strand_ids,
                                                'TILT_ANGLE(DEGREES)': tilt_angle,
                                                'TOTAL_STRAND_NUMBER': total_strand_number,
                                                'PARENT_STRAND_NUMBER': indv_strand_number,
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
                                                'BRIDGE_PAIRS': bridge_pairs,
                                                'VDW': van_der_waals,
                                                'HBOND': h_bonds,
                                                'IONIC': ionic,
                                                'SSBOND': ss_bonds,
                                                'PIPISTACK': pi_pi_stacking,
                                                'PICATION': cation_pi,
                                                'MINUS_1_POS': minus_1,
                                                'PLUS_1_POS': plus_1,
                                                'MINUS_2_POS': minus_2,
                                                'PLUS_2_POS': plus_2,
                                                'NEIGHBOURING_RESIDUES(<{}A)_FASTA'.format(
                                                    self.radius
                                                ): neighbours_fasta,
                                                'BRIDGE_PAIRS_FASTA': bridge_pairs_fasta,
                                                'VDW_FASTA': van_der_waals_fasta,
                                                'HBOND_FASTA': h_bonds_fasta,
                                                'IONIC_FASTA': ionic_fasta,
                                                'SSBOND_FASTA': ss_bonds_fasta,
                                                'PIPISTACK_FASTA': pi_pi_stacking_fasta,
                                                'PICATION_FASTA': cation_pi_fasta})
            beta_strands_df = pd.DataFrame(beta_strands_df_dict)
            beta_strands_df = beta_strands_df.dropna()
            beta_strands_df = beta_strands_df.reset_index(drop=True)
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
                                                'BRIDGE_PAIRS': bridge_pairs,
                                                'VDW': van_der_waals,
                                                'HBOND': h_bonds,
                                                'IONIC': ionic,
                                                'SSBOND': ss_bonds,
                                                'PIPISTACK': pi_pi_stacking,
                                                'PICATION': cation_pi,
                                                'MINUS_1_POS': minus_1,
                                                'PLUS_1_POS': plus_1,
                                                'MINUS_2_POS': minus_2,
                                                'PLUS_2_POS': plus_2,
                                                'NEIGHBOURING_RESIDUES(<{}A)_FASTA'.format(
                                                    self.radius
                                                ): neighbours_fasta,
                                                'BRIDGE_PAIRS_FASTA': bridge_pairs_fasta,
                                                'VDW_FASTA': van_der_waals_fasta,
                                                'HBOND_FASTA': h_bonds_fasta,
                                                'IONIC_FASTA': ionic_fasta,
                                                'SSBOND_FASTA': ss_bonds_fasta,
                                                'PIPISTACK_FASTA': pi_pi_stacking_fasta,
                                                'PICATION_FASTA': cation_pi_fasta})
            beta_strands_df = pd.DataFrame(beta_strands_df_dict)
            beta_strands_df = beta_strands_df.dropna()
            beta_strands_df = beta_strands_df.reset_index(drop=True)
            beta_strands_df.to_pickle('Beta_{}_dataframe.pkl'.format(strand_or_res))
            beta_strands_df.to_csv('Beta_{}_dataframe.csv'.format(strand_or_res))

        with open('Unprocessed_domains.txt', 'a') as unprocessed_file:
            unprocessed_file.write('\n\nERROR with strand processing:\n')
            for domain_id in set(unprocessed_list):
                unprocessed_file.write('{}\n'.format(domain_id))
