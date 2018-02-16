
import sys
import math
import itertools
import isambard
import pandas as pd
import numpy as np
import networkx as nx
from collections import OrderedDict
if __name__ == 'subroutines.naccess':
    from subroutines.run_stages import run_stages
else:
    from datagen.subroutines.run_stages import run_stages


class naccess_solv_acsblty_calcs():

    def calculate_barrel_solv_acsblty(sheets):
        # Calculates the solvent accessibility of beta-sheets in domain X to
        # identify which is the largest (and hence is assumed to form the
        # beta-barrel)
        solv_acsblty_dict = OrderedDict()

        for sheet_id in sheets:
            sheet_sub_strings = []
            with open('Beta_strands/{}.pdb'.format(sheet_id), 'r') as pdb_file:
                for line in pdb_file.readlines():
                    line_start = line[0:16]
                    line_end = line[17:]
                    new_line = line_start + ' ' + line_end
                    sheet_sub_strings.append(new_line)
            sheet_string = ''.join(sheet_sub_strings)

            print('Calculating solvent accessible surface area for {}'.format(sheet_id))
            naccess_out = isambard.external_programs.naccess.run_naccess(
                sheet_string, 'rsa', path=False, include_hetatms=True
            )
            naccess_out = naccess_out.split('\n')
            naccess_out = [line for line in naccess_out if line != '']
            # Determines the solvent accessibility of the sheet as the sum of
            # the absolute solvent accessibilities of all residues
            sheet_solv_acsblty = naccess_out[-1]
            sheet_solv_acsblty = float(sheet_solv_acsblty.split()[1])
            solv_acsblty_dict[sheet_solv_acsblty] = sheet_id

        # Selects the beta-sheet with the highest solvent accessibility
        max_solv_acsblty = max(list(solv_acsblty_dict.keys()))
        beta_barrel = solv_acsblty_dict[max_solv_acsblty]
        solv_acsblty_dict = OrderedDict()
        solv_acsblty_dict[max_solv_acsblty] = beta_barrel

        return solv_acsblty_dict

    def calculate_sandwich_solv_acsblty(sheets):
        # Calculates the solvent accessibility of beta-sheets to identify the
        # two which share the largest buried surface area and hence form the
        # beta-sandwich in domain X
        combinations = list(itertools.combinations(list(sheets.keys()), 2))
        solv_acsblty_dict = OrderedDict()
        for sheet_pair in combinations:
            # Calculates solvent accessibility of complex of both sheets in pair
            print('Calculating solvent accessible surface area for {}'.format(sheet_pair))
            sheet_sub_strings = []
            for sheet_id in sheet_pair:
                with open('Beta_strands/{}.pdb'.format(sheet_id), 'r') as pdb_file:
                    for line in pdb_file.readlines():
                        line_start = line[0:16]
                        line_end = line[17:]
                        new_line = line_start + ' ' + line_end
                        sheet_sub_strings.append(new_line)
            sheet_string = ''.join(sheet_sub_strings)

            naccess_out = isambard.external_programs.naccess.run_naccess(
                sheet_string, 'rsa', path=False, include_hetatms=True
            )
            naccess_out = naccess_out.split('\n')
            naccess_out = [line for line in naccess_out if line != '']
            # Determines the solvent accessibility of the complex as the sum of
            # the absolute solvent accessibilities of all residues
            complex_solv_acsblty = naccess_out[-1]
            complex_solv_acsblty = float(complex_solv_acsblty.split()[1])

            # Calculates solvent accessibility of individual sheets
            sum_sheet_solv_acsblty = 0
            for sheet_id in sheet_pair:
                sheet_sub_strings = []
                print('Calculating solvent accessible surface area for {}'.format(sheet_id))
                with open('Beta_strands/{}.pdb'.format(sheet_id), 'r') as pdb_file:
                    for line in pdb_file.readlines():
                        line_start = line[0:16]
                        line_end = line[17:]
                        new_line = line_start + ' ' + line_end
                        sheet_sub_strings.append(new_line)
                sheet_string = ''.join(sheet_sub_strings)
                naccess_out = isambard.external_programs.naccess.run_naccess(
                    sheet_string, 'rsa', path=False, include_hetatms=True
                )
                naccess_out = naccess_out.split('\n')
                naccess_out = [line for line in naccess_out if line != '']
                # Determines the solvent accessibility of the individual sheet as
                # the sum of the absolute solvent accessibilities of all residues
                sheet_solv_acsblty = naccess_out[-1]
                sheet_solv_acsblty = float(sheet_solv_acsblty.split()[1])
                sum_sheet_solv_acsblty += sheet_solv_acsblty

            buried_surface_area = sum_sheet_solv_acsblty - complex_solv_acsblty
            solv_acsblty_dict[buried_surface_area] = [sheet_pair[0], sheet_pair[1]]

        return solv_acsblty_dict

    def calculate_residue_solv_acsblty(domain_id):
        # Calculates solvent accessibility of individual residues (in the
        # context of the parent CATH / SCOPe domain)
        res_solv_acsblty = OrderedDict()
        sheet_sub_strings = []

        print('Calculating solvent accessible surface areas of individual '
              'residues in {}'.format(domain_id))

        with open('Entire_domains/{}.pdb'.format(domain_id), 'r') as pdb_file:
            for line in pdb_file.readlines():
                line_start = line[0:16]
                line_end = line[17:]
                new_line = line_start + ' ' + line_end
                sheet_sub_strings.append(new_line)
        sheet_string = ''.join(sheet_sub_strings)
        naccess_out = isambard.external_programs.naccess.run_naccess(
            sheet_string, 'rsa', path=False, include_hetatms=True
        )
        naccess_out = naccess_out.split('\n')
        for line in naccess_out:
            if line[0:3] in ['RES', 'HEM']:
                chain = line[8:9].strip()
                res_num = line[9:13].strip()
                ins_code = line[13:14].strip()
                chain_res_num = chain + res_num + ins_code
                res_solv_acsblty[chain_res_num] = float(line[14:22])

        return res_solv_acsblty

    def add_naccess_info_to_df(domain_id, dssp_df, sheets, code,
                               solv_acsblty_dict, res_solv_acsblty,
                               sec_struct_dfs_dict, domain_sheets_dict,
                               unprocessed_list):
        # Uses the previously calculated solvent accessibility values to
        # determine which sheet/s form/s the beta-barrel/sandwich, and then
        # prunes the dataframe (dssp_df) to remove all other sheets
        solv_acsblty_list = ['']*dssp_df.shape[0]

        # Discounts a beta-sandwich domain from further analysis if the solvent
        # accessibility calculations show that none of the retained beta-sheets
        # are in contact with one another
        if (code[0:4] in ['2.60']
                    and float(max(list(solv_acsblty_dict.keys()))) == 0.0
                ):
            unprocessed_list.append(domain_id)
            sec_struct_dfs_dict[domain_id] = None
            for sheet_id in sheets:
                domain_sheets_dict[sheet_id] = None

            return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

        else:
            sandwich = solv_acsblty_dict[max(list(solv_acsblty_dict.keys()))]

            # Removes sheets that do not form part of the sandwich / barrel
            # from further analysis
            for sheet_id in sheets:
                if sheet_id not in sandwich:
                    domain_sheets_dict[sheet_id] = None

            sheets_retained = [sheet_id.replace('{}_sheet_'.format(domain_id), '')
                               for sheet_id in sandwich]
            sub_dssp_df = dssp_df[dssp_df['SHEET_NUM'].isin(sheets_retained)]
            chain_res_num = sub_dssp_df['RES_ID'].tolist()

            for row in range(dssp_df.shape[0]):
                res_id = dssp_df['RES_ID'][row]

                if res_id not in chain_res_num:
                    dssp_df.loc[row, 'REC'] = None

                if (res_id in list(res_solv_acsblty.keys())
                        and dssp_df['ATMNAME'][row] == 'CA'
                        ):
                    solv_acsblty_list[row] = res_solv_acsblty[res_id]

        solv_acsblty_df = pd.DataFrame({'SOLV_ACSBLTY': solv_acsblty_list})
        dssp_df = pd.concat([dssp_df, solv_acsblty_df], axis=1)
        dssp_df = dssp_df[dssp_df['REC'].notnull()]
        dssp_df = dssp_df.reset_index(drop=True)
        sec_struct_dfs_dict[domain_id] = dssp_df

        return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

    def calculate_int_ext_barrel(domain_id, dssp_df, sheets,
                                 sec_struct_dfs_dict, domain_sheets_dict,
                                 unprocessed_list):
        # Determines which face of the beta-sheet forms the interior of the
        # barrel and which forms the exterior based upon whether or not the
        # Calpha atom of a residue is located closer to the central axis of the
        # barrel than the N and C atoms of that residue.

        # Checks that correct number of sheets has been retained for analysis
        if len(sheets) != 1:
            print('ERROR: more than 1 sheet retained in {} following solvent '
                  'accessibility calculation'.format(domain_id))
            sec_struct_dfs_dict[domain_id] = None
            for sheet in sheets:
                domain_sheets_dict[sheet] = None
            unprocessed_list.append(domain_id)

            return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

        # Creates dataframe of residues in barrel
        sheet_num = sheets[0].replace('{}_sheet_'.format(domain_id), '')
        sheets_df = dssp_df[dssp_df['SHEET_NUM'] == sheet_num]
        sheets_df = sheets_df.reset_index(drop=True)

        # Creates ampal object from barrel
        barrel = isambard.ampal.convert_pdb_to_ampal('Beta_strands/{}.pdb'.format(sheets[0]))

        # Finds network of interacting neighbouring strands
        networks = [network for key, network in domain_sheets_dict.items()
                    if domain_id in key]
        G = networks[0]

        nodes_dict = {}
        strands = nx.nodes(G)
        for strand in strands:
            nodes_dict[strand] = len(G.neighbors(strand))
        try:
            while min(list(nodes_dict.values())) < 2:
                for strand in strands:
                    if len(G.neighbors(strand)) < 2:
                        G.remove_node(strand)
                        del nodes_dict[strand]
                strands = nx.nodes(G)
                for strand in strands:
                    nodes_dict[strand] = len(G.neighbors(strand))
        except ValueError:
            unprocessed_list.append(domain_id)
            return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

        cycles = [strand for cycle in nx.cycle_basis(G) for strand in cycle]
        if len(set(cycles)) != len(nx.nodes(G)):
            unprocessed_list.append(domain_id)
            return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

        strands = nx.cycle_basis(G)[0]  # Finds first complete cycle

        # Finds reference axis
        count = 0
        coords_res_1 = []
        coords_res_n = []
        for strand in strands:
            count += 1

            strand_df = sheets_df[sheets_df['STRAND_NUM'] == strand]
            strand_df = strand_df.reset_index(drop=True)

            row_num = strand_df.shape[0]
            if row_num % 2 == 1:
                row_num += 1

            index_1 = (row_num / 2) - 1
            index_n = row_num / 2
            if count % 2 == 1:
                index_1 = row_num / 2
                index_n = (row_num / 2) - 1

            res_1_x = strand_df['XPOS'][index_1]
            res_1_y = strand_df['YPOS'][index_1]
            res_1_z = strand_df['ZPOS'][index_1]
            res_n_x = strand_df['XPOS'][index_n]
            res_n_y = strand_df['YPOS'][index_n]
            res_n_z = strand_df['ZPOS'][index_n]

            coords_res_1.append((res_1_x, res_1_y, res_1_z))
            coords_res_n.append((res_n_x, res_n_y, res_n_z))

        coords_res_1 = np.array(coords_res_1)
        coords_res_n = np.array(coords_res_n)
        x_coord_1 = np.sum(coords_res_1[:, 0])
        y_coord_1 = np.sum(coords_res_1[:, 1])
        z_coord_1 = np.sum(coords_res_1[:, 2])
        x_coord_n = np.sum(coords_res_n[:, 0])
        y_coord_n = np.sum(coords_res_n[:, 1])
        z_coord_n = np.sum(coords_res_n[:, 2])

        # Aligns the barrel with z = 0
        s1 = [x_coord_1, y_coord_1, z_coord_1]
        e1 = [x_coord_n, y_coord_n, z_coord_n]
        s2 = [0.0, 0.0, 0.0]
        e2 = [0.0, 0.0, 1.0]
        translation, angle, axis, point = isambard.geometry.find_transformations(
            s1, e1, s2, e2
        )
        barrel.rotate(angle, axis, point=point)  # Rotation must be performed
        # before translation
        barrel.translate(translation)

        # Generates dictionary of residues and their new xyz coordinates
        barrel_pdb_string = barrel.make_pdb().split('\n')
        xy_dict = OrderedDict()

        for line in barrel_pdb_string:
            if line[0:6].strip() in ['ATOM', 'HETATM']:
                atom_id = line[21:27].replace(' ', '') + '_' + line[12:16].strip()
                line_segments = line.split('.')
                x_coord_atom = float(line_segments[0].split()[-1] + '.' + line_segments[1][0:3])
                y_coord_atom = float(line_segments[1][3:] + '.' + line_segments[2][0:3])
                xy_coords = np.array([[x_coord_atom],
                                      [y_coord_atom]])
                xy_dict[atom_id] = xy_coords

        # Initialises records of interior / exterior facing residues
        print('Determining interior and exterior residues in {}'.format(sheets[0]))
        int_ext_list = ['']*dssp_df.shape[0]
        int_ext_dict = OrderedDict()

        # Calculates average xy coordinates of the barrel Calpha atoms
        x_sum = 0
        y_sum = 0
        count = 0
        for res_id in list(xy_dict.keys()):
            if res_id.split('_')[0] in sheets_df['RES_ID'].tolist():
                x_sum += xy_dict[res_id][0][0]
                y_sum += xy_dict[res_id][1][0]
                count += 1

        if count > 0:
            x_average = x_sum / count
            y_average = y_sum / count
            com = np.array([[x_average],
                            [y_average]])
        else:
            print('ERROR: Unable to identify xyz coordinates for interior / '
                  'exterior calculation')

        # Calculates whether a residue faces towards the interior or the
        # exterior of the barrel
        res_dict = OrderedDict()
        for index, res_id in enumerate(sheets_df['RES_ID'].tolist()):
            res_dict[res_id] = sheets_df['RESNAME'].tolist()[index]

        for res_id in list(res_dict.keys()):
            if ((not 'GLY' in res_dict[res_id])
                    and all('{}_{}'.format(res_id, x) in list(xy_dict.keys())
                            for x in ['N', 'CA', 'C', 'O', 'CB']
                            )
                ):
                # Calculates angle between C_alpha, C_beta and the centre of
                # mass
                c_alpha_x = xy_dict['{}_CA'.format(res_id)][0][0]
                c_alpha_y = xy_dict['{}_CA'.format(res_id)][1][0]

                c_beta_x = xy_dict['{}_CB'.format(res_id)][0][0]
                c_beta_y = xy_dict['{}_CB'.format(res_id)][1][0]

                c_alpha_beta_vector = np.array([[c_beta_x-c_alpha_x],
                                                [c_beta_y-c_alpha_y]])
                c_alpha_com_vector = np.array([[com[0][0]-c_alpha_x],
                                               [com[1][0]-c_alpha_y]])

                c_alpha_numerator = np.sum((c_alpha_beta_vector*c_alpha_com_vector), axis=0)[0]

                c_alpha_beta_magnitude = math.sqrt(((c_alpha_beta_vector[0][0])**2)
                                                   + ((c_alpha_beta_vector[1][0])**2))
                c_alpha_com_magnitude = math.sqrt(((c_alpha_com_vector[0][0])**2)
                                                  + ((c_alpha_com_vector[1][0])**2))
                c_alpha_denominator = c_alpha_beta_magnitude*c_alpha_com_magnitude

                cos_c_alpha_angle = (c_alpha_numerator / c_alpha_denominator)
                c_alpha_angle = math.degrees(math.acos(cos_c_alpha_angle))

                if c_alpha_angle < 90:
                    int_ext_dict[res_id] = 'interior'
                elif c_alpha_angle > 90:
                    int_ext_dict[res_id] = 'exterior'

        # Updates dataframe with solvent accessibility information
        for row in range(dssp_df.shape[0]):
            if (dssp_df['RES_ID'][row] in list(int_ext_dict.keys())
                    and dssp_df['ATMNAME'][row] == 'CA'
                    ):
                int_ext_list[row] = int_ext_dict[dssp_df['RES_ID'][row]]
        int_ext_df = pd.DataFrame({'INT_EXT': int_ext_list})
        dssp_df = pd.concat([dssp_df, int_ext_df], axis=1)
        sec_struct_dfs_dict[domain_id] = dssp_df

        return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

    """
    def calculate_int_ext_barrel(domain_id, dssp_df, sheets,
                                 sec_struct_dfs_dict, domain_sheets_dict,
                                 unprocessed_list):
        # Determines which face of the beta-sheet forms the interior of the
        # barrel and which forms the exterior. **NOTE: this function assumes
        # that adjacent residues point in opposite directions.**
        # TO DO: UPDATE TO USE BACKBONE PHI AND PSI ANGLES TO DETERMINE IF THE
        # STRAND CONTAINS A BULGE OR NOT

        # Checks that correct number of sheets has been retained for analysis
        if len(sheets) != 1:
            print('ERROR: more than 1 sheet retained in {} following solvent '
                  'accessibility calculation'.format(domain_id))
            sec_struct_dfs_dict[domain_id] = None
            for sheet in sheets:
                domain_sheets_dict[sheet] = None
            unprocessed_list.append(domain_id)

            return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

        # Initialises records of interior / exterior facing residues
        int_ext_list = ['']*dssp_df.shape[0]
        int_ext_dict = OrderedDict()

        # Calculates solvent accessibility of both faces of barrel
        sheet_sub_strings = []
        print('Calculating solvent accessible surface area for {}'.format(domain_id))

        with open('Beta_strands/{}.pdb'.format(sheets[0]), 'r') as pdb_file:
            for line in pdb_file.readlines():
                line_start = line[0:16]
                line_end = line[17:]
                new_line = line_start + ' ' + line_end
                sheet_sub_strings.append(new_line)
        sheet_string = ''.join(sheet_sub_strings)
        naccess_out = isambard.external_programs.naccess.run_naccess(
            sheet_string, 'rsa', path=False, include_hetatms=True
        )
        naccess_out = naccess_out.split('\n')
        naccess_out = [line for line in naccess_out if line != '']
        for line in naccess_out:
            if line[0:3] in ['RES', 'HEM']:
                chain = line[8:9].strip()
                res_num = line[9:13].strip()
                ins_code = line[13:14].strip()
                chain_res_num = chain + res_num + ins_code
                int_ext_dict[chain_res_num] = float(line[14:22])

        # Selects strands in barrel
        sheet_num = sheets[0].replace('{}_sheet_'.format(domain_id), '')
        sheets_df = dssp_df[dssp_df['SHEET_NUM'] == sheet_num]
        strands = [strand for strand in set(sheets_df['STRAND_NUM'].tolist())
                   if strand != '']

        # For each strand, determines which of its residues face towards the
        # interior and which face towards the exterior (assuming that
        # consecutive residues face in opposite directions)
        for strand in strands:
            strand_df = dssp_df[dssp_df['STRAND_NUM'] == strand]
            res_acsblty = OrderedDict()
            for res_id in strand_df['RES_ID'].tolist():
                if res_id != '':
                    res_acsblty[res_id] = int_ext_dict[res_id]

            # Splits the residues into opposite faces and calculates the sum of
            # their solvent accessibilities
            face_1_res = list(res_acsblty.values())[:: 2]
            face_2_res = list(res_acsblty.values())[1:: 2]
            face_1_solv_acsblty = sum(face_1_res)
            face_2_solv_acsblty = sum(face_2_res)

            # Determines which face is exterior and which is interior (the sum
            # of the solvent accessibilities of the exterior facing residues
            # will be greater than that of the interior facing residues)
            if face_1_solv_acsblty > face_2_solv_acsblty:
                ext_chain_res_num = list(res_acsblty.keys())[:: 2]
                int_chain_res_num = list(res_acsblty.keys())[1:: 2]
            elif face_1_solv_acsblty < face_2_solv_acsblty:
                ext_chain_res_num = list(res_acsblty.keys())[1:: 2]
                int_chain_res_num = list(res_acsblty.keys())[:: 2]
            else:
                print('ERROR: solvent accessibilities of two faces of '
                      '{}_strand_{} are equal - something has gone wrong with '
                      'the analysis'.format(domain_id, strand))
                sec_struct_dfs_dict[domain_id] = None
                domain_sheets_dict[sheets[0]] = None
                unprocessed_list.append(domain_id)
                return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

            # Assigns 'interior' and 'exterior' labels to the relevant res_id
            for chain_res_num in ext_chain_res_num:
                int_ext_dict[chain_res_num] = 'exterior'
            for chain_res_num in int_chain_res_num:
                int_ext_dict[chain_res_num] = 'interior'

        # Updates dataframe with solvent accessibility information
        for row in range(dssp_df.shape[0]):
            if (dssp_df['RES_ID'][row] in list(int_ext_dict.keys())
                    and dssp_df['ATMNAME'][row] == 'CA'
                ):
                int_ext_list[row] = int_ext_dict[dssp_df['RES_ID'][row]]
        int_ext_df = pd.DataFrame({'INT_EXT': int_ext_list})
        dssp_df = pd.concat([dssp_df, int_ext_df], axis=1)
        sec_struct_dfs_dict[domain_id] = dssp_df

        return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list
    """

    def calc_core_residues_sandwich(domain_id, dssp_df, sheets,
                                    sec_struct_dfs_dict, domain_sheets_dict,
                                    unprocessed_list):
        # Calculates whether each residue faces towards the interior or the
        # exterior of the sandwich from the ratio of its solvent accessibility
        # within the sandwich to its solvent accessibility within its
        # individual parent beta-sheet

        # Checks that correct number of sheets has been retained for analysis
        if len(sheets) != 2:
            print('ERROR: incorrect number of sheets retained in {} following '
                  'solvent accessibility calculation'.format(domain_id))
            sec_struct_dfs_dict[domain_id] = None
            for sheet in sheets:
                domain_sheets_dict[sheet] = None
            unprocessed_list.append(domain_id)

            return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

        # Initialises records of interior / exterior facing residues
        print('Determining interior and exterior facing residues in {}'.format(domain_id))
        core_ext_list = ['']*dssp_df.shape[0]
        core_ext_combined = OrderedDict()
        core_ext_indv = OrderedDict()

        # Calculates solvent accessibility of both sheets in pair
        print('Calculating solvent accessible surface area for residues in '
              '{}'.format(domain_id))
        sheet_sub_strings = []
        for sheet_id in sheets:
            with open('Beta_strands/{}.pdb'.format(sheet_id), 'r') as pdb_file:
                for line in pdb_file.readlines():
                    line_start = line[0:16]
                    line_end = line[17:]
                    new_line = line_start + ' ' + line_end
                    sheet_sub_strings.append(new_line)
        sheet_string = ''.join(sheet_sub_strings)
        naccess_out = isambard.external_programs.naccess.run_naccess(
            sheet_string, 'rsa', path=False, include_hetatms=True
        )
        naccess_out = naccess_out.split('\n')
        naccess_out = [line for line in naccess_out if line != '']
        for line in naccess_out:
            if line[0:3] in ['RES', 'HEM']:
                chain = line[8:9].strip()
                res_num = line[9:13].strip()
                ins_code = line[13:14].strip()
                chain_res_num = chain + res_num + ins_code
                core_ext_combined[chain_res_num] = float(line[14:22])

        # Calculates solvent accessibility of sheets in pair individually
        for sheet_id in sheets:
            sheet_sub_strings = []
            print('Calculating solvent accessible surface area for residues in '
                  '{}'.format(sheet_id))
            with open('Beta_strands/{}.pdb'.format(sheet_id), 'r') as pdb_file:
                for line in pdb_file.readlines():
                    line_start = line[0:16]
                    line_end = line[17:]
                    new_line = line_start + ' ' + line_end
                    sheet_sub_strings.append(new_line)
            sheet_string = ''.join(sheet_sub_strings)
            naccess_out = isambard.external_programs.naccess.run_naccess(
                sheet_string, 'rsa', path=False, include_hetatms=True
            )
            naccess_out = naccess_out.split('\n')
            naccess_out = [line for line in naccess_out if line != '']
            for line in naccess_out:
                if line[0:3] in ['RES', 'HEM']:
                    chain = line[8:9].strip()
                    res_num = line[9:13].strip()
                    ins_code = line[13:14].strip()
                    chain_res_num = chain + res_num + ins_code
                    core_ext_indv[chain_res_num] = float(line[14:22])

        # Determines solvent accessibility of each residue in the beta-sheet
        # assembly as compared to its individual parent beta-sheet. If the
        # solvent accessibility is lower in the assembly, the residue is
        # assigned as 'interior', if it is the same it is assigned as
        # 'exterior', otherwise an error is thrown.
        for res in list(core_ext_combined.keys()):
            solv_acsblty_combined = core_ext_combined[res]
            solv_acsblty_indv = core_ext_indv[res]

            if solv_acsblty_indv - solv_acsblty_combined > 1:
                core_ext_combined[res] = 'core'
            elif solv_acsblty_indv - solv_acsblty_combined < 1:
                core_ext_combined[res] = 'external'

        # Updates dataframe with solvent accessibility information
        for row in range(dssp_df.shape[0]):
            if (dssp_df['RES_ID'][row] in list(core_ext_combined.keys())
                    and dssp_df['ATMNAME'][row] == 'CA'
                    ):
                core_ext_list[row] = core_ext_combined[dssp_df['RES_ID'][row]]
        core_ext_df = pd.DataFrame({'CORE_OR_EXT': core_ext_list})
        dssp_df = pd.concat([dssp_df, core_ext_df], axis=1)
        sec_struct_dfs_dict[domain_id] = dssp_df

        return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list


class calculate_solvent_accessible_surface_area(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    def run_naccess(self, sec_struct_dfs_dict, domain_sheets_dict):
        # Pipeline function for running naccess to calculate solvent accessible
        # surface area
        if __name__ == 'subroutines.naccess':
            from subroutines.naccess import naccess_solv_acsblty_calcs
        else:
            from datagen.subroutines.naccess import naccess_solv_acsblty_calcs

        unprocessed_list = []
        for domain_id in list(sec_struct_dfs_dict.keys()):
            dssp_df = sec_struct_dfs_dict[domain_id]
            sheets = {key: value for key, value in domain_sheets_dict.items()
                      if domain_id in key}

            if not sheets:
                print('ERROR: No beta-sheets retained for {}'.format(domain_id))
                unprocessed_list.append(domain_id)

            if self.code[0:4] in ['2.40']:
                solv_acsblty_dict = naccess_solv_acsblty_calcs.calculate_barrel_solv_acsblty(
                    sheets
                )
                res_solv_acsblty = naccess_solv_acsblty_calcs.calculate_residue_solv_acsblty(
                    domain_id
                )
            elif self.code[0:4] in ['2.60']:
                solv_acsblty_dict = naccess_solv_acsblty_calcs.calculate_sandwich_solv_acsblty(
                    sheets
                )
                res_solv_acsblty = naccess_solv_acsblty_calcs.calculate_residue_solv_acsblty(
                    domain_id
                )
            else:
                solv_acsblty_dict = {}
                res_solv_acsblty = {}

            (sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list
             ) = naccess_solv_acsblty_calcs.add_naccess_info_to_df(
                domain_id, dssp_df, sheets, self.code, solv_acsblty_dict,
                res_solv_acsblty, sec_struct_dfs_dict, domain_sheets_dict,
                unprocessed_list
            )

        sec_struct_dfs_dict = {key: value for key, value in
                               sec_struct_dfs_dict.items() if value is not None}
        domain_sheets_dict = {key: value for key, value in
                              domain_sheets_dict.items() if value is not None}

        with open('Unprocessed_domains.txt', 'a') as unprocessed_file:
            unprocessed_file.write('\n\nError in solvent accessibility calculation:\n')
            for domain_id in set(unprocessed_list):
                unprocessed_file.write('{}\n'.format(domain_id))

        return sec_struct_dfs_dict, domain_sheets_dict

    def identify_int_ext(self, sec_struct_dfs_dict, domain_sheets_dict):
        # Pipeline script to identify interior and exterior-facing residues in
        # barrels / sandwiches
        unprocessed_list = []

        for domain_id in list(sec_struct_dfs_dict.keys()):
            dssp_df = sec_struct_dfs_dict[domain_id]
            sheets = [key for key in list(domain_sheets_dict.keys()) if
                      domain_id in key]

            if self.code[0:4] in ['2.40']:
                (sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list
                 ) = naccess_solv_acsblty_calcs.calculate_int_ext_barrel(
                    domain_id, dssp_df, sheets, sec_struct_dfs_dict,
                    domain_sheets_dict, unprocessed_list
                )
            elif self.code[0:4] in ['2.60']:
                (sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list
                 ) = naccess_solv_acsblty_calcs.calculate_int_ext_sandwich(
                    domain_id, dssp_df, sheets, sec_struct_dfs_dict,
                    domain_sheets_dict, unprocessed_list
                )

        sec_struct_dfs_dict = {key: value for key, value in
                               sec_struct_dfs_dict.items() if value is not None}
        domain_sheets_dict = {key: value for key, value in
                              domain_sheets_dict.items() if value is not None}

        with open('Unprocessed_domains.txt', 'a') as unprocessed_file:
            unprocessed_file.write('\n\nError in determination of interior and '
                                   'exterior facing residues:\n')
            for domain_id in set(unprocessed_list):
                unprocessed_file.write('{}\n'.format(domain_id))

        return sec_struct_dfs_dict
