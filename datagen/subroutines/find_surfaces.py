
import sys
import math
import itertools
import isambard
import pandas as pd
import numpy as np
import networkx as nx
from collections import OrderedDict
if __name__ == 'subroutines.find_surfaces':
    from subroutines.run_stages import run_stages
else:
    from datagen.subroutines.run_stages import run_stages


class interior_exterior_calcs():

    def find_strands_network(G):
        # Finds network of interacting neighbouring strands
        nodes_dict = {}
        strands = nx.nodes(G)

        for strand in strands:
            nodes_dict[strand] = len(G.neighbors(strand))

        # Removes edge strands to find a closed circle of interacting strands
        # (note that this approach assumes only one closed circle is present in
        # the input network)
        while min(list(nodes_dict.values())) < 2:
            for strand in strands:
                if len(G.neighbors(strand)) < 2:
                    G.remove_node(strand)
                    del nodes_dict[strand]

            strands = nx.nodes(G)
            for strand in strands:
                nodes_dict[strand] = len(G.neighbors(strand))

        # Makes list of strands in closed circle network
        strands = nx.cycle_basis(G)[0]  # Finds first complete cycle

        return strands

    def find_z_axis(strands, sheets_df):
        # Finds axis through the barrel pore as the line that passes through
        # the average xyz coordinates of the middle two residues of each
        # strand.
        count = 0
        coords_res_1 = []
        coords_res_n = []

        # Determines the indices of the central two residues in each strand.
        for strand in strands:
            count += 1

            strand_df = sheets_df[sheets_df['STRAND_NUM'] == strand]
            strand_df = strand_df.reset_index(drop=True)

            row_num = strand_df.shape[0]
            # If there is an odd number of strands, the central residue plus
            # its C-terminal residue are selected as the two central residues
            if row_num % 2 == 1:
                row_num += 1

            index_1 = (row_num / 2) - 1
            index_n = row_num / 2
            if count % 2 == 1:  # Ensures that residues closer to the
                # periplasmic side of the membrane are grouped together,
                # likewise for residues closer to the extracellular side of the
                # membrane. Note that this assumes that neighbouring strands
                # interact with one another in an antiparallel hydrogen bonding
                # arrangement (however, the code will only break if the
                # majority of strands interact with a parallel instead of a
                # parallel hydrogen bonding arrangement, which is highly unlikely)
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

        # Calculates average xyz coordinates
        coords_res_1 = np.array(coords_res_1)
        coords_res_n = np.array(coords_res_n)
        x_coord_1 = np.sum(coords_res_1[:, 0]) / coords_res_1.shape[0]
        y_coord_1 = np.sum(coords_res_1[:, 1]) / coords_res_1.shape[0]
        z_coord_1 = np.sum(coords_res_1[:, 2]) / coords_res_1.shape[0]
        x_coord_n = np.sum(coords_res_n[:, 0]) / coords_res_n.shape[0]
        y_coord_n = np.sum(coords_res_n[:, 1]) / coords_res_n.shape[0]
        z_coord_n = np.sum(coords_res_n[:, 2]) / coords_res_n.shape[0]
        xyz_coords = [x_coord_1, y_coord_1, z_coord_1, x_coord_n, y_coord_n,
                      z_coord_n]

        return xyz_coords

    def rotate_translate_barrel(domain_id, xyz_coords, sheets):
        # Aligns the barrel to z = 0 using the reference axis through the
        # barrel pore calculated in the previous step
        print('Aligning {} barrel with z = 0'.format(domain_id))

        # Creates ampal object from barrel
        barrel = isambard.ampal.convert_pdb_to_ampal(
            'Beta_strands/{}.pdb'.format(list(sheets.keys())[0])
        )

        # Aligns the barrel with z = 0
        s1 = [xyz_coords[0], xyz_coords[1], xyz_coords[2]]
        e1 = [xyz_coords[3], xyz_coords[4], xyz_coords[5]]
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

        return xy_dict

    def calc_average_coordinates(domain_id, xy_dict, sheets_df,
                                 unprocessed_list):
        # Calculates average xy coordinates of the barrel C_alpha atoms
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
            unprocessed_list.append(domain_id)

        return com, unprocessed_list

    def calc_int_ext(domain_id, sheets_df, xy_dict, com, int_ext_dict):
        # Calculates whether a residue faces towards the interior or the
        # exterior of the barrel
        print('Identifying interior and exterior residues in {}'.format(domain_id))

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
                c_alpha_beta_magnitude = math.sqrt(((c_alpha_beta_vector[0][0])**2)
                                                   + ((c_alpha_beta_vector[1][0])**2))

                c_alpha_com_vector = np.array([[com[0][0]-c_alpha_x],
                                               [com[1][0]-c_alpha_y]])
                c_alpha_com_magnitude = math.sqrt(((c_alpha_com_vector[0][0])**2)
                                                  + ((c_alpha_com_vector[1][0])**2))

                c_alpha_numerator = np.sum((c_alpha_beta_vector*c_alpha_com_vector), axis=0)[0]
                c_alpha_denominator = c_alpha_beta_magnitude*c_alpha_com_magnitude

                cos_c_alpha_angle = (c_alpha_numerator / c_alpha_denominator)
                c_alpha_angle = math.degrees(math.acos(cos_c_alpha_angle))

                if c_alpha_angle < 90:
                    int_ext_dict[res_id] = 'interior'
                elif c_alpha_angle > 90:
                    int_ext_dict[res_id] = 'exterior'

        return int_ext_dict


class barrel_int_ext():

    def int_ext_pipeline(domain_id, dssp_df, sheets, sec_struct_dfs_dict,
                         domain_sheets_dict, unprocessed_list):
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
        sheet_num = list(sheets.keys())[0].replace('{}_sheet_'.format(domain_id), '')
        sheets_df = dssp_df[dssp_df['SHEET_NUM'] == sheet_num]
        sheets_df = sheets_df.reset_index(drop=True)

        # Finds network of interacting neighbouring strands
        networks = list(sheets.values())
        G = networks[0]
        strands = interior_exterior_calcs.find_strands_network(G)

        # Finds axis through barrel pore
        xyz_coords = interior_exterior_calcs.find_z_axis(strands, sheets_df)

        # Aligns barrel with z = 0
        xy_dict = interior_exterior_calcs.rotate_translate_barrel(
            domain_id, xyz_coords, sheets
        )

        # Initialises records of interior / exterior facing residues
        int_ext_list = ['']*dssp_df.shape[0]
        int_ext_dict = OrderedDict()

        # Calculates average xy coordinates of the barrel C_alpha atoms
        com, unprocessed_list = interior_exterior_calcs.calc_average_coordinates(
            domain_id, xy_dict, sheets_df, unprocessed_list
        )

        int_ext_dict = interior_exterior_calcs.calc_int_ext(
            domain_id, sheets_df, xy_dict, com, int_ext_dict
        )

        # Updates dataframe with solvent accessibility information
        for row in range(dssp_df.shape[0]):
            res_id = dssp_df['RES_ID'][row]
            if (res_id in list(int_ext_dict.keys())
                    and dssp_df['ATMNAME'][row] == 'CA'
                    ):
                int_ext_list[row] = int_ext_dict[res_id]
        int_ext_df = pd.DataFrame({'INT_EXT': int_ext_list})
        dssp_df = pd.concat([dssp_df, int_ext_df], axis=1)
        sec_struct_dfs_dict[domain_id] = dssp_df

        return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list


class sandwich_int_ext():

    def calculate_int_ext_sandwich(domain_id, dssp_df, sheets,
                                   sec_struct_dfs_dict, domain_sheets_dict,
                                   unprocessed_list):
        if len(sheets) != 2:
            print('ERROR: incorrect number of sheets retained in {} following '
                  'solvent accessibility calculation'.format(domain_id))
            sec_struct_dfs_dict[domain_id] = None
            for sheet in sheets:
                domain_sheets_dict[sheet] = None
            unprocessed_list.append(domain_id)

            return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

        # Creates an ampal object of the beta-sandwich
        sandwich = isambard.ampal.convert_pdb_to_ampal('Beta_strands/{}.pdb'.format(domain_id))

        # Identifies pairs of edge strands
        networks = [network for key, network in domain_sheets_dict.items()
                    if domain_id in key]
        edges = []
        for G in networks:
            strands = nx.nodes(G)
            for strand in strands:
                if len(G.neighbors(strand)) == 1:
                    edges.append(strand)

        # Finds centre of mass of each edge strand
        edge_com_dict = {}
        for strand in edges:
            strand_df = dssp_df[dssp_df['STRAND_NUM'] == strand]
            x_coords = [num for num in strand_df['XPOS'].tolist() if num != '']
            x_coord = sum(x_coords) / len(x_coords)
            y_coords = [num for num in strand_df['YPOS'].tolist() if num != '']
            y_coord = sum(y_coords) / len(y_coords)
            z_coords = [num for num in strand_df['ZPOS'].tolist() if num != '']
            z_coord = sum(z_coords) / len(z_coords)
            edge_com_dict[strand] = np.array([[x_coord], [y_coord], [z_coord]])

        # Uses the distance between strand centres of mass to determine which
        # form the same edge of the sandwich
        G = networks[0]
        for H in networks[1:]:
            G = nx.compose(G, H)

        strand_1 = edges[0]
        dist_dict = {}
        for strand in edges[1:]:
            if strand_1 not in G.neighbors(strand):
                dist = math.sqrt(((edge_com_dict[strand][0][0] -
                                   edge_com_dict[strand_1][0][0])**2)
                                 + ((edge_com_dict[strand][1][0] -
                                     edge_com_dict[strand_1][1][0])**2)
                                 + ((edge_com_dict[strand][2][0] -
                                     edge_com_dict[strand_1][2][0])**2))
                dist_dict[dist] = (strand, strand_1)
        min_dist = min(list(dist_dict.keys()))
        strand_pair = dist_dict[min_dist]

        # Calculates centre of mass of C_alpha atoms across the two beta-sheets
        sheet_ids = [sheet.replace('{}_sheet_'.format(domain_id), '') for sheet
                     in sheets]
        sheets_df = dssp_df[dssp_df['SHEET_NUM'].isin(sheet_ids)]
        sheets_df = sheets_df.reset_index(drop=True)

        sheets_xyz_coords = []
        for row in range(sheets_df.shape[0]):
            x_coord = sheets_df['XPOS'][row]
            y_coord = sheets_df['YPOS'][row]
            z_coord = sheets_df['ZPOS'][row]
            sheets_xyz_coords.append((x_coord, y_coord, z_coord))

        sheets_xyz_coords = np.array(sheets_xyz_coords)
        x_coord = np.sum(sheets_xyz_coords[:, 0]) / sheets_xyz_coords.shape[0]
        y_coord = np.sum(sheets_xyz_coords[:, 1]) / sheets_xyz_coords.shape[0]
        z_coord = np.sum(sheets_xyz_coords[:, 2]) / sheets_xyz_coords.shape[0]

        com_all_coords = np.array([[x_coord],
                                   [y_coord],
                                   [z_coord]])

        # Calculates centre of mass of C_alpha atoms across the two beta-sheets
        # excluding atoms in the edge strand pair identified earlier
        sheets_xyz_coords = []
        for row in range(sheets_df.shape[0]):
            if sheets_df['STRAND_NUM'][row] not in strand_pair:
                x_coord = sheets_df['XPOS'][row]
                y_coord = sheets_df['YPOS'][row]
                z_coord = sheets_df['ZPOS'][row]
                sheets_xyz_coords.append((x_coord, y_coord, z_coord))

        sheets_xyz_coords = np.array(sheets_xyz_coords)
        x_coord = np.sum(sheets_xyz_coords[:, 0]) / sheets_xyz_coords.shape[0]
        y_coord = np.sum(sheets_xyz_coords[:, 1]) / sheets_xyz_coords.shape[0]
        z_coord = np.sum(sheets_xyz_coords[:, 2]) / sheets_xyz_coords.shape[0]

        com_sub_coords = np.array([[x_coord],
                                   [y_coord],
                                   [z_coord]])

        # Finds reference axis between the two centres of mass and aligns the
        # sandwich along this axis
        x_coord_1 = np.sum(com_all_coords[0][0])
        y_coord_1 = np.sum(com_all_coords[1][0])
        z_coord_1 = np.sum(com_all_coords[2][0])
        x_coord_n = np.sum(com_sub_coords[0][0])
        y_coord_n = np.sum(com_sub_coords[1][0])
        z_coord_n = np.sum(com_sub_coords[2][0])

        # Aligns the barrel with z = 0
        s1 = [x_coord_1, y_coord_1, z_coord_1]
        e1 = [x_coord_n, y_coord_n, z_coord_n]
        s2 = [0.0, 0.0, 0.0]
        e2 = [0.0, 0.0, 1.0]
        translation, angle, axis, point = isambard.geometry.find_transformations(
            s1, e1, s2, e2
        )
        sandwich.rotate(angle, axis, point=point)  # Rotation must be performed
        # before translation
        sandwich.translate(translation)

        # Generates dictionary of residues and their new xyz coordinates
        sandwich_pdb_string = sandwich.make_pdb().split('\n')
        xy_dict = OrderedDict()
        for line in sandwich_pdb_string:
            if line[0:6].strip() in ['ATOM', 'HETATM']:
                atom_id = line[21:27].replace(' ', '') + '_' + line[12:16].strip()
                line_segments = line.split('.')
                x_coord_atom = float(line_segments[0].split()[-1] + '.' + line_segments[1][0:3])
                y_coord_atom = float(line_segments[1][3:] + '.' + line_segments[2][0:3])
                xy_coords = np.array([[x_coord_atom],
                                      [y_coord_atom]])
                xy_dict[atom_id] = xy_coords

        # Initialises records of interior / exterior facing residues
        print('Determining interior and exterior residues in {}'.format(domain_id))
        int_ext_list = ['']*dssp_df.shape[0]
        int_ext_dict = OrderedDict()

        # Determines whether a residue is interior- or exterior-facing from the
        # angle between the centre of mass, C_alpha and C_beta
        com_coords = np.array([[np.sum(list(xy_dict.values())[0][0])],
                               [np.sum(list(xy_dict.values())[1][0])]])

        for res_id in sheets_df['RES_ID'].tolist():
            print(domain_id, res_id)
            if 'GLY' not in res_id:
                res_df = dssp_df[dssp_df['RES_ID'] == res_id]
                res_atoms = res_df['ATMNAME'].tolist()

                if all([x in res_atoms for x in ['N', 'CA', 'C', 'O', 'CB']]):
                    n_coords = xy_dict['{}_N'.format(res_id)]
                    c_coords = xy_dict['{}_C'.format(res_id)]
                    ca_coords = xy_dict['{}_CA'.format(res_id)]
                    cb_coords = xy_dict['{}_CB'.format(res_id)]

                    n_com_vector = np.array([[com_coords[0][0] - n_coords[0][0]],
                                             [com_coords[1][0] - n_coords[1][0]]])
                    n_com_magnitude = math.sqrt(((n_com_vector[0][0])**2)
                                                + ((n_com_vector[1][0])**2))

                    c_com_vector = np.array([[com_coords[0][0] - c_coords[0][0]],
                                             [com_coords[1][0] - c_coords[1][0]]])
                    c_com_magnitude = math.sqrt(((c_com_vector[0][0])**2)
                                                + ((c_com_vector[1][0])**2))

                    ca_com_vector = np.array([[com_coords[0][0] - ca_coords[0][0]],
                                              [com_coords[1][0] - ca_coords[1][0]]])
                    ca_com_magnitude = math.sqrt(((ca_com_vector[0][0])**2)
                                                 + ((ca_com_vector[1][0])**2))

                    n_ca_vector = np.array([[ca_coords[0][0] - n_coords[0][0]],
                                            [ca_coords[1][0] - n_coords[1][0]]])
                    n_ca_magnitude = math.sqrt(((n_ca_vector[0][0])**2)
                                               + ((n_ca_vector[1][0])**2))

                    c_ca_vector = np.array([[ca_coords[0][0] - c_coords[0][0]],
                                            [ca_coords[1][0] - c_coords[1][0]]])
                    c_ca_magnitude = math.sqrt(((c_ca_vector[0][0])**2)
                                               + ((c_ca_vector[1][0])**2))

                    ca_cb_vector = np.array([[cb_coords[0][0] - ca_coords[0][0]],
                                             [cb_coords[1][0] - ca_coords[1][0]]])
                    ca_cb_magnitude = math.sqrt(((ca_cb_vector[0][0])**2)
                                                + ((ca_cb_vector[1][0])**2))

                    # Calculates com-N-C_alpha angle
                    com_n_ca_numerator = np.sum((n_com_vector*n_ca_vector), axis=0)[0]
                    com_n_ca_denominator = n_com_magnitude*n_ca_magnitude
                    com_n_ca_angle = math.degrees(
                        math.acos(com_n_ca_numerator / com_n_ca_denominator)
                    )

                    # Caculates com-C-C_alpha angle
                    com_c_ca_numerator = np.sum((c_com_vector*c_ca_vector), axis=0)[0]
                    com_c_ca_denominator = c_com_magnitude*c_ca_magnitude
                    com_c_ca_angle = math.degrees(
                        math.acos(com_c_ca_numerator / com_c_ca_denominator)
                    )

                    # Selects the largest angle of the two backbone angles
                    if n_com_magnitude < c_com_magnitude:
                        backbone_angle = com_c_ca_angle

                    # Calculates com-C_alpha-C_beta angle
                    com_ca_cb_numerator = np.sum((ca_com_vector*ca_cb_vector), axis=0)[0]
                    com_ca_cb_denominator = ca_com_magnitude*ca_cb_magnitude
                    com_ca_cb_angle = math.degrees(
                        math.acos(com_ca_cb_numerator / com_ca_cb_denominator)
                    )

                    # Determines whether the residue is interior- or
                    # exterior-facing
                    if com_ca_cb_angle < 90:
                        int_ext_dict[res_id] = 'interior'
                    elif com_ca_cb_angle > 90:
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


class find_interior_exterior_surfaces(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    def identify_int_ext(self, sec_struct_dfs_dict, domain_sheets_dict):
        # Pipeline script to identify interior and exterior-facing residues in
        # barrels / sandwiches
        unprocessed_list = []

        for domain_id in list(sec_struct_dfs_dict.keys()):
            dssp_df = sec_struct_dfs_dict[domain_id]
            sheets = {key: value for key, value in domain_sheets_dict.items()
                      if domain_id in key}

            if self.code[0:4] in ['2.40']:
                (sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list
                 ) = barrel_int_ext.int_ext_pipeline(
                    domain_id, dssp_df, sheets, sec_struct_dfs_dict,
                    domain_sheets_dict, unprocessed_list
                )
            elif self.code[0:4] in ['2.60']:
                (sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list
                 ) = barrel_int_ext.int_ext_pipeline(
                    domain_id, dssp_df, sheets, sec_struct_dfs_dict,
                    domain_sheets_dict, unprocessed_list
                )

        sec_struct_dfs_dict = OrderedDict(
            {key: value for key, value in sec_struct_dfs_dict.items() if value is not None}
        )
        domain_sheets_dict = OrderedDict(
            {key: value for key, value in domain_sheets_dict.items() if value is not None}
        )

        with open('Unprocessed_domains.txt', 'a') as unprocessed_file:
            unprocessed_file.write('\n\nError in determination of interior and '
                                   'exterior facing residues:\n')
            for domain_id in set(unprocessed_list):
                unprocessed_file.write('{}\n'.format(domain_id))

        return sec_struct_dfs_dict, domain_sheets_dict
