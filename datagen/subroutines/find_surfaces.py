
import sys
import math
import copy
import itertools
import isambard
import pandas as pd
import numpy as np
import networkx as nx
from collections import OrderedDict
from shapely.geometry import Polygon, Point

if __name__ == 'subroutines.find_surfaces':
    from subroutines.run_stages import run_stages
else:
    from datagen.subroutines.run_stages import run_stages


class barrel_interior_exterior_calcs():

    def find_strands_network(G, domain_id, unprocessed_list):
        # Finds network of interacting neighbouring strands (used to remove
        # additional strands that don't form part of the main barrel (and so
        # which may disrupt z-axis alignment in subsequent steps))
        nodes_dict = {}
        strands = list(G.nodes())
        orig_strands = copy.copy(strands)

        for strand in strands:
            nodes_dict[strand] = len(list(G.neighbors(strand)))

        # If no closed circle of interacting strands, takes all strands
        # forwards for further analysis
        if len(nx.cycle_basis(G)) == 0:
            unprocessed_list.append(domain_id)
            return orig_strands, unprocessed_list

        # Removes edge strands to find a closed circle of interacting strands
        # (note that this approach assumes only one closed circle is present in
        # the input network)
        while min(list(nodes_dict.values())) < 2:
            for strand in strands:
                if len(list(G.neighbors(strand))) < 2:
                    G.remove_node(strand)
                    del nodes_dict[strand]

            strands = list(G.nodes())
            if len(strands) == 0:
                break
            else:
                for strand in strands:
                    nodes_dict[strand] = len(list(G.neighbors(strand)))

        # Makes list of strands in closed circle network
        strand_cycles = nx.cycle_basis(G)
        strands = max(strand_cycles, key=len)  # Finds longest complete cycle

        return strands, unprocessed_list

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
            'Beta_strands/{}.pdb'.format(domain_id)
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

    def calc_average_coordinates(domain_id, xy_dict, sheets_df):
        # Calculates average xy coordinates of the barrel C_alpha atoms
        x_sum = 0
        y_sum = 0
        count = 0

        for res_id in list(xy_dict.keys()):
            if res_id.split('_')[0] in sheets_df['RES_ID'].tolist():
                x_sum += xy_dict[res_id][0][0]
                y_sum += xy_dict[res_id][1][0]
                count += 1

        x_average = x_sum / count
        y_average = y_sum / count
        com = np.array([[x_average],
                        [y_average]])

        return com

    def calc_int_ext(domain_id, sheets_df, xy_dict, com, int_ext_dict):
        # Calculates whether a residue faces towards the interior or the
        # exterior of the barrel
        print('Identifying interior and exterior residues in {}'.format(domain_id))

        res_dict = OrderedDict()
        for index, res_id in enumerate(sheets_df['RES_ID'].tolist()):
            res_dict[res_id] = sheets_df['RESNAME'].tolist()[index]

        for res_id in list(res_dict.keys()):
            if (
                (not 'GLY' in res_dict[res_id])
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

                if c_alpha_angle <= 90:
                    int_ext_dict[res_id] = 'interior'
                elif c_alpha_angle > 90:
                    int_ext_dict[res_id] = 'exterior'

        return int_ext_dict


class sandwich_interior_exterior_calcs():

    def find_z_axis(sheets_df):
        # Finds axis between the two sheets of the beta-sandwiches

        # Extracts hydrogen bonded residues
        h_bonds_dict = {}
        for row in range(sheets_df.shape[0]):
            if sheets_df['SHEET?'][row] == 'E':
                res_1 = sheets_df['DSSP_NUM'][row]
                res_2 = sheets_df['H-BONDS'][row][0]
                res_3 = sheets_df['H-BONDS'][row][1]

                if res_2 != '0':
                    pair_1 = [str(min([int(res_1), int(res_2)])),
                              str(max([int(res_1), int(res_2)]))]
                    if pair_1[0] not in list(h_bonds_dict.keys()):
                        h_bonds_dict[pair_1[0]] = pair_1[1]
                if res_3 != '0':
                    pair_2 = [str(min([int(res_1), int(res_3)])),
                              str(max([int(res_1), int(res_3)]))]
                    if pair_2[0] not in list(h_bonds_dict.keys()):
                        h_bonds_dict[pair_2[0]] = pair_2[1]

        # Calculates vectors between C_alpha atoms of hydrogen bonded residues
        vectors_array = np.empty((len(list(h_bonds_dict.keys())), 3))
        for index, res_1 in enumerate(h_bonds_dict.keys()):
            res_2 = h_bonds_dict[res_1]
            xyz_1 = np.array([])
            xyz_2 = np.array([])

            for row in range(sheets_df.shape[0]):
                if sheets_df['DSSP_NUM'][row] == res_1:
                    x_1 = sheets_df['XPOS'][row]
                    y_1 = sheets_df['YPOS'][row]
                    z_1 = sheets_df['ZPOS'][row]
                    xyz_1 = np.array([[x_1],
                                      [y_1],
                                      [z_1]])
                elif sheets_df['DSSP_NUM'][row] == res_2:
                    x_2 = sheets_df['XPOS'][row]
                    y_2 = sheets_df['YPOS'][row]
                    z_2 = sheets_df['ZPOS'][row]
                    xyz_2 = np.array([[x_2],
                                      [y_2],
                                      [z_2]])
                    break

            vectors_array[index][0] = xyz_2[0][0] - xyz_1[0][0]
            vectors_array[index][1] = xyz_2[1][0] - xyz_1[1][0]
            vectors_array[index][2] = xyz_2[2][0] - xyz_1[2][0]

        vector = np.mean(vectors_array, axis=0)

        return vector

    def align_sandwich(domain_id, sheets, vector):
        # Aligns the sandwich to z = 0 using the reference axis through the
        # sandwich core calculated in the previous step
        print('Aligning {} sandwich with z = 0'.format(domain_id))

        # Creates ampal object from sandwich
        sandwich = isambard.ampal.convert_pdb_to_ampal(
            'Beta_strands/{}.pdb'.format(domain_id)
        )

        # Aligns the sandwich with z = 0
        s1 = [0.0, 0.0, 0.0]
        e1 = [vector[0], vector[1], vector[2]]
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

        x_sum = 0
        y_sum = 0
        count = 0
        for res in list(xy_dict.keys()):
            if 'CA' in res:
                x_sum += xy_dict[res][0][0]
                y_sum += xy_dict[res][1][0]
                count += 1
        x_avrg = x_sum / count
        y_avrg = y_sum / count

        with open('my_sandwich.pdb', 'w') as pdb_file:
            for line in sandwich_pdb_string:
                pdb_file.write('{}\n'.format(line))
            pdb_file.write('HETATM    1  N   DUM     1       0.000  0.000   0.000\n')
            pdb_file.write('HETATM    2  O   DUM     2       0.000  0.000   1.000\n')
            pdb_file.write('HETATM    3  C   DUM     3       {}{}0.000\n'.format(
                str(round(x_avrg, 3)).rjust(8), str(round(y_avrg, 3)).rjust(8)))

        return xy_dict

    def calc_int_ext(domain_id, dssp_df, xy_dict):
        # Determines whether a residue is interior- or exterior-facing from
        # whether its C_beta atom lies on the same side of the plane between
        # its N and C atoms as the centre of mass of the sandwich or not
        # TODO: Update from 2D to 3D (so that alignment along the z-axis is no
        # longer necessary)
        print('Determining interior and exterior residues in {}'.format(domain_id))

        int_ext_dict = {}

        x_sum = 0
        y_sum = 0
        count = 0
        for res in list(xy_dict.keys()):
            if 'CA' in res:
                x_sum += xy_dict[res][0][0]
                y_sum += xy_dict[res][1][0]
                count += 1
        x_avrg = x_sum / count
        y_avrg = y_sum / count

        res_list = list(set(dssp_df['RES_ID'].tolist()))
        for res in res_list:
            n = False
            c = False
            cb = False

            xy_sub_dict = {key: value for key, value in xy_dict.items() if res
                           in key}
            for atom in list(xy_sub_dict.keys()):
                if atom.endswith('_N'):
                    n = True
                    n_dist = math.sqrt(((xy_sub_dict[atom][0][0] - x_avrg)**2)
                                       + ((xy_sub_dict[atom][1][0] - y_avrg)**2))
                    n_x_coord = xy_sub_dict[atom][0][0]
                    n_y_coord = xy_sub_dict[atom][1][0]
                elif atom.endswith('_C'):
                    c = True
                    c_dist = math.sqrt(((xy_sub_dict[atom][0][0] - x_avrg)**2)
                                       + ((xy_sub_dict[atom][1][0] - y_avrg)**2))
                    c_x_coord = xy_sub_dict[atom][0][0]
                    c_y_coord = xy_sub_dict[atom][1][0]
                elif atom.endswith('_CB'):
                    cb = True
                    cb_x_coord = xy_sub_dict[atom][0][0]
                    cb_y_coord = xy_sub_dict[atom][1][0]

            distance_com = (((x_avrg-n_x_coord)*(c_y_coord-n_y_coord)) -
                            ((y_avrg-n_y_coord)*(c_x_coord-n_x_coord)))
            distance_cb = (((cb_x_coord-n_x_coord)*(c_y_coord-n_y_coord)) -
                           ((cb_y_coord-n_y_coord)*(c_x_coord-n_x_coord)))

            if all(x is True for x in [n, c, cb]):
            if (distance_com < 0 and distance_cb < 0) or (distance_com > 0 and distance_cb > 0):
                int_ext_dict[res] = 'interior'
            else:
                int_ext_dict[res] = 'exterior'

        return int_ext_dict


class int_ext_pipeline():

    def barrel_pipeline(domain_id, dssp_df, sheets, sec_struct_dfs_dict,
                        domain_sheets_dict, unprocessed_list):
        # Determines which face of the beta-sheet forms the interior of the
        # barrel and which forms the exterior based upon whether the angle
        # between the centre of mass of the C_alpha atoms of the barrel and the
        # residue in question is less or greater than 90 degrees. The barrel is
        # collapsed along the central axis that runs through its pore in order
        # to avoid misclassification of residues in frayed edges of the barrel.

        # Checks that correct number of sheets has been retained for analysis
        if len(list(sheets.keys())) != 1:
            print('ERROR: more than 1 sheet retained in {} following solvent '
                  'accessibility calculation'.format(domain_id))
            sec_struct_dfs_dict[domain_id] = None
            for sheet in list(sheets.keys()):
                domain_sheets_dict[sheet] = None
            unprocessed_list.append(domain_id)

            return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

        # Creates dataframe of CA atoms in barrel
        sheet_num = list(sheets.keys())[0].replace('{}_sheet_'.format(domain_id), '')
        sheets_df = dssp_df[dssp_df['SHEET_NUM'] == sheet_num]
        sheets_df = sheets_df.reset_index(drop=True)

        # Finds network of interacting neighbouring strands
        G = list(sheets.values())[0]
        strands, unprocessed_list = barrel_interior_exterior_calcs.find_strands_network(
            G, domain_id, unprocessed_list
        )

        # Checks complete circle of hydrogen bonded strands is present in domain
        if domain_id in unprocessed_list:
            print('ERROR: domain {} does not contain closed circular network '
                  'of hydrogen-bonded strands'.format(domain_id))
            sec_struct_dfs_dict[domain_id] = None
            for sheet in list(sheets.keys()):
                domain_sheets_dict[sheet] = None
            unprocessed_list.append(domain_id)

            return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

        # Finds axis through barrel pore
        xyz_coords = barrel_interior_exterior_calcs.find_z_axis(strands, sheets_df)

        # Aligns barrel with z = 0
        xy_dict = barrel_interior_exterior_calcs.rotate_translate_barrel(
            domain_id, xyz_coords, sheets
        )

        # Initialises records of interior / exterior facing residues
        int_ext_list = ['']*dssp_df.shape[0]
        int_ext_dict = OrderedDict()

        # Calculates average xy coordinates of the barrel C_alpha atoms
        com = barrel_interior_exterior_calcs.calc_average_coordinates(
            domain_id, xy_dict, sheets_df
        )

        # Calculates interior- and exterior-facing residues
        int_ext_dict = barrel_interior_exterior_calcs.calc_int_ext(
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

    def sandwich_pipeline(domain_id, dssp_df, sheets, sec_struct_dfs_dict,
                          domain_sheets_dict, unprocessed_list):
        # Calculates whether a residue faces towards the interior or the
        # exterior of the sandwich

        # Checks that correct number of sheets have been retained in the
        # beta-sandwich
        networks = [network for key, network in domain_sheets_dict.items()
                    if domain_id in key]

        if len(list(sheets.keys())) != 2 or len(networks) != 2:
            print('ERROR: incorrect number of sheets retained in {} following '
                  'solvent accessibility calculation'.format(domain_id))
            sec_struct_dfs_dict[domain_id] = None
            for sheet in list(sheets.keys()):
                domain_sheets_dict[sheet] = None
            unprocessed_list.append(domain_id)

            return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

        # Calculates z-axis between the two beta-sheets
        sheet_ids = [sheet.replace('{}_sheet_'.format(domain_id), '') for sheet
                     in list(sheets.keys()) if sheet != '']
        sheets_df = dssp_df[dssp_df['SHEET_NUM'].isin(sheet_ids)]
        sheets_df = sheets_df.reset_index(drop=True)
        vector = sandwich_interior_exterior_calcs.find_z_axis(sheets_df)

        # Aligns the beta-sandwich with z = 0
        xy_dict = sandwich_interior_exterior_calcs.align_sandwich(
            domain_id, sheets, vector
        )

        # Labels residues as interior or exterior
        int_ext_dict = sandwich_interior_exterior_calcs.calc_int_ext(
            domain_id, dssp_df, xy_dict
        )

        # Updates dataframe with solvent accessibility information
        int_ext_list = ['']*dssp_df.shape[0]
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
                 ) = int_ext_pipeline.barrel_pipeline(
                    domain_id, dssp_df, sheets, sec_struct_dfs_dict,
                    domain_sheets_dict, unprocessed_list
                )
            elif self.code[0:4] in ['2.60']:
                (sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list
                 ) = int_ext_pipeline.sandwich_pipeline(
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
