
import sys
import math
import itertools
import isambard
import pandas as pd
import numpy as np
import networkx as nx
from collections import OrderedDict
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
if __name__ == 'subroutines.find_surfaces':
    from subroutines.run_stages import run_stages
else:
    from datagen.subroutines.run_stages import run_stages


class barrel_interior_exterior_calcs():

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


class sandwich_interior_exterior_calcs():

    def find_edge_strand_pair(domain_id, dssp_df, networks, domain_sheets_dict):
        # Selects a pair of edge strands (that together form one of the ends
        # of the sandwich) as the two strands in opposite sheets whose centres
        # of mass are closest

        # Identifies edge strands in the sandwich
        edges = []
        for G in networks:
            strands = nx.nodes(G)
            for strand in strands:
                if len(G.neighbors(strand)) == 1:
                    edges.append(strand)

        # Finds centre of mass of the C_alpha atoms of each edge strand
        edge_com_dict = {}
        edge_length_dict = {}
        for strand in edges:
            strand_df = dssp_df[dssp_df['STRAND_NUM'] == strand]
            strand_df = strand_df.reset_index(drop=True)

            x_coords = [num for num in strand_df['XPOS'].tolist() if num != '']
            x_coord = sum(x_coords) / len(x_coords)
            y_coords = [num for num in strand_df['YPOS'].tolist() if num != '']
            y_coord = sum(y_coords) / len(y_coords)
            z_coords = [num for num in strand_df['ZPOS'].tolist() if num != '']
            z_coord = sum(z_coords) / len(z_coords)

            edge_com_dict[strand] = np.array([[x_coord], [y_coord], [z_coord]])
            edge_length_dict[strand_df.shape[0]] = strand

        # Combines the beta-sheet networks into a single beta-sandwich network
        G = networks[0]
        for H in networks[1:]:
            G = nx.compose(G, H)

        # Selects the longest edge strand to be one of the strands in the pair
        max_length = max(list(edge_length_dict.keys()))
        strand_1 = edge_length_dict[max_length]

        # Uses the distance between strand centres of mass to determine the
        # pair closest in space (and so are assumed to form the same edge of
        # the sandwich)
        dist_dict = {}
        for strand in edges:
            if ((strand_1 != strand) and (strand_1 not in G.neighbors(strand))):
                dist = math.sqrt(((edge_com_dict[strand][0][0] -
                                   edge_com_dict[strand_1][0][0])**2)
                                 + ((edge_com_dict[strand][1][0] -
                                     edge_com_dict[strand_1][1][0])**2)
                                 + ((edge_com_dict[strand][2][0] -
                                     edge_com_dict[strand_1][2][0])**2))
                dist_dict[dist] = (strand, strand_1)
        min_dist = min(list(dist_dict.keys()))
        strand_pair = dist_dict[min_dist]

        return strand_pair

    def find_z_axis(sheets_df, strand_pair):
        # Finds axis between the two sheets of the beta-sandwiches
        sheets_xyz_all_coords = []
        sheets_xyz_sub_coords = []

        for row in range(sheets_df.shape[0]):
            x_coord = sheets_df['XPOS'][row]
            y_coord = sheets_df['YPOS'][row]
            z_coord = sheets_df['ZPOS'][row]
            sheets_xyz_all_coords.append((x_coord, y_coord, z_coord))

            if sheets_df['STRAND_NUM'][row] not in strand_pair:
                sheets_xyz_sub_coords.append((x_coord, y_coord, z_coord))

        # Calculates centre of mass of C_alpha atoms across the two beta-sheets
        sheets_xyz_all_coords = np.array(sheets_xyz_all_coords)

        x_coord_all = np.sum(sheets_xyz_all_coords[:, 0]) / sheets_xyz_all_coords.shape[0]
        y_coord_all = np.sum(sheets_xyz_all_coords[:, 1]) / sheets_xyz_all_coords.shape[0]
        z_coord_all = np.sum(sheets_xyz_all_coords[:, 2]) / sheets_xyz_all_coords.shape[0]

        # Calculates centre of mass of C_alpha atoms across the two beta-sheets
        # excluding atoms in the edge strand pair identified earlier
        sheets_xyz_sub_coords = np.array(sheets_xyz_sub_coords)

        x_coord_sub = np.sum(sheets_xyz_sub_coords[:, 0]) / sheets_xyz_sub_coords.shape[0]
        y_coord_sub = np.sum(sheets_xyz_sub_coords[:, 1]) / sheets_xyz_sub_coords.shape[0]
        z_coord_sub = np.sum(sheets_xyz_sub_coords[:, 2]) / sheets_xyz_sub_coords.shape[0]

        xyz_coords = [x_coord_all, y_coord_all, z_coord_all, x_coord_sub,
                      y_coord_sub, z_coord_sub]

        return xyz_coords

    def align_sandwich(domain_id, xyz_coords):
        # TODO: STILL PROBLEMS WITH AXIS ALIGNMENT, NEEDS FIXING

        # Creates an ampal object of the beta-sandwich
        sandwich = isambard.ampal.convert_pdb_to_ampal('Beta_strands/{}.pdb'.format(domain_id))

        # Aligns the sandwich with z = 0
        s1 = [xyz_coords[0], xyz_coords[1], xyz_coords[2]]
        e1 = [xyz_coords[3], xyz_coords[4], xyz_coords[5]]
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

        with open('/home/ks17361/my_sandwich.pdb', 'w') as pdb_file:
            for line in sandwich_pdb_string:
                pdb_file.write('{}\n'.format(line))
            pdb_file.write('HETATM {} {} {}\n'.format(xyz_coords[0], xyz_coords[1], xyz_coords[2]))
            pdb_file.write('HETATM {} {} {}\n'.format(xyz_coords[3], xyz_coords[4], xyz_coords[5]))

        for line in sandwich_pdb_string:
            if line[0:6].strip() in ['ATOM', 'HETATM']:
                atom_id = line[21:27].replace(' ', '') + '_' + line[12:16].strip()

                line_segments = line.split('.')
                x_coord_atom = float(line_segments[0].split()[-1] + '.' + line_segments[1][0:3])
                y_coord_atom = float(line_segments[1][3:] + '.' + line_segments[2][0:3])
                xy_coords = np.array([[x_coord_atom],
                                      [y_coord_atom]])

                xy_dict[atom_id] = xy_coords

        return xy_dict

    def calc_int_ext(domain_id, dssp_df, xy_dict):
        # Determines whether a residue is interior- or exterior-facing from
        # whether its C_beta atom lies within the polygon formed from the
        # xy-coordinates of the backbone (N, C_alpha and C) atoms of the
        # sandwich
        print('Determining interior and exterior residues in {}'.format(domain_id))

        int_ext_dict = OrderedDict()

        backbone_coords = []
        for res in list(xy_dict.keys()):
            if any(x in res for x in ['N', 'CA', 'C']):
                backbone_coords.append((xy_dict[res][0][0], xy_dict[res][1][0]))
        outline = Polygon(backbone_coords)

        for res in list(xy_dict.keys()):
            if 'CB' in res:
                point = Point(xy_dict[res][0][0], xy_dict[res][1][0])
                if outline.contains(point):
                    int_ext_dict[res.split('_')[0]] = 'interior'
                elif not outline.contains(point):
                    int_ext_dict[res.split('_')[0]] = 'exterior'

        print(int_ext_dict)
        import sys
        sys.exit()

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

        print(domain_id)

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

        # Calculates interior- and exterior-facing residues
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

    def sandwich_pipeline(domain_id, dssp_df, sheets, sec_struct_dfs_dict,
                          domain_sheets_dict, unprocessed_list):
        # Calculates whether a residue faces towards the interior or the
        # exterior of the sandwich

        # Checks that correct number of sheets have been retained in the
        # beta-sandwich
        networks = [network for key, network in domain_sheets_dict.items()
                    if domain_id in key]

        if len(sheets) != 2 or len(networks) != 2:
            print('ERROR: incorrect number of sheets retained in {} following '
                  'solvent accessibility calculation'.format(domain_id))
            sec_struct_dfs_dict[domain_id] = None
            for sheet in sheets:
                domain_sheets_dict[sheet] = None
            unprocessed_list.append(domain_id)

            return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

        # Identifies pair of edge strands
        strand_pair = sandwich_interior_exterior_calcs.find_edge_strand_pair(
            domain_id, dssp_df, networks, domain_sheets_dict
        )

        # Calculates z-axis between the two beta-sheets
        sheet_ids = [sheet.replace('{}_sheet_'.format(domain_id), '') for sheet
                     in sheets if sheet != '']
        sheets_df = dssp_df[dssp_df['SHEET_NUM'].isin(sheet_ids)]
        sheets_df = sheets_df.reset_index(drop=True)
        xyz_coords = sandwich_interior_exterior_calcs.find_z_axis(
            sheets_df, strand_pair
        )

        # Aligns the beta-sandwich with z = 0
        xy_dict = sandwich_interior_exterior_calcs.align_sandwich(
            domain_id, xyz_coords
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
