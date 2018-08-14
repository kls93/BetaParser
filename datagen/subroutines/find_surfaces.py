
import math
import copy
import isambard
import pandas as pd
import numpy as np
import networkx as nx
from collections import OrderedDict

if __name__ == 'subroutines.find_surfaces':
    from subroutines.run_stages import run_stages
else:
    from datagen.subroutines.run_stages import run_stages


def align_with_z_axis(ampal_object, xyz_coords, coord_dict, xy_or_z,
                      rel_to_centre, unprocessed_list, domain_id):
    s1 = [xyz_coords[0], xyz_coords[1], xyz_coords[2]]
    e1 = [xyz_coords[3], xyz_coords[4], xyz_coords[5]]
    s2 = [0.0, 0.0, 0.0]
    e2 = [0.0, 0.0, 1.0]
    translation, angle, axis, point = isambard.geometry.find_transformations(
        s1, e1, s2, e2
    )
    ampal_object.rotate(angle, axis, point=point)  # ROTATION MUST BE PERFORMED
    # BEFORE TRANSLATION
    ampal_object.translate(translation)

    # Generates dictionary of residues and their new xyz coordinates
    object_pdb_string = ampal_object.make_pdb().split('\n')
    sub_dict = OrderedDict()

    for line in object_pdb_string:
        if line[0:6].strip() in ['ATOM', 'HETATM']:
            if len(line) > 80:
                unprocessed_list.append(domain_id)
                break
            else:
                atom_id = line[21:27].replace(' ', '') + '_' + line[12:16].strip()

                x_coord_atom = float(line[30:38])
                y_coord_atom = float(line[38:46])
                z_coord_atom = float(line[46:54])
                if xy_or_z == 'xy':
                    xy_coords = np.array([[x_coord_atom],
                                          [y_coord_atom]])
                    sub_dict[atom_id] = xy_coords
                elif xy_or_z == 'z':
                    sub_dict[atom_id] = z_coord_atom

    if rel_to_centre is True:
        avrg_z_coord = np.sum(np.array(list(sub_dict.values())))

        for atom_id in list(sub_dict.keys()):
            z_coord = sub_dict[atom_id]
            new_z_coord = abs(z_coord - avrg_z_coord)
            sub_dict[atom_id] = new_z_coord

    for atom_id in list(sub_dict.keys()):
        coord_dict[atom_id] = sub_dict[atom_id]

    return coord_dict, unprocessed_list


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

        # If no closed circle of interacting strands, the domain isn't
        # technically a barrel and so DataGen is unable to further analyse the
        # domain
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
            # its immediate C-terminal residue are selected as the two central
            # residues
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

    def rotate_translate_barrel(domain_id, xyz_coords, unprocessed_list):
        # Aligns the barrel to z-axis using the reference axis through the
        # barrel pore calculated in the previous step
        print('Aligning {} barrel with z-axis'.format(domain_id))

        # Creates AMPAL object from barrel. Protons are added (using reduce)
        # to allow DataGen to use hydrogen HA3 of a glycine residue as a proxy
        # for the C_alpha atoms of the chiral amino acids in the interior- /
        # exterior-facing calculation performed in the following step.
        barrel = isambard.external_programs.reduce.assembly_plus_protons(
            'Beta_strands/{}.pdb'.format(domain_id)
        )

        # Aligns barrel with the z-axis, returns dictionary of atom ids and
        # their corresponding x and y coordinates in the aligned barrel
        xy_dict = OrderedDict()
        xy_dict, unprocessed_list = align_with_z_axis(
            barrel, xyz_coords, xy_dict, 'xy', False, unprocessed_list,
            domain_id
        )

        return xy_dict, unprocessed_list

    def calc_average_coordinates(domain_id, xy_dict, sheets_df):
        # Calculates average xy coordinates of the barrel C_alpha atoms
        x_sum = 0
        y_sum = 0
        count = 0

        for res_id in list(xy_dict.keys()):
            if (
                res_id.split('_')[0] in sheets_df['RES_ID'].tolist()
                and res_id.split('_')[1] == 'CA'
            ):
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
        print('Determining interior and exterior residues in {}'.format(domain_id))

        res_dict = OrderedDict()
        for index, res_id in enumerate(sheets_df['RES_ID'].tolist()):
            res_dict[res_id] = sheets_df['RESNAME'].tolist()[index]

        for res_id in list(res_dict.keys()):
            if (
                (all('{}_{}'.format(res_id, x) in list(xy_dict.keys())
                     for x in ['N', 'CA', 'C', 'O', 'CB']))
                or (res_dict[res_id] == 'GLY'
                    and all('{}_{}'.format(res_id, x) in list(xy_dict.keys())
                            for x in ['N', 'CA', 'C', 'O', 'HA3']))
            ):
                # Calculates angle between C_alpha, C_beta and the centre of
                # mass
                c_alpha_x = xy_dict['{}_CA'.format(res_id)][0][0]
                c_alpha_y = xy_dict['{}_CA'.format(res_id)][1][0]

                try:
                    c_beta_x = xy_dict['{}_CB'.format(res_id)][0][0]
                    c_beta_y = xy_dict['{}_CB'.format(res_id)][1][0]
                except KeyError:
                    # HA3 is pro-S
                    c_beta_x = xy_dict['{}_HA3'.format(res_id)][0][0]
                    c_beta_y = xy_dict['{}_HA3'.format(res_id)][1][0]

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


class sandwich_strand_position_calcs():

    def find_strand_principal_components(domain_id, sheets_df):
        # Finds principal component of the C_alpha atoms of each individual
        # strand in input beta-sandwich structure
        res_ids_dict = OrderedDict()
        princ_comp_coords_dict = OrderedDict()

        # Makes ordered list of strand ids in sheet_df
        strand_ids_all = sheets_df['STRAND_NUM'].tolist()
        strand_ids = []
        for strand_id in strand_ids_all:
            if not strand_id in strand_ids:
                strand_ids.append(strand_id)

        # Calculates principal component of each strand and stores the
        # resulting coordinates in a dictionary
        for strand_id in strand_ids:
            print(
                'Calculating principal component of {} strand {}'.format(
                    domain_id, strand_id
                )
            )

            strand_df = sheets_df[sheets_df['STRAND_NUM'] == strand_id]
            strand_df = strand_df.reset_index(drop=True)

            res_ids_dict[strand_id] = strand_df['RES_ID'].tolist()

            xyz = np.zeros((strand_df.shape[0], 3))
            for row in range(strand_df.shape[0]):
                xyz[row][0] = strand_df['XPOS'][row]
                xyz[row][1] = strand_df['YPOS'][row]
                xyz[row][2] = strand_df['ZPOS'][row]

            xyz_centred = xyz - xyz.mean(axis=0)

            U, S, V = np.linalg.svd(xyz_centred, full_matrices=True)

            # NEED TO INTRODUCE CHECK TO ENSURE THAT THIS IS THE EIGENVECTOR WITH THE LARGEST CORRESPONDING EIGENVALUE!
            princ_comp_coords_centred = V[0] * np.array([[-1], [1]])
            princ_comp_coords = princ_comp_coords_centred + xyz.mean(axis=0)
            princ_comp_coords = [princ_comp_coords[0][0], princ_comp_coords[0][1],
                                 princ_comp_coords[0][2], princ_comp_coords[1][0],
                                 princ_comp_coords[1][1], princ_comp_coords[1][2]]

            princ_comp_coords_dict[strand_id] = princ_comp_coords

        return princ_comp_coords_dict

    def align_strand_to_princ_comp(domain_id, dssp_df, princ_comp_coords_dict,
                                   res_ids_dict, unprocessed_list):
        # Aligns each strand to its principal component and extracts the
        # z_coordinates of the C_alpha atoms in the aligned strand
        z_dict = OrderedDict()

        for strand_id in list(res_ids_dict.keys()):
            print(
                'Aligning {} strand {} with its principal component'.format(
                    domain_id, strand_id
                )
            )

            pdb_file_lines = []

            for row in range(dssp_df.shape[0]):
                if dssp_df['RES_ID'][row] in res_ids_dict[strand_id]:
                    pdb_file_lines.append(dssp_df['PDB_FILE_LINES'][row])

            pdb_file_lines = ('\n').join(pdb_file_lines)
            pdb_file_lines = pdb_file_lines + '\nEND\n'

            # Creates AMPAL object from sandwich PDB structure
            strand = isambard.ampal.convert_pdb_to_ampal(
                pdb_file_lines, path=False
            )

            # Aligns strand with the z-axis, returns dictionary of atom ids and
            # their corresponding x and y coordinates in the aligned strand
            z_dict, unprocessed_list = align_with_z_axis(
                strand, princ_comp_coords_dict[strand_id], z_dict, 'z', True,
                unprocessed_list, domain_id
            )

        return z_dict, unprocessed_list

    def align_sandwich_to_princ_comp(domain_id, princ_comp_coords_dict,
                                     res_ids_dict, unprocessed_list):
        # Aligns entire beta-sandwich with the principal component of the
        # longest strand, then extracts the z-coordinates of the C_alpha atoms
        # in the aligned sandwich structure
        z_dict = OrderedDict()

        strand_id = max(res_ids_dict, key=lambda x: len(res_ids_dict[x]))
        sandwich = isambard.ampal.convert_pdb_to_ampal(
            'Beta_strands/{}.pdb'.format(domain_id)
        )

        z_dict, unprocessed_list = align_with_z_axis(
            sandwich, princ_comp_coords_dict[strand_id], z_dict, 'z', True,
            unprocessed_list, domain_id
        )

        return z_dict, unprocessed_list

    def centre_z_coords(z_dict):
        # Measures z-coordinates relative to either the centre of the parent
        # strand or the centre of mass of the whole sandwich


        return z_dict

    def update_dataframe(dssp_df, z_dict_strand, z_dict_sandwich):
        # Updates dataframe with z_coordinates (both relative to the parent
        # strand and relative to the centre of mass of the entire sandwich
        # input domain)
        z_coord_parent_strand = ['']*dssp_df.shape[0]
        z_coord_parent_domain = ['']*dssp_df.shape[0]

        for row in range(dssp_df.shape[0]):
            atom_id = '{}_{}'.format(dssp_df['RES_ID'][row], dssp_df['ATMNAME'][row])
            if atom_id in list(z_dict_strand.keys()):
                z_coord_strand[row] = z_dict_strand[atom_id]
            if atom_id in list(z_dict_sandwich.keys()):
                z_coord_sandwich[row] = z_dict_sandwich[atom_id]

        z_coord_df = pd.DataFrame(OrderedDict({'strand_z_coord': z_coord_sandwich,
                                               'sandwich_z_coord': z_coord_sandwich}))
        dssp_df = pd.concat([dssp_df, z_coord_df], axis=1)

        return dssp_df


class int_ext_pipeline():

    def barrel_pipeline(domain_id, sheets, sec_struct_dfs_dict,
                        domain_sheets_dict, unprocessed_list):
        # Determines which face of the beta-sheet forms the interior of the
        # barrel and which forms the exterior based upon whether the angle
        # between the centre of mass of the C_alpha atoms of the barrel and the
        # residue in question is less or greater than 90 degrees. The barrel is
        # collapsed along the central axis that runs through its pore in order
        # to avoid misclassification of residues in the frayed edges of the
        # barrel.
        dssp_df = sec_struct_dfs_dict[domain_id]

        # Checks that correct number of sheets has been retained for analysis
        if len(list(sheets.keys())) != 1:
            print('ERROR: more than 1 sheet retained in {} following solvent '
                  'accessibility calculation'.format(domain_id))
            sec_struct_dfs_dict[domain_id] = None
            for sheet in list(sheets.keys()):
                domain_sheets_dict[sheet] = None
            unprocessed_list.append(domain_id)

            return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

        # Creates dataframe of C_alpha atoms in barrel
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

            return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

        # Finds axis through barrel pore
        xyz_coords = barrel_interior_exterior_calcs.find_z_axis(strands, sheets_df)

        # Aligns barrel with z-axis
        xy_dict, unprocessed_list = barrel_interior_exterior_calcs.rotate_translate_barrel(
            domain_id, xyz_coords, unprocessed_list
        )
        if domain_id in unprocessed_list:
            print('ERROR: Coordinates of rotated and translated barrel are > '
                  '8 characters in length')
            sec_struct_dfs_dict[domain_id] = None
            for sheet in list(sheets.keys()):
                domain_sheets_dict[sheet] = None

            return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

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

        # Updates dataframe to identify interior- / exterior-facing residues
        for row in range(dssp_df.shape[0]):
            res_id = dssp_df['RES_ID'][row]
            if (
                res_id in list(int_ext_dict.keys())
                and dssp_df['ATMNAME'][row] == 'CA'
            ):
                int_ext_list[row] = int_ext_dict[res_id]
        int_ext_df = pd.DataFrame({'INT_EXT': int_ext_list})
        dssp_df = pd.concat([dssp_df, int_ext_df], axis=1)
        sec_struct_dfs_dict[domain_id] = dssp_df

        return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

    def sandwich_pipeline(domain_id, sheets, sec_struct_dfs_dict,
                          domain_sheets_dict, unprocessed_list):
        # Calculates whether a residue faces towards the interior or the
        # exterior of the sandwich
        dssp_df = sec_struct_dfs_dict[domain_id]

        # Checks that correct number of sheets have been retained in the
        # beta-sandwich
        if len(list(sheets.keys())) != 2:
            print('ERROR: incorrect number of sheets retained in {} following '
                  'solvent accessibility calculation'.format(domain_id))
            sec_struct_dfs_dict[domain_id] = None
            for sheet in list(sheets.keys()):
                domain_sheets_dict[sheet] = None
            unprocessed_list.append(domain_id)

            return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

        # Creates dataframe of C_alpha atoms in sandwich
        sheet_ids = [sheet.replace('{}_sheet_'.format(domain_id), '') for sheet
                     in set(sheets.keys()) if sheet != '']
        sheets_df = dssp_df[dssp_df['SHEET_NUM'].isin(sheet_ids)]
        sheets_df = sheets_df.reset_index(drop=True)

        # Calculates principal component of each strand in sandwich

        # Calculates z-axis between the two beta-sheets
        vector = sandwich_interior_exterior_calcs.find_z_axis(sheets_df)

        # Aligns the beta-sandwich with z-axis
        xy_dict, unprocessed_list = sandwich_interior_exterior_calcs.align_sandwich(
            domain_id, vector, unprocessed_list
        )
        if domain_id in unprocessed_list:
            print('ERROR: Coordinates of rotated and translated sandwich are '
                  '> 8 characters in length')
            sec_struct_dfs_dict[domain_id] = None
            for sheet in list(sheets.keys()):
                domain_sheets_dict[sheet] = None

            return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

        # Labels residues as interior or exterior
        int_ext_dict = sandwich_interior_exterior_calcs.calc_int_ext(
            domain_id, dssp_df, xy_dict
        )

        # Updates dataframe to identify interior- / exterior-facing residues
        int_ext_list = ['']*dssp_df.shape[0]
        for row in range(dssp_df.shape[0]):
            if (
                dssp_df['RES_ID'][row] in list(int_ext_dict.keys())
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
            sheets = {key: value for key, value in domain_sheets_dict.items()
                      if domain_id in key}

            if self.code[0:4] in ['2.40']:
                (sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list
                 ) = int_ext_pipeline.barrel_pipeline(
                    domain_id, sheets, sec_struct_dfs_dict, domain_sheets_dict,
                    unprocessed_list
                )
            elif self.code[0:4] in ['2.60']:
                (sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list
                 ) = int_ext_pipeline.sandwich_pipeline(
                    domain_id, sheets, sec_struct_dfs_dict, domain_sheets_dict,
                    unprocessed_list
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
