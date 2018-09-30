
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


def align_with_z_axis(ampal_object, xyz_coords, xy_or_z, rel_to_centre,
                      unprocessed_list, domain_id):
    # Aligns input AMPAL object with the z_axis, and returns a dictionary of
    # coordinates from the aligned object
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

    # Generates dictionary of atoms and their new xyz coordinates
    object_pdb_string = ampal_object.make_pdb().split('\n')
    coord_dict = OrderedDict()

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
                    coord_dict[atom_id] = xy_coords
                elif xy_or_z == 'z':
                    if line[12:16].strip() == 'CA':
                        coord_dict[atom_id] = z_coord_atom

    if not domain_id in unprocessed_list:
        if rel_to_centre is True:
            z_coord_vals = list(coord_dict.values())
            avrg_z_coord = np.sum(np.array(z_coord_vals)) / len(z_coord_vals)

            for atom_id in list(coord_dict.keys()):
                z_coord = coord_dict[atom_id]
                new_z_coord = round(abs(z_coord - avrg_z_coord), 3)
                coord_dict[atom_id] = new_z_coord

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

    def find_barrel_principal_component(domain_id, barrel_strands, sheets_df):
        # Finds principal component of the C_alpha atoms of beta-strands in
        # barrel
        princ_comp_coords_dict = OrderedDict()

        # Filters sheets_df to retain Calpha atoms in strands in the main barrel
        sheets_df = sheets_df[sheets_df['STRAND_NUM'].isin(barrel_strands)]
        sheets_df = sheets_df.reset_index(drop=True)

        # Calculates principal component of strands in barrel and stores the
        # resulting coordinates in a dictionary
        print('Calculating principal component of {}'.format(domain_id)

        xyz = np.zeros((sheets_df.shape[0], 3))
        for row in range(sheets_df.shape[0]):
            xyz[row][0] = sheets_df['XPOS'][row]
            xyz[row][1] = sheets_df['YPOS'][row]
            xyz[row][2] = sheets_df['ZPOS'][row]

        xyz_centred = xyz - xyz.mean(axis=0)

        U, S, V = np.linalg.svd(xyz_centred, full_matrices=True)

        # The eigenvectors in V are ordered by the size of their
        # corresponding eigenvalues, so can simply take V[0] as max.
        princ_comp_coords_centred = V[0] * np.array([[-1], [1]])
        princ_comp_coords = princ_comp_coords_centred + xyz.mean(axis=0)
        princ_comp_coords = [princ_comp_coords[0][0], princ_comp_coords[0][1],
                             princ_comp_coords[0][2], princ_comp_coords[1][0],
                             princ_comp_coords[1][1], princ_comp_coords[1][2]]

        return sheets_df, princ_comp_coords

    def align_barrel_princ_comp_to_z(domain_id, princ_comp_coords,
                                     unprocessed_list):
        # Aligns barrel principal component with the z-axis, then extracts the
        # z-coordinates of the C_alpha atoms in the aligned sandwich structure
        print('Aligning {} barrel with z-axis'.format(domain_id))

        # Creates AMPAL object from barrel. Protons are added (using reduce)
        # to allow DataGen to use hydrogen HA3 of a glycine residue as a proxy
        # for the C_alpha atoms of the chiral amino acids in the interior- /
        # exterior-facing calculation performed in the following step.
        barrel = isambard.external_programs.reduce.assembly_plus_protons(
            'Beta_strands/{}.pdb'.format(domain_id)
        )

        xy_dict, unprocessed_list = align_with_z_axis(
            barrel, princ_comp_coords, 'xy', False, unprocessed_list, domain_id
        )

        return xy_dict, unprocessed_list

    def calc_average_coordinates(domain_id, xy_dict, sheets_df):
        # Calculates average xy coordinates of the C_alpha atoms in the strands
        # in the main barrel
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

    def calc_int_ext(domain_id, dssp_df, xy_dict, com, int_ext_dict):
        # Calculates whether residue in the barrel face towards the interior or
        # the exterior of the barrel. Note that all strands are included in
        # this analysis, not just the strands that form the main barrel.
        print('Determining interior and exterior residues in {}'.format(domain_id))

        res_dict = OrderedDict()
        for row in range(dssp_df.shape[0]):
            if dssp_df['ATMNAME'][row] == 'CA':
                res_dict[dssp_df['RES_ID'][row]] = dssp_df['RESNAME'][row]

        for res_id in list(res_dict.keys()):
            if (
                   (res_dict[res_id] != 'GLY'
                    and all('{}_{}'.format(res_id, x) in list(xy_dict.keys())
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

        # Makes ordered list of strand ids in sheet_df. (List is ordered only
        # because I prefer an ordered list of strands to be printed in the
        # terminal.)
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

            # The eigenvectors in V are ordered by the size of their
            # corresponding eigenvalues, so can simply take V[0] as max.
            princ_comp_coords_centred = V[0] * np.array([[-1], [1]])
            princ_comp_coords = princ_comp_coords_centred + xyz.mean(axis=0)
            princ_comp_coords = [princ_comp_coords[0][0], princ_comp_coords[0][1],
                                 princ_comp_coords[0][2], princ_comp_coords[1][0],
                                 princ_comp_coords[1][1], princ_comp_coords[1][2]]

            princ_comp_coords_dict[strand_id] = princ_comp_coords

        return princ_comp_coords_dict, res_ids_dict

    def align_strand_to_princ_comp(domain_id, dssp_df, princ_comp_coords_dict,
                                   res_ids_dict, unprocessed_list):
        # Aligns each strand to its principal component and extracts the
        # z_coordinates of the C_alpha atoms in the aligned strand
        for strand_id in list(res_ids_dict.keys()):
            print(
                'Aligning {} strand {} principal component with the z-axis'.format(
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
                strand, princ_comp_coords_dict[strand_id], 'z', True,
                unprocessed_list, domain_id
            )

        return z_dict, unprocessed_list

    def align_sandwich_to_princ_comp(domain_id, princ_comp_coords_dict,
                                     res_ids_dict, unprocessed_list):
        # Aligns entire beta-sandwich with the principal component of the
        # longest strand, then extracts the z-coordinates of the C_alpha atoms
        # in the aligned sandwich structure
        print('Aligning {} with the z-axis'.format(domain_id))

        strand_id = max(res_ids_dict, key=lambda x: len(res_ids_dict[x]))
        xyz = princ_comp_coords_dict[strand_id]

        # Creates AMPAL object from sandwich.
        sandwich = isambard.ampal.convert_pdb_to_ampal(
            'Beta_strands/{}.pdb'.format(domain_id)
        )

        z_dict, unprocessed_list = align_with_z_axis(
            sandwich, xyz, 'z', True, unprocessed_list, domain_id
        )

        return z_dict, unprocessed_list

    def update_dataframe(dssp_df, z_dict_strand, z_dict_sandwich):
        # Updates dataframe with z_coordinates (both relative to the parent
        # strand and relative to the centre of mass of the entire sandwich
        # input domain)
        z_coord_parent_strand = ['']*dssp_df.shape[0]
        z_coord_parent_domain = ['']*dssp_df.shape[0]

        for row in range(dssp_df.shape[0]):
            atom_id = '{}_{}'.format(dssp_df['RES_ID'][row], dssp_df['ATMNAME'][row])
            if atom_id in list(z_dict_strand.keys()):
                z_coord_parent_strand[row] = z_dict_strand[atom_id]
            if atom_id in list(z_dict_sandwich.keys()):
                z_coord_parent_domain[row] = z_dict_sandwich[atom_id]

        z_coord_df = pd.DataFrame(OrderedDict({'STRAND_Z_COORDS': z_coord_parent_strand,
                                               'SANDWICH_Z_COORDS': z_coord_parent_domain}))
        dssp_df = pd.concat([dssp_df, z_coord_df], axis=1)

        return dssp_df


class pipeline():

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

        # Creates dataframe of C_alpha atoms in barrel
        sheet_num = sheets[0].replace('{}_sheet_'.format(domain_id), '')
        sheets_df = dssp_df[dssp_df['SHEET_NUM'] == sheet_num]
        sheets_df = sheets_df.reset_index(drop=True)

        # Finds network of interacting neighbouring strands
        G = domain_sheets_dict[sheets[0]]
        strands, unprocessed_list = barrel_interior_exterior_calcs.find_strands_network(
            G, domain_id, unprocessed_list
        )

        # Checks complete circle of hydrogen bonded strands is present in domain
        if domain_id in unprocessed_list:
            print('ERROR: domain {} does not contain closed circular network '
                  'of hydrogen-bonded strands'.format(domain_id))
            sec_struct_dfs_dict[domain_id] = None
            for sheet in sheets:
                domain_sheets_dict[sheet] = None

            return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

        # Finds barrel principal component (= axis through barrel pore,
        # provided that this is axis is longer than the diameter of the barrel)
        (
        sheets_df, princ_comp_coords
        ) = barrel_interior_exterior_calcs.find_barrel_principal_component(
            domain_id, strands, sheets_df
        )

        # Aligns barrel with z-axis
        (
        xy_dict, unprocessed_list
        ) = barrel_interior_exterior_calcs.align_barrel_princ_comp_to_z(
            domain_id, princ_comp_coords, unprocessed_list
        )

        if domain_id in unprocessed_list:
            print('ERROR: Coordinates of rotated and translated barrel are > '
                  '8 characters in length')
            sec_struct_dfs_dict[domain_id] = None
            for sheet in sheets:
                domain_sheets_dict[sheet] = None

            return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

        # Calculates average xy coordinates of the barrel C_alpha atoms
        com = barrel_interior_exterior_calcs.calc_average_coordinates(
            domain_id, xy_dict, sheets_df
        )

        # Initialises records of interior / exterior facing residues
        int_ext_list = ['']*dssp_df.shape[0]
        int_ext_dict = OrderedDict()

        # Calculates interior- and exterior-facing residues
        int_ext_dict = barrel_interior_exterior_calcs.calc_int_ext(
            domain_id, dssp_df, xy_dict, com, int_ext_dict
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

        # Creates dataframe of C_alpha atoms in sandwich
        sheet_ids = [sheet.replace('{}_sheet_'.format(domain_id), '') for sheet
                     in sheets]
        sheets_df = dssp_df[dssp_df['SHEET_NUM'].isin(sheet_ids)]
        sheets_df = sheets_df.reset_index(drop=True)

        # Calculates principal component of each strand in sandwich
        (
        strand_princ_comps_dict, res_ids_dict
        ) = sandwich_strand_position_calcs.find_strand_principal_components(
            domain_id, sheets_df
        )

        # Aligns each individual strand to the z-axis via its principal
        # component
        (
        z_dict_indv_strands, unprocessed_list
        ) = sandwich_strand_position_calcs.align_strand_to_princ_comp(
            domain_id, dssp_df, strand_princ_comps_dict, res_ids_dict,
            unprocessed_list
        )

        # Aligns the entire sandwich domain to the z-axis via the principal
        # component of the longest strand
        (
        z_dict_sandwich, unprocessed_list
        ) = sandwich_strand_position_calcs.align_sandwich_to_princ_comp(
            domain_id, strand_princ_comps_dict, res_ids_dict,
            unprocessed_list
        )

        # Updates dataframe to include new z_coords (**DON'T MIX UP ORDER OF
        # STRAND AND SANDWICH DICTS**)
        dssp_df = sandwich_strand_position_calcs.update_dataframe(
            dssp_df, z_dict_indv_strands, z_dict_sandwich
        )

        if domain_id in unprocessed_list:
            print('ERROR whilst calculating principal component of sandwich')
            sec_struct_dfs_dict[domain_id] = None
            for sheet in sheets:
                domain_sheets_dict[sheet] = None
        else:
            sec_struct_dfs_dict[domain_id] = dssp_df

        return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list


class find_interior_exterior_surfaces(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    def run_pipeline(self, sec_struct_dfs_dict, domain_sheets_dict):
        # Pipeline script to identify interior and exterior-facing residues in
        # barrels / sandwiches
        unprocessed_list = []

        for domain_id in list(sec_struct_dfs_dict.keys()):
            sheets = [key for key, value in domain_sheets_dict.items() if
                      domain_id in key]

            if self.code[0:4] in ['2.40']:
                (sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list
                 ) = pipeline.barrel_pipeline(
                    domain_id, sheets, sec_struct_dfs_dict, domain_sheets_dict,
                    unprocessed_list
                )
            elif self.code[0:4] in ['2.60']:
                (sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list
                 ) = pipeline.sandwich_pipeline(
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
