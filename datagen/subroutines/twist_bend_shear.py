
import copy
import os
import isambard_dev as isambard
import networkx as nx
import numpy as np
import pandas as pd
import scipy.stats as stats
from collections import OrderedDict
if __name__ == 'subroutines.twist_bend_shear':
    from subroutines.run_stages import run_stages
else:
    from datagen.subroutines.run_stages import run_stages


def find_strand_bend(domain_id, dssp_df):
    """
    Calculates the bend angles of each residue in a strand as the angle
    between the vectors between its alpha carbon atom and the alpha carbon
    atoms of its neighbouring (+/-1) residues
    """

    print('Calculating bend for residues in {}'.format(domain_id))

    # Creates AMPAL object
    pdb = isambard.ampal.convert_pdb_to_ampal(
        'Parent_assemblies/{}.pdb'.format(domain_id)
    )

    # Calculates strand bend. All 3 residues must be in the same beta-strand.
    strand_res = []
    strand_ids = []
    for i, res_id in enumerate(dssp_df['RES_ID'].tolist()):
        if dssp_df['ATMNAME'][i] == 'CA':
            if not res_id in strand_res:
                strand_res.append(res_id)
                strand_id = dssp_df['STRAND_NUM'][i]
                strand_ids.append(strand_id)
    all_res_list = []
    for n in range(len(pdb)):
        all_res_list += list(pdb[n])
    bend_dict = {}

    for i_sub, res_id_2 in enumerate(strand_res):
        strand_id_2 = strand_ids[i_sub]

        res_1 = ''
        res_2 = ''
        res_3 = ''
        angle = ''
        for i_all, res_2 in enumerate(all_res_list):
            # Identifies AMPAL object for res_id_2
            chain_2 = res_2.ampal_parent.id
            res_num_2 = res_2.id
            ins_code_2 = res_2.insertion_code

            if res_id_2 == '{}{}{}'.format(chain_2, res_num_2, ins_code_2):
                if (i_all > 0) and (i_all < (len(all_res_list)-1)):
                    # Identifies AMPAL objects for neighbouring residues
                    res_1 = all_res_list[i_all-1]
                    chain_1 = res_1.ampal_parent.id
                    res_num_1 = res_1.id
                    ins_code_1 = res_1.insertion_code
                    res_id_1 = '{}{}{}'.format(chain_1, res_num_1, ins_code_1)

                    res_3 = all_res_list[i_all+1]
                    chain_3 = res_3.ampal_parent.id
                    res_num_3 = res_3.id
                    ins_code_3 = res_3.insertion_code
                    res_id_3 = '{}{}{}'.format(chain_3, res_num_3, ins_code_3)

                    if all(x in strand_res for x in [res_id_1, res_id_2, res_id_3]):
                        strand_id_1 = strand_ids[strand_res.index(res_id_1)]
                        strand_id_3 = strand_ids[strand_res.index(res_id_3)]

                        if (strand_id_1 == strand_id_2) and (strand_id_2 == strand_id_3):
                            calpha_2 = np.array([res_2['CA'].x, res_2['CA'].y, res_2['CA'].z])
                            calpha_1 = np.array([res_1['CA'].x, res_1['CA'].y, res_1['CA'].z])
                            calpha_3 = np.array([res_3['CA'].x, res_3['CA'].y, res_3['CA'].z])
                            vector_1 = calpha_2 - calpha_1
                            vector_2 = calpha_2 - calpha_3
                            angle = isambard.tools.geometry.angle_between_vectors(vector_1, vector_2)
                bend_dict[res_id_2] = angle
                break

    # Returns list of bend angles
    bend_angles = ['']*dssp_df.shape[0]
    for row in range(dssp_df.shape[0]):
        res_id = dssp_df['RES_ID'][row]
        if (res_id in list(bend_dict.keys())) and (dssp_df['ATMNAME'][row] == 'CA'):
            bend_angles[row] = bend_dict[res_id]

    return bend_angles


def find_strand_twist(domain_id, domain_df, strand_or_res):
    """
    Calculates strand twist for every HB and NHB pair. Residues 1 and 2 are in
    strand A, residues 3 and 4 are in strand B.
    """

    print('Calculating twist for residues in {}'.format(domain_id))

    # Creates AMPAL object
    pdb = isambard.ampal.convert_pdb_to_ampal(
        'Parent_assemblies/{}.pdb'.format(domain_id)
    )

    # Lists AMPAL objects + corresponding ids of all res in the parent assembly
    all_res_ampal = []
    all_res_ids = []
    for n in range(len(pdb)):
        all_res_ampal += list(pdb[n])
        for i, res in enumerate(pdb[n]):
            res_id = '{}{}{}'.format(
                pdb[n][i].ampal_parent.id, pdb[n][i].id, pdb[n][i].insertion_code
            )
            all_res_ids.append(res_id)

    if strand_or_res == 'res':
        twist_list_hb_pairs = ['']*domain_df.shape[0]
        angle_list_hb_pairs = ['']*domain_df.shape[0]
        twist_list_nhb_pairs = ['']*domain_df.shape[0]
        angle_list_nhb_pairs = ['']*domain_df.shape[0]
    elif strand_or_res == 'strand':
        twist_list_hb_pairs = [[] for n in range(domain_df.shape[0])]
        angle_list_hb_pairs = [[] for n in range(domain_df.shape[0])]
        twist_list_nhb_pairs = [[] for n in range(domain_df.shape[0])]
        angle_list_nhb_pairs = [[] for n in range(domain_df.shape[0])]

    # Makes sub-list of residues that form the beta-strands in the domain, plus
    # their HB and NHB pairs
    strands = []
    strands_res_full_ids = []
    hb_pairs = []
    nhb_pairs = []
    for strand in domain_df['domain_strand_ids'].tolist():
        if not strand in strands:
            strands.append(strand)
    if strand_or_res == 'res':
        strands_res_full_ids = [
            '{}_{}'.format(domain_df['domain_strand_ids'][n], domain_df['res_ids'][n])
            for n in range(domain_df.shape[0])
        ]
        hb_pairs = domain_df['hb_pairs'].tolist()
        nhb_pairs = domain_df['nhb_pairs'].tolist()
    elif strand_or_res == 'strand':
        for n in range(domain_df.shape[0]):
            sub_strands_res_full_ids = [
                '{}_{}'.format(domain_df['domain_strand_ids'][n],
                domain_df['res_ids'][n][o])
                for o in range(len(domain_df['res_ids'][n]))
            ]
            strands_res_full_ids += sub_strands_res_full_ids
            hb_pairs += domain_df['hb_pairs'][n]
            nhb_pairs += domain_df['nhb_pairs'][n]
    strands_res = [res.split('_')[-1] for res in strands_res_full_ids]
    strands_ids = ['_'.join(res.split('_')[:-1]) for res in strands_res_full_ids]

    # Calculates strand twist angles for HB and NHB pairs and adds them into the
    # input domain_df
    for res_2_index, res_2_id in enumerate(copy.deepcopy(strands_res)):
        res_2_strand_id = strands_ids[res_2_index]
        res_2_strand_index = strands.index(res_2_strand_id)
        res_2_full_id = '{}_{}'.format(res_2_strand_id, res_2_id)
        res_2 = all_res_ampal[all_res_ids.index(res_2_id)]

        for pair_index, pair_list in enumerate([hb_pairs, nhb_pairs]):
            res_1 = ''
            res_3 = ''
            res_4 = ''
            res_1_id = ''
            res_3_id = ''
            res_4_id = ''
            res_1_strand_id = ''
            res_3_strand_id = ''
            res_4_strand_id = ''
            res_1_full_id = ''
            res_3_full_id = ''
            res_4_full_id = ''

            # Identifies residue 3 (forms HB/NHB bond with residue 2)
            if pair_list[res_2_index] != []:
                res_3_id = pair_list[res_2_index][0]
                res_3 = all_res_ampal[all_res_ids.index(res_3_id)]
                try:
                    res_3_strand_id = strands_ids[strands_res.index(res_3_id)]
                    res_3_full_id = '{}_{}'.format(res_3_strand_id, res_3_id)
                except ValueError:
                    pass

            # Identifies residue 1 (residue 2 - 1 position)
            res_1_index = all_res_ids.index(res_2_id) - 1
            if 0 <= res_1_index and res_1_index < len(all_res_ids):
                res_1_id = all_res_ids[res_1_index]
                res_1 = all_res_ampal[res_1_index]
                try:
                    res_1_strand_id = strands_ids[strands_res.index(res_1_id)]
                    res_1_full_id = '{}_{}'.format(res_1_strand_id, res_1_id)
                except ValueError:
                    pass

            # Identifies residue 4 (residue 3 + 1 position)
            if res_3 != '':
                res_4_index = all_res_ids.index(res_3_id) + 1
                if 0 <= res_4_index and res_4_index < len(all_res_ids):
                    res_4_id = all_res_ids[res_4_index]
                    res_4 = all_res_ampal[res_4_index]
                    try:
                        res_4_strand_id = strands_ids[strands_res.index(res_4_id)]
                        res_4_full_id = '{}_{}'.format(res_4_strand_id, res_4_id)
                    except ValueError:
                        pass

            twist = ''
            angle = ''
            if (    all(x != '' for x in [res_1, res_2, res_3, res_4])
                and all(x in strands_res_full_ids for x in
                [res_1_full_id, res_2_full_id, res_3_full_id, res_4_full_id])
            ):
                calpha_1 = np.array([res_1['CA'].x, res_1['CA'].y, res_1['CA'].z])
                calpha_2 = np.array([res_2['CA'].x, res_2['CA'].y, res_2['CA'].z])
                calpha_3 = np.array([res_3['CA'].x, res_3['CA'].y, res_3['CA'].z])
                calpha_4 = np.array([res_4['CA'].x, res_4['CA'].y, res_4['CA'].z])

                angle = isambard.tools.geometry.dihedral(calpha_1, calpha_2, calpha_3, calpha_4)
                if abs(angle) < 90:  # Parallel beta strands
                    if angle < 0:
                        twist = 'right'
                    elif angle > 0:
                        twist = 'left'
                elif abs(angle) > 90:  # Antiparallel beta strands
                    if angle < 0:
                        twist = 'left'
                    elif angle > 0:
                        twist = 'right'

            if pair_list == hb_pairs:
                if strand_or_res == 'res':
                    twist_list_hb_pairs[res_2_index] = twist
                    angle_list_hb_pairs[res_2_index] = angle
                elif strand_or_res == 'strand':
                    twist_list_hb_pairs[res_2_strand_index].append(twist)
                    angle_list_hb_pairs[res_2_strand_index].append(angle)
            elif pair_list == nhb_pairs:
                if strand_or_res == 'res':
                    twist_list_nhb_pairs[res_2_index] = twist
                    angle_list_nhb_pairs[res_2_index] = angle
                elif strand_or_res == 'strand':
                    twist_list_nhb_pairs[res_2_strand_index].append(twist)
                    angle_list_nhb_pairs[res_2_strand_index].append(angle)

    # Appends twist angles to domain_df
    angle_df = pd.DataFrame(OrderedDict({
        'twist_hb_pairs': twist_list_hb_pairs,
        'twist_angle_hb_pairs': angle_list_hb_pairs,
        'twist_nhb_pairs': twist_list_nhb_pairs,
        'twist_angle_nhb_pairs': angle_list_nhb_pairs
    }))
    upd_domain_df = copy.deepcopy(domain_df).reset_index(drop=True)
    upd_domain_df = pd.concat([upd_domain_df, angle_df], axis=1)

    return upd_domain_df


def find_sheet_shear(domain_id, domain_df, domain_sheets_dict):
    """
    Calculates the shear number of beta-barrels.
    **First strand in calculation must not contain any beta-bulges. Currently
    does this by checking has alternating pattern of interior and exterior
    residues, and that phi and psi angles are in beta-range of Ramachandran
    space.**
    """

    print('Calculating shear for {}'.format(domain_id))

    # Determines order of interacting strands
    networks = [network for key, network in domain_sheets_dict.items()
                if domain_id in key]
    if len(networks) > 1:
        raise Exception('Only a single sheet expected in {}'.format(domain_id))
    G = networks[0]

    strand_order = nx.cycle_basis(G)
    if len(strand_order) > 1:
        raise Exception(
            'Only a single circular network expected for {}'.format(domain_id)
        )
    strand_order = strand_order[0]
    if sorted(strand_order) != sorted(list(G.nodes())):
        raise Exception(
            'Expect non-barrel strands to have been removed in {}'.format(domain_id)
        )

    # Reorganises strands to ensure first strand has an alternating pattern of
    # interior and exterior residues
    start = ''
    for strand in strand_order:
        strand_df = domain_df[
            domain_df['domain_strand_ids']=='{}_strand_{}'.format(domain_id, strand)
        ].reset_index(drop=True)
        int_ext_pattern = strand_df['int_ext'].tolist()
        even = int_ext_pattern[::2]
        odd = int_ext_pattern[1::2]
        if (len(set(even)) == 1) and (len(set(odd)) == 1) and (set(even) != set(odd)):
            start = strand
            break
    if start == '':
        raise Exception(
            'No strands without beta-bulges identified in {}'.format(domain_id)
        )
    start_index = strand_order.index(start)
    strand_order = strand_order[start_index:] + strand_order[:start_index]

    # Makes nested list, in which each sub-list records the residues in one strand
    strand_res = []
    for strand in strand_order:
        strand_df = domain_df[
            domain_df['domain_strand_ids']=='{}_strand_{}'.format(domain_id, strand)
        ].reset_index(drop=True)
        res_in_strand_n = strand_df['res_ids'].tolist()
        strand_res.append(res_in_strand_n)
    strand_res.append(copy.deepcopy(strand_res[0]))  # Shear calculation
    # measures shift in position of first strand

    # Writes array of interacting residue pairs
    coords = np.full([len(strand_res), 200], 'nan', dtype='object')
    for i, strand in enumerate(strand_res):
        if i == 0:
            for j, res_id in enumerate(strand):
                r = i
                c = 100 + j
                coords[r][c] = res_id
        elif i != 0:
            bonded_res = []
            bonded_res_coords = []

            for j, res_id in enumerate(strand):
                bp_index = domain_df['res_ids'].tolist().index(res_id)
                partner_1 = domain_df['bridge_pairs'][bp_index][0]
                partner_2 = domain_df['bridge_pairs'][bp_index][1]
                partners = [partner_1, partner_2]

                for partner in partners:
                    if (partner != '') and (partner in list(coords[i-1])):
                        r = i
                        c_vals = np.argwhere(coords[i-1] == partner)[0]
                        if len(c_vals) > 1:
                            raise Exception(
                                'Res {} listed more than once'.format(res_id)
                            )
                        c = c_vals[0]
                        coords[r][c] = res_id
                        bonded_res.append(res_id)
                        bonded_res_coords.append(c)
                        break

            # Adds residues in the strand that have so far not been added to the
            # array because they do not form bridge pairs with the preceding strand
            if len(bonded_res) == 0:
                raise Exception(
                    'Strand {} not found to interact with strand {} in domain '
                    '{}'.format(strand_order[i], strand_order[i-1], domain_id)
                )

            start_res = []
            end_res = []
            bonded_res_num = [
                int(''.join([n for n in res_id if n in '1234567890.-']))
                for res_id in bonded_res
            ]
            min_res = min(bonded_res_num[0], bonded_res_num[-1])
            max_res = max(bonded_res_num[0], bonded_res_num[-1])

            for res_id in strand:
                res_num = int(''.join([n for n in res_id if n in '1234567890.-']))
                if not res_id in bonded_res:
                    if res_num < min_res:
                        start_res.append(res_id)
                    elif res_num > max_res:
                        end_res.append(res_id)

            if all(bonded_res_coords[i] <= bonded_res_coords[i+1] for i in range(len(bonded_res_coords)-1)):
                if bonded_res_num[0] < bonded_res_num[-1]:
                    # strand = ['A10', 'A11', 'A12'], coords = [20, 21, 22]
                    # start_res = ['A8', 'A9'] => coords [18, 19]
                    # end_res = ['A13', 'A14'] => coords [23, 24]
                    min_coord = bonded_res_coords[0]
                    max_coord = bonded_res_coords[-1]
                    for j, res_id in enumerate(end_res):
                        r = i
                        c = max_coord + j + 1
                        coords[r][c] = res_id
                    for j, res_id in enumerate(start_res[::-1]):
                        r = i
                        c = min_coord - j - 1
                        coords[r][c] = res_id
                elif bonded_res_num[0] > bonded_res_num[-1]:
                    # strand = ['A12', 'A11', 'A10'], coords = [20, 21,22]
                    # start_res = ['A9', 'A8'] => coords [23, 24]
                    # end_res = ['A14', 'A13'] => coords [18, 19]
                    min_coord = bonded_res_coords[0]
                    max_coord = bonded_res_coords[-1]
                    for j, res_id in enumerate(start_res):
                        r = i
                        c = max_coord + j + 1
                        coords[r][c] = res_id
                    for j, res_id in enumerate(end_res[::-1]):
                        r = i
                        c = min_coord - j - 1
                        coords[r][c] = res_id
                else:
                    raise Exception(
                        'Strand {} in domain {} contains 1 residue'.format(strand, domain_id)
                    )
            elif all(
                bonded_res_coords[i] >= bonded_res_coords[i+1] for i in range(len(bonded_res_coords)-1)
            ):
                if bonded_res_num[0] < bonded_res_num[-1]:
                    # strand = ['A10', 'A11', 'A12'], coords = [22, 21, 20]
                    # start_res = ['A8', 'A9'] => coords [24, 23]
                    # end_res = ['A13', 'A14'] => coords [19, 18]
                    min_coord = bonded_res_coords[-1]
                    max_coord = bonded_res_coords[0]
                    for j, res_id in enumerate(start_res[::-1]):
                        r = i
                        c = max_coord + j + 1
                        coords[r][c] = res_id
                    for j, res_id in enumerate(end_res):
                        r = i
                        c = min_coord - j - 1
                        coords[r][c] = res_id
                elif bonded_res_num[0] > bonded_res_num[-1]:
                    # strand = ['A12', 'A11', 'A10'], coords = [22, 21, 20]
                    # start_res = ['A9', 'A8'] => coords [19, 18]
                    # end_res = ['A14', 'A13'] => coords [24, 23]
                    min_coord = bonded_res_coords[-1]
                    max_coord = bonded_res_coords[0]
                    for j, res_id in enumerate(end_res[::-1]):
                        r = i
                        c = max_coord + j + 1
                        coords[r][c] = res_id
                    for j, res_id in enumerate(start_res):
                        r = i
                        c = min_coord - j - 1
                        coords[r][c] = res_id
                else:
                    raise Exception(
                        'Strand {} in domain {} contains 1 residue'.format(strand, domain_id)
                    )
            else:
                raise Exception('Order of residues in strand {} domain {} not '
                                'recognised'.format(strand_order[i], domain_id))

    if not os.path.isdir('Shear_calc'):
        os.mkdir('Shear_calc')
    np.savetxt(
        'Shear_calc/{}_shear.csv'.format(domain_id), coords, fmt='%s', delimiter=','
    )

    # Calculates shear by working out difference in column number between
    # repeated residue numbers of strand X in first and last rows of the array
    keep_cols = []
    for c in range(coords.shape[1]):
        if set(coords[:,c]) == {''}:
            pass
        else:
            keep_cols.append(c)
    coords = coords[:,keep_cols]

    coord_diffs = []
    for c1, coord in np.ndenumerate(coords[0]):
        c1 = c1[0]
        if coord != '':
            c2 = np.argwhere(coords[-1] == coord)[0]
            if len(c2) > 1:
                raise Exception('Residue {} listed more than once in {}'.format(coord, domain_id))
            else:
                c2 = c2[0]
            coord_diff = np.abs(c2 - c1)
            coord_diffs.append(coord_diff)
    shear = stats.mode(np.array(coord_diffs))[0][0]
    if shear % 2 == 1:
        print('WARNING: Odd shear number ({}) recorded for {}'.format(shear, domain_id))

    # Appends shear to domain_df
    shear_df = pd.DataFrame(OrderedDict({'shear_number': [shear]*domain_df.shape[1]}))
    upd_domain_df = copy.deepcopy(domain_df).reset_index(drop=True)
    upd_domain_df = pd.concat([upd_domain_df, shear_df], axis=1)

    return upd_domain_df


class calc_twist_bend_shear(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    def find_strand_geometry(self, sec_struct_dfs_dict):
        """
        Calculates beta-strand bend and twist, plus barrel shear
        """

        shear_numbers = OrderedDict({})
        sec_struct_dfs_dict_copy = copy.deepcopy(sec_struct_dfs_dict)

        for domain_id, dssp_df in sec_struct_dfs_dict_copy.items():
            # Calculates bend and twist angles (per-residue)
            bend_angles_list = find_strand_bend(domain_id, dssp_df)
            geom_df = pd.DataFrame(OrderedDict({'BEND': bend_angles_list}))
            upd_dssp_df = pd.concat([dssp_df, geom_df], axis=1).reset_index(drop=True)
            sec_struct_dfs_dict[domain_id] = upd_dssp_df

        return sec_struct_dfs_dict
