
import pandas as pd
import networkx as nx
import numpy as np
from collections import OrderedDict
if __name__ == 'subroutines.OPM':
    from subroutines.run_stages import run_stages
else:
    from datagen.subroutines.run_stages import run_stages


class extract_barrel_info_from_OPM(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    def parse_opm(self, orig_dir):
        # Extracts strand tilt and TM information from the OPM database
        # (docs/OPM_TM_subunits.txt)
        pdb_codes = []
        chains = []
        tilt_angles = []
        tm_lists = []
        tm_segment_lists = []

        print('Creating datframe of OPM database information')

        with open('{}/docs/OPM_TM_subunits.txt'.format(orig_dir), 'r') as opm_file:
            for line in opm_file:
                line_segments = line.split('-')

                pdb_codes.append(line_segments[0][0:4])

                chain = line_segments[0][4:].strip()
                chains.append(chain)

                tilt_angle = line_segments[1].replace('Tilt:', '')
                tilt_angle = tilt_angle.replace('°', '')
                tilt_angles.append(tilt_angle.strip())

                tm_segments = '-'.join(line_segments[2:])
                tm_segments = tm_segments.replace('Segments:', '')
                tm_segments = tm_segments.split(',')
                tm_segment_lists.append(tm_segments)

                tm_residues = []
                for segment in tm_segments:
                    tm_range = ''
                    start = False
                    stop = False
                    for char in list(segment):
                        if char == ')':
                            stop = True

                        if start is True and stop is False:
                            tm_range += char

                        if char == '(':
                            start = True

                    res_min = int(tm_range.split('-')[0])
                    res_max = int(tm_range.split('-')[1])
                    tm_residues += [chain+str(num) for num in
                                    range(res_min, res_max+1)]  # No insertion
                    # code available in OPM_TM_subunits.txt
                tm_lists.append(tm_residues)

        opm_df = pd.DataFrame({'PDB_CODE': pdb_codes,
                               'CHAIN': chains,
                               'TILT_ANGLE': tilt_angles,
                               'TM_RANGE': tm_lists,
                               'SEGMENTS': tm_segment_lists})

        return opm_df

    def find_strand_tilt(self, sec_struct_dfs_dict, opm_df):
        # Determines strand tilt
        pdb_codes_list = opm_df['PDB_CODE'].tolist()
        tilt_angles = OrderedDict()

        for domain_id in list(sec_struct_dfs_dict.keys()):
            print('Calculating tilt angle for {}'.format(domain_id))

            pdb_code = domain_id[0:4]

            if pdb_code in pdb_codes_list:
                index = pdb_codes_list.index(pdb_code)
                tilt_angle = opm_df['TILT_ANGLE'][index]
                tilt_angles[domain_id] = tilt_angle
            else:
                tilt_angles[domain_id] = 'Undefined'

        return tilt_angles


class calculate_barrel_geometry(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    def find_barrel_strand_number(self, sec_struct_dfs_dict):
        # Calculates number of strands in barrel
        strand_numbers = OrderedDict()

        for domain_id, dssp_df in sec_struct_dfs_dict.items():
            print('Calculating number of strands in {}'.format(domain_id))
            strands = [strand for strand in set(dssp_df['STRAND_NUM'].tolist())
                       if strand != '']
            strand_count = len(strands)
            strand_numbers[domain_id] = strand_count

        return strand_numbers

    def find_barrel_shear_number(self, sec_struct_dfs_dict, domain_sheets_dict):
        # Calculates barrel shear number
        # CURRENTLY DOESNT WORK (ENCOUNTERS ERRORS WITH INSERTED RESIDUES ETC.)
        unprocessed_list = []
        shear_numbers = OrderedDict()

        for domain_id, dssp_df in sec_struct_dfs_dict.items():
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
                continue

            cycles = [strand for cycle in nx.cycle_basis(G) for strand in cycle]
            if len(set(cycles)) != len(nx.nodes(G)):
                unprocessed_list.append(domain_id)
                continue

            strands = nx.cycle_basis(G)[0]  # Finds first complete cycle
            res_array = np.zeros((len(strands)+1, 200))

            processed_strands = []
            strand_1 = strands[0]
            processed_strands.append(strand_1)
            strand_1_df = dssp_df[dssp_df['STRAND_NUM'] == strand_1]
            dssp_num_1 = strand_1_df['DSSP_NUM'].tolist()
            dssp_num_1 = [int(res) for res in dssp_num_1]
            dssp_num_1 = [res for res in range(min(dssp_num_1), max(dssp_num_1)+1)]
            for index, res in enumerate(dssp_num_1):
                res_array[0, 80+index] = res
            res_list_1 = res_array[0:1, ].tolist()
            res_list_1 = [res for res_sub_list in res_list_1 for res in res_sub_list]
            res_list = res_array[0:1, ].tolist()
            res_list = [res for res_sub_list in res_list for res in res_sub_list]

            current_strand = strand_1
            next_strand_df = strand_1_df
            count = 0
            diff_dict = {}
            while any(x not in processed_strands for x in G.neighbors(current_strand)):
                next_strand = G.neighbors(current_strand)[0]
                if next_strand in processed_strands:
                    next_strand = G.neighbors(current_strand)[1]
                processed_strands.append(next_strand)
                count = count + 1
                if count == 2:
                    processed_strands.remove(strand_1)

                prev_dssp_num = next_strand_df['DSSP_NUM'].tolist()
                prev_dssp_num = [int(num) for num in prev_dssp_num]

                next_strand_df = dssp_df[dssp_df['STRAND_NUM'] == next_strand]
                dssp_num = next_strand_df['DSSP_NUM'].tolist()
                dssp_num = [int(num) for num in dssp_num]
                h_bonds = next_strand_df['H-BONDS'].tolist()
                pair_index = ''
                for h_bonds_index, pair in enumerate(h_bonds):
                    if float(pair[0]) != 0.0 and float(pair[0]) in res_list:
                        pair_index = 0
                    elif float(pair[1]) != 0.0 and float(pair[1]) in res_list:
                        pair_index = 1

                dssp_num_dict = {int(pair[pair_index]): dssp_num[index] for
                                 index, pair in enumerate(h_bonds) if
                                 pair[pair_index] != '0' and int(pair[pair_index])
                                 in prev_dssp_num}
                h_bonds = [int(pair[pair_index]) for pair in h_bonds if
                           pair[pair_index] != '0' and int(pair[pair_index])
                           in prev_dssp_num]
                h_bonds = [res for res in range(min(h_bonds), max(h_bonds)+1)]

                listed_res = []
                for h_bonds_index, res in enumerate(h_bonds):
                    if (res != 0
                            and float(res) in res_list
                            and res in list(dssp_num_dict.keys())
                        ):
                        listed_res.append(int(dssp_num_dict[res]))
                        res_list_index = res_list.index(float(res))
                        res_array[count, res_list_index] = float(dssp_num_dict[res])

                lower_res = [res for res in range(min(dssp_num), min(listed_res))]
                higher_res = [res for res in range(max(listed_res)+1, max(dssp_num)+1)]

                res_list = res_array[count:count+1, ].tolist()
                res_list = [res for res_sub_list in res_list for res in res_sub_list]
                res_list_values = [num for num in res_list if num != 0.0]

                """
                diff = len(range(min(listed_res), max(listed_res)+1)) - len(listed_res)
                if diff >= 1:
                    missing_res = []
                    for num in range(min(listed_res), max(listed_res)+1):
                        if num not in listed_res:
                            missing_res.append(num)

                    for num in missing_res:
                        index_1 = res_list.index(min(listed_res, key=lambda x: abs(x-num)))
                        listed_res.remove(min(listed_res, key=lambda x: abs(x-num)))
                        index_2 = res_list.index(min(listed_res, key=lambda x: abs(x-num)))
                        indices = sorted([index_1, index_2])
                        diff_dict[diff] = indices
                """

                if len(lower_res) > 0:  # Error with overhang placement if len == 1, need to fix this
                    lower_index = res_list.index(float(min(listed_res)))
                else:
                    lower_index = 0

                if len(higher_res) > 0:
                    higher_index = res_list.index(float(max(listed_res)))
                else:
                    higher_index = 1000

                if res_list_values[0] < res_list_values[-1]:
                    if lower_index != 0:
                        for index, res in enumerate(lower_res):
                            res_list_index = (lower_index - len(lower_res)) + index
                            res_array[count, res_list_index] = res
                    if higher_index != 1000:
                        for index, res in enumerate(higher_res):
                            res_list_index = higher_index + index + 1
                            res_array[count, res_list_index] = res
                elif res_list_values[0] > res_list_values[-1]:
                    lower_res.reverse()
                    higher_res.reverse()
                    if lower_index != 0:
                        for index, res in enumerate(higher_res):
                            res_list_index = lower_index + index + 1
                            res_array[count, res_list_index] = res
                    if higher_index != 1000:
                        for index, res in enumerate(higher_res):
                            res_list_index = (higher_index - len(higher_res)) + index
                            res_array[count, res_list_index] = res

                res_list = res_array[count:count+1, ].tolist()
                res_list = [res for res_sub_list in res_list for res in res_sub_list]

                """
                for index, res in enumerate(dssp_num):
                    res_list = res_array[count:count+1, ].tolist()
                    res_list = [res for res_sub_list in res_list for res in res_sub_list]
                    if (index not in [0, len(dssp_num)-1]
                            and res not in res_list
                        ):
                        print(res)
                        print(res_list)

                        index_count = 1
                        while True:
                            try:
                                index_1 = res_list.index(float(dssp_num[index-index_count]))
                            except ValueError:
                                index_count += 1
                            else:
                                break

                        index_count = 1
                        while True:
                            try:
                                index_2 = res_list.index(float(dssp_num[index+index_count]))
                            except ValueError:
                                index_count += 1
                            else:
                                break

                        new_index = max([index_1, index_2])
                        np.insert(res_array, new_index, 0, axis=1)
                        res_array[count, new_index] = float(res)
                """

                """
                res_range = sorted([num for num in res_list if num != 0.0])
                for index, res in enumerate(res_range):
                    if (index != (len(res_range)-1)
                            and (res+1) != res_range[index+1]
                            ):
                        diff = res_range[index+1] - res
                        if diff % 2 == 0:
                            shear = shear + 1
                """

                current_strand = next_strand

            diff_dict = OrderedDict(sorted(diff_dict.items(), reverse=True))

            # Residues are now aligned correctly, but are still small errors in
            # shear number calculation due to skipped residues
            low_index_row_1 = res_list_1.index(dssp_num_1[0])
            low_index_row_n = res_list.index(dssp_num_1[0])

            if low_index_row_1 > low_index_row_n:
                shear = low_index_row_1 - low_index_row_n
                for diff in diff_dict:
                    if (low_index_row_n < diff_dict[diff][0]
                            and diff_dict[diff][1] < low_index_row_1
                        ):
                        shear += diff
                        break
            elif low_index_row_1 < low_index_row_n:
                high_index_row_1 = res_list_1.index(listed_res[-1])
                high_index_row_n = res_list.index(listed_res[-1])
                shear = high_index_row_n - high_index_row_1
                for diff in diff_dict:
                    if (high_index_row_1 < diff_dict[diff][0]
                            and diff_dict[diff][1] < high_index_row_n
                        ):
                        shear += diff
                        break

            shear_numbers[domain_id] = shear
            print(domain_id, shear)

        with open('Unprocessed_domains.txt', 'a') as unprocessed_file:
            unprocessed_file.write('\n\nBarrel formed from discontinuous sheet:\n')
            for domain_id in unprocessed_list:
                unprocessed_file.write('{}\n'.format(domain_id))

        return shear_numbers
