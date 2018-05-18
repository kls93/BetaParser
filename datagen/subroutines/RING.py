

import os
import pandas as pd
from collections import OrderedDict
if __name__ == 'subroutines.RING':
    from subroutines.run_stages import run_stages
else:
    from datagen.subroutines.run_stages import run_stages


class calculate_residue_interaction_network(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    def parse_RING_output(self, sec_struct_dfs_dict, domain_sheets_dict):
        # Appends residue interaction network information to dssp_df.
        # **TODO** Currently treats main chain / main chain, main chain / side
        # chain and side chain / side chain interactions identically, might
        # want to update this in the future.
        unprocessed_list = []

        for domain_id in list(sec_struct_dfs_dict.keys()):
            print('Identifying interacting residues for {}'.format(domain_id))

            dssp_df = sec_struct_dfs_dict[domain_id]

            interactions = {'VDW': {},
                            'HBOND': {},
                            'IONIC': {},
                            'SSBOND': {},
                            'PIPISTACK': {},
                            'PICATION': {}}

            for res_id in set(dssp_df['RES_ID'].tolist()):
                for interaction_type in list(interactions.keys()):
                    interactions[interaction_type][res_id] = []

            with open('{}{}/{}.ring'.format(
                self.ring_database, domain_id[1:3], domain_id[0:4]
            ), 'r') as ring_file:
                ring_output = ring_file.readlines()
            ring_output = ''.join(ring_output)

            if not 'NodeId' in ring_output:
                print('ERROR: {}.ring unable to be analysed by RING'.format(domain_id))
                unprocessed_list.append(domain_id)
                sec_struct_dfs_dict[domain_id] = None
                sheet_ids = [sheet_id for sheet_id in list(domain_sheets_dict.keys())
                             if sheet_id.startswith(domain_id)]
                for sheet_id in sheet_ids:
                    domain_sheets_dict[sheet_id] = None
            else:
                ring_output = ring_output.replace('\t', ' ')
                ring_output = ring_output.split('NodeId')[2]
                ring_output = ring_output.split('\n')[1:]
                ring_output = [line for line in ring_output if line != '']
                for line in ring_output:
                    aa_1 = line.split()[0].split(':')
                    aa_1 = (aa_1[0].replace('_', '') + aa_1[1].replace('_', '')
                            + aa_1[2].replace('_', ''))
                    aa_2 = line.split()[2].split(':')
                    aa_2 = (aa_2[0].replace('_', '') + aa_2[1].replace('_', '')
                            + aa_2[2].replace('_', ''))

                    interaction_type = line.split()[1].split(':')[0]
                    if interaction_type in list(interactions.keys()):
                        if (
                            aa_1 in list(interactions[interaction_type].keys())
                            and aa_2 not in interactions[interaction_type][aa_1]
                        ):
                            interactions[interaction_type][aa_1].append(aa_2)
                        if (
                            aa_2 in list(interactions[interaction_type].keys())
                            and aa_1 not in interactions[interaction_type][aa_2]
                        ):
                            interactions[interaction_type][aa_2].append(aa_1)

                ring_df_dict = OrderedDict({'VDW': ['']*dssp_df.shape[0],
                                            'HBOND': ['']*dssp_df.shape[0],
                                            'IONIC': ['']*dssp_df.shape[0],
                                            'SSBOND': ['']*dssp_df.shape[0],
                                            'PIPISTACK': ['']*dssp_df.shape[0],
                                            'PICATION': ['']*dssp_df.shape[0]})
                for row in range(dssp_df.shape[0]):
                    res_id = dssp_df['RES_ID'][row]
                    for interaction_type in list(ring_df_dict.keys()):
                        if (
                            res_id in interactions[interaction_type]
                            and dssp_df['ATMNAME'][row] == 'CA'
                        ):
                            ring_df_dict[interaction_type][row] = interactions[interaction_type][res_id]
                ring_df = pd.DataFrame(ring_df_dict)
                dssp_df = pd.concat([dssp_df, ring_df], axis=1)
                sec_struct_dfs_dict[domain_id] = dssp_df

        # Writes PDB accession codes that could not be processed to output file
        with open('Unprocessed_domains.txt', 'a') as unprocessed_file:
            unprocessed_file.write('\n\nERROR in running RING:\n')
            for domain_id in set(unprocessed_list):
                unprocessed_file.write('{}\n'.format(domain_id))

        sec_struct_dfs_dict = OrderedDict(
            {key: value for key, value in sec_struct_dfs_dict.items() if value is not None}
        )
        domain_sheets_dict = OrderedDict(
            {key: value for key, value in domain_sheets_dict.items() if value is not None}
        )

        return sec_struct_dfs_dict, domain_sheets_dict

    def identify_int_ext_sandwich(self, sec_struct_dfs_dict, domain_sheets_dict):
        # Identifies residues facing towards the interior of the sandwich as
        # the group of residues in (direct or indirect) van der Waals contact
        # with one another that has the highest average buried surface area
        unprocessed_list = []

        for domain_id in list(sec_struct_dfs_dict.keys()):
            print('Identifying interior- and exterior-facing residues in '
                  '{}'.format(domain_id))
            dssp_df = sec_struct_dfs_dict[domain_id]

            # Groups residues in van der Waals contact together
            groups = OrderedDict()
            count = 0
            for row in range(dssp_df.shape[0]):
                res_id = dssp_df['RES_ID'][row]

                in_current_group = False
                for group in list(groups.keys()):
                    if res_id in groups[group]:
                        in_current_group = True
                        groups[group] += dssp_df['VDW'][row]
                if in_current_group is False:
                    count += 1
                    if dssp_df['VDW'][row] != '':
                        groups[str(count)] = [res_id] + dssp_df['VDW'][row]  # Must be list!

            # Merges overlapping groups
            merged_group_lists = []
            for group_list in list(groups.values()):
                merged_group_lists += group_list

            while sorted(set(merged_group_lists)) != sorted(merged_group_lists):
                group_1 = groups[list(groups.keys())[0]]
                overlap = False
                for group in list(groups.keys()):
                    if (
                        group != list(groups.keys())[0]
                        and bool(set(group_1) & set(groups[group]))
                    ):
                        merged_group = group_1 + groups[group]
                        groups[group] = list(set(merged_group))
                        del groups[list(groups.keys())[0]]
                        overlap = True
                        break

                if overlap is False:
                    groups.move_to_end(group_1)

                merged_group_lists = []
                for group_list in list(groups.values()):
                    merged_group_lists += group_list

            # Checks that correct number of surfaces have been identified (2x
            # exterior, 1x interior). BE CAREFUL WITH THIS CHECK IF EXPAND CODE
            # TO LOOK AT A WIDER RANGE OF STRUCTURES!
            # PROBLEM: At the moment all residues are considered to be in van
            # der waals contact with one another => 1 surface. How to split into three?
            if len(list(groups.keys())) != 3:
                print('ERROR in determining interior- and exterior-facing '
                      'surfaces')
                unprocessed_list.append(domain_id)
                sec_struct_dfs_dict[domain_id] = None
                sheet_ids = [sheet_id for sheet_id in list(domain_sheets_dict.keys())
                             if sheet_id.startswith(domain_id)]
                for sheet_id in sheet_ids:
                    domain_sheets_dict[sheet_id] = None
            else:
                # Finds the group with the largest average buried surface area,
                # which is assumed to be the interior surface
                buried_surface_areas = {}
                dssp_ca_df = dssp_df[dssp_df['ATMNAME'] == 'CA']
                dssp_ca_df = dssp_ca_df.reset_index(drop=True)
                for group in list(groups.keys()):
                    buried_surface_area_sum = 0
                    for res_id in groups[group]:
                        buried_surface_area = dssp_ca_df['BURIED_SURFACE_AREA(%)'].tolist()[
                            dssp_ca_df['RES_ID'].tolist().index(res_id)
                        ]
                        buried_surface_area_sum += buried_surface_area
                    average_buried_surface_area = buried_surface_area_sum / len(groups[group])
                    buried_surface_areas[average_buried_surface_area] = group

                interior_group = buried_surface_areas[max(list(buried_surface_areas.keys()))]
                int_ext_dict = {}
                for group in groups:
                    for res_id in group:
                        if group == interior_group:
                            int_ext_dict[res_id] = 'interior'
                        else:
                            int_ext_dict[res_id] = 'exterior'

                # Updates dataframe with interior / exterior-facing calculation
                # results
                int_ext_list = ['']*dssp_df.shape[0]
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

        # Writes PDB accession codes that could not be processed to output file
        with open('Unprocessed_domains.txt', 'a') as unprocessed_file:
            unprocessed_file.write('\n\nERROR - failed to identify interior '
                                   'and exterior surfaces of sandwich:\n')
            for domain_id in set(unprocessed_list):
                unprocessed_file.write('{}\n'.format(domain_id))

        sec_struct_dfs_dict = OrderedDict(
            {key: value for key, value in sec_struct_dfs_dict.items() if value is not None}
        )
        domain_sheets_dict = OrderedDict(
            {key: value for key, value in domain_sheets_dict.items() if value is not None}
        )

        return sec_struct_dfs_dict, domain_sheets_dict
