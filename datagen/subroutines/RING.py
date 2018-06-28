
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

    def run_RING(self, sec_struct_dfs_dict, domain_sheets_dict):
        # Runs RING on the parent structure (biological assembly / asymmetric
        # unit, as specified by the user)
        if not os.path.isdir('ring'):
            os.mkdir('ring')
        for domain_id in list(sec_struct_dfs_dict.keys()):
            if not os.path.isfile('ring/{}.ring'.format(domain_id[0:4])):
                print('Running RING for {}'.format(domain_id[0:4]))
                os.system(
                    '/bin/ring/bin/Ring -i Parent_assemblies/{}.pdb --all -n lollipop -g 1 > ring/{}.ring'.format(
                        domain_id[0:4], domain_id[0:4]
                    )
                )

    def parse_RING_output(self, sec_struct_dfs_dict, domain_sheets_dict):
        # Appends residue interaction network information to dssp_df.
        unprocessed_list = []

        for domain_id in list(sec_struct_dfs_dict.keys()):
            print('Identifying interacting residues for {}'.format(domain_id))

            dssp_df = sec_struct_dfs_dict[domain_id]

            interactions = {'VDW_MC_MC': {},
                            'HBOND_MC_MC': {},
                            'IONIC_MC_MC': {},
                            'SSBOND_MC_MC': {},
                            'PIPISTACK_MC_MC': {},
                            'PIPISTACK_MC_MC_P': {},
                            'PIPISTACK_MC_MC_N': {},
                            'PIPISTACK_MC_MC_L': {},
                            'PIPISTACK_MC_MC_T': {},
                            'PICATION_MC_MC': {},
                            'VDW_SC_MC': {},
                            'HBOND_SC_MC': {},
                            'IONIC_SC_MC': {},
                            'SSBOND_SC_MC': {},
                            'PIPISTACK_SC_MC': {},
                            'PIPISTACK_SC_MC_P': {},
                            'PIPISTACK_SC_MC_N': {},
                            'PIPISTACK_SC_MC_L': {},
                            'PIPISTACK_SC_MC_T': {},
                            'PICATION_SC_MC': {},
                            'VDW_MC_SC': {},
                            'HBOND_MC_SC': {},
                            'IONIC_MC_SC': {},
                            'SSBOND_MC_SC': {},
                            'PIPISTACK_MC_SC': {},
                            'PIPISTACK_MC_SC_P': {},
                            'PIPISTACK_MC_SC_N': {},
                            'PIPISTACK_MC_SC_L': {},
                            'PIPISTACK_MC_SC_T': {},
                            'PICATION_MC_SC': {},
                            'VDW_SC_SC': {},
                            'HBOND_SC_SC': {},
                            'IONIC_SC_SC': {},
                            'SSBOND_SC_SC': {},
                            'PIPISTACK_SC_SC': {},
                            'PIPISTACK_SC_SC_P': {},
                            'PIPISTACK_SC_SC_N': {},
                            'PIPISTACK_SC_SC_L': {},
                            'PIPISTACK_SC_SC_T': {},
                            'PICATION_SC_SC': {}}

            for res_id in set(dssp_df['RES_ID'].tolist()):
                for interaction_type in list(interactions.keys()):
                    interactions[interaction_type][res_id] = []

            with open('ring/{}.ring'.format(domain_id[0:4]), 'r') as ring_file:
                ring_output = ring_file.readlines()
            ring_output = ''.join(ring_output)

            if not 'NodeId' in ring_output:
                print('ERROR: {} unable to be processed by RING'.format(domain_id[0:4]))
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

                    mc_or_sc = line.split()[1].split(':')[1]
                    aa_1_label = mc_or_sc.split('_')[0]
                    aa_2_label = mc_or_sc.split('_')[1]

                    interaction_type = line.split()[1].split(':')[0]
                    aa_1_interaction_type = '{}_{}_{}'.format(
                        interaction_type, aa_1_label, aa_2_label
                    )
                    aa_2_interaction_type = '{}_{}_{}'.format(
                        interaction_type, aa_2_label, aa_1_label
                    )

                    # Records interaction for both residues involved
                    if (
                        aa_1_interaction_type in list(interactions.keys())
                        and aa_2_interaction_type in list(interactions.keys())
                    ):
                        if (
                            aa_1 in list(interactions[aa_1_interaction_type].keys())
                            and aa_2 not in interactions[aa_1_interaction_type][aa_1]
                        ):
                            interactions[aa_1_interaction_type][aa_1].append(aa_2)
                        if (
                            aa_2 in list(interactions[aa_2_interaction_type].keys())
                            and aa_1 not in interactions[aa_2_interaction_type][aa_2]
                        ):
                            interactions[aa_2_interaction_type][aa_2].append(aa_1)

                    # Records the orientation of pi-pi stacking interactions
                    if interaction_type == 'PIPISTACK':
                        orientation = line.split()[-2]
                        if orientation in ['T-EF', 'T-FE']:
                            orientation = 'T'
                        aa_1_interaction_orientation = '{}_{}'.format(
                            aa_1_interaction_type, orientation
                        )
                        if aa_1 in list(interactions[aa_1_interaction_orientation].keys()):
                            interactions[aa_1_interaction_orientation][aa_1].append(aa_2)

                        aa_2_interaction_orientation = '{}_{}'.format(
                            aa_2_interaction_type, orientation
                        )
                        if aa_2 in list(interactions[aa_2_interaction_orientation].keys()):
                            interactions[aa_2_interaction_orientation][aa_2].append(aa_1)

                ring_df_dict = OrderedDict({'VDW_MC_MC': ['']*dssp_df.shape[0],
                                            'HBOND_MC_MC': ['']*dssp_df.shape[0],
                                            'IONIC_MC_MC': ['']*dssp_df.shape[0],
                                            'SSBOND_MC_MC': ['']*dssp_df.shape[0],
                                            'PIPISTACK_MC_MC': ['']*dssp_df.shape[0],
                                            'PIPISTACK_MC_MC_P': ['']*dssp_df.shape[0],
                                            'PIPISTACK_MC_MC_N': ['']*dssp_df.shape[0],
                                            'PIPISTACK_MC_MC_L': ['']*dssp_df.shape[0],
                                            'PIPISTACK_MC_MC_T': ['']*dssp_df.shape[0],
                                            'PICATION_MC_MC': ['']*dssp_df.shape[0],
                                            'VDW_SC_MC': ['']*dssp_df.shape[0],
                                            'HBOND_SC_MC': ['']*dssp_df.shape[0],
                                            'IONIC_SC_MC': ['']*dssp_df.shape[0],
                                            'SSBOND_SC_MC': ['']*dssp_df.shape[0],
                                            'PIPISTACK_SC_MC': ['']*dssp_df.shape[0],
                                            'PIPISTACK_SC_MC_P': ['']*dssp_df.shape[0],
                                            'PIPISTACK_SC_MC_N': ['']*dssp_df.shape[0],
                                            'PIPISTACK_SC_MC_L': ['']*dssp_df.shape[0],
                                            'PIPISTACK_SC_MC_T': ['']*dssp_df.shape[0],
                                            'PICATION_SC_MC': ['']*dssp_df.shape[0],
                                            'VDW_MC_SC': ['']*dssp_df.shape[0],
                                            'HBOND_MC_SC': ['']*dssp_df.shape[0],
                                            'IONIC_MC_SC': ['']*dssp_df.shape[0],
                                            'SSBOND_MC_SC': ['']*dssp_df.shape[0],
                                            'PIPISTACK_MC_SC': ['']*dssp_df.shape[0],
                                            'PIPISTACK_MC_SC_P': ['']*dssp_df.shape[0],
                                            'PIPISTACK_MC_SC_N': ['']*dssp_df.shape[0],
                                            'PIPISTACK_MC_SC_L': ['']*dssp_df.shape[0],
                                            'PIPISTACK_MC_SC_T': ['']*dssp_df.shape[0],
                                            'PICATION_MC_SC': ['']*dssp_df.shape[0],
                                            'VDW_SC_SC': ['']*dssp_df.shape[0],
                                            'HBOND_SC_SC': ['']*dssp_df.shape[0],
                                            'IONIC_SC_SC': ['']*dssp_df.shape[0],
                                            'SSBOND_SC_SC': ['']*dssp_df.shape[0],
                                            'PIPISTACK_SC_SC': ['']*dssp_df.shape[0],
                                            'PIPISTACK_SC_SC_P': ['']*dssp_df.shape[0],
                                            'PIPISTACK_SC_SC_N': ['']*dssp_df.shape[0],
                                            'PIPISTACK_SC_SC_L': ['']*dssp_df.shape[0],
                                            'PIPISTACK_SC_SC_T': ['']*dssp_df.shape[0],
                                            'PICATION_SC_SC': ['']*dssp_df.shape[0]})

                for row in range(dssp_df.shape[0]):
                    res_id = dssp_df['RES_ID'][row]
                    for interaction_type in list(ring_df_dict.keys()):
                        if (
                            res_id in list(interactions[interaction_type].keys())
                            and dssp_df['ATMNAME'][row] == 'CA'
                        ):
                            ring_df_dict[interaction_type][row] = interactions[interaction_type][res_id]
                ring_df = pd.DataFrame(ring_df_dict)
                dssp_df = pd.concat([dssp_df, ring_df], axis=1)
                sec_struct_dfs_dict[domain_id] = dssp_df

        # Writes PDB accession codes that could not be processed to output file
        with open('Unprocessed_domains.txt', 'a') as unprocessed_file:
            unprocessed_file.write('\n\nError in running RING:\n')
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
        # (via their main-chain atoms) with one another that has the highest
        # average buried surface area
        unprocessed_list = []

        for domain_id in list(sec_struct_dfs_dict.keys()):
            print('Identifying interior- and exterior-facing residues in '
                  '{}'.format(domain_id))
            dssp_df = sec_struct_dfs_dict[domain_id]
            dssp_ca_df = dssp_df[dssp_df['ATMNAME'] == 'CA']
            dssp_ca_df = dssp_ca_df.reset_index(drop=True)

            # Groups residues in van der Waals contact together
            groups = OrderedDict()
            count = 0
            res_id_list = dssp_ca_df['RES_ID'].tolist()
            for row in range(dssp_ca_df.shape[0]):
                res_id = dssp_ca_df['RES_ID'][row]

                in_current_group = False
                for group in list(groups.keys()):
                    if res_id in groups[group]:
                        in_current_group = True
                        new_group = list(set(groups[group]))
                        for new_res_id in dssp_ca_df['VDW_SC_SC'][row]:
                            if new_res_id in res_id_list and not new_res_id in new_group:
                                new_group.append(new_res_id)
                        groups[group] = new_group
                if in_current_group is False:
                    count += 1  # Generates new group key label

                    new_group = [res_id]  # Must be list!
                    for new_res_id in dssp_ca_df['VDW_SC_SC'][row]:
                        if new_res_id in res_id_list and not new_res_id in new_group:
                            new_group.append(new_res_id)
                    groups[str(count)] = new_group

            # Merges overlapping groups
            merged_group_lists = []
            for group_list in list(groups.values()):
                merged_group_lists += group_list

            while sorted(set(merged_group_lists)) != sorted(merged_group_lists):
                group_1_key = list(groups.keys())[0]
                group_1 = groups[group_1_key]
                overlap = False
                for group in list(groups.keys()):
                    if (
                        group != list(groups.keys())[0]
                        and bool(set(group_1) & set(groups[group]))
                    ):
                        merged_group = group_1 + groups[group]
                        groups[group] = list(set(merged_group))
                        del groups[group_1_key]
                        overlap = True
                        break

                if overlap is False:
                    groups.move_to_end(group_1_key)

                merged_group_lists = []
                for group_list in list(groups.values()):
                    merged_group_lists += group_list

            # Checks that correct number of surfaces have been identified (2x
            # exterior, 1x interior). BE CAREFUL WITH THIS CHECK IF EXPAND CODE
            # TO LOOK AT A WIDER RANGE OF STRUCTURES!
            # PROBLEM: At the moment all residues are considered to be in van
            # der waals contact with one another => 1 surface. How to split into three?
            if len(list(groups.keys())) < 3:
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
                num_res_in_groups = {}
                for group in list(groups.keys()):
                    num_res_in_groups[len(groups[group])] = group

                interior_group = num_res_in_groups[max(list(num_res_in_groups.keys()))]
                int_ext_dict = {}
                for group in groups:
                    for res_id in groups[group]:
                        if group == interior_group:
                            int_ext_dict[res_id] = 'interior'
                        else:
                            int_ext_dict[res_id] = 'exterior'
                print(int_ext_dict)

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
            unprocessed_file.write('\n\nFailed to identify interior and '
                                   'exterior surfaces of sandwich:\n')
            for domain_id in set(unprocessed_list):
                unprocessed_file.write('{}\n'.format(domain_id))

        sec_struct_dfs_dict = OrderedDict(
            {key: value for key, value in sec_struct_dfs_dict.items() if value is not None}
        )
        domain_sheets_dict = OrderedDict(
            {key: value for key, value in domain_sheets_dict.items() if value is not None}
        )

        return sec_struct_dfs_dict, domain_sheets_dict
