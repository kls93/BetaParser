
import os
import shutil
import pandas as pd
from collections import OrderedDict
if __name__ == 'subroutines.DSSP':
    from subroutines.run_stages import run_stages
else:
    from datagen.subroutines.run_stages import run_stages

class filter_dssp_database(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    def copy_files_from_dssp_database(self, cd_hit_domain_df):
        # Generates list of all PDB accession codes to be extracted from the
        # DSSP database
        pdb_list = cd_hit_domain_df['PDB_CODE'].tolist()

        # Copies DSSP files of listed PDB accession codes to current working
        # directory
        unprocessed_list = []
        for pdb_code in pdb_list:
            middle_characters = pdb_code[1:len(pdb_code)-1]
            cwd = os.getcwd()

            try:
                shutil.copy2(
                    '/{}/{}.dssp'.format(self.dssp_database, pdb_code),
                    '{}/DSSP_files/{}.dssp'.format(cwd, pdb_code))

                print('Copied {}.dssp'.format(pdb_code))
            except FileNotFoundError:
                unprocessed_list.append(pdb_code)

        # Writes PDB accession codes that could not be processed to output file
        with open('Unprocessed_domains.txt', 'a') as unprocessed_file:
            unprocessed_file.write('\n\nNot in DSSP database:\n')
            for pdb_code in set(unprocessed_list):
                unprocessed_file.write('{}\n'.format(pdb_code))

        # Filters cd_hit_domain_df_xyz to remove entries not in the DSSP
        # database
        dssp_domain_df = cd_hit_domain_df[~cd_hit_domain_df['PDB_CODE'].isin(unprocessed_list)]
        dssp_domain_df = dssp_domain_df.reset_index(drop=True)

        return dssp_domain_df


class beta_structure_dssp_classification(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    def extract_dssp_file_lines(self, dssp_domain_df):
        unprocessed_list = []

        # Creates dictionary of DSSP file lines for all input PDB accession codes
        dssp_file_lines = []

        for row in range(dssp_domain_df.shape[0]):
            dssp_indv_file_lines = []
            chain_num_list = dssp_domain_df['CHAIN_NUM'][row]

            with open('DSSP_files/{}.dssp'.format(dssp_domain_df['PDB_CODE'][row]), 'r') as dssp_file:
                print('Processing {}.dssp'.format(dssp_domain_df['PDB_CODE'][row]))
                for line in dssp_file:
                    if (line[11:12]+line[5:11].strip()) in chain_num_list:
                        dssp_indv_file_lines.append(line.strip('\n'))
                        index = chain_num_list.index((line[11:12]+line[5:11].strip()))
                        chain_num_list[index] = ''

            chain_num_list = [chain_num for chain_num in chain_num_list if chain_num != '']
            if len(chain_num_list) == 0:
                dssp_indv_file_lines.append('TER'.ljust(136)+'\n')
                dssp_file_lines.append(dssp_indv_file_lines)
            elif len(chain_num_list) > 0:
                unprocessed_list.append(dssp_domain_df['DOMAIN_ID'][row])

        # Writes PDB accession codes that could not be processed to output file
        with open('Unprocessed_domains.txt', 'a') as unprocessed_file:
            unprocessed_file.write('\n\nCoordinates missing from DSSP file:\n')
            for domain_id in set(unprocessed_list):
                unprocessed_file.write('{}\n'.format(domain_id))

        # Filters dssp_domain_df to remove entries whose DSSP files have
        # coordinates missing
        dssp_domain_df = dssp_domain_df[~dssp_domain_df['DOMAIN_ID'].isin(unprocessed_list)]
        dssp_domain_df = dssp_domain_df.reset_index(drop=True)

        dssp_domain_df.to_csv('Final_filtered_dataset.csv')
        dssp_domain_df.to_pickle('Final_filtered_dataset.pkl')

        dssp_residues_dict = OrderedDict()
        for index, sub_list in enumerate(dssp_file_lines):
            dssp_residues_dict[dssp_domain_df['DOMAIN_ID'][index]] = sub_list

        return dssp_residues_dict

    def get_dssp_sec_struct_df(self, dssp_residues_dict, pdb_dfs_dict):
        # Generates dataframe of relevant information in DSSP file
        dssp_dfs_dict = OrderedDict()

        for index_1, domain_id in enumerate(list(dssp_residues_dict.keys())):
            dssp_indv_file_lines = dssp_residues_dict[domain_id]
            pdb_df = pdb_dfs_dict[domain_id]

            print((index_1/len(list(dssp_residues_dict.keys())))*100)
            print('Generating dataframe summarising DSSP info for {}'.format(domain_id))

            res_num = []
            chain = []
            dssp_num = []
            sec_struct_assignment = []
            strand_number = 1
            strand_number_list = []
            sheet_number_list = []
            orientation_list = []
            bridge_pair_list = []
            phi = []
            psi = []

            for index_2, line in enumerate(dssp_indv_file_lines):
                if not line.startswith('TER'):
                    res_num.append(line[5:11].strip())
                    chain.append(line[11:12])
                    dssp_num.append(line[0:5].strip())
                    secondary_structure = line[16:17]
                    phi.append(float(line[91:97]))
                    psi.append(float(line[97:103]))
                    if secondary_structure == 'E':
                        sec_struct_assignment.append(secondary_structure)
                        strand_number_list.append(strand_number)
                        if dssp_indv_file_lines[index_2+1][16:17] != 'E':
                            strand_number = strand_number + 1
                        sheet_number_list.append(line[33:34])
                        ladder_1 = line[23:24]
                        if ladder_1 == ' ':
                            ladder_1 = ''
                        elif ladder_1.isupper():
                            ladder_1 = 'A'
                        elif ladder_1.islower():
                            ladder_1 = 'P'
                        ladder_2 = line[24:25]
                        if ladder_2 == ' ':
                            ladder_2 = ''
                        elif ladder_2.isupper():
                            ladder_2 = 'A'
                        elif ladder_2.islower():
                            ladder_2 = 'P'
                        orientation_list.append([ladder_1,ladder_2])
                        bridge_pair_list.append([line[25:29].strip(),
                                                 line[29:33].strip()])
                    elif secondary_structure != 'E':
                        sec_struct_assignment.append('')
                        strand_number_list.append('')
                        sheet_number_list.append('')
                        orientation_list.append(['', ''])
                        bridge_pair_list.append(['', ''])

            row_num = pdb_df.shape[0]

            dssp_num_extnd_df = ['']*row_num
            sec_struct_assignment_extnd_df = ['']*row_num
            strand_number_list_extnd_df = ['']*row_num
            sheet_number_list_extnd_df = ['']*row_num
            orientation_list_extnd_df = ['']*row_num
            bridge_pair_list_extnd_df = ['']*row_num
            phi_extnd_df = ['']*row_num
            psi_extnd_df = ['']*row_num

            for row in range(row_num):
                if pdb_df['ATMNAME'][row] == 'CA':
                    for index_3, value in enumerate(res_num):
                        if ((str(pdb_df['RESNUM'][row])+pdb_df['INSCODE'][row]) == res_num[index_3]
                            and pdb_df['CHAIN'][row] == chain[index_3]
                            ):
                            dssp_num_extnd_df[row] = dssp_num[index_3]
                            sec_struct_assignment_extnd_df[row] = sec_struct_assignment[index_3]
                            strand_number_list_extnd_df[row] = strand_number_list[index_3]
                            sheet_number_list_extnd_df[row] = sheet_number_list[index_3]
                            orientation_list_extnd_df[row] = orientation_list[index_3]
                            bridge_pair_list_extnd_df[row] = bridge_pair_list[index_3]
                            phi_extnd_df[row] = phi[index_3]
                            psi_extnd_df[row] = psi[index_3]
                            break

            dssp_df = pd.DataFrame({'DSSP_NUM': dssp_num_extnd_df,
                                    'SHEET?': sec_struct_assignment_extnd_df,
                                    'STRAND_NUM': strand_number_list_extnd_df,
                                    'SHEET_NUM': sheet_number_list_extnd_df,
                                    'ORIENTATION': orientation_list_extnd_df,
                                    'H-BONDS': bridge_pair_list_extnd_df,
                                    'PHI': phi_extnd_df,
                                    'PSI': psi_extnd_df})
            cols = dssp_df.columns.tolist()
            cols = ([cols[0]] + [cols[5]] + [cols[7]] + [cols[6]] + [cols[2]]
                    + [cols[1]] + [cols[3]] + [cols[4]])
            dssp_df = dssp_df[cols]

            extnd_df = pd.concat([pdb_df, dssp_df], axis=1)
            extnd_df.to_pickle('Entire_domains/{}.pkl'.format(domain_id))

            retained_chains = extnd_df[extnd_df['SHEET?']=='E']['CHAIN'].tolist()
            retained_res_num = extnd_df[extnd_df['SHEET?']=='E']['RESNUM'].tolist()
            retained_inscode = extnd_df[extnd_df['SHEET?']=='E']['INSCODE'].tolist()
            chain_res_num = [retained_chains[index]+str(retained_res_num[index])+retained_inscode[index]
                             for index, value in enumerate(retained_chains)]
            for row in range(extnd_df.shape[0]):
                if ((extnd_df['CHAIN'][row]+str(extnd_df['RESNUM'][row])+extnd_df['INSCODE'][row])
                    not in chain_res_num
                    ):
                    extnd_df.loc[row, 'REC'] = None
            filtered_extnd_df = extnd_df[extnd_df['REC'].notnull()]
            filtered_extnd_df = filtered_extnd_df.reset_index(drop=True)
            filtered_extnd_df.to_pickle('Beta_strands/{}.pkl'.format(domain_id))
            dssp_dfs_dict[domain_id] = filtered_extnd_df

        return dssp_dfs_dict

    def write_dssp_sec_struct_pdb(self, dssp_dfs_dict):
        # Writes a PDB file of the residues that DSSP classifies as forming a
        # beta-strand (secondary structure code = 'E'). Individual beta-strands
        # are separated by 'TER' cards
        for domain_id, dssp_df in dssp_dfs_dict.items():
            print('Writing PDB file of beta-strands in {}'.format(domain_id))

            strand_number_set = [strand for strand in set(dssp_df['STRAND_NUM'].tolist())
                                 if strand != '']

            with open('Beta_strands/{}.pdb'.format(domain_id), 'w') as new_pdb_file:
                for strand in strand_number_set:
                    dssp_df_strand = dssp_df[dssp_df['STRAND_NUM']==strand]
                    chain = dssp_df_strand['CHAIN'].tolist()
                    res_num = dssp_df_strand['RESNUM'].tolist()
                    inscode = dssp_df_strand['INSCODE'].tolist()
                    chain_res_num = [chain[index]+str(res_num[index])+inscode[index]
                                     for index, value in enumerate(chain)]

                    for row in range(dssp_df.shape[0]):
                        if ((dssp_df['CHAIN'][row] + str(dssp_df['RESNUM'][row])
                            + dssp_df['INSCODE'][row]) in chain_res_num
                            ):
                            new_pdb_file.write('{}\n'.format(dssp_df['FILE_LINES'][row]))

                    new_pdb_file.write('TER'.ljust(80)+'\n')
