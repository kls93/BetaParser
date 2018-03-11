
import os
import shutil
import copy
import pandas as pd
from collections import OrderedDict
if __name__ == 'subroutines.DSSP':
    from subroutines.run_stages import run_stages
else:
    from datagen.subroutines.run_stages import run_stages


class filter_dssp_database(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    def copy_files_from_dssp_database(self, cdhit_domain_df):
        # Generates list of all PDB accession codes to be extracted from the
        # DSSP database
        pdb_list = set(cdhit_domain_df['PDB_CODE'].tolist())

        # Copies DSSP files of listed PDB accession codes to current working
        # directory
        unprocessed_list = []
        for pdb_code in pdb_list:
            middle_characters = pdb_code[1:len(pdb_code)-1]
            cwd = os.getcwd()

            try:
                shutil.copy2(
                    '{}{}.dssp'.format(self.dssp_database, pdb_code),
                    '{}/DSSP_files/{}.dssp'.format(cwd, pdb_code))

                print('Copied {}.dssp'.format(pdb_code))
            except FileNotFoundError:
                unprocessed_list.append(pdb_code)

        # Writes PDB accession codes that could not be processed to output file
        with open('Unprocessed_domains.txt', 'a') as unprocessed_file:
            unprocessed_file.write('\n\nNot in DSSP database:\n')
            for pdb_code in set(unprocessed_list):
                unprocessed_file.write('{}\n'.format(pdb_code))

        # Filters cdhit_domain_df_xyz to remove entries not in the DSSP
        # database
        dssp_domain_df = cdhit_domain_df[~cdhit_domain_df['PDB_CODE'].isin(unprocessed_list)]
        dssp_domain_df = dssp_domain_df.reset_index(drop=True)

        return dssp_domain_df


class beta_structure_dssp_classification(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    def extract_dssp_file_lines(self, dssp_domain_df):
        unprocessed_list = []

        # Creates dictionary of DSSP file lines for all input PDB accession
        # codes
        dssp_residues_dict = OrderedDict()

        for row in range(dssp_domain_df.shape[0]):
            dssp_indv_file_lines = []
            chain_num_list = copy.copy(dssp_domain_df['CHAIN_NUM'][row])

            with open('DSSP_files/{}.dssp'.format(dssp_domain_df['PDB_CODE'][row]), 'r') as dssp_file:
                print('Processing {}.dssp'.format(dssp_domain_df['PDB_CODE'][row]))
                for line in dssp_file:
                    if (line[11:12].strip()+line[5:11].strip()) in chain_num_list:
                        dssp_indv_file_lines.append(line.strip('\n'))
                        index = chain_num_list.index(line[11:12].strip()+line[5:11].strip())
                        chain_num_list[index] = ''

            chain_num_list = [chain_num for chain_num in chain_num_list if chain_num != '']
            if len(chain_num_list) == 0:
                dssp_indv_file_lines.append('TER'.ljust(136)+'\n')
                dssp_residues_dict[dssp_domain_df['DOMAIN_ID'][row]] = dssp_indv_file_lines
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

        return dssp_residues_dict

    def get_dssp_sec_struct_df(self, dssp_residues_dict, all_atoms_dfs_dict):
        # Generates dataframe of relevant information in DSSP file for each
        # domain
        sec_struct_dfs_dict = OrderedDict()

        for index_1, domain_id in enumerate(list(dssp_residues_dict.keys())):
            dssp_indv_file_lines = dssp_residues_dict[domain_id]
            pdb_df = all_atoms_dfs_dict[domain_id]

            print('Generating dataframe summarising DSSP info for {}'.format(domain_id))
            print('{:0.2f}%'.format(((index_1+1)/len(list(dssp_residues_dict.keys())))*100))

            chain_res_num = []
            dssp_num = []
            sec_struct_assignment = []
            strand_number = 1
            strand_number_list = []
            sheet_number_list = []
            orientation_list = []
            bridge_pair_list = []
            lines = []

            # Extracts secondary structure information from the DSSP file lines
            for index_2, line in enumerate(dssp_indv_file_lines):
                if not line.startswith('TER'):
                    chain_res_num.append(line[11:12].strip() + line[5:11].replace(' ', ''))
                    dssp_num.append(line[0:5].strip())
                    secondary_structure = line[16:17]
                    lines.append(line.strip('\n'))
                    if secondary_structure == 'E':
                        sec_struct_assignment.append(secondary_structure)
                        strand_number_list.append(strand_number)
                        if dssp_indv_file_lines[index_2+1][16:17] != 'E':
                            strand_number = strand_number + 1
                        sheet_number_list.append(line[33:34])
                        ladder_1 = line[23:24]
                        if ladder_1 == ' ':
                            ladder_1 = ''
                        elif ladder_1.isupper():  # Uppercase => antiparallel strands
                            ladder_1 = 'A'
                        elif ladder_1.islower():  # Lowercase => parallel strands
                            ladder_1 = 'P'
                        ladder_2 = line[24:25]
                        if ladder_2 == ' ':
                            ladder_2 = ''
                        elif ladder_2.isupper():  # Uppercase => antiparallel strands
                            ladder_2 = 'A'
                        elif ladder_2.islower():  # Lowercase => parallel strands
                            ladder_2 = 'P'
                        orientation_list.append([ladder_1, ladder_2])
                        bridge_pair_list.append([line[25:29].strip(),
                                                 line[29:33].strip()])
                    elif secondary_structure != 'E':
                        sec_struct_assignment.append('')
                        strand_number_list.append('')
                        sheet_number_list.append('')
                        orientation_list.append(['', ''])
                        bridge_pair_list.append(['', ''])

            # Initialises lists of additional properties extracted from the
            # DSSP file lines to be appended to the dataframe of PDB information
            # (created by extract_beta_structure_coords.get_xyz_coords function)
            row_num = pdb_df.shape[0]
            dssp_num_extnd_df = ['']*row_num
            sec_struct_assignment_extnd_df = ['']*row_num
            strand_number_list_extnd_df = ['']*row_num
            sheet_number_list_extnd_df = ['']*row_num
            orientation_list_extnd_df = ['']*row_num
            bridge_pair_list_extnd_df = ['']*row_num
            lines_extnd_df = ['']*row_num

            # Fills lists created in the previous step
            for row in range(row_num):
                if pdb_df['ATMNAME'][row] == 'CA':  # All DSSP info
                    # (= per-residue) appended to the CA rows of the PDB info
                    # (= per-atom)
                    for index_3, value in enumerate(chain_res_num):
                        if pdb_df['RES_ID'][row] == chain_res_num[index_3]:
                            dssp_num_extnd_df[row] = dssp_num[index_3]
                            sec_struct_assignment_extnd_df[row] = sec_struct_assignment[index_3]
                            strand_number_list_extnd_df[row] = strand_number_list[index_3]
                            sheet_number_list_extnd_df[row] = sheet_number_list[index_3]
                            orientation_list_extnd_df[row] = orientation_list[index_3]
                            bridge_pair_list_extnd_df[row] = bridge_pair_list[index_3]
                            lines_extnd_df[row] = lines[index_3]
                            break

            # Appends DSSP info to dataframe of PDB info
            dssp_df = pd.DataFrame({'DSSP_FILE_LINES': lines_extnd_df,
                                    'DSSP_NUM': dssp_num_extnd_df,
                                    'SHEET?': sec_struct_assignment_extnd_df,
                                    'STRAND_NUM': strand_number_list_extnd_df,
                                    'SHEET_NUM': sheet_number_list_extnd_df,
                                    'ORIENTATION': orientation_list_extnd_df,
                                    'H-BONDS': bridge_pair_list_extnd_df})
            cols = dssp_df.columns.tolist()
            cols = ([cols[0]] + [cols[1]] + [cols[4]] + [cols[6]] + [cols[5]]
                    + [cols[3]] + [cols[2]])
            dssp_df = dssp_df[cols]

            extnd_df = pd.concat([pdb_df, dssp_df], axis=1)
            extnd_df.to_pickle('Entire_domains/{}.pkl'.format(domain_id))
            all_atoms_dfs_dict[domain_id] = extnd_df

            retained_res = extnd_df[extnd_df['SHEET?'] == 'E']['RES_ID'].tolist()
            for row in range(extnd_df.shape[0]):
                if extnd_df['RES_ID'][row] not in retained_res:
                    extnd_df.loc[row, 'REC'] = None
            filtered_extnd_df = extnd_df[extnd_df['REC'].notnull()]
            filtered_extnd_df = filtered_extnd_df.reset_index(drop=True)
            filtered_extnd_df.to_pickle('Beta_strands/{}.pkl'.format(domain_id))
            sec_struct_dfs_dict[domain_id] = filtered_extnd_df

        return all_atoms_dfs_dict, sec_struct_dfs_dict
