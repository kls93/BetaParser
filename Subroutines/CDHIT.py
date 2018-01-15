
import os
import copy
import pandas as pd

class filter_beta_structure():

    def __init__(self, run, resn, rfac, domain_df, pdb_database):
        self.run = run
        self.resn = resn
        self.rfac = rfac
        self.domain_df = domain_df
        self.pdb_database = pdb_database

    # Downloads a copy of each beta-barrel / beta-sandwich PDB structure from
    # the RCSB PDB website, and extracts its experimental method, resolution
    # and R_factor (working value) from the header information. The structures
    # are filtered to only retain those determined by X-ray diffraction to a
    # resolution of 1.6 Angstroms or higher with an R_factor (working value)
    # of 0.20 or lower.
    def resn_rfac_filter(self):
        unprocessed_list_1 = []
        unprocessed_list_2 = []
        processed_list = []
        resolution_list = []
        rfactor_list = []
        row_num = self.domain_df.shape[0]

        for row in range(row_num):
            print('Obtaining header information for {}'.format(self.domain_df['PDB_CODE'][row]))
            print('{}'.format((row/row_num)*100))

            middle_characters = self.domain_df['PDB_CODE'][row][1:3]
            cwd = os.getcwd()
            os.chdir('/{}/{}'.format(self.pdb_database, middle_characters))

            header_pdb_lines = []
            try:
                with open('{}.pdb'.format(self.domain_df['PDB_CODE'][row]), 'r') as pdb_file:
                    remark_end = False
                    for line in pdb_file:
                        if (line.replace(' ', ''))[0:7] == 'REMARK4':
                            remark_end = True

                        if remark_end is True:
                            break
                        else:
                            header_pdb_lines.append(line)
            except FileNotFoundError:
                unprocessed_list_1.append(self.domain_df['PDB_CODE'][row])

            os.chdir('{}'.format(cwd))

            resolution = 0
            rfactor = 0
            for line in header_pdb_lines:
                whitespace_remv_line = line.replace(' ', '')
                if whitespace_remv_line.startswith('EXPDTA'):
                    if not any(x in whitespace_remv_line for x in ['XRAYDIFFRACTION', 'X-RAYDIFFRACTION']):
                        unprocessed_list_2.append(self.domain_df['PDB_CODE'][row])
                        break
                elif (whitespace_remv_line.startswith('REMARK2')
                    and 'ANGSTROM' in whitespace_remv_line):
                        try:
                            resolution = float(line[23:30])
                        except:
                            ValueError
                            resolution = 0
                elif whitespace_remv_line.startswith('REMARK3RVALUE'):
                    if any(x in whitespace_remv_line for x in ['(WORKINGSET)', '(WORKINGSET,NOCUTOFF)']):
                        rfactor = whitespace_remv_line.split(':')
                        try:
                            rfactor = float(rfactor[len(rfactor)-1])
                            break
                        except:
                            ValueError
                            rfactor = 0

            if resolution == 0 or rfactor == 0:
                unprocessed_list_2.append(self.domain_df['PDB_CODE'][row])
            elif resolution <= self.resn and rfactor <= self.rfac:
                processed_list.append(self.domain_df['PDB_CODE'][row])
                resolution_list.append(resolution)
                rfactor_list.append(rfactor)

        filtered_domain_df_part_1 = self.domain_df[self.domain_df['PDB_CODE'].isin(processed_list)]
        filtered_domain_df_part_1 = filtered_domain_df_part_1.reset_index(drop=True)
        filtered_domain_df_part_2 = pd.DataFrame({'RESOLUTION': resolution_list,
                                                    'RFACTOR': rfactor_list})
        filtered_domain_df = pd.concat([filtered_domain_df_part_1, filtered_domain_df_part_2], axis=1)
        filtered_domain_df.to_pickle('CATH_{}_resn_{}_rfac_{}_pre_cd_hit.pkl'.format(self.run, self.resn, self.rfac))
        filtered_domain_df.to_csv('CATH_{}_resn_{}_rfac_{}_pre_cd_hit.csv'.format(self.run, self.resn, self.rfac))

        with open('Unprocessed_CATH_{}.txt'.format(self.run), 'w') as unprocessed_file:
            unprocessed_file.write('PDB accession code not in PDB database '
                                   '(downloaded 24/11/2017):\n')
            unprocessed_list_1 = set(unprocessed_list_1)
            for pdb in unprocessed_list_1:
                unprocessed_file.write('{}\n'.format(pdb))

            unprocessed_file.write('\n\nNot solved by X-ray diffraction OR '
                                   'failed to extract value for resolution OR '
                                   'failed to extract value for Rfactor '
                                   '(working value):\n')
            unprocessed_list_2 = set(unprocessed_list_2)
            for pdb in unprocessed_list_2:
                unprocessed_file.write('{}\n'.format(pdb))
        return filtered_domain_df

    # Generates list of beta-structure chain entries for sequence redundancy
    # filtering using the cd_hit web server
    def gen_cd_hit_list(self, filtered_domain_df):
        fasta = filtered_domain_df['DSEQS'].tolist()
        pdb_ids = filtered_domain_df['PDB_CODE'].tolist()
        pdb_chains = filtered_domain_df['CHAIN'].tolist()

        fasta_copy = copy.copy(fasta)
        fasta_repeat = []
        for index, seq in enumerate(fasta_copy):
            if seq in fasta_repeat:
                fasta[index] = None
                pdb_ids[index] = None
                pdb_chains[index] = None
            elif seq not in fasta_repeat:
                fasta_repeat.append(seq)

        fasta = [seq for seq in fasta if seq is not None]
        pdb_ids = [pdb_id for pdb_id in pdb_ids if pdb_id is not None]
        pdb_chains = [chain for chain in pdb_chains if chain is not None]

        with open(
            'CATH_{}_resn_{}_rfac_{}_domain_chain_entries_for_CD_HIT.txt'.format(
            self.run, self.resn, self.rfac), 'w'
            ) as chain_entries_file:
            for num in range(len(fasta)):
                chain_entries_file.write('>{}_{}\n'.format(pdb_ids[num], pdb_chains[num]))
                chain_entries_file.write('{}\n'.format(fasta[num]))
