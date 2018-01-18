
import os
import copy
import pandas as pd
if __name__ == 'subroutines.CDHIT':
    from subroutines.run_stages import run_stages
else:
    from datagen.subroutines.run_stages import run_stages


class filter_beta_structure(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    # Filters the all-beta structures extracted from the CATH / SCOPe database
    # to retain only those determined by X-ray diffraction of resolution and
    # Rfactor (working value) below the user-specified cutoff values.
    # (Recommended: resolution < 1.6 Angstroms, Rfactor (working value) < 0.20)
    def resn_rfac_filter(self, domain_df):
        unprocessed_list_1 = []
        unprocessed_list_2 = []
        processed_list = []
        resolution_list = []
        rfactor_list = []

        for row in range(domain_df.shape[0]):
            print('Obtaining header information for {}'.format(domain_df['PDB_CODE'][row]))
            print('{:0.2f}%'.format(((row+1)/domain_df.shape[0])*100))

            middle_characters = domain_df['PDB_CODE'][row][1:3]
            cwd = os.getcwd()
            os.chdir('{}{}'.format(self.pdb_database, middle_characters))

            header_pdb_lines = []
            try:
                with open('{}.pdb'.format(domain_df['PDB_CODE'][row]), 'r') as pdb_file:
                    remark_end = False
                    for line in pdb_file:
                        if (line.replace(' ', ''))[0:7] == 'REMARK4':
                            remark_end = True

                        if remark_end is True:
                            break
                        else:
                            header_pdb_lines.append(line)
            except FileNotFoundError:
                unprocessed_list_1.append(domain_df['PDB_CODE'][row])

            os.chdir('{}'.format(cwd))

            resolution = 0
            rfactor = 0
            for line in header_pdb_lines:
                whitespace_remv_line = line.replace(' ', '')
                if whitespace_remv_line.startswith('EXPDTA'):
                    if not any(x in whitespace_remv_line for x in ['XRAYDIFFRACTION', 'X-RAYDIFFRACTION']):
                        unprocessed_list_2.append(domain_df['PDB_CODE'][row])
                        break
                elif (whitespace_remv_line.startswith('REMARK2')
                    and 'ANGSTROM' in whitespace_remv_line):
                        try:
                            resolution = float(line[23:30])
                        except ValueError:
                            resolution = 0
                            break
                elif whitespace_remv_line.startswith('REMARK3RVALUE'):
                    if any(x in whitespace_remv_line for x in ['(WORKINGSET)', '(WORKINGSET,NOCUTOFF)']):
                        rfactor = whitespace_remv_line.split(':')
                        try:
                            rfactor = float(rfactor[len(rfactor)-1])
                            break
                        except ValueError:
                            rfactor = 0

            if resolution == 0 or rfactor == 0:
                unprocessed_list_2.append(domain_df['PDB_CODE'][row])
            elif resolution <= self.resn and rfactor <= self.rfac:
                processed_list.append(domain_df['PDB_CODE'][row])
                resolution_list.append(resolution)
                rfactor_list.append(rfactor)

        filtered_domain_df_part_1 = domain_df[domain_df['PDB_CODE'].isin(processed_list)]
        filtered_domain_df_part_1 = filtered_domain_df_part_1.reset_index(drop=True)
        filtered_domain_df_part_2 = pd.DataFrame({'RESOLUTION': resolution_list,
                                                  'RFACTOR': rfactor_list})
        filtered_domain_df = pd.concat([filtered_domain_df_part_1, filtered_domain_df_part_2], axis=1)
        filtered_domain_df.to_pickle('CDHIT_entries.pkl')
        filtered_domain_df.to_csv('CDHIT_entries.csv')

        with open('Unprocessed_domains.txt', 'w') as unprocessed_file:
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
        domain_ids = filtered_domain_df['DOMAIN_ID'].tolist()

        fasta_copy = copy.copy(fasta)
        fasta_repeat = []
        for index, seq in enumerate(fasta_copy):
            if seq in fasta_repeat:
                fasta[index] = None
                domain_ids[index] = None
            elif seq not in fasta_repeat:
                fasta_repeat.append(seq)

        fasta = [seq for seq in fasta if seq is not None]
        domain_ids = [domain_id for domain_id in domain_ids if domain_id is not None]

        with open('CDHIT_entries.txt', 'w') as chain_entries_file:
            for num in range(len(fasta)):
                chain_entries_file.write('>{}\n'.format(domain_ids[num]))
                chain_entries_file.write('{}\n'.format(fasta[num]))
