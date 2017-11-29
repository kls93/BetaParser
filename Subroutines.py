
import os
import shutil
import pandas as pd
import copy
import random
import string
import networkx as nx
import netgraph
import matplotlib.pyplot as plt
from difflib import SequenceMatcher
from collections import OrderedDict


# Defines dictionary of three- and one-letter codes of standard amino acids
amino_acids_dict = {'ALA': 'A',
                    'ARG': 'R',
                    'ASN': 'N',
                    'ASP': 'D',
                    'CYS': 'C',
                    'GLN': 'Q',
                    'GLU': 'E',
                    'GLY': 'G',
                    'HIS': 'H',
                    'ILE': 'I',
                    'LEU': 'L',
                    'LYS': 'K',
                    'MET': 'M',
                    'PHE': 'F',
                    'PRO': 'P',
                    'SER': 'S',
                    'THR': 'T',
                    'TRP': 'W',
                    'TYR': 'Y',
                    'VAL': 'V'}


class select_beta_structures():

    def __init__(self, run):
        self.run = run

    # Generates a list of the domain descriptions provided in
    # CATH_domain_description_v_4_2_0.txt
    def domain_desc_list(self):
        domains_description = []
        with open('CATH_domain_desc_v_4_2_0.txt', 'r') as domains_file:
            line_1 = False
            indv_domain_list = []
            for line in domains_file:
                if line.startswith('FORMAT'):
                    line_1 = True
                    start = True
                elif line.startswith('//'):
                    start = False

                if line_1 is True:
                    if start is True:
                        indv_domain_list.append(line)
                    elif start is False:
                        domains_description.append(''.join(indv_domain_list[0:len(indv_domain_list)]))
                        indv_domain_list = []

        domains_description = [line for line in domains_description if line != '']
        return domains_description

    # Filters the domain descriptions list for beta-structures (either
    # sandwiches or barrels depending upon the user's choice), picking out PDB
    # accession codes and sequences, whose values are stored in a dataframe.
    def domain_desc_filter(self, domains_description):
        domain_pdb_ids = []
        domain_chains = []
        domain_cathcodes = []
        domain_ids = []
        domain_dseqs = []
        domain_sseqs = []
        domain_sseqs_start_stop = []
        for domain in domains_description:
            if 'CATHCODE  {}'.format(self.run) in domain:
                dseqs_list = []
                sseqs_consec_list = []
                sseqs_list = []
                sseqs_start_stop_list = []

                domain_sublist = domain.split('\n')
                domain_sublist_copy = copy.copy(domain_sublist)

                for index, line in enumerate(domain_sublist):
                    if line.startswith('SSEQS') and domain_sublist[index+1].startswith('SSEQS'):
                        sseqs_consec_list.append(index)
                for index in sseqs_consec_list:
                    domain_sublist[index+1] = ''.join(domain_sublist_copy[index:index+2])
                    domain_sublist[index] = ''
                    domain_sublist_copy = copy.copy(domain_sublist)

                for line in domain_sublist:
                    if line.startswith('DOMAIN'):
                        domain_pdb_ids.append(line[10:14])
                        domain_ids.append(line[10:].strip())
                        chain = line[14:].strip()
                        chain = ''.join([char for char in chain if char.isalpha()])
                        domain_chains.append(chain)
                    elif line.startswith('CATHCODE'):
                        domain_cathcodes.append(line[10:])
                    elif line.startswith('DSEQS'):
                        line = line.replace('DSEQS', '')
                        line = line.replace(' ', '')
                        dseqs_list.append(line)
                    elif line.startswith('SSEQS'):
                        line = line.replace('SSEQS', '')
                        line = line.replace(' ', '')
                        sseqs_list.append(line)
                    elif line.startswith('SRANGE'):
                        line = line.replace('SRANGE', '')
                        start_stop = line.split()
                        start_stop = [item.strip() for item in start_stop if item.strip() != '']
                        sseqs_start_stop_list.append(start_stop)
                dseqs = ''.join(dseqs_list)
                domain_dseqs.append(dseqs)
                domain_sseqs.append(sseqs_list)
                domain_sseqs_start_stop.append(sseqs_start_stop_list)

        domain_df = pd.DataFrame({'PDB_CODE': domain_pdb_ids,
                                  'CHAIN': domain_chains,
                                  'CATHCODE': domain_cathcodes,
                                  'DOMAIN_ID': domain_ids,
                                  'DSEQS': domain_dseqs,
                                  'SSEQS': domain_sseqs,
                                  'SSEQS_START_STOP': domain_sseqs_start_stop})

        return domain_df


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

        with open('Unprocessed_CATH_{}_PDB_files.txt'.format(self.run), 'w') as unprocessed_file:
            unprocessed_file.write('PDB accession code not in PDB database '
                                   '(downloaded 24/11/2017):\n')
            unprocessed_list_1 = set(unprocessed_list_1)
            for pdb in unprocessed_list_1:
                unprocessed_file.write('{}\n'.format(pdb))

            unprocessed_file.write('\n\nNot solved by X-ray diffraction OR '
                                   'failed to extract value for resolution OR '
                                   'failed to extract value for Rfactor '
                                   '(working value)\n')
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


class extract_beta_structure_coords():

    def __init__(self, run, resn, rfac, pdb_database):
        self.run = run
        self.resn = resn
        self.rfac = rfac
        self.pdb_database = pdb_database

    def gen_cd_hit_dict(self, filtered_domain_df):
        # Loads the list of FASTA sequences generated by the CD-HIT web server
        fasta_list = []
        with open(
            'CATH_{}_resn_{}_rfac_{}_seqid_0.40_globalfilt_CD_HIT_output.txt'.format(
                self.run, self.resn, self.rfac
                ), 'r'
            ) as chain_entries_file:
            for seq in chain_entries_file:
                if not seq.startswith('>'):
                    fasta_list.append(seq.replace('\n', ''))

        # For each of the sequences returned by the CD-HIT web server, selects a
        # random entry from the filtered dataframe of CATH domains from which to
        # select coordinates
        df_index_list = []
        for seq in fasta_list:
            df_index_sub_list = []
            for row in range(filtered_domain_df.shape[0]):
                if seq == filtered_domain_df['DSEQS'][row]:
                    df_index_sub_list.append(row)

            rand_num = random.randint(0, len(df_index_sub_list)-1)
            index = df_index_sub_list[rand_num]
            df_index_list.append(index)

        # Filters dataframe further to retain only the domains selected in the previous
        # step
        cd_hit_domain_df = filtered_domain_df.iloc[df_index_list]
        cd_hit_domain_df = cd_hit_domain_df.reset_index(drop=True)
        return cd_hit_domain_df

    def get_xyz_coords(self, cd_hit_domain_df):
        # Extends the filtered (for resolution, R_factor (working value) and
        # sequence redundancy) dataframe to list the xyz coordinates of each
        # segment sequence (SSEQS)
        pdb_dfs_dict = OrderedDict()
        domain_residue_list = []
        unprocessed_list = []

        for row in range(cd_hit_domain_df.shape[0]):
            residue_list = []
            rec = []
            atmnum = []
            atmname = []
            conformer = []
            resname = []
            chain = []
            resnum = []
            insertion = []
            xpos = []
            ypos = []
            zpos = []
            occ = []
            bfac = []
            element = []
            charge = []

            # Downloads PDB file of the structure from the RCSB website
            print('{}'.format((row/cd_hit_domain_df.shape[0])*100))
            print('Obtaining ATOM / HETATM records for {}'.format(cd_hit_domain_df['PDB_CODE'][row]))

            middle_characters = cd_hit_domain_df['PDB_CODE'][row][1:3]
            cwd = os.getcwd()
            os.chdir('/{}/{}'.format(self.pdb_database, middle_characters))

            with open('{}.pdb'.format(cd_hit_domain_df['PDB_CODE'][row]), 'r') as pdb_file:
                pdb_file_lines = [line.strip('\n') for line in pdb_file if
                                  line[0:6].strip() in ['ATOM', 'HETATM', 'TER']]
            pdb_file_lines.append('TER'.ljust(80))

            os.chdir('{}'.format(cwd))

            # For each segment sequence in the domain, makes a list of all
            # sequences in the input PDB file that lie between the recorded
            # start and stop residue numbers and have the same chain id
            for index_1, segment in enumerate(cd_hit_domain_df['SSEQS'][row]):
                sequences = []
                indices = []

                start = cd_hit_domain_df['SSEQS_START_STOP'][row][index_1][0].replace('START=', '')
                stop = cd_hit_domain_df['SSEQS_START_STOP'][row][index_1][1].replace('STOP=', '')
                start_seq = False
                stop_seq = False
                sequence = ''
                index = []

                for index_2, line in enumerate(pdb_file_lines):
                    if index_2 != (len(pdb_file_lines)-1):
                        if line[22:27].strip() == start and line[21:22] == cd_hit_domain_df['CHAIN'][row]:
                            start_seq = True

                        if start_seq is True and stop_seq is False:
                            index.append(index_2)
                            if (line[22:27].strip() != pdb_file_lines[index_2+1][22:27].strip()
                                or pdb_file_lines[index_2+1][0:3] == 'TER'):
                                    if line[17:20].strip() in amino_acids_dict:
                                        sequence = sequence + amino_acids_dict[line[17:20].strip()]
                        elif stop_seq is True:
                            sequences.append(sequence)
                            indices.append(index)
                            sequence = ''
                            index = []
                            start_seq = False
                            stop_seq = False
                            continue

                        if (pdb_file_lines[index_2+1][0:3] == 'TER'
                            or (line[22:27].strip() == stop
                                and line[21:22] == cd_hit_domain_df['CHAIN'][row]
                                and pdb_file_lines[index_2+1][22:27].strip() != stop
                                )
                            ):
                                stop_seq = True

                sequence_identified = False
                sseqs_list = []
                for index_3, sequence in enumerate(sequences):
                    similarity = SequenceMatcher(a=segment, b=sequence).ratio()
                    if similarity > 0.95:
                        sequence_identified = True

                        pdb_file = open('CD_HIT_DSEQS/{}.pdb'.format(
                            cd_hit_domain_df['DOMAIN_ID'][row]), 'a')

                        for index_4 in indices[index_3]:
                            sseqs_list.append(pdb_file_lines[index_4][21:27].replace(' ', ''))

                            pdb_file.write('{}\n'.format(pdb_file_lines[index_4]))

                            rec.append(pdb_file_lines[index_4][0:6].strip())
                            atmnum.append(int(pdb_file_lines[index_4][6:11].strip()))
                            atmname.append(pdb_file_lines[index_4][12:16].strip())
                            conformer.append(pdb_file_lines[index_4][16:17].strip())
                            resname.append(pdb_file_lines[index_4][17:20].strip())
                            chain.append(pdb_file_lines[index_4][21:22].strip())
                            resnum.append(int(pdb_file_lines[index_4][22:26].strip()))
                            insertion.append(pdb_file_lines[index_4][26:27].strip())
                            xpos.append(float(pdb_file_lines[index_4][30:38].strip()))
                            ypos.append(float(pdb_file_lines[index_4][38:46].strip()))
                            zpos.append(float(pdb_file_lines[index_4][46:54].strip()))
                            occ.append(float(pdb_file_lines[index_4][54:60].strip()))
                            bfac.append(float(pdb_file_lines[index_4][60:66].strip()))
                            element.append(pdb_file_lines[index_4][76:78].strip())
                            charge.append(pdb_file_lines[index_4][78:80].strip())

                        pdb_file.write('TER'.ljust(80)+'\n')
                        pdb_file.close()

                        pdb_df = pd.DataFrame({'REC': rec,
                                               'ATMNUM': atmnum,
                                               'ATMNAME': atmname,
                                               'CONFORMER': conformer,
                                               'RESNAME': resname,
                                               'CHAIN': chain,
                                               'RESNUM': resnum,
                                               'INSCODE': insertion,
                                               'XPOS': xpos,
                                               'YPOS': ypos,
                                               'ZPOS': zpos,
                                               'OCC': occ,
                                               'BFAC': bfac,
                                               'ELEMENT': element,
                                               'CHARGE': charge})
                        cols = pdb_df.columns.tolist()
                        cols = ([cols[9]] + [cols[1]] + [cols[0]] + [cols[5]]
                                + [cols[10]] + [cols[3]] + [cols[11]] + [cols[7]]
                                + [cols[12]] + [cols[13]] + [cols[14]] + [cols[8]]
                                + [cols[2]] + [cols[6]] + [cols[4]])
                        pdb_df = pdb_df[cols]
                        pdb_dfs_dict[cd_hit_domain_df['DOMAIN_ID'][row]] = pdb_df

                        residue_list.extend(sseqs_list)

                        break

                if sequence_identified is False:
                    unprocessed_list.append('{}'.format(cd_hit_domain_df['DOMAIN_ID'][row]))

            chain_num_list = []
            for chain_num in residue_list:
                if chain_num not in chain_num_list:
                    chain_num_list.append(chain_num)
            domain_residue_list.append(chain_num_list)

        # Lists all segment sequences that could not be processed owing to
        # insufficient similarity between the FASTA sequence listed in the
        # CATH_domain_desc_v_4_2_0.txt and the FASTA sequence extracted from
        # the PDB file downloaded from the RCSB website in the unprocessed
        # structures file
        with open('Unprocessed_CATH_{}_PDB_files.txt'.format(self.run), 'a') as unprocessed_file:
            unprocessed_file.write('\n\nDissimilar coordinates:\n')
            for pdb_file in set(unprocessed_list):
                unprocessed_file.write('{}\n'.format(pdb_file))

        cd_hit_domain_df = cd_hit_domain_df[~cd_hit_domain_df['DOMAIN_ID'].isin(unprocessed_list)]
        cd_hit_domain_df = cd_hit_domain_df.reset_index(drop=True)

        # Appends xyz coordinates column to dataframe
        domain_residue_list = [residue_list for residue_list in
                               domain_residue_list if len(residue_list) > 0]
        domain_residue_df = pd.DataFrame({'CHAIN_NUM': domain_residue_list})
        cd_hit_domain_df = pd.concat([cd_hit_domain_df, domain_residue_df], axis=1)
        cd_hit_domain_df.to_csv('CATH_{}_resn_{}_rfac_{}_filtered.csv'.format(
            self.run, self.resn, self.rfac)
            )
        cd_hit_domain_df.to_pickle('CATH_{}_resn_{}_rfac_{}_filtered.pkl'.format(
            self.run, self.resn, self.rfac)
            )

        return cd_hit_domain_df, pdb_dfs_dict


class filter_dssp_database():

    def __init__(self, run, resn, rfac, dssp_database):
        self.run = run
        self.resn = resn
        self.rfac = rfac
        self.dssp_database = dssp_database


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
        with open('Unprocessed_CATH_{}_PDB_files.txt'.format(self.run), 'a') as unprocessed_file:
            unprocessed_file.write('\n\nNot in DSSP database:\n')
            for pdb_code in unprocessed_list:
                unprocessed_file.write('{}\n'.format(pdb_code))

        # Filters cd_hit_domain_df_xyz to remove entries not in the DSSP
        # database
        dssp_domain_df = cd_hit_domain_df[~cd_hit_domain_df['PDB_CODE'].isin(unprocessed_list)]
        dssp_domain_df = dssp_domain_df.reset_index(drop=True)

        return dssp_domain_df


class beta_structure_dssp_classification():

    def __init__(self, run, resn, rfac):
        self.run = run
        self.resn = resn
        self.rfac = rfac

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
                dssp_indv_file_lines.append('TER'.ljust(136))
                dssp_file_lines.append(dssp_indv_file_lines)
            elif len(chain_num_list) > 0:
                unprocessed_list.append(dssp_domain_df['DOMAIN_ID'][row])

        # Writes PDB accession codes that could not be processed to output file
        with open('Unprocessed_CATH_{}_PDB_files.txt'.format(self.run), 'a') as unprocessed_file:
            unprocessed_file.write('\n\n')
            unprocessed_file.write('Coordinates missing from DSSP file:\n')
            for domain_id in unprocessed_list:
                unprocessed_file.write('{}\n'.format(domain_id))

        # Filters dssp_domain_df to remove entries whose DSSP files have
        # coordinates missing
        dssp_domain_df = dssp_domain_df[~dssp_domain_df['DOMAIN_ID'].isin(unprocessed_list)]
        dssp_domain_df = dssp_domain_df.reset_index(drop=True)

        dssp_residues_dict = OrderedDict()
        for index, sub_list in enumerate(dssp_file_lines):
            dssp_residues_dict[dssp_domain_df['DOMAIN_ID'][index]] = sub_list

        return dssp_residues_dict

    def get_dssp_sec_struct_df(self, dssp_residues_dict, pdb_dfs_dict):
        # Generates dataframe of relevant information in DSSP file
        dssp_dfs_dict = OrderedDict()

        for domain_id in list(dssp_residues_dict.keys()):
            dssp_indv_file_lines = dssp_residues_dict[domain_id]
            pdb_df = pdb_dfs_dict[domain_id]

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

            for index_2, line in enumerate(dssp_indv_file_lines):
                if not line.startswith('TER'):
                    res_num.append(line[5:11].strip())
                    chain.append(line[11:12])
                    dssp_num.append(line[0:5].strip())
                    secondary_structure = line[16:17]
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

            for row in range(row_num):
                if pdb_df['ATMNAME'][row] == 'N':
                    for index, line in enumerate(res_num):
                        if ((str(pdb_df['RESNUM'][row])+pdb_df['INSCODE'][row]) == res_num[index]
                            and pdb_df['CHAIN'][row] == chain[index]
                            ):
                            dssp_num_extnd_df[row] = dssp_num[index]
                            sec_struct_assignment_extnd_df[row] = sec_struct_assignment[index]
                            strand_number_list_extnd_df[row] = strand_number_list[index]
                            sheet_number_list_extnd_df[row] = sheet_number_list[index]
                            orientation_list_extnd_df[row] = orientation_list[index]
                            bridge_pair_list_extnd_df[row] = bridge_pair_list[index]
                            break

            dssp_df = pd.DataFrame({'DSSP_NUM': dssp_num_extnd_df,
                                    'SHEET?': sec_struct_assignment_extnd_df,
                                    'STRAND_NUM': strand_number_list_extnd_df,
                                    'SHEET_NUM': sheet_number_list_extnd_df,
                                    'ORIENTATION': orientation_list_extnd_df,
                                    'H-BONDS': bridge_pair_list_extnd_df})
            cols = dssp_df.columns.tolist()
            cols = ([cols[0]] + [cols[3]] + [cols[5]] + [cols[4]] + [cols[2]]
                    + [cols[1]])
            dssp_df = dssp_df[cols]

            extnd_df = pd.concat([pdb_df, dssp_df], axis=1)
            extnd_df.to_pickle('CD_HIT_DSEQS/{}.pkl'.format(domain_id))

            retained_chains = extnd_df[extnd_df['SHEET?']=='E']['CHAIN'].tolist()
            retained_resnum = extnd_df[extnd_df['SHEET?']=='E']['RESNUM'].tolist()
            filtered_extnd_df = extnd_df[extnd_df['CHAIN'].isin(retained_chains)
                                         & extnd_df['RESNUM'].isin(retained_resnum)]
            filtered_extnd_df = filtered_extnd_df.reset_index(drop=True)
            dssp_dfs_dict[domain_id] = filtered_extnd_df

        return dssp_dfs_dict

    def write_dssp_sec_struct_pdb(self, dssp_residues_dict, dssp_dfs_dict):
        # Writes a PDB file of the residues that DSSP classifies as forming a
        # beta-strand (secondary structure code = 'E'). Individual beta-strands
        # are separated by 'TER' cards
        for domain_id in list(dssp_residues_dict.keys()):
            print('Writing PDB file of beta-strands in {}'.format(domain_id))

            dssp_df = dssp_dfs_dict[domain_id]
            strand_number_set = [strand for strand in set(dssp_df['STRAND_NUM'].tolist())
                                 if strand != '']

            with open('CD_HIT_DSEQS/{}.pdb'.format(domain_id), 'r') as pdb_file:
                pdb_file_lines = [line.strip('\n') for line in pdb_file
                                  if line[0:6].strip() in ['ATOM', 'HETATM']]

            with open('DSSP_filtered_DSEQS/{}.pdb'.format(domain_id), 'w') as new_pdb_file:
                for strand in strand_number_set:
                    dssp_df_strand = dssp_df[dssp_df['STRAND_NUM']==strand]
                    dssp_df_strand = dssp_df_strand.reset_index(drop=True)
                    for row in range(dssp_df_strand.shape[0]):
                        chain = dssp_df_strand['CHAIN'][row]
                        res_num = (str(dssp_df_strand['RESNUM'][row])+dssp_df_strand['INSCODE'][row])
                        for line in pdb_file_lines:
                            if line[21:22] == chain and line[22:27].strip() == res_num:
                                new_pdb_file.write('{}\n'.format(line))
                    new_pdb_file.write('TER'.ljust(80)+'\n')


class manipulate_beta_structure():

    def merge_sheets(self, dssp_residues_dict, dssp_dfs_dict):
        # Merges sheet names of sheets that share strands
        dssp_dfs_merged_sheets_dict = OrderedDict()

        for domain_id in list(dssp_residues_dict.keys()):
            dssp_df = dssp_dfs_dict[domain_id]
            print('Merging sheets that share beta-strands in {}'.format(domain_id))

            sheets = [sheet for sheet in set(dssp_df['SHEET_NUM'].tolist()) if sheet != '']
            sheet_strands = []
            for sheet in sheets:
                locals()['sheet_{}_strands'.format(sheet)] = []
                for strand in dssp_df[dssp_df['SHEET_NUM']==sheet]['STRAND_NUM'].tolist():
                    if strand not in locals()['sheet_{}_strands'.format(sheet)]:
                        locals()['sheet_{}_strands'.format(sheet)].append(strand)
                sheet_strands.append(locals()['sheet_{}_strands'.format(sheet)])
                del locals()['sheet_{}_strands'.format(sheet)]

            sheet_strand_numbers = [strand for strands in sheet_strands for strand in strands]
            count_1 = 0
            while len(sheet_strand_numbers) != len(set(sheet_strand_numbers)):
                strands = sheet_strands[count_1]
                sheet_strands_copy = copy.copy(sheet_strands)

                count_2 = 1
                restart = True
                while restart is True:
                    strands = sheet_strands[count_1]
                    for index in range(0, len(sheet_strands)):
                        if index == len(sheet_strands) - 1:
                            restart = False
                        if count_1 != index and strands != '':
                            intersect = set(strands).intersection(set(sheet_strands[index]))
                            if len(intersect) > 0:
                                sheet_strands_copy[count_1] = list(set(sheet_strands[count_1])|set(sheet_strands[index]))
                                sheet_strands_copy[index] = ''
                                sheet_strands = sheet_strands_copy
                                restart = True

                sheet_strand_numbers = [strand for strands in sheet_strands for strand in strands]

            sheet_strands = [strands for strands in sheet_strands if len(strands) > 2]

            for row in range(dssp_df.shape[0]):
                for strands in sheet_strands:
                    if str(dssp_df['STRAND_NUM'][row]) in strands:
                        dssp_df.loc[row, 'SHEET_NUM'] = (sheet_strands.index(strands)+1)
                        break

            dssp_dfs_merged_sheets_dict[domain_id] = dssp_df
            dssp_df.to_pickle('DSSP_filtered_DSEQS/{}.pkl'.format(domain_id))

        return dssp_dfs_merged_sheets_dict

    def identify_strand_interactions(self, dssp_residues_dict,
                                     dssp_dfs_dict):
        # Separates the beta-strands identified by DSSP into sheets, and works
        # out the strand interactions and thus loop connections both within and
        # between sheets
        domain_networks_dict = OrderedDict()

        for domain_id in list(dssp_residues_dict.keys()):
            dssp_df = dssp_dfs_dict[domain_id]

            sheets = [sheet for sheet in set(dssp_df['SHEET_NUM'].tolist()) if sheet != '']
            for sheet in sheets:
                sheet_df = dssp_df[dssp_df['SHEET_NUM']==sheet]

                strands_dict = {}
                strands = [int(strand) for strand in set(sheet_df['STRAND_NUM'].tolist())]
                for strand in strands:
                    strand_df = sheet_df[sheet_df['STRAND_NUM']==strand]

                    strand_res_num = [int(res) for res in strand_df['DSSP_NUM'].tolist()]
                    strand_min = min(strand_res_num)
                    strand_max = max(strand_res_num)

                    strands_dict[strand] = [num for num in range(strand_min, strand_max+1)]

                strand_pairs = {}
                for index, pair in enumerate(sheet_df['H-BONDS'].tolist()):
                    res_num_1 = int(sheet_df['DSSP_NUM'].tolist()[index])
                    res_num_2 = int(pair[0])
                    orientation_2 = sheet_df['ORIENTATION'].tolist()[index][0]
                    res_num_3 = int(pair[1])
                    orientation_3 = sheet_df['ORIENTATION'].tolist()[index][1]
                    for key in strands_dict:
                        if res_num_1 in strands_dict[key]:
                            strand_1 = key

                        if not res_num_2 == 0:
                            if res_num_2 in strands_dict[key]:
                                strand_2 = key
                        else:
                            strand_2 = None

                        if not res_num_3 == 0:
                            if res_num_3 in strands_dict[key]:
                                strand_3 = key
                        else:
                            strand_3 = None

                    if strand_2 is not None:
                        strand_pair = [strand_1, strand_2]
                        strand_pairs[(min(strand_pair), max(strand_pair))] = orientation_2
                    if strand_3 is not None:
                        strand_pair = [strand_1, strand_3]
                        strand_pairs[(min(strand_pair), max(strand_pair))] = orientation_3

                interacting_strands = {}
                for key in strand_pairs:
                    if key not in interacting_strands:
                        interacting_strands[key] = strand_pairs[key]

                # Creates network of the interacting strands in the sheet
                G = nx.Graph()
                for pair in interacting_strands:
                    G.add_edge(pair[0], pair[1], attr=interacting_strands[pair])
                pos = nx.circular_layout(G)
                domain_networks_dict[domain_id] = G

                # Draws plot of the network of interacting strands
                plt.clf()
                nx.draw_networkx(G, pos=pos, with_labels=True)
                nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=interacting_strands)
                plt.savefig(
                    'DSSP_filtered_DSEQS/{}_sheet_{}_network.png'.format(
                        domain_id, sheet
                        )
                    )

        return domain_networks_dict


class manipulate_beta_network():
    pass
