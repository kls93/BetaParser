
import os
import shutil
import copy
import requests
import random
import pandas as pd
from difflib import SequenceMatcher
import networkx as nx
import matplotlib.pyplot as plt

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

class gen_beta_structure_df():
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
        domain_dseqs = []
        domain_segments = []
        domain_sseqs = []
        domain_sseqs_start_stop = []
        for domain in domains_description:
            if 'CATHCODE  {}'.format(self.run) in domain:
                dseqs_list = []
                segments_list = []
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
                        chain = line[14:].strip()
                        chain = ''.join([char for char in chain if char.isalpha()])
                        domain_chains.append(chain)
                    elif line.startswith('CATHCODE'):
                        domain_cathcodes.append(line[10:])
                    elif line.startswith('DSEQS'):
                        line = line.replace('DSEQS', '')
                        line = line.replace(' ', '')
                        dseqs_list.append(line)
                    elif line.startswith('SEGMENT'):
                        segments_list.append(line[10:])
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
                domain_segments.append(segments_list)
                domain_sseqs.append(sseqs_list)
                domain_sseqs_start_stop.append(sseqs_start_stop_list)

        domain_dict = pd.DataFrame({'PDB_CODE': domain_pdb_ids,
                                    'CHAIN': domain_chains,
                                    'CATHCODE': domain_cathcodes,
                                    'DSEQS': domain_dseqs,
                                    'SEGMENT_ID': domain_segments,
                                    'SSEQS': domain_sseqs,
                                    'SSEQS_START_STOP': domain_sseqs_start_stop})

        return domain_dict


class beta_structure_df():
    def __init__(self, run, resn, rfac, domain_dict):
        self.run = run
        self.resn = resn
        self.rfac = rfac
        self.domain_dict = domain_dict

    # Downloads a copy of each beta-barrel / beta-sandwich PDB structure from
    # the RCSB PDB website, and extracts its experimental method, resolution
    # and rfactor from the header information. The structures are filtered to
    # only retain those determined by X-ray diffraction to a resolution of 1.6
    # Angstrom or higher with an Rfactor (working value) of 0.25 or lower.
    def resn_rfac_filter(self):
        unprocessed_list = []
        processed_list = []
        resolution_list = []
        rfactor_list = []
        row_num = self.domain_dict.shape[0]

        for row in range(row_num):
            print('{}'.format((row/row_num)*100))
            print('Downloading {} from the RCSB PDB website'.format(self.domain_dict['PDB_CODE'][row]))
            url = 'http://www.rcsb.org/pdb/files/{}.pdb'.format(self.domain_dict['PDB_CODE'][row].upper())
            pdb_file_lines = requests.get(url).text
            pdb_file_lines = pdb_file_lines.split('\n')

            filtered_pdb_lines = []
            remark_end = False
            for line in pdb_file_lines:
                if (line.replace(' ', ''))[0:7] == 'REMARK4':
                    remark_end = True

                if remark_end is True:
                    break
                else:
                    filtered_pdb_lines.append(line)

            resolution = 0
            rfactor = 0
            for line in filtered_pdb_lines:
                whitespace_remv_line = line.replace(' ', '')
                if whitespace_remv_line.startswith('EXPDTA'):
                    if not any(x in whitespace_remv_line for x in ['XRAY', 'X-RAY']):
                        unprocessed_list.append(self.domain_dict['PDB_CODE'][row])
                        break
                elif (whitespace_remv_line.startswith('REMARK2')
                    and 'ANGSTROMS' in whitespace_remv_line):
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
                unprocessed_list.append(self.domain_dict['PDB_CODE'][row])
            elif resolution <= self.resn and rfactor <= self.rfac:
                processed_list.append(self.domain_dict['PDB_CODE'][row])
                resolution_list.append(resolution)
                rfactor_list.append(rfactor)

        filtered_domain_dict_part_1 = self.domain_dict.loc[self.domain_dict['PDB_CODE'].isin(processed_list)]
        filtered_domain_dict_part_1 = filtered_domain_dict_part_1.reset_index(drop=True)
        filtered_domain_dict_part_2 = pd.DataFrame({'RESOLUTION': resolution_list,
                                                    'RFACTOR': rfactor_list})
        filtered_domain_dict = pd.concat([filtered_domain_dict_part_1, filtered_domain_dict_part_2], axis=1)
        filtered_domain_dict.to_pickle('CATH_{}_resn_{}_rfac_{}_pre_cd_hit.pkl'.format(self.run, self.resn, self.rfac))
        filtered_domain_dict.to_csv('CATH_{}_resn_{}_rfac_{}_pre_cd_hit.csv'.format(self.run, self.resn, self.rfac))

        with open('Unprocessed_CATH_{}_PDB_files.txt'.format(self.run), 'w') as unprocessed_file:
            unprocessed_list = set(unprocessed_list)
            for pdb in unprocessed_list:
                unprocessed_file.write('{}\n'.format(pdb))
        return filtered_domain_dict

    # Generates list of beta-structure chain entries for sequence redundancy
    # filtering using the cd_hit web server
    def gen_cd_hit_list(self, filtered_domain_dict):
        fasta = filtered_domain_dict.DSEQS.tolist()
        pdb_ids = filtered_domain_dict.PDB_CODE.tolist()
        pdb_chains = filtered_domain_dict.CHAIN.tolist()

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

        count = len(fasta)
        with open('CATH_{}_{}_{}_domain_chain_entries.txt'.format(self.run, self.resn, self.rfac), 'w') as chain_entries_file:
            for num in range(count):
                chain_entries_file.write('>{}{}\n'.format(pdb_ids[num], pdb_chains[num]))
                chain_entries_file.write('{}\n'.format(fasta[num]))


class beta_structure_coords():
    def __init__(self, run, resn, rfac):
        self.run = run
        self.resn = resn
        self.rfac = rfac

    def gen_cd_hit_dict(self, filtered_domain_dict):
        # Loads the list of FASTA sequences generated by the CD-HIT web server
        fasta_list = []
        with open('CATH_{}_resn_{}_rfac_{}_seqid_0.40_globalfilt.txt'.format(self.run, self.resn, self.rfac), 'r') as chain_entries_file:
            for seq in chain_entries_file:
                if not seq.startswith('>'):
                    fasta_list.append(seq.replace('\n', ''))

        # For each of the sequences returned by the CD-HIT web server, selects a
        # random entry from the filtered dataframe of CATH domains from which to
        # select coordinates
        df_index_list = []
        for seq in fasta_list:
            df_index_sub_list = []
            for row in range(filtered_domain_dict.shape[0]):
                if seq == filtered_domain_dict['DSEQS'][row]:
                    df_index_sub_list.append(row)

            rand_num = random.randint(0, len(df_index_sub_list)-1)
            index = df_index_sub_list[rand_num]
            df_index_list.append(index)

        # Filters dataframe further to retain only the domains selected in the previous
        # step
        cd_hit_domain_dict = filtered_domain_dict.iloc[df_index_list]
        cd_hit_domain_dict = cd_hit_domain_dict.reset_index(drop=True)
        return cd_hit_domain_dict

    def get_xyz_coords(self, cd_hit_domain_dict):
        # Extends the filtered (for resolution, Rfactor and sequence redundancy)
        # dataframe to list the xyz coordinates of each segment sequence (SSEQS)
        domain_xyz = []
        unprocessed_list = []
        count = 0

        for row in range(cd_hit_domain_dict.shape[0]):
            count = count + 1

            print('Downloading {} from the RCSB PDB website'.format(cd_hit_domain_dict['PDB_CODE'][row]))
            url = 'http://www.rcsb.org/pdb/files/{}.pdb'.format(cd_hit_domain_dict['PDB_CODE'][row].upper())
            pdb_file_lines = requests.get(url).text
            pdb_file_lines = pdb_file_lines.split('\n')
            pdb_file_lines = [line for line in pdb_file_lines
                              if line[0:6].strip() in ['ATOM', 'HETATM', 'TER']]
            pdb_file_lines.append('TER'.ljust(80))

            xyz_segments = []
            for index_1, segment in enumerate(cd_hit_domain_dict['SSEQS'][row]):
                sequences = []
                indices = []

                start = cd_hit_domain_dict['SSEQS_START_STOP'][row][index_1][0].replace('START=', '')
                stop = cd_hit_domain_dict['SSEQS_START_STOP'][row][index_1][1].replace('STOP=', '')
                start_seq = False
                stop_seq = False
                sequence = ''
                index = []

                for index_2, line in enumerate(pdb_file_lines):
                    if index_2 != (len(pdb_file_lines)-1):
                        if line[22:27].strip() == start and line[21:22] == cd_hit_domain_dict['CHAIN'][row]:
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
                                and line[21:22] == cd_hit_domain_dict['CHAIN'][row]
                                and pdb_file_lines[index_2+1][22:27].strip() != stop
                                )
                            ):
                                stop_seq = True

                xyz = []
                sequence_identified = False
                for index_3, sequence in enumerate(sequences):
                    similarity = SequenceMatcher(a=segment, b=sequence).ratio()
                    if similarity > 0.95:
                        sequence_identified = True
                        pdb_file = open('PDB_files/{}_{}.pdb'.format(cd_hit_domain_dict['PDB_CODE'][row],
                                        count), 'a')
                        for index_4 in indices[index_3]:
                            pdb_file.write('{}\n'.format(pdb_file_lines[index_4]))
                            x = float(pdb_file_lines[index_4][30:38].strip())
                            y = float(pdb_file_lines[index_4][38:46].strip())
                            z = float(pdb_file_lines[index_4][46:54].strip())
                            xyz.append([x, y, z])
                        pdb_file.write('TER'.ljust(80)+'\n')
                        pdb_file.close()
                        xyz_segments.append(xyz)
                        break

                if sequence_identified is False:
                    unprocessed_list.append('{}{}\n'.format(cd_hit_domain_dict['PDB_CODE'][row],
                                                            cd_hit_domain_dict['SEGMENT_ID'][row]))
                    xyz_segments.append(xyz)

            domain_xyz.append(xyz_segments)

        with open('Unprocessed_CATH_{}_PDB_files.txt'.format(self.run), 'a') as unprocessed_file:
            unprocessed_file.write('\n\n')
            unprocessed_file.write('Dissimilar coordinates:\n')
            for pdb_file in set(unprocessed_list):
                unprocessed_file.write('{}\n'.format(pdb_file))
            unprocessed_file.write('\n\n')

        domain_xyz = pd.DataFrame({'XYZ': domain_xyz})
        cd_hit_domain_dict_xyz = pd.concat([cd_hit_domain_dict, domain_xyz], axis=1)
        cd_hit_domain_dict_xyz.to_csv('CATH_{}_resn_{}_rfac_{}_filtered_xyz.csv'.format(self.run, self.resn, self.rfac))
        cd_hit_domain_dict_xyz.to_pickle('CATH_{}_resn_{}_rfac_{}_filtered_xyz.pkl'.format(self.run, self.resn, self.rfac))


class beta_structure_dssp_classification():
    def __init__(self, run, resn, rfac, pdb_name, pdb_code):
        self.run = run
        self.resn = resn
        self.rfac = rfac
        self.pdb_name = pdb_name
        self.pdb = pdb_code

    def extract_dssp_secondary_structure(pdb_file_lines):
        unprocessed_files = []

        pdb_res_num = []
        pdb_chains = []
        for index, line in enumerate(pdb_file_lines):
            if index != (len(pdb_file_lines)-1):
                if (line[0:3] != 'TER'
                    and line[21:27].replace(' ', '') != pdb_file_lines[index+1][21:27].replace(' ', '')
                    ):
                    pdb_res_num.append(line[22:27].replace(' ', ''))
                    pdb_chains.append(line[21:22])

        middle_characters = self.pdb[1:len(self.pdb)-1]

        # Changes directory to open DSSP file
        cwd = os.getcwd()
        os.chdir('./../../shared/structural_bioinformatics/data/{}/{}/dssp/'.format(
            middle_characters, self.pdb))

        # Generates list from DSSP file on Tombstone
        start = False
        dssp_file_lines = []
        with open('{}_1.mmol.dssp'.format(self.pdb), 'r') as dssp_file:
            for line in dssp_file:
                if start is True:
                    if (line[5:11].replace(' ', '') in pdb_res_num
                        and line[11:12].strip() in pdb_chains
                        ):
                        dssp_file_lines.append(line.strip('\n'))
                if line.strip().startswith('#'):
                    start = True
        dssp_file_lines.append('TER'.ljust(136))

        # Changes back to current working directory after DSSP file has been read
        os.chdir(cwd)

        # Checks that all required asymmetric unit chains are in the DSSP file
        # opened on Tombstone (which are of the biological assemblies rather
        # than the asymmetric unit)
        dssp_chains = [line[11:12] for line in dssp_file_lines]
        for chain in set(pdb_chains):
            if chain not in dssp_chains:
                unprocessed_files.append(self.pdb_name)


        # Generates dataframe of relevant information in DSSP file
        dssp_num = []
        secondary_structure_assignment = []
        strand_number = 1
        strand_number_list = []
        sheet_number_list = []
        orientation_list = []
        bridge_pair_list = []
        for index_1, res_num in enumerate(pdb_res_num):
            for index_2, line in enumerate(dssp_file_lines):
                if (index_2 != (len(dssp_file_lines)-1)
                    and res_num == line[5:11].replace(' ', '')
                    and pdb_chains[index_1] == line[11:12]
                    ):
                    dssp_num.append(line[0:5].strip())
                    secondary_structure = line[16:17]
                    if secondary_structure == 'E':
                        secondary_structure_assignment.append(secondary_structure)
                        strand_number_list.append(strand_number)
                        if dssp_file_lines[index_2+1][16:17] != 'E':
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
                        secondary_structure_assignment.append('')
                        strand_number_list.append('')
                        sheet_number_list.append('')
                        orientation_list.append(['', ''])
                        bridge_pair_list.append(['', ''])
                    break

        dssp_df = pd.DataFrame({'RESNUM': pdb_res_num,
                                'CHAIN': pdb_chains,
                                'DSSP_NUM': dssp_num,
                                'SHEET?': secondary_structure_assignment,
                                'STRAND_NUM': strand_number_list,
                                'SHEET_NUM': sheet_number_list,
                                'ORIENTATION': orientation_list,
                                'H-BONDS': bridge_pair_list})

        # Writes a PDB file of the residues that DSSP classifies as forming a
        # beta-strand (secondary structure code = 'E')
        with open('DSSP_PDB_files_{}/{}'.format(self.run, self.pdb_name), 'w') as pdb_file:
            dssp_df_beta_structure = dssp_df[dssp_df['SHEET?']=='E']
            for index, secondary_structure in enumerate(secondary_structure_assignment):
                if secondary_structure == 'E':
                    chain = pdb_chains[index]
                    res_num = pdb_res_num[index]
                    for line in pdb_file_lines:
                        if line[21:22] == chain and line[22:27].replace(' ', '') == res_num:
                            pdb_file.write('{}\n'.format(line))
                        elif line[0:3] == 'TER':
                            pdb_file.write('{}\n'.format(line))
                elif

        # Merges sheet names of sheets that share strands
        sheets = [sheet for sheet in set(dssp_df['SHEET_NUM'].tolist()) if sheet != '']
        sheet_strands = []
        for sheet in sheets:
            locals()['sheet_{}_strands'.format(sheet)] = [
                str(strand) for strand in set(dssp_df[dssp_df['SHEET_NUM']==sheet]['STRAND_NUM'].tolist())
                ]
            sheet_strands.append(locals()['sheet_{}_strands'.format(sheet)])
            locals()['sheet_{}_strands'.format(sheet)] = []

        sheet_strand_numbers = [strand for strands in sheet_strands for strand in strands]
        count = 0
        while len(sheet_strand_numbers) != len(set(sheet_strand_numbers)):
            sheet_strand_numbers = [strand for strands in sheet_strands for strand in strands]
            strands = sheet_strands[count]
            track = 1
            sheet_strands_copy = copy.copy(sheet_strands)
            for index in range(0, len(sheet_strands)):
                if count != index:
                    intersect = set(strands).intersection(set(sheet_strands[index]))
                    if len(intersect) > 0:
                        sheet_strands_copy[count] = list(set(sheet_strands[count])|set(sheet_strands[index]))
                        sheet_strands_copy[index] = ''
                        sheet_strands = [x for x in sheet_strands_copy if x != '']
                        break
                    elif len(intersect) == 0:
                        track = track + 1

            if track == len(sheet_strands):
                track = 0
                count = count + 1

        for row in range(dssp_df.shape[0]):
            for strands in sheet_strands:
                if str(dssp_df['STRAND_NUM'][row]) in strands:
                    dssp_df['SHEET_NUM'][row] = (sheet_strands.index(strands)+1)
                    break

        # Separates the beta-strands identified by DSSP into sheets, and works
        # out the strand interactions and thus loop connections both within and
        # between sheets
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
            pos = nx.spring_layout(G)

            # Draws plot of the network of interacting strands
            plt.clf()
            nx.draw_networkx(G, pos=pos, with_labels=True)
            nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=interacting_strands)
            plt.savefig('{}_sheet_{}_network.png'.format(self.pdb_name, sheet))

































 g
