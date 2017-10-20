
import os
import shutil
import copy
import requests
import pandas as pd
from difflib import SequenceMatcher

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

# Generates a list of the domain descriptions provided in
# CATH_domain_description_v_4_1_0.txt
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

# Filters the domain descriptions list for beta-sandwiches, picking out PDB
# accession codes and sequences (whose values are stored in a DataFrame).
domain_pdb_ids = []
domain_chains = []
domain_cathcodes = []
domain_dseqs = []
domain_segments = []
domain_sseqs = []
domain_sseqs_start_stop = []
for domain in domains_description:
    if 'CATHCODE  2.60' in domain:
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
            if line.startswith('DOMAIN'):
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

# Extends the DataFrame to list the xyz coordinates of each segment sequence (SSEQS)
domain_xyz = []
unprocessed_list = []
for row in range(0, len(domain_pdb_ids)):
    print('Downloading {} from the RCSB PDB website'.format(domain_dict['PDB_CODE'][row]))
    url = 'http://www.rcsb.org/pdb/files/{}.pdb'.format(domain_dict['PDB_CODE'][row].upper())
    pdb_file_lines = requests.get(url).text
    pdb_file_lines = pdb_file_lines.split('\n')
    pdb_file_lines = [line for line in pdb_file_lines
                      if line[0:6].strip() in ['ATOM', 'HETATM', 'TER']]

    for line in pdb_file_lines:
        whitespace_remv_line = line.replace(' ', '')
        if whitespace_remv_line.startswith('EXPDTATHEORETICAL'):
            pdb_file_lines = []
            break

    if not pdb_file_lines:
        unprocessed_list.append('{}\n'.format(domain_dict['PDB_CODE'][row]))
        break
    else:
        pdb_file_lines.append('TER'.ljust(80))

    for index_1, segment in enumerate(domain_dict['SSEQS'][row]):
        sequences = []
        indices = []

        start = domain_dict['SSEQS_START_STOP'][row][index_1][0].replace('START=', '')
        stop = domain_dict['SSEQS_START_STOP'][row][index_1][1].replace('STOP=', '')
        start_seq = False
        stop_seq = False
        sequence = ''
        index = []

        for index_2, line in enumerate(pdb_file_lines):
            if index_2 != (len(pdb_file_lines)-1):
                if line[22:27].strip() == start and line[21:22] == domain_dict['CHAIN'][row]:
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
                        and line[21:22] == domain_dict['CHAIN'][row]
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
                for index_4 in indices[index_3]:
                    x = float(pdb_file_lines[index_4][30:38].strip())
                    y = float(pdb_file_lines[index_4][38:46].strip())
                    z = float(pdb_file_lines[index_4][46:54].strip())
                    xyz.append([x, y, z])

            if sequence_identified is True:
                break

        domain_xyz.append(xyz)

        if similarity <= 0.95:
            unprocessed_list.append('{}\n'.format(domain_dict['PDB_CODE'][row]))

with open('Unprocessed_beta_sandwich_pdbs.txt', 'w') as unprocessed_file:
    for pdb_file in set(unprocessed_list):
        unprocessed_file.write('{}\n'.format(pdb_file))

domain_xyz = pd.DataFrame({'XYZ': domain_xyz})
domain_dict = pd.concat([domain_dict, domain_xyz], axis=1)
domain_dict.to_csv('CATH_2.60_domain_desc_DataFrame_v_4_1_0.csv')
domain_dict.to_pickle('CATH_2.60_domain_desc_DataFrame_v_4_1_0.pkl')
