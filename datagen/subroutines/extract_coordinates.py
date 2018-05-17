
import os
import random
import string
import shutil
import pandas as pd
from collections import OrderedDict
from difflib import SequenceMatcher
if __name__ == 'subroutines.extract_coordinates':
    from subroutines.run_stages import run_stages
    from subroutines.variables import gen_amino_acids_dict
else:
    from datagen.subroutines.run_stages import run_stages
    from datagen.subroutines.variables import gen_amino_acids_dict


class extract_beta_structure_coords(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    def gen_cdhit_dict(self, cdhit_output, filtered_domain_df):
        # Loads the list of FASTA sequences generated by the CD-HIT web server
        fasta_list = []
        with open(cdhit_output, 'r') as chain_entries_file:
            for seq in chain_entries_file:
                if not seq.startswith(('>', '#')):
                    fasta_list.append(seq.replace('\n', ''))

        # For each of the sequences returned by the CD-HIT web server, selects
        # a random entry from the filtered dataframe of CATH domains from which
        # to select coordinates
        df_index_list = []
        for seq in fasta_list:
            df_index_sub_list = []
            for row in range(filtered_domain_df.shape[0]):
                if seq == filtered_domain_df['DSEQS'][row]:
                    df_index_sub_list.append(row)

            rand_num = random.randint(0, len(df_index_sub_list)-1)
            index = df_index_sub_list[rand_num]
            df_index_list.append(index)

        # Filters dataframe further to retain only the domains selected in the
        # previous step
        cdhit_domain_df = filtered_domain_df.iloc[df_index_list]
        cdhit_domain_df = cdhit_domain_df.reset_index(drop=True)
        return cdhit_domain_df

    def copy_biological_assembly_pdb(self, cdhit_domain_df):
        # Copies the parent biological assembly PDB file of every retained
        # structure to the output files directory. NOTE I assume that the
        # asymmetric unit containing the CATH domain is listed first in model 1
        # in the biological assembly PDB file
        unprocessed_list = []

        for row in range(cdhit_domain_df.shape[0]):
            pdb_code = cdhit_domain_df['PDB_CODE'][row]
            if not os.path.isfile('Biological_assemblies/{}.pdb'.format(pdb_code)):
                try:
                    shutil.copy(
                        '{}{}/{}.pdb'.format(self.pdb_ba_database, pdb_code[1:3], pdb_code),
                        'Biological_assemblies/{}.pdb'.format(pdb_code)
                    )
                except FileNotFoundError:
                    unprocessed_list.append(cdhit_domain_df['DOMAIN_ID'][row])

        cdhit_domain_df = cdhit_domain_df[~cdhit_domain_df['DOMAIN_ID'].isin(unprocessed_list)]
        cdhit_domain_df = cdhit_domain_df.reset_index(drop=True)

        with open('Unprocessed_domains.txt', 'a') as unprocessed_file:
            unprocessed_file.write('\n\nBiological assembly file not present:\n')
            for domain_id in set(unprocessed_list):
                unprocessed_file.write('{}\n'.format(domain_id))

        return cdhit_domain_df

    def get_xyz_coords(self, cdhit_domain_df):
        # Extends the filtered (for resolution, R_factor (working value) and
        # sequence redundancy) dataframe to list the xyz coordinates of each
        # segment sequence (SSEQS)
        amino_acids_dict = gen_amino_acids_dict()
        all_atoms_dfs_dict = OrderedDict()
        domain_residue_list = []
        unprocessed_list = []

        for row in range(cdhit_domain_df.shape[0]):
            first_chain = ''
            count = 0
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
            chain_num_ins = []
            lines = []

            # Extracts ATOM / HETATM records from biological assembly PDB file
            # copied to the output files directory in the previous step
            pdb_code = cdhit_domain_df['PDB_CODE'][row]

            print('Obtaining ATOM / HETATM records for {}'.format(pdb_code))
            print('{:0.2f}%'.format(((row+1)/cdhit_domain_df.shape[0])*100))

            with open('Biological_assemblies/{}.pdb'.format(pdb_code), 'r') as pdb_file:
                pdb_file_lines = [line.strip('\n') for line in pdb_file if
                                  line[0:6].strip() in ['ATOM', 'HETATM', 'TER']]
            pdb_file_lines.append('TER'.ljust(80))
            pdb_file_lines.append('TER'.ljust(80))  # Second TER line ensures
            # that entire sequence segment is extracted if TER cards have not
            # been placed correctly in the input PDB file

            # For each segment sequence in the domain, makes a list of all
            # sequences in the input PDB file that lie between the recorded
            # start and stop residue numbers and have the same chain id
            for index_1, segment in enumerate(cdhit_domain_df['SSEQS'][row]):
                sequences = []
                indices = []

                start = cdhit_domain_df['SSEQS_START_STOP'][row][index_1][0].replace('START=', '')
                stop = cdhit_domain_df['SSEQS_START_STOP'][row][index_1][1].replace('STOP=', '')
                start_seq = False
                stop_seq = False
                sequence = ''
                index_list = []

                for index_2, line in enumerate(pdb_file_lines):
                    if index_2 != (len(pdb_file_lines)-1):
                        if (line[22: 27].strip() == start
                                and line[21:22] == cdhit_domain_df['CHAIN'][row]
                                ):
                            start_seq = True

                        if start_seq is True and stop_seq is False:
                            index_list.append(index_2)
                            if (line[22:27].strip() != pdb_file_lines[index_2+1][22:27].strip()
                                    or pdb_file_lines[index_2+1][0:3] == 'TER'
                                    ):
                                if line[17:20].strip() in amino_acids_dict:
                                    sequence = sequence + amino_acids_dict[line[17:20].strip()]
                        elif stop_seq is True:
                            sequences.append(sequence)
                            indices.append(index_list)
                            sequence = ''
                            index_list = []
                            start_seq = False
                            stop_seq = False
                            continue

                        if (pdb_file_lines[index_2+1][0:3] == 'TER'
                                or (line[22:27].strip() == stop
                                    and line[21:22] == cdhit_domain_df['CHAIN'][row]
                                    and pdb_file_lines[index_2+1][22:27].strip() != stop
                                    )
                                ):
                            stop_seq = True

                # Selects the first identified sequence from the input PDB that
                # shares greater than 95% sequence similarity with the domain
                # segment sequence in question
                count += 1
                sequence_identified = False
                for index_3, sequence in enumerate(sequences):
                    sseqs_list = []
                    similarity = SequenceMatcher(a=segment, b=sequence).ratio()
                    if similarity >= 0.95:
                        # Ensures that all selected SSEQS are in the same chain
                        if count == 1:
                            first_chain = pdb_file_lines[indices[index_3][0]][21:22].strip()
                        elif count > 1:
                            new_chain = pdb_file_lines[indices[index_3][0]][21:22].strip()
                            if new_chain != first_chain:
                                break

                        sequence_identified = True

                        for index_4 in indices[index_3]:
                            sseqs_list.append(pdb_file_lines[index_4][21:27].replace(' ', ''))

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
                            chain_num_ins.append(pdb_file_lines[index_4][21:27].replace(' ', ''))
                            # Removes alternate conformer labels
                            line_start = pdb_file_lines[index_4][:16]
                            line_end = pdb_file_lines[index_4][17:].strip('\n')
                            lines.append(line_start + ' ' + line_end)

                        residue_list.extend(sseqs_list)
                        break

                if sequence_identified is False:
                    unprocessed_list.append('{}'.format(cdhit_domain_df['DOMAIN_ID'][row]))
                    break

            # Makes a dataframe of the PDB information for the domain sequence
            # if each of its SSEQS were identified in the input PDB file
            if not cdhit_domain_df['DOMAIN_ID'][row] in unprocessed_list:
                pdb_df = pd.DataFrame(OrderedDict({'PDB_FILE_LINES': lines,
                                                   'REC': rec,
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
                                                   'CHARGE': charge,
                                                   'RES_ID': chain_num_ins}))
                all_atoms_dfs_dict[cdhit_domain_df['DOMAIN_ID'][row]] = pdb_df

            # Makes a list of residue numbers of the domain sequence if
            # successfully identified in the input PDB file
            chain_num_list = []
            for chain_num in residue_list:
                if chain_num not in chain_num_list:
                    chain_num_list.append(chain_num)
            if not cdhit_domain_df['DOMAIN_ID'][row] in unprocessed_list:
                domain_residue_list.append(chain_num_list)

        # Lists all segment sequences that could not be processed owing to
        # insufficient similarity between the FASTA sequence listed in the
        # CATH_domain_desc_v_4_2_0.txt and the FASTA sequence extracted from
        # the PDB file in the unprocessed structures file
        with open('Unprocessed_domains.txt', 'a') as unprocessed_file:
            unprocessed_file.write('\n\n{} domain unable to be identified in '
                                   'PDB file:\n'.format(self.code))
            for domain_id in set(unprocessed_list):
                unprocessed_file.write('{}\n'.format(domain_id))

        cdhit_domain_df = cdhit_domain_df[~cdhit_domain_df['DOMAIN_ID'].isin(unprocessed_list)]
        cdhit_domain_df = cdhit_domain_df.reset_index(drop=True)

        # Appends column of residue numbers to the input dataframe
        domain_residue_df = pd.DataFrame({'CHAIN_NUM': domain_residue_list})
        cdhit_domain_df = pd.concat([cdhit_domain_df, domain_residue_df], axis=1)

        return cdhit_domain_df, all_atoms_dfs_dict

    def remove_alternate_conformers(self, all_atoms_dfs_dict):
        # Retains only the most probable alternate conformers. In the case
        # where more than one conformer is equally most probable, (only) that
        # which is listed first in the input PDB file is retained. Also removes
        # all hydrogen atoms (for fair comparison of e.g. solvent accessible
        # surface area, and nearest neighbouring residues, etc., between
        # different structures).
        for domain_id in list(all_atoms_dfs_dict.keys()):
            pdb_df = all_atoms_dfs_dict[domain_id]

            alternate_conformers_chainresnum = []
            alternate_conformers_label = []
            alternate_conformers_occupancy = []

            for row in range(pdb_df.shape[0]):
                if pdb_df['CONFORMER'][row] != '':
                    alternate_conformers_chainresnum.append(pdb_df['RES_ID'][row])
                    alternate_conformers_label.append(pdb_df['CONFORMER'][row])
                    alternate_conformers_occupancy.append(float(pdb_df['OCC'][row]))

            df = pd.DataFrame({'chainresnum': alternate_conformers_chainresnum,
                               'conformer': alternate_conformers_label,
                               'occupancy': alternate_conformers_occupancy})
            df = df.drop_duplicates()
            chainresnum = df['chainresnum'].tolist()
            conformer = df['conformer'].tolist()
            occupancy = df['occupancy'].tolist()

            alternate_conformers = {}
            chainresnum_set = set(chainresnum)
            for number in chainresnum_set:
                indices = []
                a = 'A1'
                b = 'A'
                c = 0
                for index_1, value_1 in enumerate(chainresnum):
                    if value_1 == number:
                        indices.append(index_1)
                for index_2 in indices:
                    if occupancy[index_2] > c:
                        a = chainresnum[index_2]
                        b = conformer[index_2]
                        c = occupancy[index_2]
                alternate_conformers[a] = b

            for row in range(pdb_df.shape[0]):
                # Removes alternate conformers
                if (pdb_df['RES_ID'][row] in alternate_conformers
                        and pdb_df['CONFORMER'][row] not in
                        [alternate_conformers[pdb_df['RES_ID'][row]], '']
                        ):
                    pdb_df.loc[row, 'REC'] = None
                # Removes hydrogens
                if pdb_df['ELEMENT'][row] == 'H':
                    pdb_df.loc[row, 'REC'] = None

            pdb_df = pdb_df[pdb_df['REC'].notnull()]
            pdb_df = pdb_df.reset_index(drop=True)

            all_atoms_dfs_dict[domain_id] = pdb_df

        return all_atoms_dfs_dict
