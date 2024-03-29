

import copy
import os
import networkx as nx
import pandas as pd
import numpy as np
from collections import OrderedDict
if __name__ == 'subroutines.output_dataframe':
    from subroutines.run_stages import run_stages
    from subroutines.twist_bend_shear import find_strand_twist, find_sheet_shear
    from subroutines.variables import gen_amino_acids_dict
else:
    from datagen.subroutines.run_stages import run_stages
    from datagen.subroutines.twist_bend_shear import find_strand_twist, find_sheet_shear
    from datagen.subroutines.variables import gen_amino_acids_dict


class output_calcs():

    def make_res_id_to_fasta_dict(domain_id, amino_acids_dict):
        # Creates dictionary of residue ids and their corresponding 1 letter
        # amino acid codes from the input dataframe (dssp_df)

        # Extracts list of consecutive res_ids from biological assembly PDB file
        with open('Parent_assemblies/{}.pdb'.format(domain_id), 'r') as pdb_file:
            pdb_file_lines = pdb_file.readlines()
        pdb_file_lines = [line for line in pdb_file_lines if line[0:6].strip()
                          in ['ATOM', 'HETATM']]
        consec_res_id_list = []
        consec_res_name_list = []
        for line in pdb_file_lines:
            chain_res_num = line[21:27].replace(' ', '')
            res_name = line[17:20].strip()
            if not chain_res_num in consec_res_id_list:
                consec_res_id_list.append(chain_res_num)
                consec_res_name_list.append(res_name)

        # Generates res_id to FASTA code dictionary
        fasta_list = ['']*len(consec_res_id_list)
        for index, res_id in enumerate(consec_res_id_list):
            res_name = consec_res_name_list[index]
            if res_name in list(amino_acids_dict.keys()):
                fasta = amino_acids_dict[res_name]
            else:
                fasta = 'X'
            fasta_list[index] = fasta

        res_id_to_fasta_dict = dict(zip(consec_res_id_list, fasta_list))

        return res_id_to_fasta_dict, consec_res_id_list

    def determine_strand_orientation(domain_id, strand_id, strand_df, tm,
                                     discard_non_tm, opm_database,
                                     unprocessed_list, unprocessed_strands):
        # Determines orientation of input strand in domain if that domain is in
        # the OPM database. If the N-terminus is in the periplasm
        # (Z-coordinate = negative) and C-terminus is extracellular
        # (Z-coordinate = positive), the strand orientation is reversed in the
        # output dataframe.

        pdb_code = domain_id[0:4]
        res_ids_pdb = strand_df['RES_ID'].tolist()
        resname_pdb = dict(zip(strand_df['RES_ID'].tolist(), strand_df['RESNAME'].tolist()))
        conformers_pdb = dict(zip(strand_df['RES_ID'].tolist(), strand_df['CONFORMER'].tolist()))
        strand_coordinates_opm = OrderedDict()
        in_database = False
        lower_bound = 0
        upper_bound = 0
        reverse = False

        if tm:
            if os.path.isfile('{}/{}.pdb'.format(opm_database, pdb_code)):
                print('Determining strand orientation in {} strand '
                      '{}'.format(domain_id, strand_id))
                in_database = True

                with open('{}/{}.pdb'.format(opm_database, pdb_code), 'r') as opm_file:
                    for line in opm_file:
                        if line[17:20].strip() != 'DUM':  # Don't put as 'and'
                            # statement with line below (see corresponding 'else' statement)
                            if line[0:6].strip() in ['ATOM', 'HETATM']:
                                chain_res_num = line[21:27].replace(' ', '')
                                if (
                                        line[12:16].strip() == 'CA'
                                    and chain_res_num in res_ids_pdb
                                    and line[17:20].strip() == resname_pdb[chain_res_num]  # NOTE
                                    # that OPM sometimes adds in alternate conformers that were
                                    # not present in the original PDB file - consequently, I ignore
                                    # the alternate conformer label, as I am only interested in the
                                    # Calpha position (since this would not be expected to vary
                                    # greatly between conformers unless those conformers are within
                                    # a consecutive sequence of alternate conformers). This means
                                    # that the code will always take the final listed conformer's
                                    # Calpha position.
                                    # and line[16:17].strip() == conformers_pdb[chain_res_num]
                                ):
                                    strand_coordinates_opm[chain_res_num] = float(line[46:54])
                        else:
                            if float(line[46:54]) > upper_bound:
                                upper_bound = float(line[46:54])
                            elif float(line[46:54]) < lower_bound:
                                lower_bound = float(line[46:54])

                z_coordinates = list(strand_coordinates_opm.values())
                # Ensures that all residues in input domain are present in the OPM
                # structure
                if len(z_coordinates) != len(res_ids_pdb):
                    strand_coordinates_opm = OrderedDict()
                    unprocessed_list.append(domain_id)
                else:
                    # Orients strands from the extracellular environment to the
                    # periplasm
                    if z_coordinates[0] < z_coordinates[-1]:
                        reverse = True

                # Ensures that input domain is transmembrane - NOTE that one
                # non-transmembrane strand does not necessarily mean that the
                # entire domain is non-transmembrane, so do not merge this
                # statement with the one above!
                if not domain_id in unprocessed_list:
                    if not any(lower_bound < z < upper_bound for z in z_coordinates):
                        unprocessed_strands.append(strand_id)

            else:
                if discard_non_tm is True:
                    unprocessed_list.append(domain_id)

        return (in_database, reverse, res_ids_pdb, strand_coordinates_opm,
                lower_bound, upper_bound, unprocessed_list,
                unprocessed_strands)

    def determine_tm_or_ext(domain_id, strand_id, res_ids_list, in_database,
                            strand_coordinates, lower_bound, upper_bound,
                            unprocessed_list):
        # If the parent structure is in the OPM database, labels residues in an
        # input strand as either 'transmembrane', 'periplasm' or 'extracellular'
        print('Calculating positions of strands in membrane in '
              '{}'.format(domain_id))

        tm_ext_list = ['']*len(res_ids_list)

        if (
            in_database is True
            and domain_id not in unprocessed_list
        ):
            for index, res_id in enumerate(res_ids_list):
                z_coord = strand_coordinates[res_id]
                if z_coord < lower_bound:
                    tm_ext_list[index] = 'periplasm'
                elif z_coord > upper_bound:
                    tm_ext_list[index] = 'extracellular'
                else:
                    tm_ext_list[index] = 'transmembrane'

        return tm_ext_list

    def determine_strand_position_sandwich(strand_df):
        # Calculates displacement (both absolute and percentage), measured in
        # number of residues, of each residue from the centre of the strand
        strand_abs_displacement = ['']*strand_df.shape[0]
        strand_percentage_displacement = ['']*strand_df.shape[0]
        strand_centre = (strand_df.shape[0] - 1) / 2

        for row in range(strand_df.shape[0]):
            diff = abs(row - strand_centre)
            strand_abs_displacement[row] = diff

            percentage = (abs(row - strand_centre) / strand_centre) * 100
            strand_percentage_displacement[row] = round(percentage, 1)

        return strand_abs_displacement, strand_percentage_displacement

    def determine_strand_position_barrel(strand_df, tm_ext_list):
        # Calculates position of residues within the barrel TM region (position
        # is measured as a percentage of the total number of TM residues in the
        # strand, with the extracellular environment and the periplasm
        # corresponding to 0% and 100%, respectively)
        membrane_loc = ['']*strand_df.shape[0]

        if 'transmembrane' in tm_ext_list:
            max_depth = tm_ext_list.count('transmembrane')

            count = 0
            for row in range(strand_df.shape[0]):
                if tm_ext_list[row] == 'transmembrane':
                    count += 1
                    percentage = (count / max_depth) * 100
                    membrane_loc[row] = round(percentage, 1)

        return membrane_loc

    def convert_dssp_num_to_res_id(strand_df, dssp_df, dssp_to_pdb_sub_dict):
        # Converts bridge pairs extracted from DSSP file to res_ids. Bridge
        # pairs are split into hydrogen bonding and non-hydrogen bonding pairs.
        hb_list = ['']*strand_df.shape[0]
        nhb_list = ['']*strand_df.shape[0]
        bridge_pairs_list = []
        for row in range(strand_df.shape[0]):
            bridges = strand_df['BRIDGE_PAIRS'][row]
            bridges_renum = []
            for dssp_num in bridges:
                if dssp_num in list(dssp_to_pdb_sub_dict.keys()):
                    bridges_renum.append(dssp_to_pdb_sub_dict[dssp_num])
                else:
                    bridges_renum.append('')
            bridge_pairs_list.append(bridges_renum)

            hb_sub_list = []
            nhb_sub_list = []
            for index, pdb_num in enumerate(bridges_renum):
                if pdb_num != '' and strand_df['ORIENTATION'][row][index] == 'A':
                    if strand_df['HBOND_MC_MC_all'][row].count(pdb_num) == 2:
                        hb_sub_list.append(pdb_num)
                    elif not pdb_num in strand_df['HBOND_MC_MC_all'][row]:
                        nhb_sub_list.append(pdb_num)
            hb_list[row] = hb_sub_list
            nhb_list[row] = nhb_sub_list

        return hb_list, nhb_list, bridge_pairs_list

    def find_minus_plus_residues(domain_id, strand_df, dssp_df, displacement,
                                 consec_res_id_list, res_id_to_fasta_dict):
        # Generates lists of residues +/- x (where x is a user-specified number
        # of residues) from each input residue. Currently will not exclude
        # small molecule HETATM such as water if the chain_id is the same as
        # the previous residue, however it will be listed with res_id X (and so
        # excluded from further analysis anyway)
        minus_list = []
        plus_list = []
        minus_list_fasta = []
        plus_list_fasta = []

        dssp_df_res_ids = dssp_df['RES_ID'].tolist()
        strand_ids = dssp_df['STRAND_NUM'].tolist()

        for row in range(strand_df.shape[0]):
            res_id = strand_df['RES_ID'][row]
            res_id_index = consec_res_id_list.index(res_id)
            strand_id = strand_df['STRAND_NUM'][row]

            try:
                minus_res = consec_res_id_list[res_id_index - displacement]
                plus_res = consec_res_id_list[res_id_index + displacement]

                try:
                    minus_res_strand_id = strand_ids[dssp_df_res_ids.index(minus_res)]
                except ValueError:
                    minus_res_strand_id = ''

                if (
                        res_id[0:1] != minus_res[0:1]
                    or not minus_res_strand_id in [strand_id, '']
                ):  # Prevents residues in different chains and/or strands from
                    # being considered as consecutive residues
                    minus_res = ''
                    minus_res_fasta = ''
                else:
                    minus_res_fasta = res_id_to_fasta_dict[minus_res]

                try:
                    plus_res_strand_id = strand_ids[dssp_df_res_ids.index(plus_res)]
                except ValueError:
                    plus_res_strand_id = ''

                if (
                        res_id[0:1] != plus_res[0:1]
                    or not plus_res_strand_id in [strand_id, '']
                ):  # Prevents residues in different chains and/or strands from
                    # being considered as consecutive residues
                    plus_res = ''
                    plus_res_fasta = ''
                else:
                    plus_res_fasta = res_id_to_fasta_dict[plus_res]
            except IndexError:
                minus_res = ''
                plus_res = ''
                minus_res_fasta = ''
                plus_res_fasta = ''

            minus_list.append(minus_res)
            plus_list.append(plus_res)
            minus_list_fasta.append(minus_res_fasta)
            plus_list_fasta.append(plus_res_fasta)

        return minus_list, plus_list, minus_list_fasta, plus_list_fasta

    def convert_res_id_to_fasta(sub_list, res_id_to_fasta_dict, nested, chain_id):
        # Converts input list of residue ids (= chain + residue number +
        # insertion code) into a list of 1 letter amino acid codes. Separates
        # interactions by chainID.
        intra_chain_fasta_list = []
        inter_chain_fasta_list = []

        for item in sub_list:
            if nested is False:
                try:
                    fasta = res_id_to_fasta_dict[item]
                except KeyError:
                    fasta = ''

                if item[0:1] == chain_id:
                    intra_chain_fasta_list.append(fasta)
                else:
                    inter_chain_fasta_list.append(fasta)

            elif nested is True:
                intra_fasta = []
                inter_fasta = []

                for res in item:
                    try:
                        fasta_indv = res_id_to_fasta_dict[res]
                    except KeyError:
                        fasta_indv = ''

                    if res[0:1] == chain_id:
                        intra_fasta.append(fasta_indv)
                    else:
                        inter_fasta.append(fasta_indv)

                intra_chain_fasta_list.append(intra_fasta)
                inter_chain_fasta_list.append(inter_fasta)

        return intra_chain_fasta_list, inter_chain_fasta_list


def reverse_strand_lists(sub_list, reverse):
    if reverse is True:
        sub_list = sub_list[::-1]
    return sub_list


def append_to_output_lists(sub_list, long_list, res_ids_list, strand_or_res):
    if strand_or_res == 'strand':
        long_list.append(sub_list)
    elif strand_or_res == 'res':
        if isinstance(sub_list, list):
            long_list += sub_list
        else:
            long_list += [sub_list]*len(res_ids_list)

    return long_list


def reverse_and_append(residue_property, strand_df, long_list, reverse,
                       res_ids_list, strand_or_res):
    sub_list = strand_df[residue_property].tolist()
    sub_list = reverse_strand_lists(sub_list, reverse)
    long_list = append_to_output_lists(sub_list, long_list, res_ids_list, strand_or_res)

    return long_list, sub_list


def gen_interaction_lists(properties_list, interaction_type_abbrev,
                          interaction_type_long, strand_df, reverse,
                          res_ids_list, strand_or_res, res_id_to_fasta_dict,
                          chain_id):
    if interaction_type_abbrev == 'PIPISTACK':
        chain_combos = ['', '_p', '_l', '_n', '_t',
                        '_mc_mc', '_mc_mc_p', '_mc_mc_l', '_mc_mc_n', '_mc_mc_t',
                        '_mc_sc', '_mc_sc_p', '_mc_sc_l', '_mc_sc_n', '_mc_sc_t',
                        '_sc_mc', '_sc_mc_p', '_sc_mc_l', '_sc_mc_n', '_sc_mc_t',
                        '_sc_sc', '_sc_sc_p', '_sc_sc_l', '_sc_sc_n', '_sc_sc_t']
    else:
        chain_combos = ['', '_mc_mc', '_mc_sc', '_sc_mc', '_sc_sc']

    for chain_combo in chain_combos:
        (properties_list['{}{}'.format(interaction_type_long, chain_combo)],
         sub_list) = reverse_and_append(
            '{}{}'.format(interaction_type_abbrev, chain_combo.upper()),
            strand_df, properties_list['{}{}'.format(interaction_type_long, chain_combo)],
            reverse, res_ids_list, strand_or_res
        )

        (fasta_sub_list_intra, fasta_sub_list_inter
         ) = output_calcs.convert_res_id_to_fasta(
            sub_list, res_id_to_fasta_dict, True, chain_id
        )
        (properties_list['{}_fasta_intra{}'.format(interaction_type_long, chain_combo)]
         ) = append_to_output_lists(
            fasta_sub_list_intra,
            properties_list['{}_fasta_intra{}'.format(interaction_type_long, chain_combo)],
            res_ids_list, strand_or_res
        )
        (properties_list['{}_fasta_inter{}'.format(interaction_type_long, chain_combo)]
         ) = append_to_output_lists(
            fasta_sub_list_inter,
            properties_list['{}_fasta_inter{}'.format(interaction_type_long, chain_combo)],
            res_ids_list, strand_or_res
        )

    return properties_list


class gen_output(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    def identify_edge_central(self, domain_sheets_dict, sec_struct_dfs_dict):
        # Uses domain_networks_dict to identify whether strands are edge or
        # central
        for domain_id in list(sec_struct_dfs_dict.keys()):
            print('Identifying edge and central strands in {}'.format(domain_id))

            networks = [network for key, network in domain_sheets_dict.items()
                        if domain_id in key]

            G = networks[0]
            for num in range(1, len(networks)):
                H = networks[num]
                G = nx.compose(G, H)

            dssp_df = sec_struct_dfs_dict[domain_id]

            edge_or_central_dict = OrderedDict()
            df_strands = [strand for strand in set(dssp_df['STRAND_NUM'].tolist())
                          if strand != '']

            for strand_num in df_strands:
                num_of_edges = len(list(G.neighbors(strand_num)))
                if num_of_edges == 1:
                    edge_id = 'edge'
                elif num_of_edges > 1:
                    edge_id = 'central'
                edge_or_central_dict[strand_num] = edge_id

            edge_or_central_list = ['']*dssp_df.shape[0]
            for row in range(dssp_df.shape[0]):
                if (
                    dssp_df['STRAND_NUM'][row] in list(edge_or_central_dict.keys())
                    and dssp_df['ATMNAME'][row] == 'CA'
                ):
                    edge_or_central_list[row] = edge_or_central_dict[dssp_df['STRAND_NUM'][row]]
            df_edge_or_central = pd.DataFrame({'EDGE_OR_CNTRL': edge_or_central_list})
            dssp_df = pd.concat([dssp_df, df_edge_or_central], axis=1)
            sec_struct_dfs_dict[domain_id] = dssp_df

        return sec_struct_dfs_dict

    def write_beta_strand_dataframe(self, strand_or_res, sec_struct_dfs_dict,
                                    domain_sheets_dict, dssp_to_pdb_dict,
                                    tilt_angles, strand_numbers):
        # Generates dataframes of residues and strands in the retained domains.
        if __name__ == 'subroutines.output_dataframe':
            from subroutines.output_dataframe import output_calcs
        else:
            from datagen.subroutines.output_dataframe import output_calcs

        # Generates dictionary of amino acid 3 and 1 letter codes
        amino_acids_dict = gen_amino_acids_dict()

        # Initialises lists of properties to be displayed in the dataframe
        # columns
        properties_list = OrderedDict({'domain_ids': [],
                                       'domain_strand_ids': [],
                                       'sheet_number': [],
                                       'tilt_angle': [],
                                       'total_strand_number': [],
                                       'indv_strand_number': [],
                                       'res_ids': [],
                                       'res_nums': [],
                                       'edge_or_central': [],
                                       'fasta_seq': [],
                                       'z_coords': [],
                                       'strand_abs_pos': [],
                                       'strand_percentage_pos': [],
                                       'strand_z_coords': [],
                                       'sandwich_z_coords': [],
                                       'bend_angles': [],
                                       'tm_pos': [],
                                       'int_ext': [],
                                       'core_surf': [],
                                       'buried_surface_area': [],
                                       'tm_ext': [],
                                       'omega': [],
                                       'phi': [],
                                       'psi': [],
                                       'chi': [],
                                       'solv_acsblty': [],
                                       'minus_1': [],
                                       'plus_1': [],
                                       'minus_2': [],
                                       'plus_2': [],
                                       'neighbours': [],
                                       'bridge_pairs': [],
                                       'hb_pairs': [],
                                       'nhb_pairs': [],
                                       'van_der_waals': [],
                                       'h_bonds': [],
                                       'ionic': [],
                                       'ss_bonds': [],
                                       'pi_pi_stacking': [],
                                       'pi_pi_stacking_p': [],
                                       'pi_pi_stacking_n': [],
                                       'pi_pi_stacking_l': [],
                                       'pi_pi_stacking_t': [],
                                       'cation_pi': [],
                                       'van_der_waals_mc_mc': [],
                                       'h_bonds_mc_mc': [],
                                       'ionic_mc_mc': [],
                                       'ss_bonds_mc_mc': [],
                                       'pi_pi_stacking_mc_mc': [],
                                       'pi_pi_stacking_mc_mc_p': [],
                                       'pi_pi_stacking_mc_mc_n': [],
                                       'pi_pi_stacking_mc_mc_l': [],
                                       'pi_pi_stacking_mc_mc_t': [],
                                       'cation_pi_mc_mc': [],
                                       'van_der_waals_mc_sc': [],
                                       'h_bonds_mc_sc': [],
                                       'ionic_mc_sc': [],
                                       'ss_bonds_mc_sc': [],
                                       'pi_pi_stacking_mc_sc': [],
                                       'pi_pi_stacking_mc_sc_p': [],
                                       'pi_pi_stacking_mc_sc_n': [],
                                       'pi_pi_stacking_mc_sc_l': [],
                                       'pi_pi_stacking_mc_sc_t': [],
                                       'cation_pi_mc_sc': [],
                                       'van_der_waals_sc_mc': [],
                                       'h_bonds_sc_mc': [],
                                       'ionic_sc_mc': [],
                                       'ss_bonds_sc_mc': [],
                                       'pi_pi_stacking_sc_mc': [],
                                       'pi_pi_stacking_sc_mc_p': [],
                                       'pi_pi_stacking_sc_mc_n': [],
                                       'pi_pi_stacking_sc_mc_l': [],
                                       'pi_pi_stacking_sc_mc_t': [],
                                       'cation_pi_sc_mc': [],
                                       'van_der_waals_sc_sc': [],
                                       'h_bonds_sc_sc': [],
                                       'ionic_sc_sc': [],
                                       'ss_bonds_sc_sc': [],
                                       'pi_pi_stacking_sc_sc': [],
                                       'pi_pi_stacking_sc_sc_p': [],
                                       'pi_pi_stacking_sc_sc_n': [],
                                       'pi_pi_stacking_sc_sc_l': [],
                                       'pi_pi_stacking_sc_sc_t': [],
                                       'cation_pi_sc_sc': [],
                                       'minus_1_fasta': [],
                                       'plus_1_fasta': [],
                                       'minus_2_fasta': [],
                                       'plus_2_fasta': [],
                                       'neighbours_fasta_intra': [],
                                       'bridge_pairs_fasta_intra': [],
                                       'hb_pairs_fasta_intra': [],
                                       'nhb_pairs_fasta_intra': [],
                                       'van_der_waals_fasta_intra': [],
                                       'h_bonds_fasta_intra': [],
                                       'ionic_fasta_intra': [],
                                       'ss_bonds_fasta_intra': [],
                                       'pi_pi_stacking_fasta_intra': [],
                                       'pi_pi_stacking_fasta_intra_p': [],
                                       'pi_pi_stacking_fasta_intra_n': [],
                                       'pi_pi_stacking_fasta_intra_l': [],
                                       'pi_pi_stacking_fasta_intra_t': [],
                                       'cation_pi_fasta_intra': [],
                                       'van_der_waals_fasta_intra_mc_mc': [],
                                       'h_bonds_fasta_intra_mc_mc': [],
                                       'ionic_fasta_intra_mc_mc': [],
                                       'ss_bonds_fasta_intra_mc_mc': [],
                                       'pi_pi_stacking_fasta_intra_mc_mc': [],
                                       'pi_pi_stacking_fasta_intra_mc_mc_p': [],
                                       'pi_pi_stacking_fasta_intra_mc_mc_n': [],
                                       'pi_pi_stacking_fasta_intra_mc_mc_l': [],
                                       'pi_pi_stacking_fasta_intra_mc_mc_t': [],
                                       'cation_pi_fasta_intra_mc_mc': [],
                                       'van_der_waals_fasta_intra_mc_sc': [],
                                       'h_bonds_fasta_intra_mc_sc': [],
                                       'ionic_fasta_intra_mc_sc': [],
                                       'ss_bonds_fasta_intra_mc_sc': [],
                                       'pi_pi_stacking_fasta_intra_mc_sc': [],
                                       'pi_pi_stacking_fasta_intra_mc_sc_p': [],
                                       'pi_pi_stacking_fasta_intra_mc_sc_n': [],
                                       'pi_pi_stacking_fasta_intra_mc_sc_l': [],
                                       'pi_pi_stacking_fasta_intra_mc_sc_t': [],
                                       'cation_pi_fasta_intra_mc_sc': [],
                                       'van_der_waals_fasta_intra_sc_mc': [],
                                       'h_bonds_fasta_intra_sc_mc': [],
                                       'ionic_fasta_intra_sc_mc': [],
                                       'ss_bonds_fasta_intra_sc_mc': [],
                                       'pi_pi_stacking_fasta_intra_sc_mc': [],
                                       'pi_pi_stacking_fasta_intra_sc_mc_p': [],
                                       'pi_pi_stacking_fasta_intra_sc_mc_n': [],
                                       'pi_pi_stacking_fasta_intra_sc_mc_l': [],
                                       'pi_pi_stacking_fasta_intra_sc_mc_t': [],
                                       'cation_pi_fasta_intra_sc_mc': [],
                                       'van_der_waals_fasta_intra_sc_sc': [],
                                       'h_bonds_fasta_intra_sc_sc': [],
                                       'ionic_fasta_intra_sc_sc': [],
                                       'ss_bonds_fasta_intra_sc_sc': [],
                                       'pi_pi_stacking_fasta_intra_sc_sc': [],
                                       'pi_pi_stacking_fasta_intra_sc_sc_p': [],
                                       'pi_pi_stacking_fasta_intra_sc_sc_n': [],
                                       'pi_pi_stacking_fasta_intra_sc_sc_l': [],
                                       'pi_pi_stacking_fasta_intra_sc_sc_t': [],
                                       'cation_pi_fasta_intra_sc_sc': [],
                                       'neighbours_fasta_inter': [],
                                       'bridge_pairs_fasta_inter': [],
                                       'hb_pairs_fasta_inter': [],
                                       'nhb_pairs_fasta_inter': [],
                                       'van_der_waals_fasta_inter': [],
                                       'h_bonds_fasta_inter': [],
                                       'ionic_fasta_inter': [],
                                       'ss_bonds_fasta_inter': [],
                                       'pi_pi_stacking_fasta_inter': [],
                                       'pi_pi_stacking_fasta_inter_p': [],
                                       'pi_pi_stacking_fasta_inter_n': [],
                                       'pi_pi_stacking_fasta_inter_l': [],
                                       'pi_pi_stacking_fasta_inter_t': [],
                                       'cation_pi_fasta_inter': [],
                                       'van_der_waals_fasta_inter_mc_mc': [],
                                       'h_bonds_fasta_inter_mc_mc': [],
                                       'ionic_fasta_inter_mc_mc': [],
                                       'ss_bonds_fasta_inter_mc_mc': [],
                                       'pi_pi_stacking_fasta_inter_mc_mc': [],
                                       'pi_pi_stacking_fasta_inter_mc_mc_p': [],
                                       'pi_pi_stacking_fasta_inter_mc_mc_n': [],
                                       'pi_pi_stacking_fasta_inter_mc_mc_l': [],
                                       'pi_pi_stacking_fasta_inter_mc_mc_t': [],
                                       'cation_pi_fasta_inter_mc_mc': [],
                                       'van_der_waals_fasta_inter_mc_sc': [],
                                       'h_bonds_fasta_inter_mc_sc': [],
                                       'ionic_fasta_inter_mc_sc': [],
                                       'ss_bonds_fasta_inter_mc_sc': [],
                                       'pi_pi_stacking_fasta_inter_mc_sc': [],
                                       'pi_pi_stacking_fasta_inter_mc_sc_p': [],
                                       'pi_pi_stacking_fasta_inter_mc_sc_n': [],
                                       'pi_pi_stacking_fasta_inter_mc_sc_l': [],
                                       'pi_pi_stacking_fasta_inter_mc_sc_t': [],
                                       'cation_pi_fasta_inter_mc_sc': [],
                                       'van_der_waals_fasta_inter_sc_mc': [],
                                       'h_bonds_fasta_inter_sc_mc': [],
                                       'ionic_fasta_inter_sc_mc': [],
                                       'ss_bonds_fasta_inter_sc_mc': [],
                                       'pi_pi_stacking_fasta_inter_sc_mc': [],
                                       'pi_pi_stacking_fasta_inter_sc_mc_p': [],
                                       'pi_pi_stacking_fasta_inter_sc_mc_n': [],
                                       'pi_pi_stacking_fasta_inter_sc_mc_l': [],
                                       'pi_pi_stacking_fasta_inter_sc_mc_t': [],
                                       'cation_pi_fasta_inter_sc_mc': [],
                                       'van_der_waals_fasta_inter_sc_sc': [],
                                       'h_bonds_fasta_inter_sc_sc': [],
                                       'ionic_fasta_inter_sc_sc': [],
                                       'ss_bonds_fasta_inter_sc_sc': [],
                                       'pi_pi_stacking_fasta_inter_sc_sc': [],
                                       'pi_pi_stacking_fasta_inter_sc_sc_p': [],
                                       'pi_pi_stacking_fasta_inter_sc_sc_n': [],
                                       'pi_pi_stacking_fasta_inter_sc_sc_l': [],
                                       'pi_pi_stacking_fasta_inter_sc_sc_t': [],
                                       'cation_pi_fasta_inter_sc_sc': []})

        unprocessed_list = []

        # Extracts property values for each domain_id
        for domain_id in list(sec_struct_dfs_dict.keys()):
            dssp_df = sec_struct_dfs_dict[domain_id]
            dssp_to_pdb_sub_dict = dssp_to_pdb_dict[domain_id]

            # Generates dictionary of residue ids (= chain + residue number +
            # insertion code) and names
            res_id_to_fasta_dict, consec_res_id_list = output_calcs.make_res_id_to_fasta_dict(
                domain_id, amino_acids_dict
            )

            strands_list = [strand_num for strand_num in set(dssp_df['STRAND_NUM'].tolist())
                            if strand_num != '']

            unprocessed_strands = []

            for strand_num in strands_list:
                strand_df = dssp_df[dssp_df['STRAND_NUM'] == strand_num]
                strand_df = strand_df.reset_index(drop=True)
                chain_id = strand_df['CHAIN'][0]

                # If the strand's parent structure is in the OPM database,
                # determines the orientation of the strand in the OM
                if self.code[0:4] in ['2.40']:
                    tm = True
                else:
                    tm = False
                (in_database, reverse, res_ids_list, strand_coordinates,
                 lower_bound, upper_bound, unprocessed_list,
                 unprocessed_strands
                 ) = output_calcs.determine_strand_orientation(
                    domain_id, strand_num, strand_df, tm, self.discard_non_tm,
                    self.opm_database, unprocessed_list, unprocessed_strands
                )

                # Lists parent domain IDs
                properties_list['domain_ids'] = append_to_output_lists(
                    '{}'.format(domain_id), properties_list['domain_ids'],
                    res_ids_list, strand_or_res
                )

                # Gives strand a unique ID
                properties_list['domain_strand_ids'] = append_to_output_lists(
                    '{}_strand_{}'.format(domain_id, strand_num),
                    properties_list['domain_strand_ids'], res_ids_list, strand_or_res
                )

                # Lists parent sheet of each strand
                sheet = strand_df['SHEET_NUM'].tolist()[0]
                properties_list['sheet_number'] = append_to_output_lists(
                    sheet, properties_list['sheet_number'], res_ids_list,
                    strand_or_res
                )

                # Lists tilt angle of strand
                if self.code[0:4] in ['2.40']:
                    if domain_id in list(tilt_angles.keys()):
                        value = tilt_angles[domain_id]
                    else:
                        value = ''

                    properties_list['tilt_angle'] = append_to_output_lists(
                        value, properties_list['tilt_angle'], res_ids_list,
                        strand_or_res
                    )

                # Lists number of strands in barrel
                if self.code[0:4] in ['2.40']:
                    properties_list['total_strand_number'] = append_to_output_lists(
                        strand_numbers[domain_id],
                        properties_list['total_strand_number'], res_ids_list,
                        strand_or_res
                    )

                # Lists parent strand number
                properties_list['indv_strand_number'] = append_to_output_lists(
                    strand_df['STRAND_NUM'].tolist()[0],
                    properties_list['indv_strand_number'], res_ids_list,
                    strand_or_res
                )

                # Determines whether the strand is an edge or central strand
                if self.code[0:4] in ['2.60']:
                    properties_list['edge_or_central'] = append_to_output_lists(
                        strand_df['EDGE_OR_CNTRL'].tolist()[0],
                        properties_list['edge_or_central'], res_ids_list,
                        strand_or_res
                    )

                # Lists residue IDs in strand
                res_ids_list = reverse_strand_lists(res_ids_list, reverse)
                properties_list['res_ids'] = append_to_output_lists(
                    res_ids_list, properties_list['res_ids'], res_ids_list,
                    strand_or_res
                )

                # Lists residue numbers in strand
                res_nums_list = strand_df['RESNUM'].tolist()
                res_nums_list = reverse_strand_lists(res_nums_list, reverse)
                properties_list['res_nums'] = append_to_output_lists(
                    res_nums_list, properties_list['res_nums'], res_ids_list,
                    strand_or_res
                )

                # Generates FASTA sequence of strand
                res_list = strand_df['RESNAME'].tolist()
                sequence = ''
                for res in res_list:
                    if res in list(amino_acids_dict.keys()):
                        sequence += amino_acids_dict[res]
                    else:
                        sequence += 'X'
                sequence = list(sequence)

                sequence = reverse_strand_lists(sequence, reverse)
                properties_list['fasta_seq'] = append_to_output_lists(
                    sequence, properties_list['fasta_seq'], res_ids_list,
                    strand_or_res
                )

                # Generates list of strand z-coordinates
                if self.code[0:4] in ['2.40']:
                    if strand_coordinates:
                        strand_z_coords = list(strand_coordinates.values())
                    else:
                        strand_z_coords = ['']*len(res_ids_list)

                    strand_z_coords = reverse_strand_lists(strand_z_coords, reverse)
                    properties_list['z_coords'] = append_to_output_lists(
                        strand_z_coords, properties_list['z_coords'],
                        res_ids_list, strand_or_res
                    )

                # Generates list of interior and exterior facing residues in
                # the strand
                properties_list['int_ext'], int_ext_sub_list = reverse_and_append(
                    'INT_EXT', strand_df, properties_list['int_ext'], reverse,
                    res_ids_list, strand_or_res
                )

                # Generates list of strand bend values
                properties_list['bend_angles'], bend_angle_sub_list = reverse_and_append(
                    'BEND', strand_df, properties_list['bend_angles'], reverse,
                    res_ids_list, strand_or_res
                )

                # Generates list of residues that form the beta-sandwich core
                if self.code[0:4] in ['2.60']:
                    properties_list['core_surf'], core_surf_sub_list = reverse_and_append(
                        'CORE_OR_SURFACE', strand_df, properties_list['core_surf'],
                        reverse, res_ids_list, strand_or_res
                    )

                # Generates list of percentage values for the surface area of
                # each residue that is buried in the sandwich as compared with
                # the individual parent sheet
                if self.code[0:4] in ['2.60']:
                    (properties_list['buried_surface_area'], buried_surface_area_sub_list
                     ) = reverse_and_append(
                        'BURIED_SURFACE_AREA(%)', strand_df,
                        properties_list['buried_surface_area'], reverse,
                        res_ids_list, strand_or_res
                    )

                # Generates list of transmembrane, periplasmic and extracellular
                # residues in the strand - MUST COME AFTER REVERSAL (OR NOT) OF
                # RES_IDS_LIST
                if self.code[0:4] in ['2.40']:
                    tm_ext_list = output_calcs.determine_tm_or_ext(
                        domain_id, strand_num, res_ids_list, in_database,
                        strand_coordinates, lower_bound, upper_bound,
                        unprocessed_list
                    )
                    # Don't need to reverse list as will match res_ids_list
                    # order (which has already been reversed)
                    properties_list['tm_ext'] = append_to_output_lists(
                        tm_ext_list, properties_list['tm_ext'], res_ids_list,
                        strand_or_res
                    )

                # Determines positions of residues relative to centre of strand
                # in sandwiches
                if self.code[0:4] in ['2.60']:
                    (strand_abs_displacement, strand_percentage_displacement
                     ) = output_calcs.determine_strand_position_sandwich(strand_df)

                    strand_abs_displacement = reverse_strand_lists(
                        strand_abs_displacement, reverse
                    )
                    properties_list['strand_abs_pos'] = append_to_output_lists(
                        strand_abs_displacement, properties_list['strand_abs_pos'],
                        res_ids_list, strand_or_res
                    )

                    strand_percentage_displacement = reverse_strand_lists(
                        strand_percentage_displacement, reverse
                    )
                    properties_list['strand_percentage_pos'] = append_to_output_lists(
                        strand_percentage_displacement,
                        properties_list['strand_percentage_pos'], res_ids_list,
                        strand_or_res
                    )

                    strand_z_coords = strand_df['STRAND_Z_COORDS'].tolist()
                    strand_z_coords = reverse_strand_lists(
                        strand_z_coords, reverse
                    )
                    properties_list['strand_z_coords'] = append_to_output_lists(
                        strand_z_coords, properties_list['strand_z_coords'],
                        res_ids_list, strand_or_res
                    )

                    sandwich_z_coords = strand_df['SANDWICH_Z_COORDS'].tolist()
                    sandwich_z_coords = reverse_strand_lists(
                        sandwich_z_coords, reverse
                    )
                    properties_list['sandwich_z_coords'] = append_to_output_lists(
                        sandwich_z_coords, properties_list['sandwich_z_coords'],
                        res_ids_list, strand_or_res
                    )

                # Determines positions of residues in TM region of strand in
                # beta-barrel - MUST COME AFTER TM / PERIPLASM / EXTRACELLULAR
                # CALCULATION
                if self.code[0:4] in ['2.40']:
                    membrane_loc = output_calcs.determine_strand_position_barrel(
                        strand_df, tm_ext_list
                    )

                    # Don't need to reverse list as will match tm_ext_list
                    # order (which has already been reversed)
                    properties_list['tm_pos'] = append_to_output_lists(
                        membrane_loc, properties_list['tm_pos'], res_ids_list,
                        strand_or_res
                    )

                # Generates list of torsion angles along the strand
                properties_list['omega'], omega_sub_list = reverse_and_append(
                    'OMEGA', strand_df, properties_list['omega'], reverse,
                    res_ids_list, strand_or_res
                )
                properties_list['phi'], phi_sub_list = reverse_and_append(
                    'PHI', strand_df, properties_list['phi'], reverse,
                    res_ids_list, strand_or_res
                )
                properties_list['psi'], psi_sub_list = reverse_and_append(
                    'PSI', strand_df, properties_list['psi'], reverse,
                    res_ids_list, strand_or_res
                )
                properties_list['chi'], chi_sub_list = reverse_and_append(
                    'CHI', strand_df, properties_list['chi'], reverse,
                    res_ids_list, strand_or_res
                )

                # Generates list of residue solvent accessibility along the
                # strand
                properties_list['solv_acsblty'], solv_acsblty_sub_list = reverse_and_append(
                    'SOLV_ACSBLTY', strand_df, properties_list['solv_acsblty'],
                    reverse, res_ids_list, strand_or_res
                )

                # Generates list of neighbouring residues
                properties_list['neighbours'], neighbours_sub_list = reverse_and_append(
                    'NEIGHBOURS', strand_df, properties_list['neighbours'], reverse,
                    res_ids_list, strand_or_res
                )

                (neighbours_fasta_sub_list_intra, neighbours_fasta_sub_list_inter
                 ) = output_calcs.convert_res_id_to_fasta(
                    neighbours_sub_list, res_id_to_fasta_dict, True, chain_id
                )
                properties_list['neighbours_fasta_intra'] = append_to_output_lists(
                    neighbours_fasta_sub_list_intra,
                    properties_list['neighbours_fasta_intra'], res_ids_list,
                    strand_or_res
                )
                properties_list['neighbours_fasta_inter'] = append_to_output_lists(
                    neighbours_fasta_sub_list_inter,
                    properties_list['neighbours_fasta_inter'],
                    res_ids_list, strand_or_res
                )

                # Generates list of residues forming backbone hydrogen-bonds
                # with the residue in question
                (hb_pairs_list, nhb_pairs_list, bridge_pairs_list
                 ) = output_calcs.convert_dssp_num_to_res_id(
                    strand_df, dssp_df, dssp_to_pdb_sub_dict
                )
                bridge_pairs_list = reverse_strand_lists(bridge_pairs_list, reverse)
                properties_list['bridge_pairs'] = append_to_output_lists(
                    bridge_pairs_list, properties_list['bridge_pairs'],
                    res_ids_list, strand_or_res
                )
                hb_pairs_list = reverse_strand_lists(hb_pairs_list, reverse)
                properties_list['hb_pairs'] = append_to_output_lists(
                    hb_pairs_list, properties_list['hb_pairs'], res_ids_list,
                    strand_or_res
                )
                nhb_pairs_list = reverse_strand_lists(nhb_pairs_list, reverse)
                properties_list['nhb_pairs'] = append_to_output_lists(
                    nhb_pairs_list, properties_list['nhb_pairs'],
                    res_ids_list, strand_or_res
                )

                (bridge_pairs_fasta_sub_list_intra, bridge_pairs_fasta_sub_list_inter
                 ) = output_calcs.convert_res_id_to_fasta(
                    bridge_pairs_list, res_id_to_fasta_dict, True, chain_id
                )
                properties_list['bridge_pairs_fasta_intra'] = append_to_output_lists(
                    bridge_pairs_fasta_sub_list_intra,
                    properties_list['bridge_pairs_fasta_intra'],
                    res_ids_list, strand_or_res
                )
                properties_list['bridge_pairs_fasta_inter'] = append_to_output_lists(
                    bridge_pairs_fasta_sub_list_inter,
                    properties_list['bridge_pairs_fasta_inter'],
                    res_ids_list, strand_or_res
                )

                (hb_pairs_fasta_sub_list_intra, hb_pairs_fasta_sub_list_inter
                 ) = output_calcs.convert_res_id_to_fasta(
                    hb_pairs_list, res_id_to_fasta_dict, True, chain_id
                )
                properties_list['hb_pairs_fasta_intra'] = append_to_output_lists(
                    hb_pairs_fasta_sub_list_intra,
                    properties_list['hb_pairs_fasta_intra'],
                    res_ids_list, strand_or_res
                )
                properties_list['hb_pairs_fasta_inter'] = append_to_output_lists(
                    hb_pairs_fasta_sub_list_inter,
                    properties_list['hb_pairs_fasta_inter'],
                    res_ids_list, strand_or_res
                )

                (nhb_pairs_fasta_sub_list_intra, nhb_pairs_fasta_sub_list_inter
                 ) = output_calcs.convert_res_id_to_fasta(
                    nhb_pairs_list, res_id_to_fasta_dict, True, chain_id
                )
                properties_list['nhb_pairs_fasta_intra'] = append_to_output_lists(
                    nhb_pairs_fasta_sub_list_intra,
                    properties_list['nhb_pairs_fasta_intra'],
                    res_ids_list, strand_or_res
                )
                properties_list['nhb_pairs_fasta_inter'] = append_to_output_lists(
                    nhb_pairs_fasta_sub_list_inter,
                    properties_list['nhb_pairs_fasta_inter'],
                    res_ids_list, strand_or_res
                )

                # Generates list of residues in van der Waals contact
                properties_list = gen_interaction_lists(
                    properties_list, 'VDW', 'van_der_waals', strand_df,
                    reverse, res_ids_list, strand_or_res, res_id_to_fasta_dict,
                    chain_id
                )

                # Generates list of residues that form hydrogen bonds with the
                # residue in question
                properties_list = gen_interaction_lists(
                    properties_list, 'HBOND', 'h_bonds', strand_df,
                    reverse, res_ids_list, strand_or_res, res_id_to_fasta_dict,
                    chain_id
                )

                # Generates list of residues that form ionic bonds with the
                # residue in question
                properties_list = gen_interaction_lists(
                    properties_list, 'IONIC', 'ionic', strand_df,
                    reverse, res_ids_list, strand_or_res, res_id_to_fasta_dict,
                    chain_id
                )

                # Generates list of disulfide bonds
                properties_list = gen_interaction_lists(
                    properties_list, 'SSBOND', 'ss_bonds', strand_df,
                    reverse, res_ids_list, strand_or_res, res_id_to_fasta_dict,
                    chain_id
                )

                # Generates list of residues that form pi-pi stacking
                # interactions with the residue in question
                properties_list = gen_interaction_lists(
                    properties_list, 'PIPISTACK', 'pi_pi_stacking', strand_df,
                    reverse, res_ids_list, strand_or_res, res_id_to_fasta_dict,
                    chain_id
                )

                # Generates list of residues that form cation_pi interactions
                # with the residue in question
                properties_list = gen_interaction_lists(
                    properties_list, 'PICATION', 'cation_pi', strand_df,
                    reverse, res_ids_list, strand_or_res, res_id_to_fasta_dict,
                    chain_id
                )

                # Generates list of residues in +/-1 and +/2 positions to the
                # residue in question
                (minus_1_list, plus_1_list, minus_1_list_fasta, plus_1_list_fasta
                ) = output_calcs.find_minus_plus_residues(
                    domain_id, strand_df, dssp_df, 1, consec_res_id_list,
                    res_id_to_fasta_dict
                )
                (minus_2_list, plus_2_list, minus_2_list_fasta, plus_2_list_fasta
                ) = output_calcs.find_minus_plus_residues(
                    domain_id, strand_df, dssp_df, 2, consec_res_id_list,
                    res_id_to_fasta_dict
                )

                minus_1_list = reverse_strand_lists(minus_1_list, reverse)
                properties_list['minus_1'] = append_to_output_lists(
                    minus_1_list, properties_list['minus_1'], res_ids_list,
                    strand_or_res
                )
                plus_1_list = reverse_strand_lists(plus_1_list, reverse)
                properties_list['plus_1'] = append_to_output_lists(
                    plus_1_list, properties_list['plus_1'], res_ids_list,
                    strand_or_res
                )
                minus_1_list_fasta = reverse_strand_lists(minus_1_list_fasta, reverse)
                properties_list['minus_1_fasta'] = append_to_output_lists(
                    minus_1_list_fasta, properties_list['minus_1_fasta'],
                    res_ids_list, strand_or_res
                )
                plus_1_list_fasta = reverse_strand_lists(plus_1_list_fasta, reverse)
                properties_list['plus_1_fasta'] = append_to_output_lists(
                    plus_1_list_fasta, properties_list['plus_1_fasta'],
                    res_ids_list, strand_or_res
                )
                minus_2_list = reverse_strand_lists(minus_2_list, reverse)
                properties_list['minus_2'] = append_to_output_lists(
                    minus_2_list, properties_list['minus_2'], res_ids_list,
                    strand_or_res
                )
                plus_2_list = reverse_strand_lists(plus_2_list, reverse)
                properties_list['plus_2'] = append_to_output_lists(
                    plus_2_list, properties_list['plus_2'], res_ids_list,
                    strand_or_res
                )
                minus_2_list_fasta = reverse_strand_lists(minus_2_list_fasta, reverse)
                properties_list['minus_2_fasta'] = append_to_output_lists(
                    minus_2_list_fasta, properties_list['minus_2_fasta'],
                    res_ids_list, strand_or_res
                )
                plus_2_list_fasta = reverse_strand_lists(plus_2_list_fasta, reverse)
                properties_list['plus_2_fasta'] = append_to_output_lists(
                    plus_2_list_fasta, properties_list['plus_2_fasta'],
                    res_ids_list, strand_or_res
                )

            # Ensures that domain is transmembrane
            if self.discard_non_tm and sorted(set(strands_list)) == sorted(set(unprocessed_strands)):
                unprocessed_list.append(domain_id)

        # Removes unnecessary columns from dataset
        if self.code[0:4] in ['2.40']:
            unwanted_columns = ['sheet_number', 'edge_or_central',
                                'strand_abs_pos', 'strand_percentage_pos',
                                'strand_z_coords', 'sandwich_z_coords',
                                'core_surf', 'buried_surface_area']
        elif self.code[0:4] in ['2.60']:
            unwanted_columns = ['tilt_angle', 'total_strand_number',
                                'z_coords', 'tm_pos', 'tm_ext']
        for col in unwanted_columns:
            del properties_list[col]

        # Converts properties_list dict into dataframe
        for index, domain_id in enumerate(properties_list['domain_ids']):
            if domain_id in unprocessed_list:
                for property in list(properties_list.keys()):
                    properties_list[property][index] = np.nan
        beta_strands_df = pd.DataFrame(properties_list)
        beta_strands_df = beta_strands_df.dropna()
        beta_strands_df = beta_strands_df.reset_index(drop=True)
        for domain_id in unprocessed_list:
            for key in list(copy.deepcopy(sec_struct_dfs_dict).keys()):
                if domain_id == key:
                    del sec_struct_dfs_dict[key]
            for key in list(copy.deepcopy(domain_sheets_dict).keys()):
                if key.startswith(domain_id):
                    del domain_sheets_dict[key]

        # Calculates strand twist angles for HB and NHB pairs
        domain_dfs = []
        for domain_id in list(sec_struct_dfs_dict.keys()):
            domain_df = copy.deepcopy(
                beta_strands_df[beta_strands_df['domain_ids']==domain_id]
            ).reset_index(drop=True)
            if domain_df.empty:
                raise Exception('No data recorded for {}'.format(domain_id))
            domain_df = find_strand_twist(domain_id, domain_df, strand_or_res)
            domain_dfs.append(domain_df)
        beta_strands_df = pd.concat(domain_dfs, axis=0).reset_index(drop=True)

        # Calculates barrel shear number
        if (self.code[0:4] in ['2.40']) and (strand_or_res == 'res'):
            domain_dfs = []
            for domain_id in list(sec_struct_dfs_dict.keys()):

                domain_df = copy.deepcopy(
                    beta_strands_df[beta_strands_df['domain_ids']==domain_id]
                ).reset_index(drop=True)
                if domain_df.empty:
                    raise Exception('No data recorded for {}'.format(domain_id))
                domain_df = find_sheet_shear(
                    domain_id, domain_df, domain_sheets_dict
                )
                domain_dfs.append(domain_df)
            beta_strands_df = pd.concat(domain_dfs, axis=0).reset_index(drop=True)

        # Saves output dataframe
        beta_strands_df.to_pickle('Beta_{}_dataframe.pkl'.format(strand_or_res))
        beta_strands_df.to_csv('Beta_{}_dataframe.csv'.format(strand_or_res),
                               index=False)

        with open('Unprocessed_domains.txt', 'a') as unprocessed_file:
            unprocessed_file.write('\n\nError with strand processing:\n')
            for domain_id in set(unprocessed_list):
                unprocessed_file.write('{}\n'.format(domain_id))
