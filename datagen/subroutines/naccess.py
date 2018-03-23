
import sys
import math
import itertools
import isambard
import pandas as pd
import numpy as np
import networkx as nx
from collections import OrderedDict
if __name__ == 'subroutines.naccess':
    from subroutines.run_stages import run_stages
else:
    from datagen.subroutines.run_stages import run_stages


class naccess_solv_acsblty_calcs():

    def calculate_barrel_solv_acsblty(sheets):
        # Calculates the solvent accessibility of beta-sheets in domain X to
        # identify which is the largest (and hence is assumed to form the
        # beta-barrel)
        solv_acsblty_dict = OrderedDict()

        for sheet_id in sheets:
            sheet_sub_strings = []
            with open('Beta_strands/{}.pdb'.format(sheet_id), 'r') as pdb_file:
                for line in pdb_file.readlines():
                    line_start = line[0:16]
                    line_end = line[17:]
                    new_line = line_start + ' ' + line_end  # Removes alternate
                    # conformer labels to prevent problems with running naccess
                    # (for some reason naccess encounters an error if alternate
                    # conformers are removed without also removing the
                    # conformer id of the retained conformer)
                    sheet_sub_strings.append(new_line)
            sheet_string = ''.join(sheet_sub_strings)

            # Runs NACCESS to calculate solvent accessible surface area of
            # every residue in domain
            print('Calculating solvent accessible surface area for {}'.format(sheet_id))
            naccess_out = isambard.external_programs.naccess.run_naccess(
                sheet_string, 'rsa', path=False, include_hetatms=True
            )
            naccess_out = naccess_out.split('\n')
            naccess_out = [line for line in naccess_out if line != '']
            # Determines the solvent accessibility of the sheet as the sum of
            # the absolute solvent accessibilities of all residues
            sheet_solv_acsblty = naccess_out[-1]
            sheet_solv_acsblty = float(sheet_solv_acsblty.split()[1])
            solv_acsblty_dict[sheet_solv_acsblty] = sheet_id

        # Selects the beta-sheet with the highest solvent accessibility
        max_solv_acsblty = max(list(solv_acsblty_dict.keys()))
        beta_barrel = solv_acsblty_dict[max_solv_acsblty]
        solv_acsblty_dict = OrderedDict()
        solv_acsblty_dict[max_solv_acsblty] = [beta_barrel]  # Must be a list!

        return solv_acsblty_dict

    def calculate_sandwich_solv_acsblty(sheets):
        # Calculates the solvent accessibility of beta-sheets to identify the
        # two which share the largest buried surface area and hence form the
        # beta-sandwich in domain X

        # Generates all possible combinations of beta-sheets in the domain
        combinations = list(itertools.combinations(list(sheets.keys()), 2))
        solv_acsblty_dict = OrderedDict()
        for sheet_pair in combinations:
            # Calculates solvent accessibility of complex of both sheets in pair
            print('Calculating solvent accessible surface area for {}'.format(sheet_pair))
            sheet_sub_strings = []
            for sheet_id in sheet_pair:
                with open('Beta_strands/{}.pdb'.format(sheet_id), 'r') as pdb_file:
                    for line in pdb_file.readlines():
                        line_start = line[0:16]
                        line_end = line[17:]
                        new_line = line_start + ' ' + line_end  # Removes
                        # alternate conformer labels to prevent problems with
                        # running naccess (for some reason naccess encounters
                        # an error if alternate conformers are removed without
                        # also removing the conformer id of the retained
                        # conformer)
                        sheet_sub_strings.append(new_line)
            sheet_string = ''.join(sheet_sub_strings)

            # Runs NACCESS to calculate solvent accessible surface area of
            # every residue in sheet pair complex
            naccess_out = isambard.external_programs.naccess.run_naccess(
                sheet_string, 'rsa', path=False, include_hetatms=True
            )
            naccess_out = naccess_out.split('\n')
            naccess_out = [line for line in naccess_out if line != '']
            # Determines the solvent accessibility of the complex as the sum of
            # the absolute solvent accessibilities of all residues
            complex_solv_acsblty = naccess_out[-1]
            complex_solv_acsblty = float(complex_solv_acsblty.split()[1])

            # Calculates solvent accessibility of individual sheets
            sum_sheet_solv_acsblty = 0
            for sheet_id in sheet_pair:
                sheet_sub_strings = []
                print('Calculating solvent accessible surface area for {}'.format(sheet_id))
                with open('Beta_strands/{}.pdb'.format(sheet_id), 'r') as pdb_file:
                    for line in pdb_file.readlines():
                        line_start = line[0:16]
                        line_end = line[17:]
                        new_line = line_start + ' ' + line_end  # Removes
                        # alternate conformer labels to prevent problems with
                        # running naccess (for some reason naccess encounters
                        # an error if alternate conformers are removed without
                        # also removing the conformer id of the retained
                        # conformer)
                        sheet_sub_strings.append(new_line)
                sheet_string = ''.join(sheet_sub_strings)

                # Runs NACCESS to calculate solvent accessible surface area of
                # every residue in individual sheet
                naccess_out = isambard.external_programs.naccess.run_naccess(
                    sheet_string, 'rsa', path=False, include_hetatms=True
                )
                naccess_out = naccess_out.split('\n')
                naccess_out = [line for line in naccess_out if line != '']
                # Determines the solvent accessibility of the individual sheet as
                # the sum of the absolute solvent accessibilities of all residues
                sheet_solv_acsblty = naccess_out[-1]
                sheet_solv_acsblty = float(sheet_solv_acsblty.split()[1])
                sum_sheet_solv_acsblty += sheet_solv_acsblty

            # Calculates absolute surface area buried in the sheet complex from
            # the difference between the summed solvent accessibility of the
            # individual sheets and the solvent accessibility of the sheet pair
            # complex
            buried_surface_area = sum_sheet_solv_acsblty - complex_solv_acsblty
            solv_acsblty_dict[buried_surface_area] = [sheet_pair[0], sheet_pair[1]]

        return solv_acsblty_dict

    def calculate_residue_solv_acsblty(domain_id):
        # Calculates solvent accessibility of individual residues (in the
        # context of the parent CATH / SCOPe domain)
        res_solv_acsblty = OrderedDict()
        sheet_sub_strings = []

        print('Calculating solvent accessible surface areas of individual '
              'residues in {}'.format(domain_id))

        with open('Entire_domains/{}.pdb'.format(domain_id), 'r') as pdb_file:
            for line in pdb_file.readlines():
                line_start = line[0:16]
                line_end = line[17:]
                new_line = line_start + ' ' + line_end  # Removes alternate
                # conformer labels to prevent problems with running naccess
                # (for some reason naccess encounters an error if alternate
                # conformers are removed without also removing the
                # conformer id of the retained conformer)
                sheet_sub_strings.append(new_line)
        sheet_string = ''.join(sheet_sub_strings)

        # Runs NACCESS to calculate solvent accessible surface area of every
        # residue in domain
        naccess_out = isambard.external_programs.naccess.run_naccess(
            sheet_string, 'rsa', path=False, include_hetatms=True
        )
        naccess_out = naccess_out.split('\n')
        for line in naccess_out:
            if line[0:3] in ['RES', 'HEM']:
                chain = line[8:9].strip()
                res_num = line[9:13].strip()
                ins_code = line[13:14].strip()
                chain_res_num = chain + res_num + ins_code
                res_solv_acsblty[chain_res_num] = float(line[14:22])

        return res_solv_acsblty

    def add_naccess_info_to_df(domain_id, dssp_df, sheets, code,
                               solv_acsblty_dict, res_solv_acsblty,
                               sec_struct_dfs_dict, domain_sheets_dict,
                               unprocessed_list):
        # Uses the previously calculated solvent accessibility values to
        # determine which sheet/s form/s the beta-barrel/sandwich, and then
        # prunes the dataframe (dssp_df) to remove all other sheets

        # Initialises list of residue solvent accessibility values for
        # incorporation into the domain dataframe (dssp_df)
        solv_acsblty_list = ['']*dssp_df.shape[0]

        # Discounts a beta-sandwich domain from further analysis if the solvent
        # accessibility calculations show that none of the retained beta-sheets
        # are in contact with one another
        if (code[0:4] in ['2.60']
                    and max(list(solv_acsblty_dict.keys())) == 0.0
                ):
            unprocessed_list.append(domain_id)
            sec_struct_dfs_dict[domain_id] = None
            for sheet_id in sheets:
                domain_sheets_dict[sheet_id] = None

            return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

        else:
            sandwich = solv_acsblty_dict[max(list(solv_acsblty_dict.keys()))]

            # Removes sheets that do not form part of the sandwich / barrel
            # from further analysis
            for sheet_id in sheets:
                if not sheet_id in sandwich:
                    domain_sheets_dict[sheet_id] = None

            sheets_retained = [sheet_id.replace('{}_sheet_'.format(domain_id), '')
                               for sheet_id in sandwich]
            sub_dssp_df = dssp_df[dssp_df['SHEET_NUM'].isin(sheets_retained)]
            chain_res_num = sub_dssp_df['RES_ID'].tolist()

            for row in range(dssp_df.shape[0]):
                res_id = dssp_df['RES_ID'][row]

                if not res_id in chain_res_num:
                    dssp_df.loc[row, 'REC'] = None

                if (res_id in list(res_solv_acsblty.keys())
                        and dssp_df['ATMNAME'][row] == 'CA'
                        ):
                    solv_acsblty_list[row] = res_solv_acsblty[res_id]

        solv_acsblty_df = pd.DataFrame({'SOLV_ACSBLTY': solv_acsblty_list})
        dssp_df = pd.concat([dssp_df, solv_acsblty_df], axis=1)
        dssp_df = dssp_df[dssp_df['REC'].notnull()]
        dssp_df = dssp_df.reset_index(drop=True)
        sec_struct_dfs_dict[domain_id] = dssp_df

        return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

    def write_domain_pdb(domain_id, dssp_df, sheets):
        # Writes a PDB file of the residues that form the beta-barrel /
        # -sandwich as determined from sheet surface solvent accessibility
        # calculations performed in the previous steps. Individual beta-strands
        # are separated by 'TER' cards
        print('Writing PDB file of beta-strands in {}'.format(domain_id))

        strand_number_set = [strand for strand in set(dssp_df['STRAND_NUM'].tolist())
                             if strand != '']

        with open('Beta_strands/{}.pdb'.format(domain_id), 'w') as new_pdb_file:
            for strand in strand_number_set:
                dssp_df_strand = dssp_df[dssp_df['STRAND_NUM'] == strand]
                chain_res_num = dssp_df_strand['RES_ID'].tolist()

                for row in range(dssp_df.shape[0]):
                    if dssp_df['RES_ID'][row] in chain_res_num:
                        new_pdb_file.write('{}\n'.format(dssp_df['PDB_FILE_LINES'][row]))
                new_pdb_file.write('TER'.ljust(80)+'\n')

    def calc_core_residues_sandwich(domain_id, dssp_df, sheets,
                                    sec_struct_dfs_dict, domain_sheets_dict,
                                    unprocessed_list):
        # Calculates whether each residue faces towards the interior or the
        # exterior of the sandwich from the ratio of its solvent accessibility
        # within the sandwich to its solvent accessibility within its
        # individual parent beta-sheet

        # Checks that correct number of sheets has been retained for analysis
        if len(sheets) != 2:
            print('ERROR: incorrect number of sheets retained in {} following '
                  'solvent accessibility calculation'.format(domain_id))
            sec_struct_dfs_dict[domain_id] = None
            for sheet in sheets:
                domain_sheets_dict[sheet] = None
            unprocessed_list.append(domain_id)

            return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

        print('Determining interior and exterior facing residues in {}'.format(domain_id))

        # Initialises records of core / surface residues
        core_ext_list = ['']*dssp_df.shape[0]
        buried_surface_area_list = ['']*dssp_df.shape[0]
        core_ext_combined = OrderedDict()
        core_ext_indv = OrderedDict()
        buried_surface_area_dict = OrderedDict()

        # Calculates solvent accessibility of both sheets in pair
        print('Calculating solvent accessible surface area for residues in '
              '{}'.format(domain_id))
        sheet_sub_strings = []
        for sheet_id in sheets:
            with open('Beta_strands/{}.pdb'.format(sheet_id), 'r') as pdb_file:
                for line in pdb_file.readlines():
                    line_start = line[0:16]
                    line_end = line[17:]
                    new_line = line_start + ' ' + line_end  # Removes alternate
                    # conformer labels to prevent problems with running naccess
                    # (for some reason naccess encounters an error if alternate
                    # conformers are removed without also removing the
                    # conformer id of the retained conformer)
                    sheet_sub_strings.append(new_line)
        sheet_string = ''.join(sheet_sub_strings)

        # Runs NACCESS to calculate solvent accessible surface area of every
        # residue in beta-sheet pair
        naccess_out = isambard.external_programs.naccess.run_naccess(
            sheet_string, 'rsa', path=False, include_hetatms=True
        )
        naccess_out = naccess_out.split('\n')
        naccess_out = [line for line in naccess_out if line != '']
        for line in naccess_out:
            if line[0:3] in ['RES', 'HEM']:
                chain = line[8:9].strip()
                res_num = line[9:13].strip()
                ins_code = line[13:14].strip()
                chain_res_num = chain + res_num + ins_code
                core_ext_combined[chain_res_num] = float(line[14:22])

        # Calculates solvent accessibility of sheets in pair individually
        for sheet_id in sheets:
            sheet_sub_strings = []
            print('Calculating solvent accessible surface area for residues in '
                  '{}'.format(sheet_id))
            with open('Beta_strands/{}.pdb'.format(sheet_id), 'r') as pdb_file:
                for line in pdb_file.readlines():
                    line_start = line[0:16]
                    line_end = line[17:]
                    new_line = line_start + ' ' + line_end  # Removes alternate
                    # conformer labels to prevent problems with running naccess
                    # (for some reason naccess encounters an error if alternate
                    # conformers are removed without also removing the
                    # conformer id of the retained conformer)
                    sheet_sub_strings.append(new_line)
            sheet_string = ''.join(sheet_sub_strings)

            # Runs NACCESS to calculate solvent accessible surface area of every
            # residue in domain
            naccess_out = isambard.external_programs.naccess.run_naccess(
                sheet_string, 'rsa', path=False, include_hetatms=True
            )
            naccess_out = naccess_out.split('\n')
            naccess_out = [line for line in naccess_out if line != '']
            for line in naccess_out:
                if line[0:3] in ['RES', 'HEM']:
                    chain = line[8:9].strip()
                    res_num = line[9:13].strip()
                    ins_code = line[13:14].strip()
                    chain_res_num = chain + res_num + ins_code
                    core_ext_indv[chain_res_num] = float(line[14:22])

        # Determines solvent accessibility of each residue in the beta-sheet
        # assembly as compared to its individual parent beta-sheet. If the
        # solvent accessibility is reduced by >= 20% in the assembly, the
        # residue is classified as 'core', else it is classified as 'surface'.
        for res in list(core_ext_combined.keys()):
            solv_acsblty_combined = core_ext_combined[res]
            solv_acsblty_indv = core_ext_indv[res]

            if solv_acsblty_indv != 0 and solv_acsblty_combined == 0:
                buried_surface_area = solv_acsblty_indv
                buried_surface_area_dict[res] = buried_surface_area
                if buried_surface_area >= 5:
                    # 5 is an arbitrarily selected value that I have found to
                    # discriminate well between core and surface residues
                    core_ext_combined[res] = 'core'
                elif buried_surface_area < 5:
                    core_ext_combined[res] = 'surface'

            elif solv_acsblty_indv != 0 and solv_acsblty_combined != 0:
                buried_surface_area = round((((solv_acsblty_indv - solv_acsblty_combined)
                                              / solv_acsblty_indv) * 100), 1)
                buried_surface_area_dict[res] = buried_surface_area
                if buried_surface_area >= 20:
                    # 20% is an arbitrarily selected value that I have found to
                    # discriminate well between core and surface residues
                    core_ext_combined[res] = 'core'
                elif buried_surface_area < 20:
                    core_ext_combined[res] = 'surface'

        # Updates dataframe with solvent accessibility information
        for row in range(dssp_df.shape[0]):
            res_id = dssp_df['RES_ID'][row]
            if (res_id in list(core_ext_combined.keys())
                    and dssp_df['ATMNAME'][row] == 'CA'
                    ):
                core_ext_list[row] = core_ext_combined[res_id]
            if (res_id in list(buried_surface_area_dict.keys())
                    and dssp_df['ATMNAME'][row] == 'CA'
                    ):
                buried_surface_area_list[row] = buried_surface_area_dict[res_id]
        core_ext_df = pd.DataFrame({'CORE_OR_SURFACE': core_ext_list,
                                    'BURIED_SURFACE_AREA(%)': buried_surface_area_list})
        dssp_df = pd.concat([dssp_df, core_ext_df], axis=1)
        sec_struct_dfs_dict[domain_id] = dssp_df

        return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list


class calculate_solvent_accessible_surface_area(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    def calc_sasa(self, sec_struct_dfs_dict, domain_sheets_dict):
        # Pipeline function for running naccess to calculate solvent accessible
        # surface area

        unprocessed_list_1 = []
        unprocessed_list_2 = []

        for domain_id in list(sec_struct_dfs_dict.keys()):
            dssp_df = sec_struct_dfs_dict[domain_id]
            sheets = {key: value for key, value in domain_sheets_dict.items()
                      if domain_id in key}

            if (
                (not sheets)
                or (self.code[0:4] in ['2.60'] and len(list(sheets.keys())) < 2)
            ):
                print('ERROR: No / only 1 beta-sheet/s retained for {}'.format(domain_id))
                unprocessed_list_1.append(domain_id)
                continue

            if self.code[0:4] in ['2.40']:
                solv_acsblty_dict = naccess_solv_acsblty_calcs.calculate_barrel_solv_acsblty(
                    sheets
                )
                res_solv_acsblty = naccess_solv_acsblty_calcs.calculate_residue_solv_acsblty(
                    domain_id
                )
            elif self.code[0:4] in ['2.60']:
                solv_acsblty_dict = naccess_solv_acsblty_calcs.calculate_sandwich_solv_acsblty(
                    sheets
                )
                res_solv_acsblty = naccess_solv_acsblty_calcs.calculate_residue_solv_acsblty(
                    domain_id
                )
            else:
                solv_acsblty_dict = {}
                res_solv_acsblty = {}

            (sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list_2
             ) = naccess_solv_acsblty_calcs.add_naccess_info_to_df(
                domain_id, dssp_df, sheets, self.code, solv_acsblty_dict,
                res_solv_acsblty, sec_struct_dfs_dict, domain_sheets_dict,
                unprocessed_list_2
            )
            naccess_solv_acsblty_calcs.write_domain_pdb(domain_id, dssp_df, sheets)

        sec_struct_dfs_dict = OrderedDict(
            {key: value for key, value in sec_struct_dfs_dict.items() if value is not None}
        )
        domain_sheets_dict = OrderedDict(
            {key: value for key, value in domain_sheets_dict.items() if value is not None}
        )

        with open('Unprocessed_domains.txt', 'a') as unprocessed_file:
            unprocessed_file.write('\n\nError in network generation - no / '
                                   'only 1 sheet/s retained for solvent '
                                   'accessibility calculations:\n')
            for domain_id in set(unprocessed_list_1):
                unprocessed_file.write('{}\n'.format(domain_id))
            unprocessed_file.write('\n\nError in solvent accessibility calculation:\n')
            for domain_id in set(unprocessed_list_2):
                unprocessed_file.write('{}\n'.format(domain_id))

        return sec_struct_dfs_dict, domain_sheets_dict

    def identify_core_surface(self, sec_struct_dfs_dict, domain_sheets_dict):
        # Pipeline function to identify residues that form part of the
        # (hydrophobic) core of a beta-sandwich domain

        unprocessed_list = []

        for domain_id in list(sec_struct_dfs_dict.keys()):
            dssp_df = sec_struct_dfs_dict[domain_id]
            sheets = [key for key in list(domain_sheets_dict.keys()) if
                      domain_id in key]

            (sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list
             ) = naccess_solv_acsblty_calcs.calc_core_residues_sandwich(
                domain_id, dssp_df, sheets, sec_struct_dfs_dict,
                domain_sheets_dict, unprocessed_list
            )

        sec_struct_dfs_dict = OrderedDict(
            {key: value for key, value in sec_struct_dfs_dict.items() if value is not None}
        )
        domain_sheets_dict = OrderedDict(
            {key: value for key, value in domain_sheets_dict.items() if value is not None}
        )

        with open('Unprocessed_domains.txt', 'a') as unprocessed_file:
            unprocessed_file.write('\n\nError in finding core and external '
                                   'residues:\n')
            for domain_id in set(unprocessed_list):
                unprocessed_file.write('{}\n'.format(domain_id))

        return sec_struct_dfs_dict, domain_sheets_dict
