
import sys
import itertools
import pickle
import pandas as pd
import isambard.external_programs.naccess as naccess
from collections import OrderedDict
if __name__ == 'subroutines.naccess':
    from subroutines.run_stages import run_stages
else:
    from datagen.subroutines.run_stages import run_stages

class naccess_solv_acsblty_calcs():

    def calculate_sandwich_solv_acsblty(sheets):
        combinations = list(itertools.combinations(list(sheets.keys()), 2))
        solv_acsblty_dict = OrderedDict()
        for sheet_pair in combinations:
            # Calculates solvent accessibility of both sheets in pair
            print('Calculating solvent accessible surface area for {}'.format(sheet_pair))
            sheet_string = ''
            for sheet_id in sheet_pair:
                with open('Beta_strands/{}.pdb'.format(sheet_id), 'r') as pdb_file:
                    sheet_string += ''.join(pdb_file.readlines())

            naccess_out = naccess.run_naccess(sheet_string, 'rsa', path=False,
                                              include_hetatms=True)
            naccess_out = naccess_out.split('\n')
            naccess_out = [line for line in naccess_out if line != '']
            domain_solv_acsblty = naccess_out[-1]
            domain_solv_acsblty = float(domain_solv_acsblty.split()[1])

            # Calculates solvent accessibility of individual sheets
            sum_sheet_solv_acsblty = 0
            for sheet_id in sheet_pair:
                print('Calculating solvent accessible surface area for {}'.format(sheet_id))
                with open('Beta_strands/{}.pdb'.format(sheet_id), 'r') as pdb_file:
                    sheet_string = ''.join(pdb_file.readlines())
                naccess_out = naccess.run_naccess(sheet_string, 'rsa',
                                                  path=False,
                                                  include_hetatms=True)
                naccess_out = naccess_out.split('\n')
                naccess_out = [line for line in naccess_out if line != '']
                sheet_solv_acsblty = naccess_out[-1]
                sheet_solv_acsblty = float(sheet_solv_acsblty.split()[1])
                sum_sheet_solv_acsblty += sheet_solv_acsblty

            solv_acsblty_ratio = domain_solv_acsblty / sum_sheet_solv_acsblty
            solv_acsblty_dict[solv_acsblty_ratio] = [sheet_pair[0], sheet_pair[1]]

        return solv_acsblty_dict

    def calculate_barrel_solv_acsblty(sheets):
        solv_acsblty_dict = OrderedDict()

        for sheet_id in sheets:
            with open('Beta_strands/{}.pdb'.format(sheet_id), 'r') as pdb_file:
                sheet_string = ''.join(pdb_file.readlines())

            print('Calculating solvent accessible surface area for {}'.format(sheet_id))
            naccess_out = naccess.run_naccess(sheet_string, 'rsa', path=False,
                                              include_hetatms=True)
            naccess_out = naccess_out.split('\n')
            naccess_out = [line for line in naccess_out if line != '']
            domain_solv_acsblty = naccess_out[-1]
            domain_solv_acsblty = float(domain_solv_acsblty.split()[1])
            solv_acsblty_dict[domain_solv_acsblty] = sheet_id

        max_solv_acsblty = max(list(solv_acsblty_dict.keys()))
        beta_barrel = solv_acsblty_dict[max_solv_acsblty]
        solv_acsblty_dict = OrderedDict()
        solv_acsblty_dict[max_solv_acsblty] = beta_barrel

        return solv_acsblty_dict

    def calculate_residue_solv_acsblty(domain_id):
        # Calculates solvent accessibility of individual residues (in the
        # context of the parent CATH / SCOPe domain)
        res_solv_acsblty = OrderedDict()
        with open('Entire_domains/{}.pdb'.format(domain_id), 'r') as pdb_file:
            sheet_string = ''.join(pdb_file.readlines())

        print('Calculating solvent accessible surface areas of individual '
              'residues in {}'.format(domain_id))
        naccess_out = naccess.run_naccess(sheet_string, 'rsa', path=False,
                                          include_hetatms=True)
        naccess_out = naccess_out.split('\n')
        for line in naccess_out:
            if line[0:3] == 'RES':
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
        solv_acsblty_list = ['']*dssp_df.shape[0]

        # Discounts 'sandwich' from further analysis if the solvent
        # accessiblity calculations show that none of the retained beta-sheets
        # are in contact with one another
        if (code[0:4] not in ['2.40']
            and min(list(solv_acsblty_dict.keys())) >= 1
            ):
            unprocessed_list.append(domain_id)
            sec_struct_dfs_dict[domain_id] = None
            for sheet_id in sheets:
                domain_sheets_dict[sheet_id] = None

            return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

        else:
            sandwich = solv_acsblty_dict[min(list(solv_acsblty_dict.keys()))]

            # Removes sheets outside of the sandwich / barrel from further
            # analysis
            for sheet_id in sheets:
                if sheet_id not in sandwich:
                    domain_sheets_dict[sheet_id] = None

            sheets_retained = [sheet_id.replace('{}_sheet_'.format(domain_id), '')
                               for sheet_id in sandwich]
            sub_dssp_df = dssp_df[dssp_df['SHEET_NUM'].isin(sheets_retained)]

            chain_res_num = sub_dssp_df['RES_ID'].tolist()

            for row in range(dssp_df.shape[0]):
                res_id = dssp_df['RES_ID'][row]

                if res_id not in chain_res_num:
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

class calculate_solvent_accessible_surface_area(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    def run_naccess(self, sec_struct_dfs_dict, domain_sheets_dict):
        # Pipeline function for running naccess to calculate solvent accessible
        # surface area
        if __name__ == 'subroutines.naccess':
            from subroutines.naccess import naccess_solv_acsblty_calcs
        else:
            from datagen.subroutines.naccess import naccess_solv_acsblty_calcs

        unprocessed_list = []
        for domain_id in list(sec_struct_dfs_dict.keys()):
            dssp_df = sec_struct_dfs_dict[domain_id]
            sheets = {key: value for key, value in domain_sheets_dict.items()
                      if domain_id in key}

            if not sheets:
                print('ERROR: No beta-sheets retained for {}'.format(domain_id))

            if self.code[0:4] in ['2.40']:
                solv_acsblty_dict = naccess_solv_acsblty_calcs.calculate_barrel_solv_acsblty(
                    sheets
                    )
                res_solv_acsblty = naccess_solv_acsblty_calcs.calculate_residue_solv_acsblty(domain_id)
            elif self.code[0:4] in ['2.60']:
                solv_acsblty_dict = naccess_solv_acsblty_calcs.calculate_sandwich_solv_acsblty(
                    sheets
                    )
                res_solv_acsblty = naccess_solv_acsblty_calcs.calculate_residue_solv_acsblty(domain_id)
            else:
                solv_acsblty_dict = {}
                res_solv_acsblty = {}

            (sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list
            ) = naccess_solv_acsblty_calcs.add_naccess_info_to_df(
                domain_id, dssp_df, sheets, self.code, solv_acsblty_dict,
                res_solv_acsblty, sec_struct_dfs_dict, domain_sheets_dict,
                unprocessed_list
                )

        sec_struct_dfs_dict = {key: value for key, value in
                               sec_struct_dfs_dict.items() if value is not None}
        domain_sheets_dict = {key: value for key, value in
                              domain_sheets_dict.items() if value is not None}

        with open('Unprocessed_domains.txt', 'a') as unprocessed_file:
            unprocessed_file.write('\n\nError in solvent accessibility calculation:\n')
            for domain_id in set(unprocessed_list):
                unprocessed_file.write('{}\n'.format(domain_id))

        return sec_struct_dfs_dict, domain_sheets_dict

    def identify_int_ext(self, sec_struct_dfs_dict, domain_sheets_dict):
        # Calculates the sum of the solvent accessibilities of every other
        # residue, and from this determines which face towards the interior
        # of the sandwich and which towards the exterior
        # import isambard.external_programs.naccess as naccess
        unprocessed_list = []

        for domain_id in list(sec_struct_dfs_dict.keys()):
            dssp_df = sec_struct_dfs_dict[domain_id]
            sheets = [key for key in list(domain_sheets_dict.keys()) if
                      domain_id in key]
            int_ext_list = ['']*dssp_df.shape[0]
            int_ext_dict = OrderedDict()

            if self.code[0:4] in ['2.40'] and len(sheets) > 1:
                print('ERROR: more than 1 sheet retained following solvent '
                      'accessibility calculation')
                sec_struct_dfs_dict[domain_id] = None
                unprocessed_list.append(domain_id)
                continue
            elif self.code[0:4] in ['2.60'] and len(sheets) > 2:
                print('ERROR: more than 2 sheets retained following solvent '
                      'accessibility calculation')
                sec_struct_dfs_dict[domain_id] = None
                unprocessed_list.append(domain_id)
                continue

            # Calculates solvent accessibility of both sheets in pair
            print('Calculating solvent accessible surface area for {}'.format(domain_id))
            sheet_string = ''
            for sheet_id in sheets:
                with open('Beta_strands/{}.pdb'.format(sheet_id), 'r') as pdb_file:
                    sheet_string += ''.join(pdb_file.readlines())
            naccess_out = naccess.run_naccess(sheet_string, 'rsa', path=False,
                                              include_hetatms=True)
            naccess_out = naccess_out.split('\n')
            naccess_out = [line for line in naccess_out if line != '']
            for line in naccess_out:
                if line[0:3] == 'RES':
                    chain = line[8:9].strip()
                    res_num = line[9:13].strip()
                    ins_code = line[13:14].strip()
                    chain_res_num = chain + res_num + ins_code
                    int_ext_dict[chain_res_num] = float(line[14:22])

            # Selects every other residue in each beta-sheet
            for sheet_id in sheets:
                sheet_id = sheet_id.replace('{}_sheet_'.format(domain_id), '')
                sheet_df = dssp_df[dssp_df['SHEET_NUM']==sheet_id]
                sheet_res = sheet_df['RES_ID'].tolist()
                res_acsblty = [value for key, value in int_ext_dict.items() if
                               key in sheet_res]

                face_1_res = res_acsblty[::2]
                face_2_res = res_acsblty[1::2]
                face_1_solv_acsblty = sum(face_1_res)
                face_2_solv_acsblty = sum(face_2_res)

                if face_1_solv_acsblty > face_2_solv_acsblty:
                    ext_chain_res_num = face_1_res
                    int_chain_res_num = face_2_res
                elif face_1_solv_acsblty < face_2_solv_acsblty:
                    ext_chain_res_num = face_2_res
                    int_chain_res_num = face_1_res
                else:
                    print('ERROR: solvent accessibilities of two faces of {} '
                          'are equal - something has gone wrong with the '
                          'analysis'.format(sheet_id))
                    sec_struct_dfs_dict[domain_id] = None
                    unprocessed_list.append(domain_id)
                    break

                for chain_res_num in ext_chain_res_num:
                    int_ext_dict[chain_res_num] = 'exterior'
                for chain_res_num in int_chain_res_num:
                    int_ext_dict[chain_res_num] = 'interior'

            for row in range(dssp_df.shape[0]):
                chain_res_num = dssp_df['RES_ID'][row]
                if (chain_res_num in list(int_ext_dict.keys())
                    and dssp_df['ATMNAME'][row] == 'CA'
                    ):
                    int_ext_list[row] = int_ext_dict[chain_res_num]
            int_ext_df = pd.DataFrame({'INT_EXT': int_ext_list})
            dssp_df = pd.concat([dssp_df, int_ext_df], axis=1)
            sec_struct_dfs_dict[domain_id] = dssp_df

        with open('Unprocessed_domains.txt', 'a') as unprocessed_file:
            unprocessed_file.write('\n\nError in determination of interior and '
                                   'exterior facing residues:\n')
            for domain_id in set(unprocessed_list):
                unprocessed_file.write('{}\n'.format(domain_id))

        sec_struct_dfs_dict = {key: value for key, value in
                               sec_struct_dfs_dict.items() if value is not None}

        return sec_struct_dfs_dict
