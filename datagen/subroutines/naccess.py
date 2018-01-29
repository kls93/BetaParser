
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

    def calculate_sandwich_solv_acsblty(sec_struct_df, sheets):
        combinations = list(itertools.combinations(list(sheets.keys()), 2))
        solv_acsblty_dict = OrderedDict()
        for sheet_pair in combinations:
            # Calculates solvent accessibility of both sheets in pair
            print('Calculating solvent accessible surface area for {}'.format(domain_id))
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

    def calculate_barrel_solv_acsblty(sec_struct_df, sheets):
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

    def add_naccess_info_to_df(domain_id, sec_struct_df, sheets, code,
                               solv_acsblty_dict, res_solv_acsblty,
                               sec_struct_dfs_dict, domain_sheets_dict,
                               unprocessed_list):
        # Checks that ratio of solvent accessible surface area of the two
        # sheets in complex is not greater than the total surface
        # area of the individual sheets
        res_id_list = ['']*sec_struct_df.shape[0]
        solv_acsblty_list = ['']*sec_struct_df.shape[0]

        if solv_acsblty_dict:
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

                for sheet_id in sheets:
                    if sheet_id not in sandwich:
                        domain_sheets_dict[sheet_id] = None

                sheets_retained = [sheet_id.replace('{}_sheet_'.format(domain_id), '')
                                   for sheet_id in sandwich]
                sub_sec_struct_df = sec_struct_df[sec_struct_df['SHEET_NUM'].isin(sheets_retained)]
                chain = sub_sec_struct_df['CHAIN'].tolist()
                res_num = sub_sec_struct_df['RESNUM'].tolist()
                inscode = sub_sec_struct_df['INSCODE'].tolist()
                chain_res_num_list = [chain[index]+str(res_num[index])+inscode[index]
                                      for index, value in enumerate(chain)]
                for row in range(sec_struct_df.shape[0]):
                    chain_res_num = (sec_struct_df['CHAIN'][row]
                                     +str(sec_struct_df['RESNUM'][row])
                                     +sec_struct_df['INSCODE'][row])
                    if chain_res_num not in chain_res_num_list:
                        sec_struct_df.loc[row, 'REC'] = None
                    elif (chain_res_num in chain_res_num_list and
                        sec_struct_df['ATMNAME'][row] == 'CA'
                        ):
                        res_id_list[row] = chain_res_num

                    if (chain_res_num in list(res_solv_acsblty.keys())
                        and sec_struct_df['ATMNAME'][row] == 'CA'
                        ):
                        solv_acsblty_list[row] = res_solv_acsblty[chain_res_num]

        solv_acsblty_df = pd.DataFrame({'SOLV_ACSBLTY': solv_acsblty_list})
        res_id_df = pd.DataFrame({'RES_ID': res_id_list})
        sec_struct_df = pd.concat([sec_struct_df, solv_acsblty_df,
                                  res_id_df], axis=1)
        sec_struct_df = sec_struct_df[sec_struct_df['REC'].notnull()]
        sec_struct_df = sec_struct_df.reset_index(drop=True)
        sec_struct_dfs_dict[domain_id] = sec_struct_df

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

        beta_structure = naccess_solv_acsblty_calcs()

        unprocessed_list = []
        for domain_id in list(sec_struct_dfs_dict.keys()):
            sec_struct_df = sec_struct_dfs_dict[domain_id]
            sheets = {key: value for key, value in domain_sheets_dict.items()
                      if domain_id in key}

            if self.code[0:4] in ['2.40']:
                solv_acsblty_dict = naccess_solv_acsblty_calcs.calculate_barrel_solv_acsblty(
                    sec_struct_df, sheets
                    )
                res_solv_acsblty = naccess_solv_acsblty_calcs.calculate_residue_solv_acsblty(domain_id)
            elif self.code[0:4] in ['2.60']:
                solv_acsblty_dict = naccess_solv_acsblty_calcs.calculate_sandwich_solv_acsblty(
                    sec_struct_df, sheets
                    )
                res_solv_acsblty = naccess_solv_acsblty_calcs.calculate_residue_solv_acsblty(domain_id)
            else:
                solv_acsblty_dict = {}
                res_solv_acsblty = {}

            (sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list
            ) = naccess_solv_acsblty_calcs.add_naccess_info_to_df(
                domain_id, sec_struct_df, sheets, self.code, solv_acsblty_dict,
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
        import isambard.external_programs.naccess as naccess

        for domain_id in list(sec_struct_dfs_dict.keys()):
            sec_struct_df = sec_struct_dfs_dict[domain_id]
            sheets = [key for key in list(domain_sheets_dict.keys()) if
                      domain_id in key]
            int_ext_list = ['']*sec_struct_df.shape[0]
            int_ext_dict = OrderedDict()
            if len(sheets) > 2:
                sys.exit('Error - more than 2 sheets retained following '
                         'solvent accessibility calculation')
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
            face_1_res = [value for value in list(int_ext_dict.values())[::2]]
            face_2_res = [value for value in list(int_ext_dict.values())[1::2]]
            face_1_solv_acsblty = sum(face_1_res)
            face_2_solv_acsblty = sum(face_2_res)
            if face_1_solv_acsblty > face_2_solv_acsblty:
                ext_chain_res_num = list(int_ext_dict.keys())[::2]
                int_chain_res_num = list(int_ext_dict.keys())[1::2]
            elif face_1_solv_acsblty < face_2_solv_acsblty:
                ext_chain_res_num = list(int_ext_dict.keys())[1::2]
                int_chain_res_num = list(int_ext_dict.keys())[::2]
            else:
                sys.exit('Solvent accessibilities of two faces of {} are equal - '
                         'something has gone wrong with the analysis'.format(sheet_id))
            for chain_res_num in ext_chain_res_num:
                int_ext_dict[chain_res_num] = 'exterior'
            for chain_res_num in int_chain_res_num:
                int_ext_dict[chain_res_num] = 'interior'

            for row in range(sec_struct_df.shape[0]):
                chain_res_num = (sec_struct_df['CHAIN'][row]
                                 +str(sec_struct_df['RESNUM'][row])
                                 +sec_struct_df['INSCODE'][row])
                if (chain_res_num in list(int_ext_dict.keys())
                    and sec_struct_df['ATMNAME'][row] == 'CA'
                    ):
                    int_ext_list[row] = int_ext_dict[chain_res_num]
            int_ext_df = pd.DataFrame({'INT_EXT': int_ext_list})
            sec_struct_df = pd.concat([sec_struct_df, int_ext_df], axis=1)
            sec_struct_df = sec_struct_df.reset_index(drop=True)
            sec_struct_dfs_dict[domain_id] = sec_struct_df

        return sec_struct_dfs_dict
