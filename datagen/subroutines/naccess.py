
import itertools
import pickle
import pandas as pd
from collections import OrderedDict
if __name__ == 'subroutines.naccess':
    from subroutines.run_stages import run_stages
else:
    from datagen.subroutines.run_stages import run_stages

class calculate_solvent_accessible_surface_area(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    def run_naccess(self, sec_struct_dfs_dict, domain_sheets_dict):
        import isambard.external_programs.naccess as naccess

        unprocessed_list = []

        for domain_id in list(sec_struct_dfs_dict.keys()):
            sec_struct_df = sec_struct_dfs_dict[domain_id]
            sheets = {key: value for key, value in domain_sheets_dict.items()
                      if domain_id in key}
            combinations = list(itertools.combinations(list(sheets.keys()), 2))
            solv_acsblty_list = ['']*sec_struct_df.shape[0]
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

            # Calculates solvent accessibility of individual residues (in the
            # context of the parent CATH / SCOPe domain)
            res_solv_acsblty = {}
            with open('Entire_domains/{}.pdb'.format(domain_id), 'r') as pdb_file:
                sheet_string = ''.join(pdb_file.readlines())
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

            # Checks that ratio of solvent accessible surface area of the two
            # sheets in complex is not greater than 90% of the total surface
            # area of the individual sheets
            if min(list(solv_acsblty_dict.keys())) > 1:
                unprocessed_list.append(domain_id)
                sec_struct_dfs_dict[domain_id] = None
                for sheet_id in sheets:
                    domain_sheets_dict[sheet_id] = None
            else:
                sandwich = solv_acsblty_dict[min(list(solv_acsblty_dict.keys()))]

                for sheet_id in sheets:
                    if sheet_id not in sandwich:
                        domain_sheets_dict[sheet_id] = None

                sheets_retained = [sheet_id.strip('{}_sheet_'.format(domain_id))
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
                    if (chain_res_num in list(res_solv_acsblty.keys())
                        and sec_struct_df['ATMNAME'][row] == 'CA'
                        ):
                        solv_acsblty_list[row] = res_solv_acsblty[chain_res_num]
                solv_acsblty_df = pd.DataFrame({'SOLV_ACSBLTY': solv_acsblty_list})
                sec_struct_df = pd.concat([sec_struct_df, solv_acsblty_df], axis=1)
                sec_struct_df = sec_struct_df[sec_struct_df['REC'].notnull()]
                sec_struct_df = sec_struct_df.reset_index(drop=True)
                sec_struct_dfs_dict[domain_id] = sec_struct_df

        sec_struct_dfs_dict = {key: value for key, value in
                               sec_struct_dfs_dict.items() if value is not None}
        domain_sheets_dict = {key: value for key, value in
                              domain_sheets_dict.items() if value is not None}

        return sec_struct_dfs_dict, domain_sheets_dict



    def identify_int_ext(self):
        # Uses solvent accessibility to identify the interior and exterior face
        # of each strand
        return sec_struct_dfs_dict
