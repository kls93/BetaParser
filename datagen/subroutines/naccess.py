
import itertools
import pickle
import pandas as pd
from collections import OrderedDict

class calculate_solvent_accessible_surface_area():

    def __init__(self, run, resn, rfac, dssp_dfs_dict, domain_networks_dict,
                 domain_sheets_dict):
        self.run = run
        self.resn = resn
        self.rfac = rfac
        self.dssp_dfs_dict = dssp_dfs_dict
        self.sheets_dict = domain_sheets_dict

    def run_naccess(self):
        import isambard_dev.external_programs.naccess as naccess

        unprocessed_list = []

        for domain_id in list(self.dssp_dfs_dict.keys()):
            sheets = [key for key in list(self.sheets_dict.keys())
                      if domain_id in key]
            combinations = list(itertools.combinations(sheets, 2))
            solv_acsblty_dict = OrderedDict()
            for sheet_pair in combinations:
                # Calculates solvent accessibility of both sheets in pair
                print('Calculating solvent accessible surface area for {}'.format(domain_id))
                sheet_string = ''
                for sheet_id in sheet_pair:
                    with open('DSSP_filtered_DSEQS/{}.pdb'.format(sheet_id), 'r') as pdb_file:
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
                    print('Calculating solvent accessible surface area for {} sheet '
                          '{}'.format(domain_id, sheet_id))
                    naccess_out = naccess.run_naccess('DSSP_filtered_DSEQS/{}.pdb'.format(sheet_id),
                                                      'rsa', include_hetatms=True)
                    naccess_out = naccess_out.split('\n')
                    naccess_out = [line for line in naccess_out if line != '']
                    sheet_solv_acsblty = naccess_out[-1]
                    sheet_solv_acsblty = float(sheet_solv_acsblty.split()[1])
                    sum_sheet_solv_acsblty += sheet_solv_acsblty

                solv_acsblty_ratio = domain_solv_acsblty / sum_sheet_solv_acsblty
                solv_acsblty_dict[solv_acsblty_ratio] = [sheet_pair[0], sheet_pair[1]]

            # Checks that ratio of solvent accessible surface area of the two
            # sheets in complex is not greater than 90% of the total surface
            # area of the individual sheets
            if min(list(solv_acsblty_dict.keys())) > 0.9:
                unprocessed_list.append(domain_id)
                self.dssp_dfs_dict[domain_id] = None
                for sheet_id in sheets:
                    self.sheets_dict[sheet_id] = None
            else:
                sandwich = solv_acsblty_dict[min(list(solv_acsblty_dict.keys()))]

                for sheet_id in sheets:
                    if sheet_id not in sandwich:
                        self.sheets_dict[sheet_id] = None

                dssp_df = self.dssp_dfs_dict[domain_id]
                sheets_retained = [sheet_id.strip('{}_sheet_'.format(domain_id))
                                   for sheet_id in sandwich]
                sub_dssp_df = dssp_df[dssp_df['SHEET_NUM'].isin(sheets_retained)]
                chain = sub_dssp_df['CHAIN'].tolist()
                res_num = sub_dssp_df['RESNUM'].tolist()
                inscode = sub_dssp_df['INSCODE'].tolist()
                chain_res_num = [chain[index]+str(res_num[index])+inscode[index]
                                 for index, value in enumerate(chain)]
                for row in range(dssp_df.shape[0]):
                    if ((dssp_df['CHAIN'][row]+str(dssp_df['RESNUM'][row])+dssp_df['INSCODE'][row])
                        not in chain_res_num
                        ):
                        dssp_df.loc[row, 'REC'] = None
                dssp_df = dssp_df[dssp_df['REC'].notnull()]
                dssp_df = dssp_df.reset_index(drop=True)
                self.dssp_dfs_dict[domain_id] = dssp_df

        fltrd_dssp_dfs_dict = {key: value for key, value in self.dssp_dfs_dict.items()
                               if value is not None}
        fltrd_sheets_dict = {key: value for key, value in self.sheets_dict.items()
                             if value is not None}

        with open(
            'CATH_{}_resn_{}_rfac_{}_filtered_domain_networks_dict.pkl'.format(
                self.run, self.resn, self.rfac
                ), 'wb'
            ) as pickle_file:
            pickle.dump((fltrd_dssp_dfs_dict, fltrd_sheets_dict), pickle_file)
