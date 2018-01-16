
import pickle
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from collections import OrderedDict

class manipulate_beta_structure():

    def __init__(self, run, resn, rfac):
        self.run = run
        self.resn = resn
        self.rfac = rfac

    def identify_strand_interactions(self, dssp_dfs_dict):
        # Separates the beta-strands identified by DSSP into sheets, and works
        # out the strand interactions and thus loop connections both within and
        # between sheets
        unprocessed_list = []
        domain_sheets_dict = OrderedDict()
        domain_networks_dict = OrderedDict()

        for domain_id, dssp_df in dssp_dfs_dict.items():
            edge_labels = {}

            sheets = [sheet for sheet in set(dssp_df['SHEET_NUM'].tolist()) if sheet != '']
            for sheet in sheets:
                sheet_df = dssp_df[dssp_df['SHEET_NUM']==sheet]
                sheet_df = sheet_df.reset_index(drop=True)

                strands_dict = {}
                strands = [int(strand) for strand in set(sheet_df['STRAND_NUM'].tolist())]
                for strand in strands:
                    strand_df = sheet_df[sheet_df['STRAND_NUM']==strand]

                    strand_res_num = [int(res) for res in strand_df['DSSP_NUM'].tolist()]
                    strand_min = min(strand_res_num)
                    strand_max = max(strand_res_num)

                    strands_dict[strand] = [num for num in range(strand_min, strand_max+1)]

                strand_pairs = {}
                for index_1, pair in enumerate(sheet_df['H-BONDS'].tolist()):
                    res_num_1 = int(sheet_df['DSSP_NUM'].tolist()[index_1])
                    res_num_2 = int(pair[0])
                    orientation_2 = sheet_df['ORIENTATION'].tolist()[index_1][0]
                    res_num_3 = int(pair[1])
                    orientation_3 = sheet_df['ORIENTATION'].tolist()[index_1][1]
                    strand_1 = None
                    strand_2 = None
                    strand_3 = None
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

                    if strand_1 is not None and strand_2 is not None:
                        strand_pair = [strand_1, strand_2]
                        strand_pairs[(min(strand_pair), max(strand_pair))] = orientation_2
                    if strand_1 is not None and strand_3 is not None:
                        strand_pair = [strand_1, strand_3]
                        strand_pairs[(min(strand_pair), max(strand_pair))] = orientation_3

                # Writes pdb file of individual sheets
                with open('DSSP_filtered_DSEQS/{}_sheet_{}.pdb'.format(domain_id, sheet), 'w') as sheet_pdb_file:
                    for strand in strands:
                        strand_df = sheet_df[sheet_df['STRAND_NUM']==strand]
                        chain = strand_df['CHAIN'].tolist()
                        res_num = strand_df['RESNUM'].tolist()
                        inscode = strand_df['INSCODE'].tolist()
                        chain_res_num = [chain[index]+str(res_num[index])+inscode[index]
                                         for index, value in enumerate(chain)]
                        for row in range(dssp_df.shape[0]):
                            if ((dssp_df['CHAIN'][row] + str(dssp_df['RESNUM'][row])
                                + dssp_df['INSCODE'][row]) in chain_res_num
                                ):
                                sheet_pdb_file.write('{}\n'.format(dssp_df['FILE_LINES'][row]))
                        sheet_pdb_file.write('TER'.ljust(80)+'\n')

                # Creates network of the interacting strands in the sheet
                interacting_strands = {}
                for key in strand_pairs:
                    if key not in interacting_strands:
                        interacting_strands[key] = strand_pairs[key]

                G = nx.Graph()
                for pair in interacting_strands:
                    G.add_edge(pair[0], pair[1], attr=interacting_strands[pair])

                # Discards sheets with fewer than 3 beta-strands
                if nx.number_of_nodes(G) > 2:
                    domain_sheets_dict['{}_sheet_{}'.format(domain_id, sheet)] = G
                    for key in strand_pairs:
                        if key not in edge_labels:
                            edge_labels[key] = strand_pairs[key]

            # Draws plot of the network of interacting strands
            print('Plotting network of interacting beta-strands in {}'.format(domain_id))
            interacting_sheets = {key: value for key, value in domain_sheets_dict.items()
                                  if domain_id in key}
            if len(interacting_sheets) < 2:
                unprocessed_list.append(domain_id)
            else:
                sheet_1 = list(interacting_sheets.keys())[0]
                sheets_2_to_n = list(interacting_sheets.keys())[1:]
                G = interacting_sheets[sheet_1]
                for sheet_n in sheets_2_to_n:
                    H = interacting_sheets[sheet_n]
                    G = nx.compose(G, H)
                domain_networks_dict[domain_id] = G

                pos = nx.circular_layout(G)
                plt.clf()
                nx.draw_networkx(G, pos=pos, with_labels=True, node_shape='v')
                nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels)
                plt.savefig(
                    'DSSP_filtered_DSEQS/{}_network.png'.format(domain_id)
                    )

        with open('Unprocessed_CATH_{}.txt'.format(self.run), 'a') as unprocessed_file:
            unprocessed_file.write('\n\nFewer than 2 beta-sheets of > 2 '
                                   'strands identified in the domain:\n')
            for domain_id in set(unprocessed_list):
                unprocessed_file.write('{}\n'.format(domain_id))

        dssp_dfs_dict = {key: value for key, value in dssp_dfs_dict.items() if key not
                         in unprocessed_list}

        with open(
            'CATH_{}_resn_{}_rfac_{}_domain_networks_dict.pkl'.format(
                self.run, self.resn, self.rfac
                ), 'wb'
            ) as pickle_file:
            pickle.dump((dssp_dfs_dict, domain_networks_dict,
                         domain_sheets_dict), pickle_file)
