
import pickle
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from collections import OrderedDict
if __name__ == 'subroutines.generate_network':
    from subroutines.run_stages import run_stages
else:
    from datagen.subroutines.run_stages import run_stages


class manipulate_beta_structure(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    def identify_strand_interactions(self, sec_struct_dfs_dict):
        # Separates the beta-strands identified by DSSP into sheets, and works
        # out the strand interactions and thus loop connections both within and
        # between sheets
        domain_sheets_dict = OrderedDict()

        for domain_id, dssp_df in sec_struct_dfs_dict.items():
            edge_labels = {}
            domain_sheets = {}

            sheets = [sheet for sheet in set(dssp_df['SHEET_NUM'].tolist()) if sheet != '']
            for sheet in sheets:
                sheet_df = dssp_df[dssp_df['SHEET_NUM']==sheet]
                sheet_df = sheet_df.reset_index(drop=True)

                # Creates a dictionary of each of the strand labels plus the
                # residue ranges they span in an individual beta-sheet
                strands_dict = {}
                strands = [int(strand) for strand in set(sheet_df['STRAND_NUM'].tolist())]
                for strand in strands:
                    strand_df = sheet_df[sheet_df['STRAND_NUM']==strand]

                    strand_res_num = [int(res) for res in strand_df['DSSP_NUM'].tolist()]
                    strand_min = min(strand_res_num)
                    strand_max = max(strand_res_num)

                    strands_dict[strand] = [num for num in range(strand_min, strand_max+1)]

                # Identifies pairs of hydrogen-bonded strands plus their
                # relative orientation (antiparallel or parallel) from the
                # strand numbers of the bridge-paired residues listed in the
                # DSSP file info
                strand_pairs = {}
                for index_1, pair in enumerate(sheet_df['H-BONDS'].tolist()):
                    res_num_1 = int(sheet_df['DSSP_NUM'].tolist()[index_1])
                    res_num_2 = int(pair[0])
                    orientation_2 = sheet_df['ORIENTATION'].tolist()[index_1][0]
                    res_num_3 = int(pair[1])
                    orientation_3 = sheet_df['ORIENTATION'].tolist()[index_1][1]
                    strand_1, strand_2, strand_3 = None, None, None

                    """
                    # Commented out to avoid warnings for beta-bridge pairs
                    # (DSSP code = 'B' rather than 'E') - might want to include
                    # these in my analysis in the future (at which point the
                    # warning given by this section woul be useful)
                    #
                    # Checks that all bridge pair residues are classified
                    # (according to DSSP) as being located within beta-strands
                    for res_num in (res_num_1, res_num_2, res_num_3):
                        if (not res_num in [num for num_range in
                            list(strands_dict.values()) for num in num
                            _range]+[0]
                            ):
                            print('Warning: {} unable to be located in any of '
                                  'the beta-strands listed in the DSSP '
                                  'classification.'
                                  '\nCheck the H-BONDS (bridge pairs) '
                                  'formed by residue {} (DSSP ID)'.format(
                                  res_num, res_num))
                    """

                    # Identifies corresponding strand numbers of bridge pairs
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
                        strand_pair = (min(strand_pair), max(strand_pair))
                        if strand_pair not in strand_pairs:
                            strand_pairs[strand_pair] = orientation_2
                    if strand_1 is not None and strand_3 is not None:
                        strand_pair = [strand_1, strand_3]
                        strand_pair = (min(strand_pair), max(strand_pair))
                        if strand_pair not in strand_pairs:
                            strand_pairs[strand_pair] = orientation_3

                # Writes pdb file of individual sheets
                print('Writing PDB file of beta-strands in {} sheet {}'.format(
                      domain_id, sheet))
                with open('Beta_strands/{}_sheet_{}.pdb'.format(domain_id, sheet), 'w') as sheet_pdb_file:
                    for strand in strands:
                        strand_df = sheet_df[sheet_df['STRAND_NUM']==strand]
                        chain = strand_df['CHAIN'].tolist()
                        res_num = strand_df['RESNUM'].tolist()
                        inscode = strand_df['INSCODE'].tolist()
                        chain_res_num = [chain[index]+str(res_num[index])
                                         +inscode[index] for index, value in
                                         enumerate(chain)]
                        for row in range(dssp_df.shape[0]):
                            if ((dssp_df['CHAIN'][row] + str(dssp_df['RESNUM'][row])
                                + dssp_df['INSCODE'][row]) in chain_res_num
                                ):
                                sheet_pdb_file.write('{}\n'.format(dssp_df['PDB_FILE_LINES'][row]))
                        sheet_pdb_file.write('TER'.ljust(80)+'\n')

                # Creates network of the interacting strands in the sheet
                G = nx.Graph()
                for pair in strand_pairs:
                    G.add_edge(pair[0], pair[1], attr=strand_pairs[pair])

                # Discards sheets with fewer than 3 beta-strands
                if nx.number_of_nodes(G) > 2:
                    domain_sheets['{}_sheet_{}'.format(domain_id, sheet)] = G
                    domain_sheets_dict['{}_sheet_{}'.format(domain_id, sheet)] = G
                    # Adds strand pairs in reatined sheets to edge_labels
                    # dictionary
                    for key in strand_pairs:
                        if key not in edge_labels:
                            edge_labels[key] = strand_pairs[key]

            # Draws plot of the network of interacting strands (all retained
            # sheets in a domain are plotted on the same graph)
            print('Plotting network of interacting beta-strands in {}'.format(domain_id))
            sheet_1 = list(domain_sheets.keys())[0]
            sheets_2_to_n = list(domain_sheets.keys())[1:]
            G = domain_sheets[sheet_1]
            for sheet_n in sheets_2_to_n:
                H = domain_sheets[sheet_n]
                G = nx.compose(G, H)

            pos = nx.circular_layout(G)
            plt.clf()
            nx.draw_networkx(G, pos=pos, with_labels=True, node_shape='v')
            nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels)
            plt.savefig(
                'Beta_strands/{}_network.png'.format(domain_id)
                )

        # Pickles dictionaries of: 1) dataframes of atoms in residues of the
        # selected secondary structure type (currently restricted to
        # beta-strands only), 2) individual sheet networks and 3) individual
        # domain networks. These can then be ported into my dockerised version
        # of ISAMBARD, in which I can run naccess to calculate residue and
        # sheet solvent accessibilities to work out which sheets interact with
        # one another
        with open('Domain_networks_dict.pkl', 'wb') as pickle_file:
            pickle.dump((sec_struct_dfs_dict, domain_sheets_dict), pickle_file)

        return domain_sheets_dict
