
import copy
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from collections import OrderedDict
if __name__ == 'subroutines.generate_network':
    from subroutines.run_stages import run_stages
else:
    from datagen.subroutines.run_stages import run_stages


class network_calcs():

    def create_strands_dict(sheet_df):
        # Creates a dictionary of each of the strand labels plus the residue
        # ranges they span in an individual beta-sheet
        strands_dict = {}
        strands = [int(strand) for strand in set(sheet_df['STRAND_NUM'].tolist())]
        for strand in strands:
            strand_df = sheet_df[sheet_df['STRAND_NUM'] == strand]

            strand_res_num = [int(res) for res in strand_df['DSSP_NUM'].tolist()]
            strand_min = min(strand_res_num)
            strand_max = max(strand_res_num)

            strands_dict[strand] = [num for num in range(strand_min, strand_max+1)]

        return strands_dict

    def identify_strand_interactions(sheet_df, strands_dict):
        # Identifies pairs of hydrogen-bonded strands plus their relative
        # orientation (antiparallel or parallel) from the strand numbers of the
        # bridge-paired residues listed in the DSSP file info
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
            for key in list(strands_dict.keys()):
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

        return strand_pairs

    def create_network(domain_id, sheet, strand_pairs, edge_labels,
                       domain_sheets, domain_sheets_dict):
        # Creates network of the interacting strands in the sheet
        G = nx.Graph()
        for pair in strand_pairs:
            G.add_edge(pair[0], pair[1], attr=strand_pairs[pair])

        # Discards sheets with fewer than 3 beta-strands
        if nx.number_of_nodes(G) > 2:
            domain_sheets['{}_sheet_{}'.format(domain_id, sheet)] = G
            domain_sheets_dict['{}_sheet_{}'.format(domain_id, sheet)] = G
            # Adds strand pairs in retained sheets to edge_labels
            # dictionary
            for key in strand_pairs:
                if key not in edge_labels:
                    edge_labels[key] = strand_pairs[key]

        return edge_labels, domain_sheets, domain_sheets_dict

    def write_network_pdb(domain_id, dssp_df, domain_sheets):
        # Writes pdb file of individual sheets
        for sheet in list(domain_sheets.keys()):
            sheet = sheet.replace('{}_sheet_'.format(domain_id), '')
            print('Writing PDB file of beta-strands in {} sheet {}'.format(
                  domain_id, sheet))

            with open(
                'Beta_strands/{}_sheet_{}.pdb'.format(domain_id, sheet), 'w'
            ) as sheet_pdb_file:
                sheet_df = dssp_df[dssp_df['SHEET_NUM'] == sheet]
                sheet_df = sheet_df.reset_index(drop=True)

                strands = sorted([int(strand) for strand in set(sheet_df['STRAND_NUM'].tolist())])
                for strand in strands:
                    strand_df = sheet_df[sheet_df['STRAND_NUM'] == strand]
                    chain_res_num = strand_df['RES_ID'].tolist()

                    for row in range(dssp_df.shape[0]):
                        if dssp_df['RES_ID'][row] in chain_res_num:
                            sheet_pdb_file.write('{}\n'.format(dssp_df['PDB_FILE_LINES'][row]))
                    sheet_pdb_file.write('TER'.ljust(80)+'\n')

    def draw_network(domain_id, domain_sheets, edge_labels):
        # Draws plot of the network of interacting strands (all retained sheets
        # in a domain are plotted on the same graph)
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


class calculate_beta_network(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    # SPLIT DOWN INTO SMALLER FUNCTIONS TO AVOID INDENTATION BUGS ETC.!
    def generate_network(self, sec_struct_dfs_dict):
        # Separates the beta-strands identified by DSSP into sheets, and works
        # out the strand interactions and thus loop connections both within and
        # between sheets
        if __name__ == 'subroutines.generate_network':
            from subroutines.generate_network import network_calcs
        else:
            from datagen.subroutines.generate_network import network_calcs

        sec_struct_dfs_dict_copy = copy.copy(sec_struct_dfs_dict)
        domain_sheets_dict = OrderedDict()

        for domain_id, dssp_df in sec_struct_dfs_dict_copy.items():
            edge_labels = {}
            domain_sheets = {}

            sheets = [sheet for sheet in set(dssp_df['SHEET_NUM'].tolist()) if sheet != '']
            for sheet in sheets:
                sheet_df = dssp_df[dssp_df['SHEET_NUM'] == sheet]
                sheet_df = sheet_df.reset_index(drop=True)

                strands_dict = network_calcs.create_strands_dict(sheet_df)
                strand_pairs = network_calcs.identify_strand_interactions(
                    sheet_df, strands_dict
                )
                edge_labels, domain_sheets, domain_sheets_dict = network_calcs.create_network(
                    domain_id, sheet, strand_pairs, edge_labels, domain_sheets,
                    domain_sheets_dict
                )

            if domain_sheets:
                network_calcs.write_network_pdb(domain_id, dssp_df, domain_sheets)
                network_calcs.draw_network(domain_id, domain_sheets, edge_labels)
            else:
                sec_struct_dfs_dict[domain_id] = None

        sec_struct_dfs_dict = {key: value for key, value in
                               sec_struct_dfs_dict.items() if value is not None}

        return domain_sheets_dict, sec_struct_dfs_dict
