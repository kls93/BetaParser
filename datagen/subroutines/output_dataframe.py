
import networkx as nx
import pandas as pd
if __name__ == 'subroutines.output_dataframe':
    from subroutines.run_stages import run_stages
    from subroutines.variables import gen_amino_acids_dict
else:
    from datagen.subroutines.run_stages import run_stages
    from datagen.subroutines.variables import gen_amino_acids_dict

class gen_output(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    def identify_edge_central(self, domain_sheets_dict, sec_struct_dfs_dict):
        # Uses domain_networks_dict to identify whether strands are edge or
        # central
        for domain_id in list(sec_struct_dfs_dict.keys()):
            networks = [network for key, network in domain_sheets_dict.items()
                        if domain_id in key]
            G = networks[0]
            for num in range(1, len(networks)):
                H = networks[1]
                G = nx.compose(G, H)
            sec_struct_df = sec_struct_dfs_dict[domain_id]

            edge_or_central = {}
            strands = set(sec_struct_df['STRAND_NUM'].tolist())
            for strand_num in strands:
                if strand_num in G:
                    num_of_edges = len(G.neighbors(strand_num))
                    if num_of_edges == 1:
                        edge_id = 'edge'
                    elif num_of_edges > 1:
                        edge_id = 'central'
                    edge_or_central[strand_num] = edge_id


            df_strands = sec_struct_df['STRAND_NUM'].tolist()
            df_edge_or_central = ['']*len(df_strands)
            for index, value in enumerate(df_strands):
                if value in edge_or_central:
                    df_edge_or_central[index] = edge_or_central[value]
            df_edge_or_central = pd.DataFrame({'EDGE_OR_CNTRL': df_edge_or_central})
            sec_struct_df = pd.concat([sec_struct_df, df_edge_or_central], axis=1)
            sec_struct_dfs_dict[domain_id] = sec_struct_df

        return sec_struct_dfs_dict

    def write_beta_strand_dataframe(self, sec_struct_dfs_dict):
        # Generates dataframe of beta-strands
        amino_acids_dict = gen_amino_acids_dict()

        domain_strand_ids = []
        edge_or_central = []
        residues = []
        int_ext = []
        bckbn_phi_psi = []
        solv_ascblty = []

        for domain_id, sec_struct_df in sec_struct_dfs_dict.items():
            strands_list = set(sec_struct_df['STRAND_NUM'].tolist())
            strands_list = [strand_num for strand_num in strands_list if
                            strand_num != '']

            for strand_num in strands_list:
                strand_df = sec_struct_df[sec_struct_df['STRAND_NUM']==strand_num]
                stranddf = strand_df.reset_index(drop=True)

                # Gives strand a unique ID
                domain_strand_ids.append('{}_strand_{}'.format(domain_id, strand_num))

                # Determines whether the strand is an edge or central strand
                edge_or_central.append(set(strand_df['EDGE_OR_CNTRL'].tolist()))

                # Generates FASTA sequence of strand
                residues.append([amino_acids_dict[residue] for residue in
                                strand_df['RESNAME'].tolist()])

                # Generates list of interior and exterior facing residues in
                # the strand
                # int_ext.append(strand_df['INT_OR_EXT'].tolist())

                # Generates list of phi and psi angles along the strand
                phi = strand_df['PHI'].tolist()
                psi = strand_df['PSI'].tolist()
                bckbn_geom = [[phi[index], psi[index]] for index, value in
                              enumerate(phi)]
                bckbn_phi_psi.append(bckbn_geom)

                # Generates list of residue solvent accessibility along the
                # strand
                solv_ascblty.append(strand_df['SOLV_ACSBLTY'].tolist())

        beta_strands_df = pd.DataFrame({'STRAND_ID': domain_strand_ids,
                                        'EDGE_OR_CNTRL': edge_or_central,
                                        'RESIDUES': residues,
                                        'BCKBN_GEOM': bckbn_phi_psi,
                                        'SOLV_ACSBLTY': solv_ascblty})
        beta_strands_df.to_pickle('Beta_strands_dataframe.pkl')
        beta_strands_df.to_csv('Beta_strands_dataframe.csv')

        """
        int_indices = [index for index, value in enumerate(int_ext) if
                       value == 'interior']
        int_domain_strand_ids = [value for index, value in enumerate(domain_strand_ids)
                                 if index in int_indices]
        int_edge_or_central = [value for index, value in enumerate(edge_or_central)
                               if index in int_indices]
        int_residues = [value for index, value in enumerate(residues)
                        if index in int_indices]
        interior = [value for index, value in enumerate(int_ext)
                    if index in int_indices]
        int_bckbn_phi_psi = [value for index, value in enumerate(bckbn_phi_psi)
                             if index in int_indices]
        int_solv_ascblty = [value for index, value in enumerate(solv_ascblty)
                            if index in int_indices]
        int_beta_strands_df = pd.DataFrame({'STRAND_ID': int_domain_strand_ids,
                                            'EDGE_OR_CNTRL': int_edge_or_central,
                                            'RESIDUES': int_residues,
                                            'INT_EXT': interior,
                                            'BCKBN_GEOM': int_bckbn_phi_psi,
                                            'SOLV_ACSBLTY': int_solv_ascblty})
        int_beta_strands_df.to_pickle('Interior_beta_strands_dataframe.pkl')
        int_beta_strands_df.to_csv('Interior_beta_strands_dataframe.csv')

        ext_indices = [index for index, value in enumerate(int_ext) if
                       value == 'exterior']
        ext_domain_strand_ids = [value for index, value in enumerate(domain_strand_ids)
                                 if index in ext_indices]
        ext_edge_or_central = [value for index, value in enumerate(edge_or_central)
                               if index in ext_indices]
        ext_residues = [value for index, value in enumerate(residues)
                        if index in ext_indices]
        exterior = [value for index, value in enumerate(int_ext)
                    if index in ext_indices]
        ext_bckbn_phi_psi = [value for index, value in enumerate(bckbn_phi_psi)
                             if index in ext_indices]
        ext_solv_ascblty = [value for index, value in enumerate(solv_ascblty)
                            if index in ext_indices]
        ext_beta_strands_df = pd.DataFrame({'STRAND_ID': ext_domain_strand_ids,
                                            'EDGE_OR_CNTRL': ext_edge_or_central,
                                            'RESIDUES': ext_residues,
                                            'INT_EXT': exterior,
                                            'BCKBN_GEOM': ext_bckbn_phi_psi,
                                            'SOLV_ACSBLTY': ext_solv_ascblty})
        ext_beta_strands_df.to_pickle('Exterior_beta_strands_dataframe.pkl')
        ext_beta_strands_df.to_csv('Exterior_beta_strands_dataframe.csv')
        """
