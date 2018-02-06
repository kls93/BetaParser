
import os
import sys
import shutil
import pickle
import pandas as pd
from collections import OrderedDict

prompt = '> '


class run_stages():

    def __init__(self, run_parameters):
        self.run_parameters = run_parameters
        self.code = self.run_parameters['id']
        self.pdb_database = self.run_parameters['pdbdatabase']
        self.dssp_database = self.run_parameters['dsspdatabase']
        self.resn = float(self.run_parameters['resolution'])
        self.rfac = float(self.run_parameters['rfactor'])

    def run_stage_1_cath(self, orig_dir):
        # Runs stage 1 of the DataGen pipeline, extracting sequences of the
        # structural domain of interest from the CATH database
        if __name__ == 'subroutines.run_stages':
            from subroutines.CATH import (
                gen_domain_desc_list, domain_desc_filter
                )
            from subroutines.CDHIT import filter_beta_structure
        else:
            from datagen.subroutines.CATH import (
                gen_domain_desc_list, domain_desc_filter
                )
            from datagen.subroutines.CDHIT import filter_beta_structure

        # Generates a list of the domain descriptions provided in
        # CATH_domain_description_v_4_2_0.txt. Then filters the domain
        # descriptions list for beta-structures (the type dependent upon the
        # earlier user input), picking out PDB accession codes and sequences
        # (whose values are stored in the 'domain_df' dataframe).
        domains_desc = gen_domain_desc_list(orig_dir)
        domain_df = domain_desc_filter(self.code, domains_desc)

        # Filters the domain_df for X-ray structures with resolution < 1.6
        # Angstroms (to allow distinction of hydrogen bonds) and R_factor
        # (working value) < 0.20. Writes a file listing all PDB ids that meet
        # these criteria suitable for uploading to the CD-HIT web server.
        beta_structure = filter_beta_structure(self.run_parameters)
        filtered_domain_df = beta_structure.resn_rfac_filter(domain_df)
        beta_structure.gen_cdhit_list(filtered_domain_df)

    def run_stage_1_scope(self, orig_dir):
        return

    def run_stage_2(self, cdhit_entries, cdhit_output):
        if __name__ == 'subroutines.run_stages':
            from subroutines.extract_coordinates import extract_beta_structure_coords
            from subroutines.DSSP import (
                filter_dssp_database, beta_structure_dssp_classification
                )
            from subroutines.generate_network import calculate_beta_network
        else:
            from datagen.subroutines.extract_coordinates import extract_beta_structure_coords
            from datagen.subroutines.DSSP import (
                filter_dssp_database, beta_structure_dssp_classification
                )
            from datagen.subroutines.generate_network import calculate_beta_network

        # Loads the dataframe generated in previous steps
        filtered_domain_df = pd.read_pickle(cdhit_entries)

        # Obtains xyz coordinates for the sequences output from CD-HIT
        if os.path.isdir('Entire_domains'):
            shutil.rmtree('Entire_domains')
        os.mkdir('Entire_domains')

        beta_structure = extract_beta_structure_coords(self.run_parameters)
        cdhit_domain_df = beta_structure.gen_cdhit_dict(
            cdhit_output, filtered_domain_df
            )
        cdhit_domain_df, all_atoms_dfs_dict = beta_structure.get_xyz_coords(
            cdhit_domain_df
            )

        # Copies required DSSP files from database (on hard drive) to local
        # machine
        if os.path.isdir('DSSP_files'):
            shutil.rmtree('DSSP_files')
        os.mkdir('DSSP_files')

        filtered_files = filter_dssp_database(self.run_parameters)
        dssp_domain_df = filtered_files.copy_files_from_dssp_database(cdhit_domain_df)

        # Extracts beta-strands (as classified by DSSP) from the beta-structure
        # domains
        beta_structure = beta_structure_dssp_classification(self.run_parameters)
        dssp_residues_dict = beta_structure.extract_dssp_file_lines(dssp_domain_df)

        shutil.rmtree('DSSP_files')

        if os.path.isdir('Beta_strands'):
            shutil.rmtree('Beta_strands')
        os.mkdir('Beta_strands')

        all_atoms_dfs_dict, sec_struct_dfs_dict = beta_structure.get_dssp_sec_struct_df(
            dssp_residues_dict, all_atoms_dfs_dict
            )
        beta_structure.write_dssp_sec_struct_pdb(sec_struct_dfs_dict)
        del dssp_residues_dict  # To save memory and reduce the number of
        # variables considered

        # Combines the beta-strands into sheets and translates the identified
        # beta-strand interactions into a network
        beta_structure = calculate_beta_network(self.run_parameters)
        domain_sheets_dict, sec_struct_dfs_dict = beta_structure.generate_network(
            sec_struct_dfs_dict
            )

        # Pickles variables required for running stage 3 (which currently has
        # to run within ISAMBARD, hence the division of stages 2-4)
        with open('Input_ISAMBARD_variables.pkl', 'wb') as pickle_file:
            pickle.dump((sec_struct_dfs_dict, domain_sheets_dict), pickle_file)

    def run_stage_3(self):
        # To be run within ISAMBARD
        if __name__ == 'subroutines.run_stages':
            from subroutines.naccess import calculate_solvent_accessible_surface_area
        else:
            from datagen.subroutines.naccess import calculate_solvent_accessible_surface_area

        with open('Input_ISAMBARD_variables.pkl', 'rb') as pickle_file:
            sec_struct_dfs_dict, domain_sheets_dict = pickle.load(pickle_file)

        beta_structure = calculate_solvent_accessible_surface_area(self.run_parameters)
        sec_struct_dfs_dict, domain_sheets_dict = beta_structure.run_naccess(
            sec_struct_dfs_dict, domain_sheets_dict
            )
        sec_struct_dfs_dict = beta_structure.identify_int_ext(
            sec_struct_dfs_dict, domain_sheets_dict
            )

        with open('Output_ISAMBARD_variables.pkl', 'wb') as pickle_file:
            pickle.dump((sec_struct_dfs_dict, domain_sheets_dict), pickle_file)

    def run_stage_4(self, orig_dir, opm_database):
        if __name__ == 'subroutines.run_stages':
            from subroutines.OPM import (
                extract_strand_tilt_and_TM_regions, calculate_barrel_geometry
                )
            from subroutines.output_dataframe import gen_output
        else:
            from datagen.subroutines.OPM import (
                extract_strand_tilt_and_TM_regions, calculate_barrel_geometry
                )
            from datagen.subroutines.output_dataframe import gen_output

        with open('Output_ISAMBARD_variables.pkl', 'rb') as pickle_file:
            sec_struct_dfs_dict, domain_sheets_dict = pickle.load(pickle_file)

        output = gen_output(self.run_parameters)
        if self.code[0:4] in ['2.40']:
            beta_structure = extract_strand_tilt_and_TM_regions(self.run_parameters)
            opm_df = beta_structure.parse_opm(orig_dir)
            tilt_angles = beta_structure.find_strand_tilt(
                sec_struct_dfs_dict, opm_df
                )
            barrel_structure = calculate_barrel_geometry(self.run_parameters)
            strand_numbers = barrel_structure.find_barrel_strand_number(
                sec_struct_dfs_dict
                )
            shear_numbers = barrel_structure.find_barrel_shear_number(
                sec_struct_dfs_dict, domain_sheets_dict
                )
            del opm_df  # To save memory and reduce the number of variables
            # considered
        else:
            tilt_angles = OrderedDict()
            strand_numbers = OrderedDict()
            shear_numbers = OrderedDict()
            sec_struct_dfs_dict = output.identify_edge_central(
                domain_sheets_dict, sec_struct_dfs_dict
                )

        output.write_beta_strand_dataframe(
            sec_struct_dfs_dict, opm_database, tilt_angles, strand_numbers,
            shear_numbers
            )
