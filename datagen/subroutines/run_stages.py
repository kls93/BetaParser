
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
        self.pdb_au_database = self.run_parameters['pdbaudatabase']
        self.pdb_ba_database = self.run_parameters['pdbbadatabase']
        self.dssp_database = self.run_parameters['dsspdatabase']
        self.opm_database = self.run_parameters['opmdatabase']
        self.ring_database = self.run_parameters['ringdatabase']
        self.resn = float(self.run_parameters['resolution'])
        self.rfac = float(self.run_parameters['rfactor'])
        self.suffix = self.run_parameters['suffix']
        self.discard_non_tm = self.run_parameters['discardnontm']

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
        domain_df = domain_desc_filter(self.code, domains_desc, self.discard_non_tm)

        # Filters the domain_df for X-ray structures below user-specified
        # resolution and R_factor (working value) cutoffs. Writes a file
        # listing all PDB ids that meet these criteria suitable for uploading
        # to the CD-HIT web server.
        beta_structure = filter_beta_structure(self.run_parameters)
        filtered_domain_df = beta_structure.resn_rfac_filter(domain_df)
        beta_structure.gen_cdhit_list(filtered_domain_df)

    def run_stage_1_scope(self, orig_dir):
        return

    def run_stage_2(self, cdhit_entries, cdhit_output):
        # To be run on local machine within ISAMBARD (in docker container).
        if __name__ == 'subroutines.run_stages':
            from subroutines.extract_coordinates import extract_beta_structure_coords
            from subroutines.DSSP import beta_structure_dssp_classification
            from subroutines.generate_network import calculate_beta_network
        else:
            from datagen.subroutines.extract_coordinates import extract_beta_structure_coords
            from datagen.subroutines.DSSP import beta_structure_dssp_classification
            from datagen.subroutines.generate_network import calculate_beta_network

        # Loads the dataframe generated in previous steps
        filtered_domain_df = pd.read_pickle(cdhit_entries)

        # Obtains xyz coordinates for the sequences output from CD-HIT
        beta_structure = extract_beta_structure_coords(self.run_parameters)
        cdhit_domain_df = beta_structure.gen_cdhit_dict(
            cdhit_output, filtered_domain_df
        )

        if os.path.isdir('Parent_assemblies'):
            shutil.rmtree('Parent_assemblies')
        os.mkdir('Parent_assemblies')
        cdhit_domain_df = beta_structure.copy_parent_assembly_pdb(
            cdhit_domain_df, self.suffix
        )
        cdhit_domain_df, all_atoms_dfs_dict = beta_structure.get_xyz_coords(
            cdhit_domain_df
        )

        all_atoms_dfs_dict = beta_structure.remove_alternate_conformers(
            all_atoms_dfs_dict
        )

        # Extracts beta-strands (as classified by DSSP) from the beta-structure
        # domains
        beta_structure = beta_structure_dssp_classification(self.run_parameters)
        dssp_residues_dict, dssp_domain_df = beta_structure.extract_dssp_file_lines(
            cdhit_domain_df
        )

        if os.path.isdir('Beta_strands'):
            shutil.rmtree('Beta_strands')
        os.mkdir('Beta_strands')
        (all_atoms_dfs_dict, sec_struct_dfs_dict, dssp_to_pdb_dict
         ) = beta_structure.get_dssp_sec_struct_df(
            dssp_residues_dict, all_atoms_dfs_dict
        )
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
            pickle.dump((sec_struct_dfs_dict, domain_sheets_dict, dssp_to_pdb_dict
                         ), pickle_file)

    def run_stage_3(self, radius):
        # To be run within ISAMBARD. **NOTE** solvent accessibility
        # calculations must be run first.
        if __name__ == 'subroutines.run_stages':
            from subroutines.naccess import calculate_solvent_accessible_surface_area
            from subroutines.find_surfaces import find_interior_exterior_surfaces
            from subroutines.dihedral_angles import dihedral_angles
            from subroutines.neighbouring_residues import nearest_neighbours
        else:
            from datagen.subroutines.naccess import calculate_solvent_accessible_surface_area
            from datagen.subroutines.find_surfaces import find_interior_exterior_surfaces
            from datagen.subroutines.dihedral_angles import dihedral_angles
            from datagen.subroutines.neighbouring_residues import nearest_neighbours

        # Loads pickled variables generated from running stage 2.
        with open('Input_ISAMBARD_variables.pkl', 'rb') as pickle_file:
            (sec_struct_dfs_dict, domain_sheets_dict, dssp_to_pdb_dict
             ) = pickle.load(pickle_file)

        # Runs NACCESS to calculate the solvent accessible surface area of
        # each individual residue
        beta_structure = calculate_solvent_accessible_surface_area(self.run_parameters)
        sec_struct_dfs_dict, domain_sheets_dict = beta_structure.calc_sasa(
            sec_struct_dfs_dict, domain_sheets_dict
        )
        # For beta-sandwiches, calculates the surface area buried by each
        # individual residue within the sandwich core - this value is then used
        # to classify the residue as 'core' or 'surface'
        if self.code[0:4] in ['2.60']:
            sec_struct_dfs_dict, domain_sheets_dict = beta_structure.identify_core_surface(
                sec_struct_dfs_dict, domain_sheets_dict
            )

        # For each individual residue, calculates whether it faces in towards
        # or outwards from the centre of the beta_barrel/sheet
        if self.code[0:4] in ['2.40']:
            beta_structure = find_interior_exterior_surfaces(self.run_parameters)
            sec_struct_dfs_dict, domain_sheets_dict = beta_structure.identify_int_ext(
                sec_struct_dfs_dict, domain_sheets_dict
            )

        # Calculates the backbone and side chain dihedral angles of each
        # individual residue
        beta_structure = dihedral_angles(self.run_parameters)
        sec_struct_dfs_dict = beta_structure.calc_dihedral_angles(sec_struct_dfs_dict)

        # For each residue, determines the RES_IDs of all residues with at
        # least one atom within a user-specified radius (measured in Angstroms)
        # of its C_alpha atom
        beta_structure = nearest_neighbours(self.run_parameters, radius)
        sec_struct_dfs_dict = beta_structure.calculate_nearest_neighbours(
            sec_struct_dfs_dict
        )

        # Pickles variables required for running stage 4
        with open('Output_ISAMBARD_variables.pkl', 'wb') as pickle_file:
            pickle.dump((sec_struct_dfs_dict, domain_sheets_dict,
                         dssp_to_pdb_dict, radius), pickle_file)

    def run_stage_4(self, orig_dir):
        # To be run on local machine (in docker container).
        if __name__ == 'subroutines.run_stages':
            from subroutines.RING import calculate_residue_interaction_network
            from subroutines.OPM import (
                extract_barrel_info_from_OPM, calculate_barrel_geometry
            )
            from subroutines.output_dataframe import gen_output
        else:
            from datagen.subroutines.RING import calculate_residue_interaction_network
            from datagen.subroutines.OPM import (
                extract_barrel_info_from_OPM, calculate_barrel_geometry
            )
            from datagen.subroutines.output_dataframe import gen_output

        with open('Output_ISAMBARD_variables.pkl', 'rb') as pickle_file:
            (sec_struct_dfs_dict, domain_sheets_dict, dssp_to_pdb_dict, radius
             ) = pickle.load(pickle_file)

        # Calculates residue interaction network of the amino acids in the
        # selected beta-strands of the CATH domain in the context of its parent
        # biological assembly
        res_int_network = calculate_residue_interaction_network(self.run_parameters)
        res_int_network.run_RING(sec_struct_dfs_dict, domain_sheets_dict)
        sec_struct_dfs_dict, domain_sheets_dict = res_int_network.parse_RING_output(
            sec_struct_dfs_dict, domain_sheets_dict
        )

        if self.code[0:4] in ['2.60']:
            sec_struct_dfs_dict, domain_sheets_dict = res_int_network.identify_int_ext_sandwich(
                sec_struct_dfs_dict, domain_sheets_dict
            )

        # For beta-barrel domains, calculates the strand number, plus, if the
        # barrel is in the OPM database, extracts its average strand tilt
        # number from the database
        # TODO: Fix shear number calculation
        tilt_angles = OrderedDict()
        strand_numbers = OrderedDict()
        shear_numbers = OrderedDict()
        if self.code[0:4] in ['2.40']:
            beta_structure = extract_barrel_info_from_OPM(self.run_parameters)

            opm_df = beta_structure.parse_opm(orig_dir)
            tilt_angles = beta_structure.find_strand_tilt(
                sec_struct_dfs_dict, opm_df
            )
            barrel_structure = calculate_barrel_geometry(self.run_parameters)
            strand_numbers = barrel_structure.find_barrel_strand_number(
                sec_struct_dfs_dict
            )
            del opm_df  # To save memory and reduce the number of variables
            # considered

        # Classifies each beta-strand as either 'edge' or 'central' based upon
        # the number of strands with which it forms backbone hydrogen-bonding
        # interactions
        output = gen_output(self.run_parameters, radius)
        if self.code[0:4] in ['2.60']:
            sec_struct_dfs_dict = output.identify_edge_central(
                domain_sheets_dict, sec_struct_dfs_dict
            )

        # Writes a csv file of the beta-barrel/sandwich dataset organised such
        # that each row in the file represents an individual beta-strand.
        output.write_beta_strand_dataframe(
            'strand', sec_struct_dfs_dict, dssp_to_pdb_dict, tilt_angles,
            strand_numbers, shear_numbers
        )
        # Writes a csv file of the beta-barrel/sandwich dataset organised such
        # that each row in the file represents an individual residue.
        output.write_beta_strand_dataframe(
            'res', sec_struct_dfs_dict, dssp_to_pdb_dict, tilt_angles,
            strand_numbers, shear_numbers
        )
