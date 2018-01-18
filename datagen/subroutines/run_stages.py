
import os
import sys
import shutil
import pandas as pd

prompt = '> '

def gen_run_parameters(args):
    # Determines which stage of analysis the user wants to perform - currently
    # the analysis cannot be run in its entirety owing to the need to run the
    # data through programs (the CDHIT web server and naccess) that I cannot
    # run on my local machine
    print('Specify analysis stage:')
    stage = ''
    while stage not in ['1', '2', '3', '4']:
        stage = input(prompt).lower()
        if not stage in ['1', '2', '3', '4']:
            print('Analysis stage not recognised')
        else:
            break

    # Sets as many run parameters as possible from the input file (if provided)
    run_parameters = {}
    if vars(args)['input_file']:
        try:
            with open('/{}'.format(vars(args)['input_file'].strip('/')), 'r') as input_file:
                for line in input_file:
                    key = line.split(':')[0].replace(' ', '').lower()
                    value = line.split(':')[1].replace('\n', '').strip()
                    if key in ['workingdirectory', 'pdbdatabase', 'dsspdatabase']:
                        value = value.replace('\\', '/')  # For windows file paths
                        value = '/{}/'.format(value.strip('/'))
                    else:
                        value = value.replace(' ', '').lower()
                    run_parameters[key] = value
        except FileNotFoundError:
            sys.exit('Absolute path to input file not recognised')

    # Requires user input if the structural database (CATH or SCOPe) is not
    # specified in the input file / is not recognised
    if 'structuredatabase' in run_parameters:
        if run_parameters['structuredatabase'].upper() not in ['CATH', 'SCOP']:
            print('DataGen can currently only parse the CATH and SCOP databases\n'
                  '- please select one of these databases to continue')
            run_parameters.pop('structuredatabase')
    if not 'structuredatabase' in run_parameters:
        print('CATH or SCOP database?')
        database = ''
        while database not in ['CATH', 'SCOP']:
            database = input(prompt).upper()
            if not database in ['CATH', 'SCOP']:
                print('DataGen can currently only parse the CATH and '
                      'SCOP databases\n'
                      '- please select one of these databases to continue')
            else:
                run_parameters['structuredatabase'] = database
                break
    # Ensures database name is uppercase (required for consistent directory naming)
    run_parameters['structuredatabase'] = run_parameters['structuredatabase'].upper()

    # Requires user input if the (all-beta) structural domain the user wishes
    # to analyse is not specified in the input file / is not recognised
    if 'id' in run_parameters:
        if (run_parameters['structuredatabase'] == 'CATH'
            and not run_parameters['id'].startswith('2')
            ):
            print('DataGen is currently only suitable for generation and '
                  'analysis of\nall-beta structures')
            run_parameters.pop('id')
        elif (run_parameters['structuredatabase'] == 'SCOP'
            and not run_parameters['id'].startswith('b')
            ):
            print('DataGen is currently only suitable for generation and '
                  'analysis of\nall-beta structures')
            run_parameters.pop('id')
    if not 'id' in run_parameters:
        if run_parameters['structuredatabase'] == 'CATH':
            print('Specify CATHCODE:')
            run = ''
            while not run.startswith('2'):
                run = input(prompt).lower()
                if run.startswith('2'):
                    run_parameters['id'] = run
                    break
                else:
                    print('DataGen is currently only suitable for '
                          'generation and analysis of\nall-beta structures')
        elif run_parameters['structuredatabase'] == 'SCOP':
            print('Specify SCOP code:')
            run = ''
            while not run.startswith('b'):
                run = input(prompt).lower()
                if run.startswith('b'):
                    run_parameters['id'] = run
                    break
                else:
                    print('DataGen is currently only suitable for '
                          'generation and analysis of\nall-beta structures')

    # Requires user input if the absolute file path of the working directory is
    # not specified in the input file / is not recognised
    if 'workingdirectory' in run_parameters:
        if not os.path.isdir(run_parameters['workingdirectory']):
            print('Specified working directory not recognised')
            run_parameters.pop('workingdirectory')
    if not 'workingdirectory' in run_parameters:
        print('Specify absolute file path of working directory:')
        directory = ''
        while not os.path.isdir(directory):
            directory = '/{}/'.format(input(prompt).strip('/'))
            if not os.path.isdir(directory):
                print('Specified working directory not recognised')
            else:
                run_parameters['workingdirectory'] = directory
                break

    # Requires user input if the absolute file path of the (locally saved) PDB
    # database is not specified in the input file / is not recognised
    if 'pdbdatabase' in run_parameters:
        if not os.path.isdir(run_parameters['pdbdatabase']):
            print('Specified directory for PDB database not recognised')
            run_parameters.pop('pdbdatabase')
    if not 'pdbdatabase' in run_parameters:
        print('Specify absolute file path of PDB database:')
        pdb_database = ''
        while not os.path.isdir(pdb_database):
            pdb_database = '/{}/'.format(input(prompt).strip('/'))
            if not os.path.isdir(pdb_database):
                print('Specified directory for PDB database not recognised')
            else:
                run_parameters['pdbdatabase'] = pdb_database
                break

    # Requires user input if the absolute file path of the (locally saved) DSSP
    # database is not specified in the input file / is not recognised
    if 'dsspdatabase' in run_parameters:
        if not os.path.isdir(run_parameters['dsspdatabase']):
            print('Specified directory for DSSP database not recognised')
            run_parameters.pop('dsspdatabase')
    if not 'dsspdatabase' in run_parameters:
        print('Specify absolute file path of DSSP database:')
        dssp_database = ''
        while not os.path.isdir(dssp_database):
            dssp_database = '/{}/'.format(input(prompt).strip('/'))
            if not os.path.isdir(dssp_database):
                print('Specified directory for DSSP database not recognised')
            else:
                run_parameters['dsspdatabase'] = dssp_database
                break

    # Requires user input if the resolution threshold for the dataset to be
    # generated is not specified in the input file / is not recognised
    if 'resolution' in run_parameters:
        try:
            resn = float(run_parameters['resolution'])
            if resn <= 0:
                print('Specified resolution cutoff must be greater than 0')
                run_parameters.pop('resolution')
        except ValueError:
            print('Specified resolution cutoff must be a number')
            run_parameters.pop('resolution')
    if not 'resolution' in run_parameters:
        print('Select resolution cutoff:')
        resn = 0
        while resn == 0:
            resn = input(prompt)
            try:
                resn = float(resn)
                if resn <= 0:
                    print('Specified resolution cutoff must be greater than 0')
                    resn = 0
                else:
                    run_parameters['resolution'] = resn
                    break
            except ValueError:
                print('Specified resolution cutoff must be a number')
                resn = 0

    # Requires user input if the R_factor (working value) threshold for the
    # dataset to be generated is not specified in the input file / is not
    # recognised
    if 'rfactor' in run_parameters:
        try:
            rfac = float(run_parameters['rfactor'])
            if rfac < 0 or rfac > 1:
                print('Specified Rfactor (working value) cutoff must be between 0 and 1')
                run_parameters.pop('rfactor')
        except ValueError:
            print('Specified Rfactor (working value) cutoff must be a number')
            run_parameters.pop('rfactor')
    if not 'rfactor' in run_parameters:
        print('Select Rfactor (working value) cutoff:')
        rfac = 0
        while rfac == 0:
            rfac = input(prompt)
            try:
                rfac = float(rfac)
                if rfac < 0 or rfac > 1:
                    rfac = 0
                    print('Specified Rfactor (working value) cutoff must be '
                          'between 0 and 1')
                else:
                    run_parameters['rfactor'] = rfac
                    break
            except ValueError:
                rfac = 0
                print('Specified Rfactor (working value) cutoff must be a '
                      'number')

    # Creates and / or sets the output directory as the current working
    # directory
    os.chdir('{}'.format(run_parameters['workingdirectory']))
    dir_name = '{}_{}_resn_{}_rfac_{}'.format(
        run_parameters['structuredatabase'], run_parameters['id'],
        run_parameters['resolution'], run_parameters['rfactor']
        )
    if stage == '1':
        if os.path.isdir(dir_name):
            shutil.rmtree(dir_name)
        os.mkdir(dir_name)
        os.chdir(dir_name)
    else:
        os.chdir(dir_name)

    # Writes run parameters to a txt file
    with open('Run_parameters_stage_{}.txt'.format(stage), 'w') as parameters_file:
        parameters_file.write('Structure database: {}\n'.format(run_parameters['structuredatabase']) +
                              'ID: {}\n'.format(run_parameters['id']) +
                              'Working directory: {}\n'.format(run_parameters['workingdirectory']) +
                              'PDB database: {}\n'.format(run_parameters['pdbdatabase']) +
                              'DSSP database: {}\n'.format(run_parameters['dsspdatabase']) +
                              'Resolution: {}\n'.format(run_parameters['resolution']) +
                              'Rfactor: {}\n'.format(run_parameters['rfactor']))

    return stage, run_parameters


def find_cd_hit_input(stage, args):
    # Locates input file of CDHIT filtered FASTA sequences required for stage
    # 2 of the analysis pipeline
    if vars(args)['sequences']:
        files = ['/{}'.format(input_file.strip('/')) for input_file in
                 vars(args)['sequences']]
        cdhit_entries = ''
        cdhit_output = ''
        for input_file in files:
            if input_file[-4:] == '.pkl':
                cdhit_entries = input_file
            elif input_file[-4:] == '.txt':
                cdhit_output = input_file
        if cdhit_entries != '' and not os.path.isfile(cdhit_entries):
            print('Absolute path to CDHIT input pkl file not recognised')
            cdhit_entries = ''
        if cdhit_output != '' and not os.path.isfile(cdhit_output):
            print('Absolute file path to CDHIT output txt file not recognised')
            cdhit_output = ''
    else:
        cdhit_entries = ''
        cdhit_output = ''

    while not os.path.isfile(cdhit_entries):
        print('Specify absolute file path of input pkl file of FASTA '
              'sequences fed into CDHIT')
        cdhit_entries = '/{}'.format(input(prompt).strip('/'))
        if not os.path.isfile(cdhit_entries):
            print('Specified file path not recognised')
        else:
            break

    while not os.path.isfile(cdhit_output):
        print('Specify absolute file path of txt file of filtered FASTA '
             'sequences output from CDHIT')
        cdhit_output = '/{}'.format(input(prompt).strip('/'))
        if not os.path.isfile(cdhit_output):
            print('Specified file path not recognised')
        else:
            break

    return cdhit_entries, cdhit_output


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
        # these criteria suitable for uploading to the cd_hit web server.
        beta_structure = filter_beta_structure(self.run_parameters)
        filtered_domain_df = beta_structure.resn_rfac_filter(domain_df)
        beta_structure.gen_cd_hit_list(filtered_domain_df)


    def run_stage_1_scope(self, orig_dir):
        return


    def run_stage_2(self, cdhit_entries, cdhit_output):
        if __name__ == 'subroutines.run_stages':
            from subroutines.extract_coordinates import extract_beta_structure_coords
            from subroutines.DSSP import (
                filter_dssp_database, beta_structure_dssp_classification
                )
            from subroutines.generate_network import manipulate_beta_structure
        else:
            from datagen.subroutines.extract_coordinates import extract_beta_structure_coords
            from datagen.subroutines.DSSP import (
                filter_dssp_database, beta_structure_dssp_classification
                )
            from datagen.subroutines.generate_network import manipulate_beta_structure

        # Loads the dataframe generated in previous steps
        filtered_domain_df = pd.read_pickle(cdhit_entries)

        # Obtains xyz coordinates for the sequences output from CD-HIT
        if os.path.isdir('Entire_domains'):
            shutil.rmtree('Entire_domains')
        os.mkdir('Entire_domains')

        beta_structure = extract_beta_structure_coords(self.run_parameters)
        cd_hit_domain_df = beta_structure.gen_cd_hit_dict(
            cdhit_output, filtered_domain_df
            )
        cd_hit_domain_df, all_atoms_dfs_dict = beta_structure.get_xyz_coords(
            cd_hit_domain_df
            )

        # Copies required DSSP files from database (on hard drive) to local
        # machine
        if os.path.isdir('DSSP_files'):
            shutil.rmtree('DSSP_files')
        os.mkdir('DSSP_files')

        filtered_files = filter_dssp_database(self.run_parameters)
        dssp_domain_df = filtered_files.copy_files_from_dssp_database(cd_hit_domain_df)

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

        # Combines the beta-strands into sheets and translates the identified
        # beta-strand interactions into a network
        beta_structure = manipulate_beta_structure(self.run_parameters)
        beta_structure.identify_strand_interactions(sec_struct_dfs_dict)


    def run_stage_3(self):
        with open(
            'CATH_{}_resn_{}_rfac_{}_domain_networks_dict.pkl'.format(
                run, resn, rfac
                ), 'rb'
            ) as pickle_file:
            (sec_struct_dfs_dict, domain_networks_dict, domain_sheets_dict) = pickle.load(pickle_file)

        beta_structure = calculate_solvent_accessible_surface_area(
            run=run, resn=resn, rfac=rfac, sec_struct_dfs_dict=sec_struct_dfs_dict,
            domain_sheets_dict=domain_sheets_dict
            )
        beta_structure.run_naccess()
        return


    def run_stage_4(self):
        return
