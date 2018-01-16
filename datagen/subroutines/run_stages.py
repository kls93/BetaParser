
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

    # Sets run parameters
    run_parameters = {}
    if vars(args)['input_file']:
        try:
            with open('{}'.format(vars(args)['input_file']), 'r') as input_file:
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
    # specified in the input file
    if not 'structuredatabase' in run_parameters:
        print('CATH or SCOPe database?')
        database = ''
        while database not in ['cath', 'scope']:
            database = input(prompt).lower()
            if not database in ['cath', 'scope']:
                print('DataGen can currently only parse the CATH and '
                      'SCOPe databases\n'
                      '- please select one of these databases to continue')
            else:
                run_parameters['structuredatabase'] = database
                break

    # Requires user input if the (all-beta) structural domain the user wishes
    # to analyse is not specified in the input file
    if not 'id' in run_parameters:
        if run_parameters['structuredatabase'] == 'cath':
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
        elif run_parameters['structuredatabase'] == 'scope':
            print('Specify SCOPe code:')
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
    # not specified in the input file
    if not 'workingdirectory' in run_parameters:
        print('Specify absolute file path of working directory:')
        directory = ''
        while not os.path.isdir(directory):
            directory = '/{}/'.format(input(prompt).strip('/'))
            if not os.path.isdir(directory):
                print('Specified directory not recognised')
            else:
                run_parameters['workingdirectory'] = directory
                break

    # Requires user input if the absolute file path of the (locally saved) PDB
    # database is not specified in the input file
    if not 'pdbdatabase' in run_parameters:
        print('Specify absolute file path of PDB database:')
        pdb_database = ''
        while not os.path.isdir(pdb_database):
            pdb_database = '/{}/'.format(input(prompt).strip('/'))
            if not os.path.isdir(pdb_database):
                print('Specified directory not recognised')
            else:
                run_parameters['pdbdatabase'] = pdb_database
                break

    # Requires user input if the absolute file path of the (locally saved) DSSP
    # database is not specified in the input file
    if not 'dsspdatabase' in run_parameters:
        print('Specify absolute file path of DSSP database:')
        dssp_database = ''
        while not os.path.isdir(dssp_database):
            dssp_database = '/{}/'.format(input(prompt).strip('/'))
            if not os.path.isdir(dssp_database):
                print('Specified directory not recognised')
            else:
                run_parameters['dsspdatabase'] = dssp_database
                break

    # Requires user input if the resolution threshold for the dataset to be
    # generated is not specified in the input file
    if not 'resolution' in run_parameters:
        print('Select resolution cutoff:')
        resn = 0
        while resn == 0:
            resn = input(prompt)
            try:
                resn = float(resn)
                if resn <= 0:
                    resn = 0
                    print('Specified resolution cutoff must be greater than 0')
                else:
                    run_parameters['resolution'] = resn
                    break
            except ValueError:
                resn = 0
                print('Specified resolution cutoff must be a number')

    # Requires user input if the R_factor (working value) threshold for the
    # dataset to be generated is not specified in the input file
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

    # Creates output directory
    if stage == '1':
        os.chdir('{}'.format(run_parameters['workingdirectory']))
        dir_name = '{}_{}_resn_{}_rfac_{}'.format(
            run_parameters['structuredatabase'], run_parameters['id'],
            run_parameters['resolution'], run_parameters['rfactor']
            )
        if os.path.isdir(dir_name):
            shutil.rmtree(dir_name)
        os.mkdir(dir_name)
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


class run_stages():

    def __init__(self, run_parameters):
        self.code_parameters = run_parameters
        self.code = self.code_parameters['id']
        self.pdb_database = self.code_parameters['pdbdatabase']
        self.dssp_database = self.code_parameters['dsspdatabase']
        self.resn = float(self.code_parameters['resolution'])
        self.rfac = float(self.code_parameters['rfactor'])

    def run_stage_1_cath(self, orig_dir):
        # Runs stage 1 of the DataGen pipeline, extracting sequences of the
        # structural domain of interest from the CATH database
        from subroutines.CATH import gen_domain_desc_list, domain_desc_filter
        from subroutines.CDHIT import filter_beta_structure

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
        beta_structure = filter_beta_structure(self.resn, self.rfac, domain_df,
                                               self.pdb_database)
        filtered_domain_df = beta_structure.resn_rfac_filter()
        beta_structure.gen_cd_hit_list(filtered_domain_df)


    def run_stage_1_scope(self, orig_dir):
        return


    def run_stage_2(self):
        from subroutines.extract_coordinates import extract_beta_structure_coords
        from subroutines.DSSP import (filter_dssp_database,
                                      beta_structure_dssp_classification)
        from subroutines.generate_network import manipulate_beta_structure

        # Loads the dataframe generated in previous steps
        filtered_domain_df = pd.read_pickle(
            'Filtered_datasets_pre_cd_hit.pkl'.format(run, resn, rfac)
            )

        # Obtains xyz coordinates for the sequences output from CD-HIT
        beta_structure = extract_beta_structure_coords(self.code, self.resn,
                                                       self.rfac,
                                                       self.pdb_database)
        cd_hit_domain_df = beta_structure.gen_cd_hit_dict(filtered_domain_df)

        if os.path.isdir('CD_HIT_DSEQS'):
            shutil.rmtree('CD_HIT_DSEQS')
        os.mkdir('CD_HIT_DSEQS')

        cd_hit_domain_df, pdb_dfs_dict = beta_structure.get_xyz_coords(
            cd_hit_domain_df
            )

        # Copies required DSSP files from database (on hard drive) to local machine
        if os.path.isdir('DSSP_files'):
            shutil.rmtree('DSSP_files')
        os.mkdir('DSSP_files')

        filtered_files = filter_dssp_database(self.code, self.resn, self.rfac,
                                              self.dssp_database)
        dssp_domain_df = filtered_files.copy_files_from_dssp_database(
            cd_hit_domain_df
            )

        # Extracts beta-strands (as classified by DSSP) from the beta-structure domains
        beta_structure = beta_structure_dssp_classification(
            self.code, self.resn, self.rfac
            )
        dssp_residues_dict = beta_structure.extract_dssp_file_lines(
            self.dssp_domain_df
            )

        shutil.rmtree('DSSP_files')

        if os.path.isdir('DSSP_filtered_DSEQS'):
            shutil.rmtree('DSSP_filtered_DSEQS')
        os.mkdir('DSSP_filtered_DSEQS')

        dssp_dfs_dict = beta_structure.get_dssp_sec_struct_df(
            dssp_residues_dict, pdb_dfs_dict
            )
        beta_structure.write_dssp_sec_struct_pdb(dssp_dfs_dict)

        # Combines the beta-strands into sheets and translates the identified
        # beta-strand interactions into a network
        beta_structure = manipulate_beta_structure(self.code, self.resn, self.rfac)
        beta_structure.identify_strand_interactions(dssp_dfs_dict)


    def run_stage_3(self):
        return


    def run_stage_4(self):
        return
