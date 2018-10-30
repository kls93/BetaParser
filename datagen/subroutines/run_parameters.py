
import os
import sys
import shutil
import pandas as pd
from collections import OrderedDict

prompt = '> '


def gen_run_parameters(args):
    # Sets as many run parameters as possible from the input file (if provided)
    run_parameters = OrderedDict()
    if vars(args)['input_file'] is not None:
        try:
            with open(vars(args)['input_file'], 'r') as input_file:
                for line in input_file:
                    key = line.split(':')[0].replace(' ', '').lower()
                    value = line.split(':')[1].replace('\n', '').strip()
                    if key in ['workingdirectory', 'pdbaudatabase',
                               'pdbbadatabase', 'dsspdatabase', 'opmdatabase',
                               'ringdatabase', 'cdhitsequencefiles',
                               'dataframes']:  # Only include file paths in this list!
                        value = value.replace('\\', '/')  # For windows file paths
                        if key == 'cdhitsequencefiles':
                            try:
                                files = [file_name.strip() for file_name in
                                         value.split(',')]
                                value = ['/{}'.format(files[0].strip('/')),
                                         '/{}'.format(files[1].strip('/'))]
                            except:
                                value = []
                        elif key == 'dataframes':
                            value = '/{}'.format(value.strip('/'))
                        else:
                            value = '/{}/'.format(value.strip('/'))
                    elif key == 'betadesigner':
                        pass
                    else:
                        value = value.replace(' ', '').lower()
                    run_parameters[key] = value
        except FileNotFoundError:
            sys.exit('Absolute path to input file not recognised')

    # Determines if DataGen is being run within the BetaDesigner program.
    # ** THIS PARAMETER MUST BE SET FIRST! **
    if vars(args)['betadesigner'] is True:
        run_parameters['betadesigner'] = True
    else:
        run_parameters['betadesigner'] = False

    # Requires user input if the analysis stage is not specified in the input
    # file / is not recognised. Note that currently the analysis cannot be run
    # in its entirety owing to the need to run the data through programs (the
    # CDHIT web server and naccess) that I cannot run on my local machine.
    # ** THIS PARAMETER MUST BE SET SECOND! **
    if 'stage' in run_parameters:
        stage = run_parameters['stage']
        if not stage in ['1', '2', '3', '4']:
            print('Analysis stage not recognised')
            run_parameters.pop('stage')
    if not 'stage' in run_parameters:
        print('Specify analysis stage:')
        stage = ''
        while stage not in ['1', '2', '3', '4']:
            stage = input(prompt).lower()
            if not stage in ['1', '2', '3', '4']:
                print('Analysis stage not recognised')
            else:
                run_parameters['stage'] = stage
                break
    if run_parameters['betadesigner'] is True:
        if run_parameters['stage'] == '1':
            sys.exit('Only stages 2-4 are compatible with the BetaDesigner '
                     'software')

    # Requires user input if the structural database (CATH or SCOPe) is not
    # specified in the input file / is not recognised
    if run_parameters['betadesigner'] is False:
        if 'structuredatabase' in run_parameters:
            run_parameters['structuredatabase'] = run_parameters['structuredatabase'].upper()
            if run_parameters['structuredatabase'] not in ['CATH', 'SCOP']:
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
    else:
        run_parameters['structuredatabase'] = 'CATH'

    # Requires user input if the (all-beta) structural domain the user wishes
    # to analyse is not specified in the input file / is not recognised
    if run_parameters['betadesigner'] is False:
        ids_dict = {'CATH': '2',
                    'SCOP': 'b'}
        if 'id' in run_parameters:
                if (
                    type(run_parameters['id']) == list
                        and not
                        all(run_parameters['id'][index].startswith(ids_dict[run_parameters['structuredatabase']])
                            for index, code in enumerate(run_parameters['id']))
                    ) or (
                    type(run_parameters['id']) != list
                        and
                        run_parameters['id'][0] != ids_dict[run_parameters['structuredatabase']]
                ):
                    print('DataGen is currently only suitable for generation and '
                          'analysis of\nall-beta structures')
                    run_parameters.pop('id')
        if not 'id' in run_parameters:
            print('Specify list of CATHCODEs:')
            run = []
            while (
                       len(run) == 0
                    or any(run[index][0] != ids_dict[run_parameters['structuredatabase']] for index, code in enumerate(run))
            ):
                run = input(prompt).lower()
                run = run.replace(' ', '')
                run = run.replace('[', '')
                run = run.replace(']', '')
                run = [cathcode for cathcode in run.split(',')]
                if (
                        len(run) != 0
                    and all(run[index].startswith(ids_dict[run_parameters['structuredatabase']])
                            for index, code in enumerate(run))
                ):
                    run_parameters['id'] = run
                    break
                else:
                    print('DataGen is currently only suitable for '
                          'generation and analysis of\nall-beta structures')
        # Joins list of CATH / SCOP database codes together
        if type(run_parameters['id']) == list:
            run_parameters['id'] = '_'.join(run_parameters['id'])

    # Converts 'barrel'/'sandwich' labelling used in BetaDesigner program into
    # CATHCODEs that can be interpreted by DataGen
    else:
        if 'id' in run_parameters:
            if run_parameters['id'].lower() == 'barrel':
                id = '2.40'
                run_parameters['id'] = id
            elif run_parameters['id'].lower() == 'sandwich':
                id = '2.60'
                run_parameters['id'] = id

            if not run_parameters['id'] in ['2.40', '2.60']:
                print('Structure type not recognised')
                run_parameters.pop('id')

        if not 'id' in run_parameters:
            print('Beta-barrel or beta-sandwich?')
            id = ''
            while not id in ['2.40', '2.60']:
                id = input(prompt)
                if id.lower() == 'barrel':
                    id = '2.40'
                elif id.lower() == 'sandwich':
                    id = '2.60'

                if id in ['2.40', '2.60']:
                    run_parameters['id'] = id
                    break
                else:
                     print('Structure type not recognised')

    # Requires user to specify whether they want to analyse each input domain
    # in the context of the parent asymmetric unit or biological assembly. By
    # default is set to the asymmetric unit when BetaDesigner is run in the
    # context of the BetaDesigner program.
    if run_parameters['betadesigner'] is False:
        if 'auorba' in run_parameters:
            if not run_parameters['auorba'].lower() in ['au', 'ba']:
                print('DataGen can only analyse input domains in the context of '
                      'either the asymmetric unit or the biological assembly')
                run_parameters.pop('auorba')
        if not 'auorba' in run_parameters:
            print('Select whether to analyse each input domain in the context of '
                  'either the parent asymmetric unit (\'au\') or the parent '
                  'biological assembly (\'ba\')')
            auorba = ''
            while not auorba in ['au', 'ba']:
                auorba = input(prompt).lower()
                if auorba in ['au', 'ba']:
                    run_parameters['auorba'] = auorba
                    break
                else:
                    print('DataGen can only analyse input domains in the context of '
                          'either the asymmetric unit or the biological assembly')
    else:
        run_parameters['auorba'] = 'au'

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
    # database (asymmetric units) is not specified in the input file / is not
    # recognised
    if run_parameters['betadesigner'] is False:
        if 'pdbaudatabase' in run_parameters:
            if not os.path.isdir(run_parameters['pdbaudatabase']):
                print('Specified directory for PDB asymmetric unit database not recognised')
                run_parameters.pop('pdbaudatabase')
        if not 'pdbaudatabase' in run_parameters:
            print('Specify absolute file path of PDB asymmetric unit database:')
            pdb_au_database = ''
            while not os.path.isdir(pdb_au_database):
                pdb_au_database = '/{}/'.format(input(prompt).strip('/'))
                if not os.path.isdir(pdb_au_database):
                    print('Specified directory for asymmetric unit PDB database not recognised')
                else:
                    run_parameters['pdbaudatabase'] = pdb_au_database
                    break
    else:
        run_parameters['pdbaudatabase'] = ''

    # Requires user input if the absolute file path of the (locally saved) PDB
    # database (biological assemblies) is not specified in the input file / is
    # not recognised and the user has previously specified that they want to
    # analyse each input domain in the context of its parent biological assembly
    if run_parameters['betadesigner'] is False:
        if run_parameters['auorba'] == 'ba':
            if 'pdbbadatabase' in run_parameters:
                if not os.path.isdir(run_parameters['pdbbadatabase']):
                    print('Specified directory for PDB biological assembly database not recognised')
                    run_parameters.pop('pdbbadatabase')
            if not 'pdbbadatabase' in run_parameters:
                print('Specify absolute file path of PDB biological assembly database:')
                pdb_ba_database = ''
                while not os.path.isdir(pdb_ba_database):
                    pdb_ba_database = '/{}/'.format(input(prompt).strip('/'))
                    if not os.path.isdir(pdb_ba_database):
                        print('Specified directory for biological assembly PDB database not recognised')
                    else:
                        run_parameters['pdbbadatabase'] = pdb_ba_database
                        break
        else:
            run_parameters['pdbbadatabase'] = run_parameters['pdbaudatabase']
    else:
        run_parameters['pdbbadatabase'] = ''

    # Requires user input if the absolute file path of the (locally saved) DSSP
    # database is not specified in the input file / is not recognised
    if run_parameters['betadesigner'] is False:
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
    else:
        run_parameters['dsspdatabase'] = ''

    # Requires user input if the absolute file path of the (locally saved) OPM
    # database is not specified in the input file / is not recognised
    if 'opmdatabase' in run_parameters:
        if not os.path.isdir(run_parameters['opmdatabase']):
            print('Specified directory for OPM database not recognised')
            run_parameters.pop('opmdatabase')
    if not 'opmdatabase' in run_parameters:
        if (
            (run_parameters['betadesigner'] is False)
            or
            (run_parameters['betadesigner'] is True and run_parameters['id'] == '2.40')
        ):
            print('Specify absolute file path of OPM database:')
            opm_database = ''
            while not os.path.isdir(opm_database):
                opm_database = '/{}/'.format(input(prompt).strip('/'))
                if not os.path.isdir(opm_database):
                    print('Specified directory for OPM database not recognised')
                else:
                    run_parameters['opmdatabase'] = opm_database
                    break
        else:
            run_parameters['betadesigner'] = ''

    # Requires user input if the absolute file path of the (locally saved) RING
    # database is not specified in the input file / is not recognised
    if run_parameters['betadesigner'] is False:
        if 'ringdatabase' in run_parameters:
            if not os.path.isdir(run_parameters['ringdatabase']):
                print('Specified directory for RING database not recognised')
                run_parameters.pop('ringdatabase')
        if not 'ringdatabase' in run_parameters:
            print('Specify absolute file path of RING database:')
            ring_database = ''
            while not os.path.isdir(ring_database):
                ring_database = '/{}/'.format(input(prompt).strip('/'))
                if not os.path.isdir(ring_database):
                    print('Specified directory for RING database not recognised')
                else:
                    run_parameters['ringdatabase'] = ring_database
                    break
    else:
        run_parameters['ringdatabase'] = ''

    # Locates input file of CDHIT filtered FASTA sequences required for stage
    # 2 of the analysis pipeline
    if run_parameters['stage'] == '2' and run_parameters['betadesigner'] is False:
        cdhit_entries = ''
        cdhit_output = ''
        if 'cdhitsequencefiles' in run_parameters:
            files = run_parameters['cdhitsequencefiles']
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

        while not os.path.isfile(cdhit_entries) or not cdhit_entries.endswith('.pkl'):
            print('Specify absolute file path of input pkl file of FASTA '
                  'sequences fed into CDHIT')
            cdhit_entries = '/{}'.format(input(prompt).strip('/'))
            if not os.path.isfile(cdhit_entries):
                print('Specified file path not recognised')
            elif cdhit_entries[-4:] != '.pkl':
                print('Specified file is not a pkl file')
            else:
                break

        while not os.path.isfile(cdhit_output) or not cdhit_output.endswith('.txt'):
            print('Specify absolute file path of txt file of filtered FASTA '
                  'sequences output from CDHIT')
            cdhit_output = '/{}'.format(input(prompt).strip('/'))
            if not os.path.isfile(cdhit_output):
                print('Specified file path not recognised')
            elif cdhit_output[-4:] != '.txt':
                print('Specified file is not a txt file')
            else:
                break

        run_parameters['cdhitsequencefiles'] = {'cdhit_entries': cdhit_entries,
                                                'cdhit_output': cdhit_output}
    else:
        run_parameters['cdhitsequencefiles'] = ''

    # Requires user input if the absolute file path of the (locally saved) RING
    # database is not specified in the input file / is not recognised
    if run_parameters['stage'] == '2' and run_parameters['betadesigner'] is True:
        if 'dataframes' in run_parameters:
            if not os.path.isfile(run_parameters['dataframes']):
                print('Specified absolute file path of pkl file of '
                      'BetaDesigner output not recognised')
                run_parameters.pop('dataframes')
        if not 'dataframes' in run_parameters:
            print('Specify absolute file path of pkl file output from BetaDesigner:')
            bd_output = ''
            while not os.path.isfile(bd_output):
                bd_output = '/{}'.format(input(prompt).strip('/'))
                if not os.path.isfile(bd_output):
                    print('Specified absolute file path of pkl file of '
                          'BetaDesigner output not recognised')
                else:
                    run_parameters['dataframes'] = bd_output
                    break
    else:
        run_parameters['dataframes'] = ''

    # Requires user input if the resolution threshold for the dataset to be
    # generated is not specified in the input file / is not recognised
    if run_parameters['betadesigner'] is False:
        if 'resolution' in run_parameters:
            try:
                resn = float(run_parameters['resolution'])
                if resn <= 0:
                    print('Specified resolution cutoff must be greater than 0')
                    run_parameters.pop('resolution')
                else:
                    run_parameters['resolution'] = resn
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
    else:
        run_parameters['resolution'] = ''

    # Requires user input if the R_factor (working value) threshold for the
    # dataset to be generated is not specified in the input file / is not
    # recognised
    if run_parameters['betadesigner'] is False:
        if 'rfactor' in run_parameters:
            try:
                rfac = float(run_parameters['rfactor'])
                if rfac < 0 or rfac > 1:
                    print('Specified Rfactor (working value) cutoff must be between 0 and 1')
                    run_parameters.pop('rfactor')
                else:
                    run_parameters['rfactor'] = rfac
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
                    if rfac <= 0 or rfac > 1:
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
    else:
        run_parameters['rfactor'] = ''

    # Determines radius of sphere for location of nearest neighbours
    if run_parameters['stage'] == '3':
        if 'radius' in run_parameters:
            radius = run_parameters['radius']
            try:
                radius = float(radius)
                if radius <= 0:
                    print('Specified radius must be greater than 0')
                    run_parameters.pop('radius')
                else:
                    run_parameters['radius'] = radius
            except ValueError:
                print('Specified radius must be a number')
                run_parameters.pop('radius')
        if not 'radius' in run_parameters:
            radius = 0
            print('Select radius of sphere for identification of neighbouring '
                  'residues:')
            while radius == 0:
                radius = input(prompt)
                try:
                    radius = float(radius)
                    if radius <= 0:
                        print('Specified radius must be greater than 0')
                        radius = 0
                    else:
                        run_parameters['radius'] = radius
                        break
                except ValueError:
                    print('Specified radius must be a number')
                    radius = 0
    else:
        run_parameters['radius'] = ''

    # Requires user input if a suffix for the PDB files in the biological
    # assembly database is not specified in the input file
    if run_parameters['betadesigner'] is False:
        if run_parameters['auorba'] == 'ba':
            if not 'suffix' in run_parameters:
                print('Specify suffix of biological assembly PDB files:')
                suffix = input(prompt)
                run_parameters['suffix'] = suffix
        else:
            run_parameters['suffix'] = ''
    else:
        run_parameters['suffix'] = ''

    # Determines whether or not the user wants to keep only transmembrane
    # structures.
    if (
        run_parameters['betadesigner'] is False
        and run_parameters['id'][0:4] in ['2.40', 'b']
    ):
        if 'discardnontm' in run_parameters:
            try:
                run_parameters['discardnontm'] = run_parameters['discardnontm'].lower()
                if not run_parameters['discardnontm'] in ['yes', 'no', 'y', 'n', 'true', 'false']:
                    print('Discard TM structures selection not recognised')
                    run_parameters.pop('discardnontm')
            except SyntaxError:
                print('Discard TM structures selection not recognised')
                run_parameters.pop('discardnontm')

        if not 'discardnontm' in run_parameters:
            print('Discard non-TM structures?')
            discard_non_tm = ''
            while not discard_non_tm in ['yes', 'no', 'y', 'n', 'true', 'false']:
                discard_non_tm = input(prompt).lower()
                if discard_non_tm in ['yes', 'no', 'y', 'n', 'true', 'false']:
                    run_parameters['discardnontm'] = discard_non_tm
                    break
                else:
                    print('Input not recognised - please enter "yes" or "no"')

        if run_parameters['discardnontm'] in ['yes', 'y', 'true']:
            run_parameters['discardnontm'] = True
        elif run_parameters['discardnontm'] in ['no', 'n', 'false']:
            run_parameters['discardnontm'] = False
    else:
        run_parameters['discardnontm'] = False

    # Creates and / or sets the output directory as the current working
    # directory
    if run_parameters['betadesigner'] is True:
        dir_name = run_parameters['workingdirectory']
    else:
        dir_name = '{}{}_{}_resn_{}_rfac_{}_{}/'.format(
            run_parameters['workingdirectory'],
            run_parameters['structuredatabase'], run_parameters['id'],
            run_parameters['resolution'], run_parameters['rfactor'],
            run_parameters['auorba']
        )

    if run_parameters['stage'] == '1':
        if os.path.isdir(dir_name):
            shutil.rmtree(dir_name)
        os.mkdir(dir_name)
        os.chdir(dir_name)
    else:
        os.chdir(dir_name)

    # Writes run parameters to a txt file
    with open('Run_parameters_stage_{}.txt'.format(run_parameters['stage']), 'w') as parameters_file:
        # Don't change the parameter names listed in Run_parameters.txt without
        # also changing these parameter names in the main body of code.
        # NOTE doesn't list BetaDesigner parameter since it is set with a
        # command line flag rather than in the program input file.
        if run_parameters['betadesigner'] is True:
            run_parameters['structuredatabase'] = ''

        for key, value in run_parameters.items():
            if value != '':
                parameters_file.write('{}: {}\n'.format(key, value))
            elif value == '':
                if key == 'suffix' and run_parameters['betadesigner'] is False:
                    parameters_file.write('{}: {}\n'.format(key, value))

    return run_parameters
