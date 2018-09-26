
import sys
import argparse

# Script to run the DataGen pipeline to extract a dataset of PDB structures of
# an (all-beta) structural domain defined in either the CATH or the SCOPe
# database. A summary of the structural characteristics of each structure in
# the dataset (for e.g. further bioinformatics analysis) is also generated.


def main():
    if __name__ == '__main__':
        from subroutines.run_parameters import (
            gen_run_parameters
        )
        from subroutines.run_stages import run_stages
    else:
        from datagen.subroutines.run_parameters import (
            gen_run_parameters
        )
        from datagen.subroutines.run_stages import run_stages
    orig_dir = sys.path[0].split('datagen')[0]

    # Reads in command line inputs
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', help='OPTIONAL: Specifies the '
                        'absolute file path of an input file listing run '
                        'parameters')
    parser.add_argument('--betadesigner', action='store_true', help='OPTIONAL: '
                        'Used to specify that the program is being run within '
                        'the BetaDesigner program')
    args = parser.parse_args()

    # Extracts run parameters and initialises run_stages object
    run_parameters = gen_run_parameters(args)
    stage = run_parameters['stage']
    analysis = run_stages(run_parameters)

    # Generates txt file of sequences appropriately formatted to be fed into
    # the CDHIT web server
    if stage in ['1']:
        if run_parameters['betadesigner'] is True:
            print('Analysis stage 1 cannot be run within BetaDesigner')
        else:
            # Determines whether the user wants to extract structures from the
            # CATH or from the SCOPe structural database
            if run_parameters['structuredatabase'] == 'CATH':
                analysis.run_stage_1_cath(orig_dir)
            elif run_parameters['structuredatabase'] == 'SCOP':
                analysis.run_stage_1_scope(orig_dir)

    # Extracts PDB structures and structural information for each of the
    # sequences listed in a CDHIT output txt file
    elif stage in ['2']:
        analysis.run_stage_2()

    # Runs naccess upon each structure to calculate the solvent accessible
    # surface area of its beta-sheets and thus identify those which interact
    elif stage in ['3']:
        analysis.run_stage_3()

    # Summarises the structural characteristics of the dataset in an output
    # dataframe / csv file
    elif stage in ['4']:
        analysis.run_stage_4(orig_dir)


# Calls 'main' function if datagen.py is run as a script
if __name__ == '__main__':
    main()
