
import os
import argparse

# Script to run the DataGen pipeline to extract a dataset of PDB structures of
# an (all-beta) structural domain defined in either the CATH or the SCOPe
# database. A summary of the structural characteristics of each structure in
# the dataset (for e.g. further bioinformatics analysis) is also generated.

def main():
    from subroutines.run_stages import gen_run_parameters, run_stages
    orig_dir = os.getcwd()

    # Reads in command line inputs
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--input_file', help='Specifies an input file '
                        'listing run parameters')
    args = parser.parse_args()

    # Initialises run_stages object and extracts run parameters
    stage, run_parameters = gen_run_parameters(args)
    analysis = run_stages(run_parameters)

    # Generates txt file of sequences appropriately formatted to be fed into
    # the CDHIT web server
    if stage in ['1']:
        # Determines whether the user wants to extract structures from the CATH
        # or from the SCOPe structural database
        if run_parameters['structuredatabase'] == 'cath':
            analysis.run_stage_1_cath(orig_dir)
        elif run_parameters['structuredatabase'] == 'scope':
            analysis.run_stage_1_scope(orig_dir)

    # Extracts PDB structures and structural information for each of the
    # sequences listed in a CDHIT output txt file
    elif stage in ['2']:
        analysis.run_stage_2()
    # Runs naccess upon each structure to calculate the solvent accessible
    # surface area of its beta-sheets and thus identify those which interact
    elif stage in ['3']:
        analysis.run_stage_3()
    # Analyses the summary of structural characteristics of the dataset via
    # random forest machine learning
    elif stage in ['4']:
        analysis.run_stage_4()

# Calls 'main' function if datagen.py is run as a script
if __name__ == '__main__':
    main()
