
import sys
import argparse

# Script to run the DataGen pipeline to extract a dataset of PDB structures of
# an (all-beta) structural domain defined in either the CATH or the SCOPe
# database. A summary of the structural characteristics of each structure in
# the dataset (for e.g. further bioinformatics analysis) is also generated.


def main():
    if __name__ == '__main__':
        from subroutines.run_parameters import (
            gen_run_parameters, determine_if_discard_tm, find_cdhit_input,
            find_radius
        )
        from subroutines.run_stages import run_stages
    else:
        from datagen.subroutines.run_parameters import (
            gen_run_parameters, determine_if_discard_tm, find_cdhit_input,
            find_radius
        )
        from datagen.subroutines.run_stages import run_stages
    orig_dir = sys.path[0].split('datagen')[0]

    # Reads in command line inputs
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', help='OPTIONAL: Specifies the '
                        'absolute file path of an input file listing run '
                        'parameters')
    parser.add_argument('--tm', help='OPTIONAL: Specifies whether or not to '
                        'retain only TM structures')
    parser.add_argument('-s', '--sequences', nargs='+', help='OPTIONAL: '
                        'Specifies the absolute file path of an input file of '
                        'CDHIT filtered FASTA sequences required for stage 2 '
                        'of the analysis pipeline')
    parser.add_argument('--radius', help='OPTIONAL: Specifies the radius of '
                        'the sphere drawn around each atom when calculating '
                        'the identities of its nearest neighbouring residues')
    args = parser.parse_args()

    # Extracts run parameters and initialises run_stages object
    stage, run_parameters = gen_run_parameters(args)
    analysis = run_stages(run_parameters)

    # Generates txt file of sequences appropriately formatted to be fed into
    # the CDHIT web server
    if stage in ['1']:
        # Determines whether the user wants to extract structures from the CATH
        # or from the SCOPe structural database
        discard_tm = determine_if_discard_tm(args, run_parameters)
        if run_parameters['structuredatabase'] == 'CATH':
            analysis.run_stage_1_cath(orig_dir, discard_tm)
        elif run_parameters['structuredatabase'] == 'SCOP':
            analysis.run_stage_1_scope(orig_dir, discard_tm)

    # Extracts PDB structures and structural information for each of the
    # sequences listed in a CDHIT output txt file
    elif stage in ['2']:
        cdhit_entries, cdhit_output = find_cdhit_input(args)
        analysis.run_stage_2(cdhit_entries, cdhit_output)

    # Runs naccess upon each structure to calculate the solvent accessible
    # surface area of its beta-sheets and thus identify those which interact
    elif stage in ['3']:
        radius = find_radius(args)
        analysis.run_stage_3(radius)

    # Summarises the structural characteristics of the dataset in an output
    # dataframe / csv file
    elif stage in ['4']:
        analysis.run_stage_4(orig_dir)


# Calls 'main' function if datagen.py is run as a script
if __name__ == '__main__':
    main()
