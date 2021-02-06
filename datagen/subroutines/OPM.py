
import pandas as pd
import networkx as nx
import numpy as np
from collections import OrderedDict
if __name__ == 'subroutines.OPM':
    from subroutines.run_stages import run_stages
else:
    from datagen.subroutines.run_stages import run_stages


class extract_barrel_info_from_OPM(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    def parse_opm(self, orig_dir):
        # Extracts strand tilt and TM information from the OPM database
        # (docs/OPM_TM_subunits.txt)
        pdb_codes = []
        chains = []
        tilt_angles = []
        tm_lists = []
        tm_segment_lists = []

        print('Creating datframe of OPM database information')

        with open('{}/docs/OPM_TM_subunits.txt'.format(orig_dir), 'r') as opm_file:
            for line in opm_file:
                line_segments = line.split('-')

                pdb_codes.append(line_segments[0][0:4])

                chain = line_segments[0][4:].strip()
                chains.append(chain)

                tilt_angle = line_segments[1].replace('Tilt:', '')
                tilt_angle = tilt_angle.replace('°', '')
                tilt_angles.append(tilt_angle.strip())

                tm_segments = '-'.join(line_segments[2:])
                tm_segments = tm_segments.replace('Segments:', '')
                tm_segments = tm_segments.split(',')
                tm_segment_lists.append(tm_segments)

                tm_residues = []
                for segment in tm_segments:
                    tm_range = ''
                    start = False
                    stop = False
                    for char in list(segment):
                        if char == ')':
                            stop = True

                        if start is True and stop is False:
                            tm_range += char

                        if char == '(':
                            start = True

                    res_min = int(tm_range.split('-')[0])
                    res_max = int(tm_range.split('-')[1])
                    tm_residues += [chain+str(num) for num in
                                    range(res_min, res_max+1)]  # No insertion
                    # code available in OPM_TM_subunits.txt
                tm_lists.append(tm_residues)

        opm_df = pd.DataFrame(OrderedDict({'PDB_CODE': pdb_codes,
                                           'CHAIN': chains,
                                           'TILT_ANGLE': tilt_angles,
                                           'TM_RANGE': tm_lists,
                                           'SEGMENTS': tm_segment_lists}))

        return opm_df

    def find_strand_tilt(self, sec_struct_dfs_dict, opm_df):
        # Determines strand tilt
        pdb_codes_list = opm_df['PDB_CODE'].tolist()
        tilt_angles = OrderedDict()

        for domain_id in list(sec_struct_dfs_dict.keys()):
            print('Calculating tilt angle for {}'.format(domain_id))

            pdb_code = domain_id[0:4]

            if pdb_code in pdb_codes_list:
                index = pdb_codes_list.index(pdb_code)
                tilt_angle = opm_df['TILT_ANGLE'][index]
                tilt_angles[domain_id] = tilt_angle
            else:
                tilt_angles[domain_id] = 'Undefined'

        return tilt_angles


class calculate_barrel_geometry(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    def find_barrel_strand_number(self, sec_struct_dfs_dict):
        # Calculates number of strands in barrel
        strand_numbers = OrderedDict()

        for domain_id, dssp_df in sec_struct_dfs_dict.items():
            print('Calculating number of strands in {}'.format(domain_id))
            strands = [strand for strand in set(dssp_df['STRAND_NUM'].tolist())
                       if strand != '']
            strand_count = len(strands)
            strand_numbers[domain_id] = strand_count

        return strand_numbers
