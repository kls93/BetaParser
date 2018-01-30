
import pandas as pd
from collections import OrderedDict
if __name__ == 'subroutines.OPM':
    from subroutines.run_stages import run_stages
else:
    from datagen.subroutines.run_stages import run_stages

class extract_strand_tilt_and_TM_regions(run_stages):

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
                    tm_residues = tm_residues + [chain+str(num) for num in
                                                 range(res_min, res_max+1)]
                tm_lists.append(tm_residues)

        opm_df = pd.DataFrame({'PDB_CODE': pdb_codes,
                               'CHAIN': chains,
                               'TILT_ANGLE': tilt_angles,
                               'TM_RANGE': tm_lists,
                               'SEGMENTS': tm_segment_lists})

        return opm_df

    def find_strand_tilt(self, sec_struct_dfs_dict, opm_df):
        # Determines strand tilt
        pdb_codes_list = opm_df['PDB_CODE'].tolist()
        tilt_angles = OrderedDict()

        for domain_id in list(sec_struct_dfs_dict.keys()):
            pdb_code = domain_id[0:4]

            if pdb_code in pdb_codes_list:
                index = pdb_codes_list.index(pdb_code)
                tilt_angle = opm_df['TILT_ANGLE'][index]
                tilt_angles[domain_id] = tilt_angle

            else:
                tilt_angles[domain_id] = 'Undefined'

        return tilt_angles

    def find_barrel_strand_and_shear_number(self, sec_struct_dfs_dict, opm_df):
        strand_numbers = OrderedDict()
        shear_numbers = OrderedDict()

        for row in range(opm_df.shape[0]):
            tm_segments = opm_df['SEGMENTS'][row]
            strand_numbers[opm_df['PDB_CODE'][row]] = len(tm_segments)

        # starting_strand = min(tm_segments, key=len)

        return strand_numbers, shear_numbers
