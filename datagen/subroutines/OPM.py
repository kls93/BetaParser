
import pandas as pd
import scipy.stats
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


class calculate_barrel_geometry(run_stages):

    def __init__(self, run_parameters):
        run_stages.__init__(self, run_parameters)

    def find_barrel_strand_number(self, sec_struct_dfs_dict):
        strand_numbers = OrderedDict()

        for domain_id, sec_struct_df in sec_struct_dfs_dict.items():
            strands = [int(strand) for strand in
                       set(sec_struct_df['STRAND_NUM'].tolist()) if strand != '']
            strand_count = max(strands)
            strand_numbers[domain_id] = strand_count

        return strand_numbers

    def find_barrel_shear_number(self, sec_struct_dfs_dict):
        shear_numbers = OrderedDict()

        for domain_id, sec_struct_df in sec_struct_dfs_dict.items():
            print(domain_id)  # DELETE ME
            shear_estimates = []

            strand_1_df = sec_struct_df[sec_struct_df['STRAND_NUM']==1]
            strand_1_res = [int(num) for num in strand_1_df['DSSP_NUM'].tolist()]

            dssp_nums = sec_struct_df['DSSP_NUM'].tolist()

            for residue in strand_1_res:
                print(residue)  # DELETE ME
                index = dssp_nums.index(str(residue))
                bp = int(sec_struct_df['H-BONDS'][index][0])

                if bp != 0:
                    old_bp = residue
                    count = 0
                    strand_num = 0
                    shear = 0

                    while strand_num != 1:
                        count = count + 1
                        print(old_bp, bp)
                        index = dssp_nums.index(str(bp))
                        strand_num = sec_struct_df['STRAND_NUM'][index]
                        print(strand_num)  # DELETE ME

                        bp_a = int(sec_struct_df['H-BONDS'][index][0])
                        print(bp_a)  # DELETE ME
                        bp_b = int(sec_struct_df['H-BONDS'][index][1])
                        print(bp_b)  # DELETE ME
                        if bp_a == old_bp:
                            new_bp = bp_b
                            if count == 1:
                                ref_bp = bp_a
                        elif bp_b == old_bp:
                            new_bp = bp_a
                            if count == 1:
                                ref_bp = bp_b
                        else:
                            print('ERROR: Strand pairs do not match')
                            import sys  # DELETE ME
                            sys.exit()
                            break

                        print(new_bp)  # DELETE ME
                        if new_bp != 0:
                            old_bp = bp  # These two lines must stay in order!
                            bp = new_bp  # These two lines must stay in order!
                        elif new_bp == 0:
                            shear = shear + 2
                            end = True
                            try:
                                dssp_nums.index(str(bp-2))
                            except ValueError:
                                end = False
                                old_bp = old_bp - 2
                                bp = bp + 2

                            try:
                                dssp_nums.index(str(bp+2))
                            except ValueError:
                                end = False
                                old_bp = old_bp + 2
                                bp = bp - 2

                            if end is True:
                                diff = abs(bp - ref_bp)
                                print(bp)
                                print(ref_bp)
                                print(diff)
                                shear = shear + diff
                                break

                print('SHEAR_NUMBER')
                print(shear)  # DELETE ME
                import sys  # DELETE ME
                sys.exit()
                shear_estimates.append(shear)

            shear = scipy.stats.mode(shear_estimates)
            shear_numbers[domain_id] = shear

        return shear_numbers
