
import
if __name__ == 'subroutines.OPM':
    from subroutines.run_stages import run_stages
else:
    from datagen.subroutines.run_stages import run_stages

class extract_strand_tilt_and_TM_regions(run_stages):

    def __init__(self, run_paramters):
        run_stages.__init__(self, run_parameters)

    def parse_opm(self, orig_dir):
        # Extracts strand tilt and TM information from the OPM database
        # (docs/OPM_TM_subunits.txt)
        pdb_codes = []
        chains = []
        tilt_angles = []
        tm_residues = []
        with open('{}/docs/OPM_TM_subunits.txt'.format(orig_dir), 'r') as opm_file:
            for line in opm_file:
                line_segments = line.split('-')
                pdb_codes.append(line_segments[0][0:4])
                chain = line_segments[0][4:].strip()
                chains.append(chain)
                tilt_angle = line_segments[1].replace('Tilt:', '')
                tilt_angle = tilt_angle.replace('°', '')
                tilt_angles.append(tilt_angle.strip())
                tm_list = []
                tm_segments = ''.join(line_segments[2:])
                tm_segments = tm_segments.replace('Segments:', '')
                tm_segments = tm_segments.split(',')
                for segment in tm_segments:
                    start = False
                    stop = True
                    for char in list(segment):
                        if char == '(':
                            start = True
                        elif char == ')':
                            stop = False
                        if start is True and stop is False:
                            tm_list.append(char)
                            break
                    res_min = int(tm_list.split('-')[0])
                    res_max = int(tm_list.split('-')[1])
                    tm_range = [chain+str(num) for num in range(res_min, res_max+1)]
                    for residue in tm_range:
                        tm_residues.append(residue)

        opm_df = pd.DataFrame({'PDB_CODE': pdb_code,
                               'CHAIN': chain,
                               'TILT_ANGLE': tilt_angles,
                               'TM_RANGE':})

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
                # unprocessed_list

        return tilt_angles

    def find_strand_TM_regions(self, sec_struct_dfs_dict, opm_df):
        # Determines strand TM regions
        for domain_id in list(sec_struct_dfs_dict.keys()):
            sec_struct_df = sec_struct_dfs_dict[domain_id]
            pdb_code = domain_id[0:4]
            tm_or_ext = ['']* sec_struct_df.shape[0]

            if pdb_code in pdb_codes_list:
                index = pdb_codes_list.index(pdb_code)
                tm_strands = opm_df['TM_RANGE'][index]

                for row in range(sec_struct_df.shape[0]):
                    chain_res_num = (sec_struct_df['CHAIN'][row]
                                     + str(sec_struct_df['RESNUM'][row])
                                     + sec_struct_df['INSCODE'][row])
                    if sec_struct_df['ATMNAME'][row] == 'CA':
                        if chain_res_num in tm_strands:
                            tm_or_ext[row] = 'transmembrane'
                        else:
                            tm_or_ext[row] = 'external'
            else:
                # unprocessed_list

            tm_or_ext_df = pd.DataFrame({'TM_OR_EXT': tm_or_ext})
            sec_struct_df = pd.concat([sec_struct_df, tm_or_ext], axis=1)
            sec_struct_dfs_dict[domain_id] = sec_struct_df

        return sec_struct_dfs_dict
