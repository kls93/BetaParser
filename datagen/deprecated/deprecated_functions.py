
    def merge_sheets(self, dssp_residues_dict, dssp_dfs_dict):
        # Merges sheet names of sheets that share strands
        dssp_dfs_merged_sheets_dict = OrderedDict()

        for domain_id in list(dssp_residues_dict.keys()):
            dssp_df = dssp_dfs_dict[domain_id]
            print('Merging sheets that share beta-strands in {}'.format(domain_id))

            sheets = [sheet for sheet in set(dssp_df['SHEET_NUM'].tolist()) if sheet != '']
            sheet_strands = []
            for sheet in sheets:
                locals()['sheet_{}_strands'.format(sheet)] = []
                for strand in dssp_df[dssp_df['SHEET_NUM'] == sheet]['STRAND_NUM'].tolist():
                    if strand not in locals()['sheet_{}_strands'.format(sheet)]:
                        locals()['sheet_{}_strands'.format(sheet)].append(strand)
                sheet_strands.append(locals()['sheet_{}_strands'.format(sheet)])
                del locals()['sheet_{}_strands'.format(sheet)]

            sheet_strand_numbers = [strand for strands in sheet_strands for strand in strands]
            count_1 = 0
            while len(sheet_strand_numbers) != len(set(sheet_strand_numbers)):
                strands = sheet_strands[count_1]
                sheet_strands_copy = copy.copy(sheet_strands)

                count_2 = 1
                restart = True
                while restart is True:
                    strands = sheet_strands[count_1]
                    for index in range(0, len(sheet_strands)):
                        if index == len(sheet_strands) - 1:
                            restart = False
                        if count_1 != index and strands != '':
                            intersect = set(strands).intersection(set(sheet_strands[index]))
                            if len(intersect) > 0:
                                sheet_strands_copy[count_1] = list(
                                    set(sheet_strands[count_1]) | set(sheet_strands[index]))
                                sheet_strands_copy[index] = ''
                                sheet_strands = sheet_strands_copy
                                restart = True

                sheet_strand_numbers = [strand for strands in sheet_strands for strand in strands]

            sheet_strands = [strands for strands in sheet_strands if len(strands) > 2]

            for row in range(dssp_df.shape[0]):
                for strands in sheet_strands:
                    if str(dssp_df['STRAND_NUM'][row]) in strands:
                        dssp_df.loc[row, 'SHEET_NUM'] = (sheet_strands.index(strands)+1)
                        break

            dssp_dfs_merged_sheets_dict[domain_id] = dssp_df
            dssp_df.to_pickle('DSSP_filtered_DSEQS/{}.pkl'.format(domain_id))

        return dssp_dfs_merged_sheets_dict


def calc_int_ext_sandwich():
    com_coords = np.array([[np.sum(list(xy_dict.values())[0][0])],
                           [np.sum(list(xy_dict.values())[1][0])]])

        for res_id in sheets_df['RES_ID'].tolist():
            print(domain_id, res_id)
            if 'GLY' not in res_id:
                res_df = dssp_df[dssp_df['RES_ID'] == res_id]
                res_atoms = res_df['ATMNAME'].tolist()

                if all([x in res_atoms for x in ['N', 'CA', 'C', 'O', 'CB']]):
                    n_coords = xy_dict['{}_N'.format(res_id)]
                    c_coords = xy_dict['{}_C'.format(res_id)]
                    ca_coords = xy_dict['{}_CA'.format(res_id)]
                    cb_coords = xy_dict['{}_CB'.format(res_id)]

                    n_com_vector = np.array([[com_coords[0][0] - n_coords[0][0]],
                                             [com_coords[1][0] - n_coords[1][0]]])
                    n_com_magnitude = math.sqrt(((n_com_vector[0][0])**2)
                                                + ((n_com_vector[1][0])**2))

                    c_com_vector = np.array([[com_coords[0][0] - c_coords[0][0]],
                                             [com_coords[1][0] - c_coords[1][0]]])
                    c_com_magnitude = math.sqrt(((c_com_vector[0][0])**2)
                                                + ((c_com_vector[1][0])**2))

                    ca_com_vector = np.array([[com_coords[0][0] - ca_coords[0][0]],
                                              [com_coords[1][0] - ca_coords[1][0]]])
                    ca_com_magnitude = math.sqrt(((ca_com_vector[0][0])**2)
                                                 + ((ca_com_vector[1][0])**2))

                    n_ca_vector = np.array([[ca_coords[0][0] - n_coords[0][0]],
                                            [ca_coords[1][0] - n_coords[1][0]]])
                    n_ca_magnitude = math.sqrt(((n_ca_vector[0][0])**2)
                                               + ((n_ca_vector[1][0])**2))

                    c_ca_vector = np.array([[ca_coords[0][0] - c_coords[0][0]],
                                            [ca_coords[1][0] - c_coords[1][0]]])
                    c_ca_magnitude = math.sqrt(((c_ca_vector[0][0])**2)
                                               + ((c_ca_vector[1][0])**2))

                    ca_cb_vector = np.array([[cb_coords[0][0] - ca_coords[0][0]],
                                             [cb_coords[1][0] - ca_coords[1][0]]])
                    ca_cb_magnitude = math.sqrt(((ca_cb_vector[0][0])**2)
                                                + ((ca_cb_vector[1][0])**2))

                    # Calculates com-N-C_alpha angle
                    com_n_ca_numerator = np.sum((n_com_vector*n_ca_vector), axis=0)[0]
                    com_n_ca_denominator = n_com_magnitude*n_ca_magnitude
                    com_n_ca_angle = math.degrees(
                        math.acos(com_n_ca_numerator / com_n_ca_denominator)
                    )

                    # Caculates com-C-C_alpha angle
                    com_c_ca_numerator = np.sum((c_com_vector*c_ca_vector), axis=0)[0]
                    com_c_ca_denominator = c_com_magnitude*c_ca_magnitude
                    com_c_ca_angle = math.degrees(
                        math.acos(com_c_ca_numerator / com_c_ca_denominator)
                    )

                    # Selects the largest angle of the two backbone angles
                    if n_com_magnitude < c_com_magnitude:
                        backbone_angle = com_c_ca_angle

                    # Calculates com-C_alpha-C_beta angle
                    com_ca_cb_numerator = np.sum((ca_com_vector*ca_cb_vector), axis=0)[0]
                    com_ca_cb_denominator = ca_com_magnitude*ca_cb_magnitude
                    com_ca_cb_angle = math.degrees(
                        math.acos(com_ca_cb_numerator / com_ca_cb_denominator)
                    )

                    # Determines whether the residue is interior- or
                    # exterior-facing
                    if com_ca_cb_angle < 90:
                        int_ext_dict[res_id] = 'interior'
                    elif com_ca_cb_angle > 90:
                        int_ext_dict[res_id] = 'exterior'


def find_barrel_shear_number(self, sec_struct_dfs_dict):
    for domain_id, sec_struct_df in sec_struct_dfs_dict.items():
        shear_estimates = []

        sec_struct_df = sec_struct_df[sec_struct_df['SHEET?'] == 'E']
        sec_struct_df = sec_struct_df.reset_index(drop=True)

        min_strand_num = min([num for num in set(sec_struct_df['STRAND_NUM'].tolist())
                              if num != ''])
        max_strand_num = max([num for num in set(sec_struct_df['STRAND_NUM'].tolist())
                              if num != ''])
        for row in range(sec_struct_df.shape[0]):
            if sec_struct_df['STRAND_NUM'][row] in [min_strand_num, max_strand_num]:
                res_1 = sec_struct_df['DSSP_NUM'][row]
                print(res_1)  # DELETE ME
                h_bonded_res_2 = int(sec_struct_df['H-BONDS'][row][0])
                h_bonded_res_3 = int(sec_struct_df['H-BONDS'][row][1])

                h_bonded_res_2_a = None
                h_bonded_res_2_b = None
                if h_bonded_res_2 != 0:
                    row_res_2 = (sec_struct_df['DSSP_NUM'].tolist()).index(str(h_bonded_res_2))
                    h_bonded_res_2_a = int(sec_struct_df['H-BONDS'][row_res_2][0])
                    h_bonded_res_2_b = int(sec_struct_df['H-BONDS'][row_res_2][1])

                h_bonded_res_3_a = None
                h_bonded_res_3_b = None
                if h_bonded_res_3 != 0:
                    row_res_3 = (sec_struct_df['DSSP_NUM'].tolist()).index(str(h_bonded_res_3))
                    h_bonded_res_3_a = int(sec_struct_df['H-BONDS'][row_res_3][0])
                    h_bonded_res_3_b = int(sec_struct_df['H-BONDS'][row_res_3][1])

                shear_pairs = [[h_bonded_res_2_b, h_bonded_res_3_a],
                               [h_bonded_res_2_a, h_bonded_res_3_b]]
                if (not any(x in [None, 0] for x in [dssp_num for pair in
                                                     shear_pairs for dssp_num in pair])
                    ):
                    for pair in shear_pairs:
                        print(pair)  # DELETE ME
                        if int(res_1) in pair:
                            print(pair[0], pair[1])  # DELETE ME
                            shear_estimate = abs(pair[0] - pair[1])
                            print(shear_estimate)  # DELETE ME
                            shear_estimates.append(shear_estimate)

        print(shear_estimates)  # DELETE ME
        import sys
        sys.exit()

        shear = scipy.stats.mode(shear_estimates)
        shear_numbers[domain_id] = shear

    import sys
    sys.exit()  # DELETE ME
    print(shear_numbers)  # DELETE ME


def find_barrel_shear_number():
    for domain_id, sec_struct_df in sec_struct_dfs_dict.items():
        print(domain_id)  # DELETE ME
        shear_estimates = []

        strand_1_df = sec_struct_df[sec_struct_df['STRAND_NUM'] == 4]
        strand_1_res = [int(num) for num in strand_1_df['DSSP_NUM'].tolist()]

        dssp_nums = sec_struct_df['DSSP_NUM'].tolist()

        for residue in strand_1_res:
            print(residue)  # DELETE ME
            index = dssp_nums.index(str(residue))
            bp = int(sec_struct_df['H-BONDS'][index][0])
            ref_bp = int(sec_struct_df['H-BONDS'][index][1])

            if 0 not in [bp, ref_bp]:
                old_bp = residue
                count = 0
                strand_num = 0
                shear = 0

                while strand_num != 4:
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
                    elif bp_b == old_bp:
                        new_bp = bp_a
                    else:
                        print('Failed to calculate shear number')
                        break

                    print(new_bp)  # DELETE ME
                    index = dssp_nums.index(str(bp))
                    strand_num = sec_struct_df['STRAND_NUM'][index]
                    if strand_num == 4:
                        diff = abs(old_bp - ref_bp)
                        shear = shear + diff
                        break

                    if new_bp != 0:
                        old_bp = bp  # These two lines must stay in order!
                        bp = new_bp  # These two lines must stay in order!
                        index = dssp_nums.index(str(bp))
                        strand_num = sec_struct_df['STRAND_NUM'][index]
                    elif new_bp == 0:
                        shear = shear + 2
                        try:
                            dssp_nums.index(str(bp-2))
                        except ValueError:
                            old_bp = old_bp - 2
                            bp = bp + 2

                        try:
                            dssp_nums.index(str(bp+2))
                        except ValueError:
                            end = False
                            old_bp = old_bp + 2
                            bp = bp - 2

                print('SHEAR_NUMBER')
                print(shear)  # DELETE ME
                import sys  # DELETE ME
                sys.exit()
                shear_estimates.append(shear)

        shear = scipy.stats.mode(shear_estimates)
        shear_numbers[domain_id] = shear


def calculate_alpha_angle(sheets_df, barrel):
    # Calculates average angle between beta-strands and central axis
    alpha_angles = []
    for strand in list(set(sheets_df['STRAND_NUM'].tolist())):
        strand_df = sheets_df[sheets_df['STRAND_NUM'] == strand]
        strand_df = strand_df.reset_index(drop=True)
        max_row = strand_df.shape[0] - 1

        x_1 = barrel['{}'.format(strand_df['CHAIN'][0])]['{}'.format(
            strand_df['RESNUM'][0])]['CA'].x
        y_1 = barrel['{}'.format(strand_df['CHAIN'][0])]['{}'.format(
            strand_df['RESNUM'][0])]['CA'].y
        z_1 = barrel['{}'.format(strand_df['CHAIN'][0])]['{}'.format(
            strand_df['RESNUM'][0])]['CA'].z
        x_n = barrel['{}'.format(strand_df['CHAIN'][max_row])]['{}'.format(
            strand_df['RESNUM'][max_row])]['CA'].x
        y_n = barrel['{}'.format(strand_df['CHAIN'][max_row])]['{}'.format(
            strand_df['RESNUM'][max_row])]['CA'].y
        z_n = barrel['{}'.format(strand_df['CHAIN'][max_row])]['{}'.format(
            strand_df['RESNUM'][max_row])]['CA'].z

        vector_1 = np.array([[x_n-x_1],
                             [y_n-y_1],
                             [z_n-z_1]])
        vector_2 = np.array([[0], [0], [1]])
        angle = isambard.tools.geometry.angle_between_vectors(
            vector_1, vector_2, radians=True
        )
        alpha_angles.append(angle)
    alpha = np.mean(alpha_angles)

    return alpha


def calc_shear_number(domain_id, sheets_df, alpha):
    print('Calculating {} barrel shear number'.format(domain_id))

    n = len(set(sheets_df['STRAND_NUM'].tolist()))
    shear = math.tan(alpha)*n*(4.4/3.3)
    shear = int(math.ceil(float(abs(shear)/2.0))*2.0)

    print(set(sheets_df['STRAND_NUM'].tolist()))
    print(n)
    print(alpha)
    print(shear)
