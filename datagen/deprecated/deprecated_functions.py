
    def find_z_axis(strands, sheets_df):
        # Finds axis through the barrel pore as the line that passes through
        # the average xyz coordinates of the middle two residues of each
        # strand.
        count = 0
        coords_res_1 = []
        coords_res_n = []

        # Determines the indices of the central two residues in each strand.
        for strand in strands:
            count += 1

            strand_df = sheets_df[sheets_df['STRAND_NUM'] == strand]
            strand_df = strand_df.reset_index(drop=True)

            row_num = strand_df.shape[0]
            # If there is an odd number of residues, the central residue plus
            # its immediate C-terminal residue are selected as the two central
            # residues
            if row_num % 2 == 1:
                row_num += 1

            index_1 = (row_num / 2) - 1
            index_n = row_num / 2
            if count % 2 == 1:  # Ensures that residues closer to the
                # periplasmic side of the membrane are grouped together,
                # likewise for residues closer to the extracellular side of the
                # membrane. Note that this assumes that neighbouring strands
                # interact with one another in an antiparallel hydrogen bonding
                # arrangement (however, the code will only break if the
                # majority of strands interact with a parallel instead of a
                # parallel hydrogen bonding arrangement, which is highly unlikely)
                index_1 = row_num / 2
                index_n = (row_num / 2) - 1

            res_1_x = strand_df['XPOS'][index_1]
            res_1_y = strand_df['YPOS'][index_1]
            res_1_z = strand_df['ZPOS'][index_1]
            res_n_x = strand_df['XPOS'][index_n]
            res_n_y = strand_df['YPOS'][index_n]
            res_n_z = strand_df['ZPOS'][index_n]

            coords_res_1.append((res_1_x, res_1_y, res_1_z))
            coords_res_n.append((res_n_x, res_n_y, res_n_z))

        # Calculates average xyz coordinates
        coords_res_1 = np.array(coords_res_1)
        coords_res_n = np.array(coords_res_n)
        x_coord_1 = np.sum(coords_res_1[:, 0]) / coords_res_1.shape[0]
        y_coord_1 = np.sum(coords_res_1[:, 1]) / coords_res_1.shape[0]
        z_coord_1 = np.sum(coords_res_1[:, 2]) / coords_res_1.shape[0]
        x_coord_n = np.sum(coords_res_n[:, 0]) / coords_res_n.shape[0]
        y_coord_n = np.sum(coords_res_n[:, 1]) / coords_res_n.shape[0]
        z_coord_n = np.sum(coords_res_n[:, 2]) / coords_res_n.shape[0]
        xyz_coords = [x_coord_1, y_coord_1, z_coord_1, x_coord_n, y_coord_n,
                      z_coord_n]

        return xyz_coords


class sandwich_interior_exterior_calcs():

    def find_z_axis(sheets_df):
        # Finds axis between the two sheets of the beta-sandwiches as the
        # average vector between its hydrogen bonded residues
        """
        # Generates dataframe of longest strand and its longest neighbouring
        # strand
        strands = set(sheets_df['STRAND_NUM'].tolist())
        strand_len = [0]*max(strands)

        for row in range(sheets_df.shape[0]):
            if sheets_df['ATMNAME'][row] == 'CA':
                strand_len[(sheets_df['STRAND_NUM'][row]) - 1] += 1

        longest_strand = strand_len.index(max(strand_len)) + 1
        sheets_df = sheets_df[~sheets_df['STRAND_NUM'].isin([longest_strand])]
        sheets_df = sheets_df.reset_index(drop=True)
        sheets_df = sheets_df.iloc[::5, :]
        sheets_df = sheets_df.reset_index(drop=True)
        """

        # Extracts hydrogen bonded residues
        h_bonds_dict = {}
        for row in range(sheets_df.shape[0]):
            if sheets_df['SHEET?'][row] == 'E':
                res_1 = sheets_df['DSSP_NUM'][row]
                res_2 = sheets_df['BRIDGE_PAIRS'][row][0]
                res_3 = sheets_df['BRIDGE_PAIRS'][row][1]

                if res_2 != '0':
                    pair_1 = [str(min([int(res_1), int(res_2)])),
                              str(max([int(res_1), int(res_2)]))]
                    if pair_1[0] not in list(h_bonds_dict.keys()):
                        h_bonds_dict[pair_1[0]] = pair_1[1]
                if res_3 != '0':
                    pair_2 = [str(min([int(res_1), int(res_3)])),
                              str(max([int(res_1), int(res_3)]))]
                    if pair_2[0] not in list(h_bonds_dict.keys()):
                        h_bonds_dict[pair_2[0]] = pair_2[1]

        # Calculates vectors between C_alpha atoms of hydrogen bonded residues
        vectors_array = np.empty((len(list(h_bonds_dict.keys())), 3))
        for index, res_1 in enumerate(list(h_bonds_dict.keys())):
            res_2 = h_bonds_dict[res_1]
            xyz_1 = np.array([])
            xyz_2 = np.array([])

            for row in range(sheets_df.shape[0]):
                if sheets_df['DSSP_NUM'][row] == res_1:
                    x_1 = sheets_df['XPOS'][row]
                    y_1 = sheets_df['YPOS'][row]
                    z_1 = sheets_df['ZPOS'][row]
                    xyz_1 = np.array([[x_1],
                                      [y_1],
                                      [z_1]])
                elif sheets_df['DSSP_NUM'][row] == res_2:
                    x_2 = sheets_df['XPOS'][row]
                    y_2 = sheets_df['YPOS'][row]
                    z_2 = sheets_df['ZPOS'][row]
                    xyz_2 = np.array([[x_2],
                                      [y_2],
                                      [z_2]])
                    break

            # Excludes beta-bridge pairs from the z-axis calculation (since
            # these residues will be included in the DSSP hydrogen-bonding
            # interactions, but not in the dataframe of residues in the
            # beta-sheets)
            if xyz_1.size and xyz_2.size:
                # TODO: Ensure that only consider vectors pointing in the same
                # direction (select central vector between two largest strands,
                # then if subsequent vectors are > 90 apart then ignore?)
                vectors_array[index][0] = xyz_2[0][0] - xyz_1[0][0]
                vectors_array[index][1] = xyz_2[1][0] - xyz_1[1][0]
                vectors_array[index][2] = xyz_2[2][0] - xyz_1[2][0]

        vector = np.mean(vectors_array, axis=0)

        return vector

    def align_sandwich(domain_id, vector, unprocessed_list):
        # Aligns the sandwich to z-axis using the reference axis through the
        # sandwich core calculated in the previous step
        print('Aligning {} sandwich with z-axis'.format(domain_id))

        # Creates AMPAL object from sandwich. Protons are added (using reduce)
        # to allow DataGen to use hydrogen HA3 of a glycine residue as a proxy
        # for the C_alpha atoms of the chiral amino acids in the interior- /
        # exterior-facing calculation performed in the following step.
        sandwich = isambard.external_programs.reduce.assembly_plus_protons(
            'Beta_strands/{}.pdb'.format(domain_id)
        )

        # Aligns the sandwich with z-axis
        s1 = [0.0, 0.0, 0.0]
        e1 = [vector[0], vector[1], vector[2]]
        s2 = [0.0, 0.0, 0.0]
        e2 = [0.0, 0.0, 1.0]
        translation, angle, axis, point = isambard.geometry.find_transformations(
            s1, e1, s2, e2
        )
        sandwich.rotate(angle, axis, point=point)  # ROTATION MUST BE PERFORMED
        # BEFORE TRANSLATION
        sandwich.translate(translation)

        # Generates dictionary of residues and their new xyz coordinates
        sandwich_pdb_string = sandwich.make_pdb().split('\n')
        xy_dict = OrderedDict()

        for line in sandwich_pdb_string:
            if line[0:6].strip() in ['ATOM', 'HETATM']:
                if len(line) > 80:
                    unprocessed_list.append(domain_id)
                    break
                else:
                    atom_id = line[21:27].replace(' ', '') + '_' + line[12:16].strip()

                    x_coord_atom = float(line[30:38])
                    y_coord_atom = float(line[38:46])
                    xy_coords = np.array([[x_coord_atom],
                                          [y_coord_atom]])
                    xy_dict[atom_id] = xy_coords

        return xy_dict, unprocessed_list

    def calc_int_ext(domain_id, dssp_df, xy_dict):
        # Determines whether a residue is interior- or exterior-facing, based
        # upon whether or not its C_beta atom lies on the same side of the
        # plane between its N and C atoms as the centre of mass of the sandwich

        print('Determining interior and exterior residues in {}'.format(domain_id))

        int_ext_dict = {}

        # Calculates centre of mass of C_alpha atoms in sandwich
        x_sum = 0
        y_sum = 0
        count = 0
        for res in list(xy_dict.keys()):
            if 'CA' in res:
                x_sum += xy_dict[res][0][0]
                y_sum += xy_dict[res][1][0]
                count += 1
        x_avrg = x_sum / count
        y_avrg = y_sum / count

        # For each residue, calculates whether the centre of mass and its
        # C-beta atom lie on the same side of the plane formed between the
        # residue's N and C atoms (= 'interior'), or on different sides
        # (= 'exterior') (using the equation d=(x−x1)(y2−y1)−(y−y1)(x2−x1),
        # where (x1, y1) and (x2, y2) describe the plane and (x, y) the point
        # of interest (see https://math.stackexchange.com/questions/274712/
        # calculate-on-which-side-of-a-straight-line-is-a-given-point-located
        # ?noredirect=1&lq=1))
        res_ids_dict = {}
        for row in range(dssp_df.shape[0]):
            if dssp_df['RES_ID'][row] not in res_ids_dict:
                res_ids_dict[dssp_df['RES_ID'][row]] = dssp_df['RESNAME'][row]

        for res in list(set(dssp_df['RES_ID'].tolist())):
            n = False
            c = False
            cb = False

            xy_sub_dict = {key: value for key, value in xy_dict.items() if res
                           in key}
            for atom in list(xy_sub_dict.keys()):
                if atom.endswith('_N'):
                    n = True
                    n_dist = math.sqrt(((xy_sub_dict[atom][0][0] - x_avrg)**2)
                                       + ((xy_sub_dict[atom][1][0] - y_avrg)**2))
                    n_x_coord = xy_sub_dict[atom][0][0]
                    n_y_coord = xy_sub_dict[atom][1][0]
                elif atom.endswith('_C'):
                    c = True
                    c_dist = math.sqrt(((xy_sub_dict[atom][0][0] - x_avrg)**2)
                                       + ((xy_sub_dict[atom][1][0] - y_avrg)**2))
                    c_x_coord = xy_sub_dict[atom][0][0]
                    c_y_coord = xy_sub_dict[atom][1][0]
                elif (
                    (atom.endswith('_CB'))
                    or
                    (res_ids_dict[atom.split('_')[0]] == 'GLY' and atom.endswith('HA3'))
                ):
                    cb = True
                    cb_x_coord = xy_sub_dict[atom][0][0]
                    cb_y_coord = xy_sub_dict[atom][1][0]

            if all(x is True for x in [n, c, cb]):
                distance_com = (((x_avrg-n_x_coord)*(c_y_coord-n_y_coord)) -
                                ((y_avrg-n_y_coord)*(c_x_coord-n_x_coord)))
                distance_cb = (((cb_x_coord-n_x_coord)*(c_y_coord-n_y_coord)) -
                               ((cb_y_coord-n_y_coord)*(c_x_coord-n_x_coord)))

                if (
                    (distance_com < 0 and distance_cb < 0)
                    or
                    (distance_com > 0 and distance_cb > 0)
                ):
                    int_ext_dict[res] = 'interior'  # 2 points lie on the same side of the plane
                else:
                    int_ext_dict[res] = 'exterior'  # 2 points lie either side of the plane

        return int_ext_dict


def merge_models_biological_assembly_pdbs(pdb_file_lines):
    # Combines multiple models in biological assembly PDB files into a single
    # model (which requires updating the chain ids and atom numbers)
    alphabet = string.ascii_uppercase + string.ascii_lowercase + '0123456789'
    new_pdb_file_lines = []

    pdb_file_lines = pdb_file_lines.replace('ENDMDL', 'MODEL')
    sections = pdb_file_lines.split('MODEL')
    sections = [section for section in sections if section.strip() != '\n']

    for index, section in enumerate(sections):
        lines = section.split('\n')
        if index == 0:
            new_pdb_file_lines += lines
        elif index == len(sections) - 1 and 'END' in section:
            new_pdb_file_lines += lines
        else:
            model = int(lines[0])
            max_atom_num = int(lines[-1][6:11])
            lines = lines[1:]
            if model == 1:
                chain_ids = []
                for line in lines:
                    if line[0:3] != 'TER':
                        new_pdb_file_lines.append(line)
                    else:
                        new_pdb_file_lines.append('TER'.ljust(80))

                    chain_id = line[21:22]
                    if chain_id in alphabet:
                        alphabet = alphabet.replace(chain_id, '')
                    if not chain_id in chain_ids:
                        chain_ids.append(chain_id)

            else:
                new_chain_ids = []
                for line in lines:
                    if line[0:3] == 'TER':
                        new_pdb_file_lines.append('TER'.ljust(80))
                    else:
                        max_atom_num += 1
                        new_atom_num = str(max_atom_num)
                        if len(new_atom_num) > 5:
                            new_atom_num = new_atom_num[-5:]
                        new_atom_num = new_atom_num.rjust(5)

                        new_chain_id = alphabet[chain_ids.index(line[21:22])]
                        if not new_chain_id in new_chain_ids:
                            new_chain_ids.append(new_chain_id)

                        new_line = (line[0:6] + new_atom_num + line[11:21] +
                                    new_chain_id + line[22:])
                        new_pdb_file_lines.append(new_line)

                for chain_id in new_chain_ids:
                    alphabet = alphabet.replace(chain_id, '')

    return new_pdb_file_lines

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
