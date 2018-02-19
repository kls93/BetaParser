
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


def calculate_int_ext_barrel(domain_id, dssp_df, sheets,
                             sec_struct_dfs_dict, domain_sheets_dict,
                             unprocessed_list):
    # Determines which face of the beta-sheet forms the interior of the
    # barrel and which forms the exterior. **NOTE: this function assumes
    # that adjacent residues point in opposite directions.**
    # TO DO: UPDATE TO USE BACKBONE PHI AND PSI ANGLES TO DETERMINE IF THE
    # STRAND CONTAINS A BULGE OR NOT

    # Checks that correct number of sheets has been retained for analysis
    if len(sheets) != 1:
        print('ERROR: more than 1 sheet retained in {} following solvent '
              'accessibility calculation'.format(domain_id))
        sec_struct_dfs_dict[domain_id] = None
        for sheet in sheets:
            domain_sheets_dict[sheet] = None
        unprocessed_list.append(domain_id)

        return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

    # Initialises records of interior / exterior facing residues
    int_ext_list = ['']*dssp_df.shape[0]
    int_ext_dict = OrderedDict()

    # Calculates solvent accessibility of both faces of barrel
    sheet_sub_strings = []
    print('Calculating solvent accessible surface area for {}'.format(domain_id))

    with open('Beta_strands/{}.pdb'.format(sheets[0]), 'r') as pdb_file:
        for line in pdb_file.readlines():
            line_start = line[0:16]
            line_end = line[17:]
            new_line = line_start + ' ' + line_end
            sheet_sub_strings.append(new_line)
    sheet_string = ''.join(sheet_sub_strings)
    naccess_out = isambard.external_programs.naccess.run_naccess(
        sheet_string, 'rsa', path=False, include_hetatms=True
    )
    naccess_out = naccess_out.split('\n')
    naccess_out = [line for line in naccess_out if line != '']
    for line in naccess_out:
        if line[0:3] in ['RES', 'HEM']:
            chain = line[8:9].strip()
            res_num = line[9:13].strip()
            ins_code = line[13:14].strip()
            chain_res_num = chain + res_num + ins_code
            int_ext_dict[chain_res_num] = float(line[14:22])

    # Selects strands in barrel
    sheet_num = sheets[0].replace('{}_sheet_'.format(domain_id), '')
    sheets_df = dssp_df[dssp_df['SHEET_NUM'] == sheet_num]
    strands = [strand for strand in set(sheets_df['STRAND_NUM'].tolist())
               if strand != '']

    # For each strand, determines which of its residues face towards the
    # interior and which face towards the exterior (assuming that
    # consecutive residues face in opposite directions)
    for strand in strands:
        strand_df = dssp_df[dssp_df['STRAND_NUM'] == strand]
        res_acsblty = OrderedDict()
        for res_id in strand_df['RES_ID'].tolist():
            if res_id != '':
                res_acsblty[res_id] = int_ext_dict[res_id]

        # Splits the residues into opposite faces and calculates the sum of
        # their solvent accessibilities
        face_1_res = list(res_acsblty.values())[:: 2]
        face_2_res = list(res_acsblty.values())[1:: 2]
        face_1_solv_acsblty = sum(face_1_res)
        face_2_solv_acsblty = sum(face_2_res)

        # Determines which face is exterior and which is interior (the sum
        # of the solvent accessibilities of the exterior facing residues
        # will be greater than that of the interior facing residues)
        if face_1_solv_acsblty > face_2_solv_acsblty:
            ext_chain_res_num = list(res_acsblty.keys())[:: 2]
            int_chain_res_num = list(res_acsblty.keys())[1:: 2]
        elif face_1_solv_acsblty < face_2_solv_acsblty:
            ext_chain_res_num = list(res_acsblty.keys())[1:: 2]
            int_chain_res_num = list(res_acsblty.keys())[:: 2]
        else:
            print('ERROR: solvent accessibilities of two faces of '
                  '{}_strand_{} are equal - something has gone wrong with '
                  'the analysis'.format(domain_id, strand))
            sec_struct_dfs_dict[domain_id] = None
            domain_sheets_dict[sheets[0]] = None
            unprocessed_list.append(domain_id)
            return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list

        # Assigns 'interior' and 'exterior' labels to the relevant res_id
        for chain_res_num in ext_chain_res_num:
            int_ext_dict[chain_res_num] = 'exterior'
        for chain_res_num in int_chain_res_num:
            int_ext_dict[chain_res_num] = 'interior'

    # Updates dataframe with solvent accessibility information
    for row in range(dssp_df.shape[0]):
        if (dssp_df['RES_ID'][row] in list(int_ext_dict.keys())
                and dssp_df['ATMNAME'][row] == 'CA'
                ):
            int_ext_list[row] = int_ext_dict[dssp_df['RES_ID'][row]]
    int_ext_df = pd.DataFrame({'INT_EXT': int_ext_list})
    dssp_df = pd.concat([dssp_df, int_ext_df], axis=1)
    sec_struct_dfs_dict[domain_id] = dssp_df

    return sec_struct_dfs_dict, domain_sheets_dict, unprocessed_list


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
