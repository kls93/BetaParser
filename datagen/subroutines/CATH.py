
import os
import pandas as pd
from collections import OrderedDict

if __name__ == 'subroutines.CATH':
    from subroutines.variables import gen_tm_pdb_codes_list
else:
    from datagen.subroutines.variables import gen_tm_pdb_codes_list


def gen_domain_desc_list(orig_dir):
    # Generates a list of the domain descriptions provided in
    # CATH_domain_description_v_4_2_0.txt
    domains_desc = []
    with open(
        '{}/docs/CATH_domain_desc_v_4_2_0.txt'.format(orig_dir), 'r'
    ) as domains_file:
        line_1 = False
        indv_domain_list = []
        for line in domains_file:
            if line.startswith('FORMAT'):
                line_1 = True
                start = True
            elif line.startswith('//'):
                start = False

            if line_1 is True:
                if start is True:
                    indv_domain_list.append(line)
                elif start is False:
                    domains_desc.append(''.join(indv_domain_list))
                    indv_domain_list = []

    domains_desc = [line for line in domains_desc if line != '']
    return domains_desc


def domain_desc_filter(cathcode, domains_desc, discard_non_tm):
    # Filters the domain descriptions list for beta-structures (either
    # sandwiches or barrels depending upon the user's choice), picking out PDB
    # accession codes and sequences, whose values are stored in a dataframe.

    tm_pdb_codes = gen_tm_pdb_codes_list()

    domain_pdb_ids = []
    domain_ids = []
    domain_chains = []
    domain_cathcodes = []
    domain_dseqs = []
    domain_sseqs = []
    domain_sseqs_start_stop = []
    for domain in domains_desc:
        # Prevents cathcodes at the same level of the hierarchy with
        # overlapping codes (e.g. 2.60.40.10 and 2.60.40.1090) from being
        # mistaken for one another
        for code in cathcode.split('_'):
            if (
                    (code.count('.') == 3 and 'CATHCODE  {}\n'.format(code) in domain)
                    or
                    (code.count('.') < 3 and 'CATHCODE  {}.'.format(code) in domain)
            ):
                # Discards structures whose PDB codes are not in the OPM
                # database if the user has set discard_non_tm to True
                domain_sublist = domain.split('\n')
                pdb_code = [line for line in domain_sublist if
                            line.startswith('DOMAIN')][0][10:14]
                if (
                        (discard_non_tm is True and pdb_code in tm_pdb_codes)
                        or
                        (discard_non_tm is False)
                ):
                    dseqs_list = []
                    sseqs_consec_list = []
                    sseqs_list = []
                    sseqs_start_stop_list = []

                    for index, line in enumerate(domain_sublist):
                        if line.startswith('SSEQS') and domain_sublist[index+1].startswith('SSEQS'):
                            sseqs_consec_list.append(index)
                    for index in sseqs_consec_list:
                        domain_sublist[index+1] = ''.join(domain_sublist[index:index+2])
                        domain_sublist[index] = ''

                    for line in domain_sublist:
                        if line.startswith('DOMAIN'):
                            domain_pdb_ids.append(line[10:14])
                            domain_ids.append(line[10:].strip())
                            chain = line[14:].strip()
                            chain = ''.join([char for char in chain if char.isalpha()])
                            domain_chains.append(chain)
                        elif line.startswith('CATHCODE'):
                            domain_cathcodes.append(line[10:])
                        elif line.startswith('DSEQS'):
                            line = line.replace('DSEQS', '')
                            line = line.replace(' ', '')
                            dseqs_list.append(line)
                        elif line.startswith('SSEQS'):
                            line = line.replace('SSEQS', '')
                            line = line.replace(' ', '')
                            sseqs_list.append(line)
                        elif line.startswith('SRANGE'):
                            line = line.replace('SRANGE', '')
                            start_stop = line.split()
                            start_stop = [item.strip() for item in start_stop if item.strip() != '']
                            sseqs_start_stop_list.append(start_stop)
                    dseqs = ''.join(dseqs_list)
                    domain_dseqs.append(dseqs)
                    domain_sseqs.append(sseqs_list)
                    domain_sseqs_start_stop.append(sseqs_start_stop_list)

    domain_df = pd.DataFrame(OrderedDict({'PDB_CODE': domain_pdb_ids,
                                          'DOMAIN_ID': domain_ids,
                                          'CHAIN': domain_chains,
                                          'CATHCODE': domain_cathcodes,
                                          'DSEQS': domain_dseqs,
                                          'SSEQS': domain_sseqs,
                                          'SSEQS_START_STOP': domain_sseqs_start_stop}))

    return domain_df
