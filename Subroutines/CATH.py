
import copy
import pandas as pd

class select_beta_structures():

    def __init__(self, run):
        self.run = run

    # Generates a list of the domain descriptions provided in
    # CATH_domain_description_v_4_2_0.txt
    def domain_desc_list(self):
        domains_description = []
        with open('CATH_domain_desc_v_4_2_0.txt', 'r') as domains_file:
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
                        domains_description.append(''.join(indv_domain_list[0:len(indv_domain_list)]))
                        indv_domain_list = []

        domains_description = [line for line in domains_description if line != '']
        return domains_description

    # Filters the domain descriptions list for beta-structures (either
    # sandwiches or barrels depending upon the user's choice), picking out PDB
    # accession codes and sequences, whose values are stored in a dataframe.
    def domain_desc_filter(self, domains_description):
        domain_pdb_ids = []
        domain_chains = []
        domain_cathcodes = []
        domain_ids = []
        domain_dseqs = []
        domain_sseqs = []
        domain_sseqs_start_stop = []
        for domain in domains_description:
            if 'CATHCODE  {}'.format(self.run) in domain:
                dseqs_list = []
                sseqs_consec_list = []
                sseqs_list = []
                sseqs_start_stop_list = []

                domain_sublist = domain.split('\n')
                domain_sublist_copy = copy.copy(domain_sublist)

                for index, line in enumerate(domain_sublist):
                    if line.startswith('SSEQS') and domain_sublist[index+1].startswith('SSEQS'):
                        sseqs_consec_list.append(index)
                for index in sseqs_consec_list:
                    domain_sublist[index+1] = ''.join(domain_sublist_copy[index:index+2])
                    domain_sublist[index] = ''
                    domain_sublist_copy = copy.copy(domain_sublist)

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

        domain_df = pd.DataFrame({'PDB_CODE': domain_pdb_ids,
                                  'CHAIN': domain_chains,
                                  'CATHCODE': domain_cathcodes,
                                  'DOMAIN_ID': domain_ids,
                                  'DSEQS': domain_dseqs,
                                  'SSEQS': domain_sseqs,
                                  'SSEQS_START_STOP': domain_sseqs_start_stop})

        return domain_df
