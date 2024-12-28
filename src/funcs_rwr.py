import pandas as pd
import networkx as nx

from multiprocessing import Pool
from funcs_shared import *


def run_rwr_into_df_for_save(cur_input):

    # Parse input
    cur_g, cur_gene = cur_input

    try:
        rwr_res = nx.pagerank(cur_g, personalization={cur_gene: 1})
        rwr_res_df = pd_from_dict_keyasindex(rwr_res, [cur_gene], key_as_index=True)
        
        return [cur_gene, rwr_res_df]

    except ZeroDivisionError:

        if ";;" in cur_gene:
            individual_seeds = cur_gene.split(';;')
            individual_seeds_dict = {indseed: 1 for indseed in individual_seeds}

            try:
                rwr_res = nx.pagerank(cur_g, personalization=individual_seeds_dict)
                rwr_res_df = pd_from_dict_keyasindex(rwr_res, [cur_gene], key_as_index=True)
                rwr_res_df = rwr_res_df/len(individual_seeds)
                
                return [cur_gene, rwr_res_df]

            except ZeroDivisionError:
                return [cur_gene, pd.DataFrame()]
        else:
            return [cur_gene, pd.DataFrame()]


def prep_all_rwr_results(
        background_network_dir, cur_res_dir, cur_res_dir_rwr, parallel_num):

    # Prep background network
    background_network_df = \
        pd.read_csv(background_network_dir, sep='|', header=True)
    background_network = \
        nx.from_pandas_edgelist(background_network_df, 'Gene_1_ENSG', 'Gene_2_ENSG')

    # Prep all nodes as seeds
    seed_list = list(background_network.nodes())
    seed_list.sort()
    all_nodes_list = seed_list.copy()

    #   Split seed list for reference to recall the result from each seed later
    seed_list_pieces = split_to_size_n(seed_list, 100)
    reference_df_listolist = list()
    reference_df_cols = ['dict_num', 'seed_list']

    for slpind, slp in enumerate(seed_list_pieces):

        slpind_str = str(slpind).zfill(8)
        reference_df_listolist.append([slpind_str, ','.join(slp)])

        slp_rwr_dict_dir = f"{cur_res_dir_rwr}/RWR_result_dict_{slpind_str}.pickle"

        if not check_if_file_exists(slp_rwr_dict_dir):
            
            report_time(f"Getting {slpind}-th set of RWR results for {background_network}")
            slp_rwr_dict = dict()

            p = Pool(parallel_num)
            rwr_input_listolist = [[background_network, sl] for sl in slp]
            slp_rwr_res_dfs_listolist = \
                p.map(run_rwr_into_df_for_save, rwr_input_listolist)

            for slp_wdf in slp_rwr_res_dfs_listolist:
                slp_rwr_dict[slp_wdf[0]] = slp_wdf[1]
        
            save_as_pickle(slp_rwr_dict, slp_rwr_dict_dir)

        else:
            report_time(f"Getting {slpind}-th set of RWR results for {background_network} already done.")

    # Save the reference table
    reference_df = pd_from_listolist(reference_df_listolist, colnames=reference_df_cols)
    pd_to_csv(reference_df, f'{cur_res_dir}/reference.csv')

    return all_nodes_list