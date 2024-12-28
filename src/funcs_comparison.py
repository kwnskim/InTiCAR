import sys
import random
import pandas as pd
import networkx as nx

from multiprocessing import Pool
from funcs_shared import *


def prep_ref_index_dict(rwr_res_ref_df, giv_list):

    rwr_res_ref_df_parsed = rwr_res_ref_df.copy()

    # Split into seed list separated by ','
    rwr_res_ref_df_parsed['seed_list'] = \
        rwr_res_ref_df_parsed['seed_list'].apply(lambda x: x.split(','))
    
    # Find which of the elements in each row belongs to the list of interest
    rwr_res_ref_df_parsed['seed_list_prepped'] = \
        rwr_res_ref_df_parsed['seed_list'].apply(
            lambda x: [ele for ele in x if ele in giv_list])

    # Parse only the files with > 0 element
    rwr_res_ref_df_parsed['seed_list_prepped_len'] = \
        rwr_res_ref_df_parsed['seed_list_prepped'].map(len)
    
    len_mask_over_0 = rwr_res_ref_df_parsed['seed_list_prepped_len'] != 0

    rwr_res_ref_df_parsed = \
        rwr_res_ref_df_parsed[len_mask_over_0]

    return rwr_res_ref_df_parsed


def collect_rwr_from_giv_list(
        rwr_res_ref_df, giv_list, rwr_file_temp):

    # Acquire which dictionary's index has the genes in the given list
    rwr_res_ref_df_parsed = \
        prep_ref_index_dict(rwr_res_ref_df, giv_list)
    
    # If no dataframe file contains the required list as seed, stop.
    if rwr_res_ref_df_parsed.empty:
        sys.exit(f'Stopped, since no seed genes can be found in the current network')

    index_dict = \
        dict(zip(rwr_res_ref_df_parsed['dict_num'], 
                rwr_res_ref_df_parsed['seed_list_prepped']))
    
    # Collect the RWR results for the giv_list
    cur_giv_list_rwr_dict = dict()

    for cur_num in [*index_dict]:
    
        # Import the current number's RWR results
        cur_rwr_file_dir = rwr_file_temp.replace("##cur_num##", cur_num)
        cur_rwr_dict = import_pickle(cur_rwr_file_dir)

        # Which RWR results do we need?
        cur_non_seed_random = index_dict[cur_num]
        for cnsr in cur_non_seed_random:
            cur_giv_list_rwr_dict[cnsr] = cur_rwr_dict[cnsr]

    return cur_giv_list_rwr_dict


def organize_rwr_df_from_rwr_dict(rwr_dict):

    giv_seeds = [*rwr_dict]
    rwr_df_listolist = list()

    for gs in giv_seeds:
        rwr_df_listolist.append(rwr_dict[gs])

    return multiple_horizontal_concat(rwr_df_listolist)


def acquire_and_save_overall_rwr_df(
        itc_rwr_seed_save_dir, rwr_res_all_ref_df, seed_list, rwr_file_temp):
    
    # Check if already organized for this seed list:
    rwr_res_df_dir = f'{itc_rwr_seed_save_dir}/ITC_RWR.csv'

    if not check_if_file_exists(rwr_res_df_dir):
        
        report_time(f"Prepping RWR for the ITCs")

        #   Organize the results from the required seeds
        cur_seed_rwr_dict = \
            collect_rwr_from_giv_list(
                rwr_res_all_ref_df, seed_list, rwr_file_temp)

        #   Save the RWR results as dataframe
        cur_seed_rwr_df = \
            organize_rwr_df_from_rwr_dict(cur_seed_rwr_dict)
        pd_to_csv(cur_seed_rwr_df, rwr_res_df_dir, inbool=True)

    else:
        report_time(f"Importing RWR from previously calculated ITCs")
        cur_seed_rwr_df = csv_to_pd(rwr_res_df_dir, indcol=0)

    return cur_seed_rwr_df



def prep_random_rwr_df(input_list):

    cur_ind, random_pool_list, len_itcs_list, \
        rwr_res_all_ref_df, rand_rwr_seed_save_raw_dir, rwr_file_temp = input_list

    report_time(f'  Working on {cur_ind}-th random gene samples')
    cur_rand_ind = str(cur_ind).zfill(4)

    # Randomly select seeds
    random_seed_list = \
        random.sample(random_pool_list, len_itcs_list)

    #   Organize the results from the randomly selected seeds
    cur_random_rwr_res_df_dir = \
        f'{rand_rwr_seed_save_raw_dir}/{cur_rand_ind}_RWR.csv'
    
    cur_random_rwr_dict = \
        collect_rwr_from_giv_list(
            rwr_res_all_ref_df, random_seed_list, rwr_file_temp)

    #   Save the RWR results as dataframe
    cur_random_rwr_df = \
        organize_rwr_df_from_rwr_dict(cur_random_rwr_dict)
    pd_to_csv(cur_random_rwr_df, cur_random_rwr_res_df_dir, inbool=True)

    return


def acquire_and_save_overall_rwr_df_random(
        rand_rwr_seed_save_raw_dir, random_pool_list, len_itcs_list, rwr_file_temp, rwr_res_all_ref_df, parallel_num):

    # Check if 100 random samples have already been collected:
    random_rwr_file_list = \
        acquire_file_list_in_dir(rand_rwr_seed_save_raw_dir)

    if len(random_rwr_file_list) != 100: # 1 files per each of 100 sets

        report_time(f'Prepping RWR for 100 sets of random genes')        
        
        input_listolist = \
            [[cur_ind, random_pool_list, len_itcs_list,
                rwr_res_all_ref_df, rand_rwr_seed_save_raw_dir, rwr_file_temp] 
                    for cur_ind in range(100)]

        pool = Pool(parallel_num)
        pool.map(prep_random_rwr_df, input_listolist)
    
    else:
        report_time(f'Prepped RWR for 100 sets of random genes already done')

    return 


def get_disease_related_genes(disgenes_df, giv_dis):

    # Get the disease genes from the parsed dataframe
    giv_dis_mask = disgenes_df['Disease'] == giv_dis
    giv_dis_df = disgenes_df[giv_dis_mask]
    giv_dis_genes = giv_dis_df['Genes'].values.tolist()[0]

    return giv_dis_genes


def get_rwr_for_giv_genes(rwr_df, giv_ensg):

    cur_rwr_df = rwr_df.copy()

    # RWR for dataframes
    giv_ensg_mask = cur_rwr_df.index.isin(giv_ensg)
    rwr_df_giv_ensg_df = cur_rwr_df[giv_ensg_mask]

    return rwr_df_giv_ensg_df


def acquire_rwr_dis_df_vals_overthres(itcs_rwr_df_dis, rwr_base_med):

    # Get the values over the minimum threshold for meaningful relation
    itcs_rwr_df_dis_overthres_mask = \
        itcs_rwr_df_dis >= rwr_base_med
    itcs_rwr_df_dis_overthres = \
        itcs_rwr_df_dis[itcs_rwr_df_dis_overthres_mask]

    return itcs_rwr_df_dis_overthres


def save_rwr_dis_df_vals_overthres(
        dis, itcs_rwr_df, disgenes_list, rwr_base_med,
        itc_rwr_seed_save_raw_dir_disgenes,
        itc_rwr_seed_save_raw_dir_disgenes_overthres):
    
    dis_rwr_df_dir = \
        f'{itc_rwr_seed_save_raw_dir_disgenes}/{dis}_RWR.csv'
    dis_rwr_df_dir_overquant = \
        f'{itc_rwr_seed_save_raw_dir_disgenes_overthres}/{dis}_RWR.csv'

    if not check_if_file_exists(dis_rwr_df_dir_overquant):
        
        #   RWR values for the disease genes from the target seeds
        itcs_rwr_df_dis = \
            get_rwr_for_giv_genes(itcs_rwr_df, disgenes_list)
        
        #   Save for further use:
        pd_to_csv(itcs_rwr_df_dis, dis_rwr_df_dir, inbool=True)

        # Parse out the values overthres RWR
        itcs_rwr_df_dis_overthres = \
            acquire_rwr_dis_df_vals_overthres(itcs_rwr_df_dis, rwr_base_med)
        
        #   Save for further use:
        pd_to_csv(itcs_rwr_df_dis_overthres, dis_rwr_df_dir_overquant, inbool=True)

    else:
        #   Acquire the previously saved dataframe
        itcs_rwr_df_dis_overthres = csv_to_pd(dis_rwr_df_dir_overquant, indcol=0)

    return itcs_rwr_df_dis_overthres


def acquire_random_for_cur_dis(input_list):

    disgenes_list, cur_ind, rwr_quant_thres, \
    rand_rwr_seed_save_raw_dir, \
    rand_rwr_seed_save_raw_dir_disgenes_dis, \
    rand_rwr_seed_save_raw_dir_disgenes_overthres_dis = input_list

    cur_rand_ind = str(cur_ind).zfill(4)

    rand_rwr_df_dir_overquant = \
        f'{rand_rwr_seed_save_raw_dir_disgenes_overthres_dis}/{cur_rand_ind}_RWR.csv'

    if not check_if_file_exists(rand_rwr_df_dir_overquant):

        # Import the result for RWR for this round of random
        cur_random_rwr_res_df_dir = \
            f'{rand_rwr_seed_save_raw_dir}/{cur_rand_ind}_RWR.csv'
        cur_random_rwr_df = csv_to_pd(cur_random_rwr_res_df_dir, indcol=0)

        cur_random_rwr_df_dis_overthres = \
            save_rwr_dis_df_vals_overthres(
                cur_rand_ind, cur_random_rwr_df, disgenes_list, rwr_quant_thres,
                rand_rwr_seed_save_raw_dir_disgenes_dis,
                rand_rwr_seed_save_raw_dir_disgenes_overthres_dis)

    else:
        #   Acquire the previously saved dataframe
        cur_random_rwr_df_dis_overthres = csv_to_pd(rand_rwr_df_dir_overquant, indcol=0)
    
    return cur_random_rwr_df_dis_overthres


def compare_rwr_at_dg_from_itcs_rgs(
        cur_res_dir, cur_res_dir_rwr, cur_res_dir_comparison, 
        itcs_list, all_nodes_list, diseases_list, disgenes_df, 
        rwr_base_med, parallel_num):

    """ Prep dirs """
    #   RWR dataframe collection
    #       Target
    itc_rwr_seed_save_dir = f"{cur_res_dir_comparison}/ITCs_Seeds_RWR"
    create_dir_if_absent(itc_rwr_seed_save_dir)

    itc_rwr_seed_save_raw_dir_disgenes = f"{itc_rwr_seed_save_dir}/DiseaseGenes"
    create_dir_if_absent(itc_rwr_seed_save_raw_dir_disgenes)

    itc_rwr_seed_save_raw_dir_disgenes_overthres = f"{itc_rwr_seed_save_dir}/DiseaseGenes_Overthres"
    create_dir_if_absent(itc_rwr_seed_save_raw_dir_disgenes_overthres)

    #       Non-target random
    rand_rwr_seed_save_dir = f"{cur_res_dir_comparison}/RandomGenes_Seeds_RWR"
    create_dir_if_absent(rand_rwr_seed_save_dir)

    rand_rwr_seed_save_raw_dir = f"{rand_rwr_seed_save_dir}/Raw"
    create_dir_if_absent(rand_rwr_seed_save_raw_dir)

    rand_rwr_seed_save_raw_dir_disgenes = f"{rand_rwr_seed_save_dir}/DiseaseGenes"
    create_dir_if_absent(rand_rwr_seed_save_raw_dir_disgenes)

    rand_rwr_seed_save_raw_dir_disgenes_overthres = f"{rand_rwr_seed_save_dir}/DiseaseGenes_Overthres"
    create_dir_if_absent(rand_rwr_seed_save_raw_dir_disgenes_overthres)

    #   RWR comparison collection
    comparison_collection_dir = f"{cur_res_dir_comparison}/RWR_Comparison_Collection"
    create_dir_if_absent(comparison_collection_dir)

    """ Prep RWR from all acquired from A """
    rwr_file_temp = f'{cur_res_dir_rwr}/RWR_result_dict_##cur_num##.pickle'
    rwr_res_all_ref_dir = f"{cur_res_dir}/reference.csv"
    rwr_res_all_ref_df = csv_to_pd(rwr_res_all_ref_dir, dtp=str)
    
    """ Prep results from required seed_list """
    itcs_rwr_df = \
        acquire_and_save_overall_rwr_df(
            cur_res_dir, rwr_res_all_ref_df, itcs_list, rwr_file_temp)

    """ Prep results from 100 sets of random genes """
    random_pool_list = list(set(all_nodes_list) - set(itcs_list))
    random_pool_list.sort()
    len_itcs_list = len(itcs_list)

    acquire_and_save_overall_rwr_df_random(
        rand_rwr_seed_save_raw_dir, random_pool_list, len_itcs_list, 
            rwr_file_temp, rwr_res_all_ref_df, parallel_num)
    
    """ Work on each disease to acquire RWR values over threshold at disease genes, and save """
    len_diseases_list = len(diseases_list)
    for disind, dis in enumerate(diseases_list):

        if disind % 20 == 0:
            report_time(f'  Currently working on {disind}-th disease {dis} (Out of {len_diseases_list})')

        #   Acquire disease gene list
        disgenes_list = \
            get_disease_related_genes(disgenes_df, dis)
        
        #       Check if none of the disease genes exist in the network
        check_disgenes = set(all_nodes_list).intersection(set(disgenes_list))

        if len(check_disgenes) == 0:
            report_time(f'  No disease genes in the current network for {dis}')
            continue

        # Prep the RWR tables for the usage later:
        itcs_rwr_df_dis_overthres = \
            save_rwr_dis_df_vals_overthres(
                dis, itcs_rwr_df, disgenes_list, rwr_base_med,
                itc_rwr_seed_save_raw_dir_disgenes,
                itc_rwr_seed_save_raw_dir_disgenes_overthres)
        
        # Repeat for the random cases
        rand_rwr_seed_save_raw_dir_disgenes_dis = \
            f"{rand_rwr_seed_save_raw_dir_disgenes}/{dis}"
        create_dir_if_absent(rand_rwr_seed_save_raw_dir_disgenes_dis)

        rand_rwr_seed_save_raw_dir_disgenes_overthres_dis = \
            f"{rand_rwr_seed_save_raw_dir_disgenes_overthres}/{dis}"
        create_dir_if_absent(rand_rwr_seed_save_raw_dir_disgenes_overthres_dis)

        input_listolist = \
            [[disgenes_list, cur_ind, rwr_base_med,
                    rand_rwr_seed_save_raw_dir, 
                    rand_rwr_seed_save_raw_dir_disgenes_dis,
                    rand_rwr_seed_save_raw_dir_disgenes_overthres_dis] 
                        for cur_ind in range(100)]

        pool = Pool(50)
        random_rwr_df_dis_overthres_list = \
            pool.map(acquire_random_for_cur_dis, input_listolist)

    return