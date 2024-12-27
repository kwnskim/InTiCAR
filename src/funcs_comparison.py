import sys
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
        cur_res_dir, rwr_res_all_ref_df, seed_list, rwr_file_temp):
    
    # Check if already organized for this seed list:
    rwr_res_df_dir = f'{cur_res_dir}/ITC_RWR.csv'

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


def 


def compare_rwr_at_dg_from_itcs_rgs(
        cur_res_dir, cur_res_dir_rwr, cur_res_dir_comparison, itcs_list):

    """ Prep dirs """
    #   RWR dataframe collection
    #       Target
    itc_rwr_seed_save_dir = f"{cur_res_dir_comparison}/ITCs_Seeds_RWR"
    create_dir_if_absent(itc_rwr_seed_save_dir)

    itc_rwr_seed_save_raw_dir = f"{itc_rwr_seed_save_dir}/Raw"
    create_dir_if_absent(itc_rwr_seed_save_raw_dir)

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
    
    """Prep results from required seed_list"""
    cur_seed_rwr_df = \
        acquire_and_save_overall_rwr_df(
            cur_res_dir, rwr_res_all_ref_df, itcs_list, rwr_file_temp)
    len_seed_list = len(cur_seed_rwr_df.columns.tolist())



    return