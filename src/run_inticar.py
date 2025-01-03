#!/usr/bin/env python

""" 
Prep the paths and libraries
"""
import os
import sys
import argparse
import warnings
warnings.filterwarnings('ignore')

from funcs_shared import *
from funcs_rwr import *
from funcs_dis_spe_itcs import *

os.chdir('.')
sys.path.append('../')


""" Main function """
def run_inticar(**kwargs):

    """ Prep kwargs """
    background_network_dir = kwargs['background_network']
    genes_of_interest_list_dir = kwargs['disease_genes_of_interest']
    disgenes_df_dir = kwargs['disease_genes_full_collection']
    zthres = kwargs['modified_z_threshold']
    parallel_num = kwargs['parallel_num']

    """ Prep variables """
    background_network_name = background_network_dir.rsplit('/', 1)[1].split('.')[0]
    
    genes_of_interest_name = genes_of_interest_list_dir.rsplit('/', 1)[1].split('.')[0]
    genes_of_interest_list = acquire_list_from_file(genes_of_interest_list_dir)
    
    itcs_list_dir = '../prep_files/ITCs_list.csv'
    itcs_list = acquire_list_from_file(itcs_list_dir)
    
    disgenes_df = csv_to_pd(disgenes_df_dir)
    diseases_list, disgenes_df = prep_diseases_n_disgenes(disgenes_df)
    
    # Add the user's query
    diseases_list.append(genes_of_interest_name)
    queried_df = \
        pd.DataFrame([[genes_of_interest_name, genes_of_interest_list]], 
                        columns=['Disease', 'Genes'])
    disgenes_df = multiple_vertical_concat([disgenes_df, queried_df])
    
    rwr_base_med = 2.2883992816268055e-05

    """ Prep dirs """
    res_upper_dir = "../results"
    create_dir_if_absent(res_upper_dir)

    cur_res_dir = f'{res_upper_dir}/{background_network_name}___{genes_of_interest_name}'
    create_dir_if_absent(cur_res_dir)

    """ Run RWRs for all genes as seeds & save for efficient analyses """
    cur_res_dir_rwr = f'{cur_res_dir}/A_rwr_collection'
    create_dir_if_absent(cur_res_dir_rwr)

    all_nodes_list, itcs_list = \
        prep_all_rwr_results(background_network_dir, background_network_name, 
                                cur_res_dir, cur_res_dir_rwr, itcs_list, parallel_num)

    """ Compare RWR values at the diseases genes from ITCs and Random genes """
    cur_res_dir_comparison = f'{cur_res_dir}/B_comparison'
    create_dir_if_absent(cur_res_dir_comparison)

    disgenes_rwr_norm_df = \
        compare_rwr_at_dg_from_itcs_rgs(
            cur_res_dir, cur_res_dir_rwr, cur_res_dir_comparison, 
            itcs_list, all_nodes_list, diseases_list, disgenes_df,
            rwr_base_med, parallel_num)

    """ Find disease-specific ITCs """
    cur_res_dir_dis_spe_itcs = f'{cur_res_dir}/C_disease_specific_ITCs'
    create_dir_if_absent(cur_res_dir_dis_spe_itcs)
    
    disease_specific_itcs_dict = \
        acquire_disease_specific_itcs(
            cur_res_dir_dis_spe_itcs, itcs_list, diseases_list, 
            disgenes_rwr_norm_df, zthres, parallel_num)
        
    """ Return the ITCs for the queried genes """
    genes_of_interest_itcs = \
        disease_specific_itcs_dict[genes_of_interest_name]
    
    pd_to_csv(genes_of_interest_itcs, 
                f'{cur_res_dir}/ITCs_for_{genes_of_interest_name}.csv')
    
    return 


""" Prep executable code """
if __name__ == '__main__':

    # Argument parsing
    parser = argparse.ArgumentParser(
                description="Run InTiCAR to acquire specific inter-tissue communicators for a given set of genes from a given in-silico biological network. \
                                Please make sure to check the instruction on the repository before running.",
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-b", "--background_network", 
                            help="Provide the directory for the background network",
                            type=str, default='../prep_files/example_network.csv')
    
    parser.add_argument("-g", "--disease_genes_of_interest", 
                            help="Provide the directory with the list of genes related to the disease of interest",
                            type=str, default='../prep_files/example_GOI.csv')

    parser.add_argument("-d", "--disease_genes_full_collection", 
                            help="Provide the directory with collection of genes for all available diseases",
                            type=str, default='../prep_files/DiseaseGenes_parsed.csv')

    parser.add_argument("-t", "--modified_z_threshold", 
                        help="Provide the threshold for modified z-score",
                        type=int, default=5)
    
    parser.add_argument("-p", "--parallel_num", 
                            help="Provide a number of cores to use for parallel processing",
                            type=int, default=50)

    args = parser.parse_args()
    args_vars = vars(args)

    # Pass on the arguments into the function
    run_inticar(**args_vars)