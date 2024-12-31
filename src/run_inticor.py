#!/usr/bin/env python

""" 
Prep the paths and libraries
"""
import os
import sys
import argparse

from funcs_shared import *
from funcs_rwr import *
from funcs_dis_spe_itcs import *

os.chdir('.')
sys.path.append('../')


""" Main function """
def run_inticor(**kwargs):

    """ Prep kwargs """
    background_network_dir = kwargs['background_network']
    genes_of_interest_list_dir = kwargs['disease_genes_of_interest']
    disgenes_df_dir = kwargs['disease_genes_full_collection']
    parallel_num = kwargs['parallel_num']

    """ Prep variables """
    background_network_name = background_network_dir.rsplit('/', 1)[1].split('.')[0]
    genes_of_interest_name = genes_of_interest_list_dir.rsplit('/', 1)[1].split('.')[0]
    genes_of_interest_list = acquire_list_from_file(genes_of_interest_list_dir)
    itcs_list = acquire_list_from_file(itcs_list_dir); itcs_list.sort()
    disgenes_df = csv_to_pd(disgenes_df_dir)
    diseases_list, disgenes_df = prep_diseases_n_disgenes(disgenes_df)
    rwr_base_med = 2.2883992816268055e-05

    """ Prep dirs """
    res_upper_dir = "../results"
    create_dir_if_absent(res_upper_dir)

    cur_res_dir = f'{res_upper_dir}/{background_network_name}_{genes_of_interest_name}'
    create_dir_if_absent(cur_res_dir)

    """ Run RWRs for all genes as seeds & save for efficient analyses """
    cur_res_dir_rwr = f'{cur_res_dir}/A_rwr_collection'
    create_dir_if_absent(cur_res_dir_rwr)

    all_nodes_list = \
        prep_all_rwr_results(background_network_dir, cur_res_dir, cur_res_dir_rwr, parallel_num)

    """ Compare RWR values at the diseases genes from ITCs and Random genes """
    cur_res_dir_comparison = f'{cur_res_dir}/B_comparison'
    create_dir_if_absent(cur_res_dir_comparison)

    compare_rwr_at_dg_from_itcs_rgs(
        cur_res_dir, cur_res_dir_rwr, cur_res_dir_comparison, 
        itcs_list, all_nodes_list, diseases_list, disgenes_df,
        rwr_base_med, parallel_num)

    return 


""" Prep executable code """
if __name__ == '__main__':

    # Argument parsing
    parser = argparse.ArgumentParser(
                description="Run InTiCoR to acquire specific inter-tissue communicators for a given set of genes from a given in-silico biological network. \
                                Please make sure to check the instruction on the repository before running.",
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-b", "--background_network", 
                            help="Provide the directory for the background network",
                            type=str, default='../prep_files/example_network.csv')
    
    parser.add_argument("-g", "--disease_genes_of_interest", 
                            help="Provide the directory with the list of genes related to the disease of interest",
                            type=str, default='../prep_files/example_genes_of_interest.csv')

    parser.add_argument("-d", "--disease_genes_full_collection", 
                            help="Provide the directory with collection of genes for all available diseases",
                            type=str, default='../prep_files/DiseaseGenes_parsed.csv')

    parser.add_argument("-p", "--parallel_num", 
                            help="Provide a number of cores to use for parallel processing",
                            type=int, default=50)

    args = parser.parse_args()
    args_vars = vars(args)

    # Pass on the arguments into the function
    run_inticor(**args_vars)