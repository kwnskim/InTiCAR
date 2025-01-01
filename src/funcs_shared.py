import os
import pickle as pk
import pandas as pd
import numpy as np

from pathlib import Path
from datetime import datetime


def acquire_list_from_file(file_necessary):
    rf_open = open(file_necessary, 'rt', encoding='UTF8')
    rf = rf_open.readlines()
    return [r.strip('\n').strip() for r in rf]


def report_time(str):
    now = datetime.now()
    current_time = now.strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{current_time}] {str} \n")
    return


def check_if_file_exists(file_dir):
    return Path(file_dir).exists()


def check_if_dir_exists(req_dir):
    return os.path.exists(req_dir)


def create_dir_if_absent(req_dir):
    is_present = check_if_dir_exists(req_dir)
    if not is_present:
        os.makedirs(req_dir)
    return 


def save_as_pickle(data, file_name):
    with open(file_name, 'wb') as f:
        pk.dump(data, f, pk.HIGHEST_PROTOCOL)
    return


def import_pickle(file_necessary):

    with open(file_necessary, 'rb') as f:
        content = pk.load(f)

    return content


def split_to_size_n(lst, n):
    return [lst[i:i+n] for i in range(0, len(lst), n)]


def pd_from_listolist(listolist, colnames):
    return pd.DataFrame(listolist, columns = colnames)


def pd_to_csv(given_csv, given_file_name, inbool=False, sepa='|', heada=True):
    return given_csv.to_csv(
        path_or_buf=given_file_name, sep=sepa, index=inbool, header=heada)


def csv_to_pd(given_file_name, sepa="|", dtp=None, indcol=None, heada=0):
    return pd.read_csv(given_file_name, sep=sepa, dtype=dtp, index_col=indcol, header=heada)


def pd_from_dict_keyasindex(giv_dict, colnames, key_as_index=True):

    if key_as_index:
        return pd.DataFrame.from_dict(giv_dict, orient='index', columns=colnames)
    
    else:
        giv_df = pd.DataFrame.from_dict(giv_dict, orient='index')
        giv_df.reset_index(inplace=True)
        giv_df.columns = colnames

        return giv_df
    

def multiple_horizontal_concat(df_list):
    return pd.concat(df_list, axis=1)


def multiple_vertical_concat(df_list, igind=False):
    return pd.concat(df_list, axis=0, ignore_index=igind)


def acquire_file_list_in_dir(cur_dir, full_dir=True):

    if full_dir:
        file_path_list = \
            [os.path.join(cur_dir, ct) 
                for ct in os.listdir(cur_dir) 
                    if os.path.isfile(os.path.join(cur_dir, ct))]

        return file_path_list
    
    else: 
        file_name_list = \
            [ct for ct in os.listdir(cur_dir)
                if os.path.isfile(os.path.join(cur_dir, ct))]
        
        return file_name_list
    

def prep_diseases_n_disgenes(disgenes_df):

    # Split the genes into list before moving on :)
    disgenes_df['Genes'] = disgenes_df['Genes'].apply(lambda x: x.split(';;'))
    diseases_list = disgenes_df['Disease'].values.tolist()

    return diseases_list, disgenes_df


def modified_zscore(data, consistency_correction=1.4826):
    
    median = np.median(data)
    deviation_from_med = np.array(data) - median
    mad = np.median(np.abs(deviation_from_med))

    if mad != 0:
        mod_zscore = deviation_from_med/(consistency_correction*mad)
    else: 
        meanad = np.mean(np.abs(deviation_from_med))
        mod_zscore = deviation_from_med/(1.253314*meanad)

    return mod_zscore