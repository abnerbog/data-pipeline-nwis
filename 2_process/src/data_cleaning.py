# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 12:15:06 2023

@author: AbnerBogan
"""

import os 
import pickle

def load_timeseries_info(in_file):
    """

    Load in timeseries data and metadata.

    Parameters
    ----------
    in_file : string
        Pickled Python file containing timeseries data and metadata.

    Returns
    -------
    timeseries_info : list
        Loaded timeseries data and metadata.

    """
    
    # load file via unpickling
    with open(in_file, 'rb') as file:
        timeseries_info = pickle.load(file)
        
    return timeseries_info

def remove_missing_data(df):
    """
    
    Remove rows from dataframe where measured value is missing.

    Parameters
    ----------
    df : pandas dataframe
        Dataframe containing daily timeseries data.

    Returns
    -------
    df_cleaned : pandas dataframe
        Dataframe containing cleaned daily timeseries data.

    """
    
    # create mask to identify rows with -999999 or nan in data column
    data_column = df.columns[-2]
    mask = (df[data_column] == -999999) | df[data_column].isna()
    
    # use mask to drop rows from dataframe
    df_cleaned = df[~mask]
    
    return df_cleaned

def save_cleaned_timeseries_info(cleaned_timeseries_info,output_file):
    """
    
    Save cleaned timeseries data and metadata as a pickled Python object.

    Parameters
    ----------
    cleaned_timeseries_info : list
        Cleaned timeeseries data and metadata (dataframe and dictionary respectively)
        from NWIS query.
        
    out_file : string
        Output file containing cleaned timeseries data and metadata.

    Returns
    -------
    None.

    """    
    
    # # create temporary folder if it doesn't already exist
    # output_dir = os.path.dirname(output_file)
    # if not os.path.exists(output_dir):
    #     os.makedirs(output_dir)

    # save to file via pickling
    with open(output_file, 'wb') as file:
        pickle.dump(cleaned_timeseries_info, file)

def main(in_file,out_file):
    """
    
    Main function.

    Parameters
    ----------
    in_file : string
        Pickled python file containing timeseries data and metadata.

    out_file : string
        Pickled python file containing cleaned timeseries data and metadata.

    Returns
    -------
    None.
    
    """
   
    # retrieve timeseries data and metadata
    timeser_info = load_timeseries_info(in_file)
    # create loop to update list of timeseries data and metadata
    N = len(timeser_info)
    for i in range(N):
        # extract timeseries data
        df = timeser_info[i]['data']
        
        # clean timeseries data and update dataframe with cleaned data
        df = remove_missing_data(df)
        timeser_info[i]['data'] = df
        
    # save timeseries data and metadata to file 
    save_cleaned_timeseries_info(timeser_info,out_file)


if __name__ == '__main__':
    
################################# INPUTS ######################################
    
    # search from 1_fetch folder
    input_file = snakemake.input['in_file']
    # 'cleaned_timeseries_info.pkl'
    output_file = snakemake.output['out_file']

###############################################################################

    main(input_file,output_file)

