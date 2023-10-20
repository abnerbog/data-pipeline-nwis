# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 11:53:20 2023

@author: AbnerBogan
"""

import pickle 
import numpy as np

def load_cleaned_timeseries_info(in_file):
    """

    Load in cleaned timeseries data and metadata.

    Parameters
    ----------
    in_file : string
        Pickled Python file containing cleaned timeseries data and metadata.

    Returns
    -------
    timeseries_info : list
        Loaded timeseries data (cleaned) and metadata.

    """
    
    # load file via unpickling
    with open(in_file, 'rb') as file:
        timeseries_info = pickle.load(file)
        
    return timeseries_info

def get_number_sites(df,data_column):
    """
    
    Get number of unique sites reporting timeseries data over time period

    Parameters
    ----------
    df : pandas dataframe
        Dataframe containing cleaned daily timeseries data.
    data_column : string
        Column header containing statistical value of interest (e.g. 00010_Mean).

    Returns
    -------
    number_sites : int
        Number of unique sites from timeseries data.

    """
    
    number_sites = df[data_column].nunique()
    
    return number_sites
    
def get_data_average(df,data_column):
    """
    
    Get average value of statistical value of interest (e.g. 00010_Mean) over 
    time period.

    Parameters
    ----------
    df : pandas dataframe
        Dataframe containing cleaned daily timeseries data.
    data_column : string
        Column header containing statistical value of interest (e.g. 00010_Mean).
        
    Returns
    -------
    average_temp : float
        Averaged value of statistical value of interest from timeseries data.

    """

    average_temp = np.average(df[data_column])
    
    return average_temp

def save_processed_timeseries_info(processed_timeseries_info,out_file):
    """
    
    Save processed timeseries data and metadata as a pickled Python object.

    Parameters
    ----------
    processed_timeseries_info : list
        Processed and cleaned timeeseries data and metadata (dataframe and 
        dictionary respectively) from NWIS query.
        
    out_file : string
        Output file containing processed and cleaned timeseries data and metadata.

    Returns
    -------
    None.

    """    
    
    # pickle (serialize) the data to a file
    with open(out_file, 'wb') as file:
        pickle.dump(processed_timeseries_info, file)

def main(in_file,out_file):
    """

    Parameters
    ----------
    in_file : TYPE
        DESCRIPTION.
    out_file : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    # retrieve cleaned timeseries data and metadata
    cleaned_timeseries_info = load_cleaned_timeseries_info(in_file)
    
    # create loop to update list of cleaned timeseries data and metadata
    N = len(cleaned_timeseries_info)
    for i in range(N):
        df = cleaned_timeseries_info[i]['data']
        
        # retrieve column header name containing stastical value of interest (e.g. 00010_Mean)
        if i == 0:
            data_column = df.columns[-2]
        
        # get unique number of sites and average value 
        number_sites = get_number_sites(df,data_column)
        average_val = get_data_average(df, data_column)
        # update dataframe and dictionary with processed value and additional metadata
        cleaned_timeseries_info[i]['data'] = average_val
        cleaned_timeseries_info[i]['metadata']['stat'] = data_column
        cleaned_timeseries_info[i]['metadata']['number sites'] = number_sites
        
    save_processed_timeseries_info(cleaned_timeseries_info,out_file)


if __name__ == '__main__':

################################# INPUTS ######################################

    # in_file = 'cleaned_timeseries_info.pkl'
    input_file = snakemake.input['in_file']
    
    # out_file = 'processed_timeseries_info.pkl
    output_file = snakemake.output['out_file']

###############################################################################

    main(input_file,output_file)

        

    
    