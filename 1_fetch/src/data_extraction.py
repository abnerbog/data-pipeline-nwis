# -*- coding: utf-8 -*-
"""
Created on Sun Oct 15 16:12:54 2023

@author: AbnerBogan
"""

import pickle
import dataretrieval.nwis as nwis

def load_hu_metadata(in_file):
    """

    Load in hydrologic unit (HU) metadata.

    Parameters
    ----------
    in_file : string
        Pickled Python file containing HU metadata.

    Returns
    -------
    hu_metadata : list
        Loaded HU metadata.

    """
    
    # load file via unpickling
    with open(in_file, 'rb') as file:
        huc_metadata = pickle.load(file)
        
    return huc_metadata

def get_timeseries_data(hu_metadata,start,end,param,stat):
    """

    Parameters
    ----------
    hu_metadata : list
        HU metadata.
    start : string
        Start date for obtaining records of data (must be in format YYYY-MM-DD).
    end : string
        End date for obtaining records of data (must be in format YYYY-MM-DD).
    param : string
        USGS parameter code (5-digits) to identify the value measured and the 
        units of measure.
    stat : string
        USGS parameter code (5-digits) to report a specific statistic of a given
        dataset.

    Returns
    -------
    df : pandas dataframe
        Dataframe containing daily timeseries data (using get_dv function).

    """
    
    # extract huc
    huc = hu_metadata['huc']

    # get daily stats of a given parameter in a given time range 
    df,_ = nwis.get_dv(huc=huc,start=start,end=end,statCd=stat,parameterCd=param)
    
    return df

def save_timeseries_info(timeseries_info,output_file):
    """
    
    Save timeseries data and metadata as a pickled Python object.

    Parameters
    ----------
    timeseries_info : list
        Timeseries data and metadata (dataframe and dictionary respectively)
        from NWIS query.
        
    out_file : string
        Output file name containing timeseries data and metadata.

    Returns
    -------
    None.

    """

    # pickle (serialize) the data to a file
    with open(output_file, 'wb') as file:
        pickle.dump(timeseries_info, file)
  
def main(start,end,param,stat,in_file,out_file):
    """
    
    Main function.

    Parameters
    ----------
    start : string
        Start date for obtaining records of data (must be in format YYYY-MM-DD).
    end : string
        End date for obtaining records of data (must be in format YYYY-MM-DD).
    param : string
        USGS parameter code (5-digits) to identify the value measured and the 
        units of measure.
    stat : string
        USGS parameter code (5-digits) to report a specific statistic of a given
        dataset.
    out_file : string
        Output file containing timeseries data and metadata.

    Returns
    -------
    None.

    """
    
    # initialize list of timeseries data and metadata
    timeser_info = []
    # retrieve hu metadata
    hu_metadata = load_hu_metadata(in_file)
    
    # create loop to populate list of timeseries data and metadata
    N = len(hu_metadata)
    for i in range(N):
        # extract timeseries data as dataframe
        df = get_timeseries_data(hu_metadata[i],start,end,param,stat)
        # populate list with dictionary containing timeseries data and metadata
        timeser_info.append({'data':df,'metadata':hu_metadata[i]})
    
    # save timeseries data and metadata to file
    save_timeseries_info(timeser_info,out_file)
    

if __name__ == '__main__':

################################# INPUTS ###################################### 

    parameter_code = str(snakemake.params['param_code'])

    statistic_code = str(snakemake.params['stat_code'])

    start_time = str(snakemake.params['start'])

    end_time = str(snakemake.params['end'])
    
    input_file = snakemake.input['in_file']

    output_file = snakemake.output['out_file'] 


###############################################################################
  
    main(start_time,end_time,parameter_code,statistic_code,input_file,output_file)  
    
    

