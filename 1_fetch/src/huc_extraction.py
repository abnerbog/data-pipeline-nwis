# -*- coding: utf-8 -*-
"""
Created on Sun Oct 15 15:42:03 2023

@author: AbnerBogan
"""

import os
import pickle
import pandas as pd
from pygeohydro import WBD

def load_huc(csv_file):
    """
    
    Load the hydrologic unit codes (HUC) of interest from .csv file.
    
    Parameters
    ----------
    csv_file : string
        .csv file containing hydrologic unit code (HUC) information and huc 
        type (e.g., huc2, huc4)

    Returns
    -------
    huc_vals : list
        HUC values expressed as a list of strings.
        
    huc_type: string
        HUC type (e.g., huc2, huc4). Note this value is taken from the column
        header in csv_file containing HUC values.

    """

    # read csv file containing huc values
    df = pd.read_csv(csv_file,dtype=str)
    # extract column name describing huc type
    huc_type = df.columns[0]
    # read in list of huc values
    huc = df[huc_type].tolist()
        
    return huc,huc_type


def get_hu_metadata(code,code_type):
    """

    Get hydrologic unit (HU) metadata.

    Parameters
    ----------
    code : string
        Individual HUC value.
    code_type : string
        HUC type (e.g., huc2, huc4).

    Returns
    -------
    name : string
        Name of hydrologic unit (e.g., Mid Atlantic Region)
    boundary : multipolygon object
        Polygon coordinates to identify official delineations for a given 
        hydrologic unit (sourced from Watershed Boundary Dataset).
    """
    
    # get hydrologic unit (hu) metadata 
    wbd = WBD(code_type)
    hu = wbd.byids(code_type,code)
    # extract name of hydrologic unit from geodataframe
    name = hu['name'][0] 
    # extract boundary multipolygon object from geodataframe
    boundary = hu['geometry'][0]


    return name,boundary

def save_hu_metadata(metadata,output_file):
    """
    
    Save hydrologic unit metadata as a pickled Python object.
    
    Parameters
    ----------
    metadata : list
        Metadata (name, huc and boundary) associated with a given hydrologic 
        unit.
        
    out_file : string
        Output file containing HU metadata.

    Returns
    -------
    None.

    """
    
    # save to file via pickling
    with open(output_file, 'wb') as file:
        pickle.dump(metadata,file)
        
    
def main(in_file,out_file):
    """
    
    Main function.

    Parameters
    ----------
    in_file : string
        .csv file containing hydrologic unit code (HUC) information and huc 
        type (e.g., huc2, huc4)

    out_file : string
        Output file containing hydrologic unit (HU) metadata

    Returns
    -------
    None.

    """
    
    # retrieve the hydrologic unit codes (huc) from input file
    huc,huc_type = load_huc(in_file)

    # initialize list of hu metadata
    hu_metadata = []
    
    # create loop to populate list of hu metadata
    N = len(huc)
    for i in range(N): 
        # retrieve the name and boundary of a given hu
        name, boundary = get_hu_metadata(huc[i],huc_type)
        # populate list with dictionary containing hu metadata
        hu_metadata.append({'name':name,'huc':huc[i],'polygon':boundary})   
      
    # save hu metadata to file
    save_hu_metadata(hu_metadata,out_file)

    
if __name__ == '__main__':
    
################################# INPUTS ######################################
    
    
    # 'huc_inputs.csv'
    huc_file = snakemake.input['csv_file']
    
    # 'hu_metadata.pkl'
    hu_metadata_file = snakemake.output['out_file'] 
    
    # print(hu_metadata_file)

###############################################################################

    main(huc_file,hu_metadata_file)