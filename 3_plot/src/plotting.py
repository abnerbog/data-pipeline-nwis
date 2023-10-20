# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 16:39:49 2023

@author: AbnerBogan
"""

import pickle
import folium
import pandas as pd
import geopandas as gpd
from branca.colormap import LinearColormap

def load_processed_timeseries_info(in_file):
    """

    Load in processed timeseries data and metadata.

    Parameters
    ----------
    in_file : string
        Pickled Python file containing processed and cleaned timeseries data and metadata.

    Returns
    -------
    timeseries_info : list
        Loaded timeseries data (processed and cleaned) and metadata.

    """
    
    # load file via unpickling
    with open(in_file, 'rb') as file:
        timeseries_info = pickle.load(file)
        
    return timeseries_info

def parse_timeseries_info(timeseries_info):
    """
    
    Parse dictionary of timeseries data and metadata.

    Parameters
    ----------
    timeseries_info : list
        Processed and cleaned timeeseries data and metadata (dataframe and 
        dictionary respectively) from NWIS query.
        
    Returns
    -------
    polygon_list : list
        List of multipolygon objects for given hydrologic units.
    value_list : TYPE
        List of statistical value of for given hydrologic units and associated query.
    name_list : TYPE
        List of region names for given hydrologic units. 
    param_stat_list : TYPE
        List of descriptions of statistical value of interest.
    site_list : TYPE
        List of unique sites for given hydrologic units and associated query.

    """

    # initialize lists containing data and metadata
    polygon_list = []
    value_list = []
    name_list = []
    param_stat_list = []
    site_list = []
    
    # create loop to populate the list
    N = len(timeseries_info)
    for i in range(N):
        # boundary multipolygon object from geodataframe
        polygon_list.append(timeseries_info[i]['metadata']['polygon'])
        # statistical value of interest (rounded to 1 decimal point)
        value_list.append(round(timeseries_info[i]['data'],1))
        # name of hydrologic unit (e.g. Mid Atlantic)
        name_list.append(timeseries_info[i]['metadata']['name'][:-7])
        # description of statistical value of interest (e.g. 00010_Mean) 
        param_stat_list.append(timeseries_info[i]['metadata']['stat'])
        # number of unique sites in given query
        site_list.append(timeseries_info[i]['metadata']['number sites']) 
        
    return polygon_list, value_list, name_list, param_stat_list, site_list

def create_geodataframe(polygon_list, value_list, name_list, param_stat_list, site_list):
    """

    Parameters
    ----------
    polygon_list : list
        List of multipolygon objects for given hydrologic units.
    value_list : TYPE
        List of statistical value of for given hydrologic units and associated query.
    name_list : TYPE
        List of region names for given hydrologic units. 
    param_stat_list : TYPE
        List of descriptions of statistical value of interest.
    site_list : TYPE
        List of unique sites for given hydrologic units and associated query.

    Returns
    -------
    geodf : geo dataframe
        Geodataframe containing all data and metadata for a given geographical region.

    """
    
    # create geodataframe from geoseries and the data
    data = {'region': name_list, 'value': value_list, 'value description':param_stat_list,'number of sites':site_list}

    # create a geoseries from polygons
    geometry_series = gpd.GeoSeries(polygon_list)
    # update geodataframe with geoseries
    geometry_df = gpd.GeoDataFrame(data, geometry=geometry_series, crs='EPSG:4326')

    return geometry_df      

def get_style(feature, colormap):
    """

    Parameters
    ----------
    feature : dict
        GeoJSON feature containing properties for style and other features.
    colormap : LinearColormap
        LinearColormap used to map data values to colors.

    Returns
    -------
    dict_value : dict
        Dictionary specifying the style for the GeoJSON feature, including fillColor, color, and weight.

    """
    
    data_value = feature['properties']['value']
    # Check if the data_value is null or empty
    if data_value is None or pd.isna(data_value):
        # if value is empty, set a default color (e.g., gray)
        return {
            'fillColor': 'gray',
            'color': 'black',
            'weight': 2
        }
    else:
        # if the value is valid, set the fill color based on the data value
        return {
            'fillColor': colormap(data_value),
            'color': 'black',
            'weight': 2
        }

def create_plot(in_file,out_file):
    """
    
    Create a choropleth plot.

    Parameters
    ----------
    in_file : string
        Pickled python file containing cleaned and processed timeseries data and metadata.
    out_file : string
        Choropleth plot name.

    Returns
    -------
    None.

    """
    # Retrieve processed timeseries data and metadata
    timeser_info = load_processed_timeseries_info(in_file)

    # Parse processed timeseries data and metadata
    polygons, values, names, param_stats, sites = parse_timeseries_info(timeser_info)

    # Create geo dataframe
    geo_df = create_geodataframe(polygons, values, names, param_stats, sites)

    # Create folium map object (location starts around US)
    m = folium.Map(location=[39.8283, -98.5795], zoom_start=4)

    # Define a colormap to map data values to colors
    colormap = LinearColormap(
        colors=['green', 'yellow', 'red'],  # Define the color range
        vmin=geo_df['value'].min(), # Minimum data value
        vmax=geo_df['value'].max()  # Maximum data value
    )

    # Add the geodataframe as geojson objects with the modified style function
    folium.GeoJson(
        geo_df,
        name='Choropleth',
        style_function=lambda feature, colormap=colormap: get_style(feature, colormap),  # Use the separate style function
        highlight_function=lambda x: {'weight': 3, 'color': 'red'},
        tooltip=folium.features.GeoJsonTooltip(fields=['region', 'value', 'value description', 'number of sites'], labels=True, sticky=True)
    ).add_to(m)

    # Add a legend for the choropleth
    colormap.add_to(m)

    # save plot
    m.save(out_file)


if __name__ == '__main__':
    
################################# INPUTS ######################################

    # in_file = 'processed_timeseries_info.pkl
    input_file = snakemake.input['in_file']
    
    # out_file = 'choropleth_plot.html
    output_file = snakemake.output['out_file']

###############################################################################

    create_plot(input_file,output_file)