# -*- coding: utf-8 -*-
"""

This module contains a collection of functions for extracting and processing
CellProfiler outputs.

Written by Vladislav Belyy (UCSF)

"""
import numpy as np
import pandas as pd
import h5py


def get_data_cp_hdf5 (file, obj, data_fields=[]):
    """Extract data from a CellProfiler HDF5 file as a Pandas dataframe.
    Returns data from all available fields unless specific ones are specified
    
    Args:
        file (str): Full path to HDF5 file written by Cell Profiler.
        obj (str): Name of Cell Profiler object for which the HDF5 file 
            contains valid measurements.
        data_fields (list of str): Names of fields for which to get data. If 
            left empty, function returns all available fields. Defaults to an 
            empty tuple.
    Returns:
        df (Pandas dataframe): all requested data from obj_name.
    """

    with h5py.File(file, 'r') as f:
        measurement_keys = list(f['Measurements'].keys())

        if len(measurement_keys) > 1:
            raise Exception("More than one dataset in file!"
                            "Not sure which one to use.")
        else:
            key = measurement_keys[0]

        if not data_fields: # Specific fields not requested; grab all of them
            data_fields = list(f['Measurements'][key][obj])
            
        # Extract data and save as a dictionary
        data_dict= {} # dictionary to store parsed data
        #prev_imagenumber = [] # Check if images are listed in the same order
        prev_data_length = None
        for field in data_fields:
            data = f['Measurements'][key][obj][field]['data']
            index = np.array(f['Measurements'][key][obj][field]['index'])
            
            # Build data array based on sorted indices
            index.sort(axis=0)
            data_trunc = []
            for row in index:
                data_trunc.extend(data[row[1]:row[2]])
            
            # Check whether data length matches that of the previous field
            data_length = np.sum(index[:,2] - index[:,1])
            if not prev_data_length: prev_data_length = data_length
            elif prev_data_length != data_length:
                data_trunc = np.full(prev_data_length, np.nan)
                print("Note: ", field, "has been replaced with NaNs due to a"
                      "length mismatch")
            
            data_dict.update({field : data_trunc})
          
          # Convert data to dataframe and return
        df = pd.DataFrame(data_dict)
    return df

def get_data_cp_csv (file, data_fields=None):
    """Extract data from a CellProfiler CSV file as a Pandas dataframe.
    Returns data from all available fields unless specific ones are specified.
    
    Args:
        file (str): Full path to HDF5 file written by Cell Profiler.
        data_fields (list of str): Names of fields for which to get data. If 
            left empty, function returns all available fields. Defaults to
            None.
    Returns:
        df (Pandas dataframe): all requested data from obj_name.
    """


    df = pd.read_csv(file, usecols=data_fields)
    return df


def get_relationship_cp_hdf5 (file, obj1, obj2, modules=[]):
    """Extract relationship information from a Cell Profiler HDF5 output file 
    as a Pandas dataframe.
    
    Args:
        file (str): Full path to HDF5 file written by Cell Profiler.
        obj1 (str): Name of "parent" Cell Profiler object.
        obj2 (str): Name of "child" Cell Profiler object.
        modules (list of str): List of modules (sub-folders of "Relationship")
            that this function is allowed to search through as it looks for
            obj1 and obj2. If left empty, function looks through all available 
            modules and uses the first one that contains obj1 and obj2. 
            Defaults to an empty list.
            
    Returns:
        df (pandas dataframe): relationship data, as a Pandas df.       
    """
    
    
    with h5py.File(file, 'r') as f:
        measurement_keys = list(f['Measurements'].keys())

        if len(measurement_keys) > 1:
            raise Exception("More than one dataset in file!"
                            "Not sure which one to use.")
        else:
            key = measurement_keys[0]
            rel_path = '/Measurements/'+key+'/Relationship/'

        if not modules: # Module not specified; try all of them
            modules = list(f[rel_path])
        
        # Look for expected parent-child relationship in allowed modules
        module_found = False
        for mod in modules:
           full_path = rel_path+mod+'/Parent/'+obj1+'/'+obj2

           if full_path in f:
               module_found = True
               break
        if not module_found:
            raise Exception("Relationship pair not found!")
                   
        # Extract data and save as a dictionary
        datasets = list(f[full_path])
        
        data_dict= {} # dictionary to store parsed data
        for d in datasets:
            data = f[full_path + '/' + d]
            data_dict.update({d : data})
                
        # Convert data to dataframe and return
        df = pd.DataFrame(data_dict)
        
    return df

def add_image_prop_to_cells (cells, images, prop):
    """Add a column for an image-specific property to the 'cells' dataframe.
    
    This function is useful for plotting the distribution of individual cell 
    measurements against properties that are common to the entire image, such
    as timepoint, type of treatment, etc. The 'cells' dataframe is appended
    with a new column that contains the value of the specified property 'prop'
    for each cell. 'prop' must be a valid field in the 'images' dataframe.
    
    Args:
        cells (pandas dataframe): Collection of cell-by-cell properties. Must 
            contain the 'ImageNumber' column.
        images (pandas dataframe): Collection of image properties. Must contain
            an 'ImageNumber' column as well as a column matching the argument
            'prop'
        prop (str): Name of the image-wide property to be added to the 'cells'
            dataframe. Must be a valid column name in 'images'.

            
    Returns:
        None (the dataframe 'cells' gets modified in place).      
    """
    
    if prop in cells.columns:
        raise Warning('Column named '+prop+' already exists in cells!')
    
    property_values = []
    
    for index, cell in cells.iterrows():
        curr_value = images.loc[images['ImageNumber'] == cell['ImageNumber'], 
                                prop].item()

        property_values.append(curr_value)
    
    cells[prop] = property_values
    
    return None

def add_child_prop_to_cells (cells, children, prop, rel_col, result_name,
                             statistic='mean'):
    """Add a column for a statistic of child objects to the 'cells' dataframe.
    
    Finds all child objects for each cell in 'cells', calculates the desired
    statistic of the chosen property ('prop') for that cell, and adds the 
    statistic to the 'cells' dataframe as a new column. A use case would be to
    calculate the mean intensity of all clusters in each cell.
    
    Args:
        cells (pandas dataframe): Collection of cell-by-cell properties.
        children (pandas dataframe): Collection of properties of the child
            objects. Must contain a column with a name that matches 'prop'.
        prop (str): Name of the property to be added to the 'cells'
            dataframe.
        rel_col (str): Name of the column in 'children' that contains the
            object ID of the parent.
        result_name (str): Name of the new 'cells' column holding the 
            calculated statistic.
        statistic (str): type of statistic to calculate for the chilren. Valid
            values are 'mean', 'median', and 'sum'. Defaults to 'mean'.

    Returns:
        None (the dataframe 'cells' gets modified in place).      
    """
      
    
    # Define available statistic functions
    def mean(x):
        return np.mean(x)
    def median(x):
        return np.median(x)
    def sum_x(x):
        return np.sum(x)
    stat_fun = {'mean': mean,
                'median' : median,
                'sum': sum_x}
    
    property_values = []
    for index, cell in cells.iterrows():
        
        # Find indices of children that match both ImageNumber and ObjectNumber
        # of the current cell
        children_in_image = children['ImageNumber'] == cell['ImageNumber']
        children_with_cell_id = children[rel_col] == cell['ObjectNumber']
        children_in_cell = children_in_image & children_with_cell_id
        
        # Extract the desired property from children of the current cell
        prop_values = children.loc[children_in_cell, prop]
        
        if prop_values.empty:
            curr_value = np.nan     
        else:
            curr_value = stat_fun[statistic](prop_values)
        
        property_values.append(curr_value)
        
    cells[result_name] = property_values
    
    return None

def bootstrap_cell_prop (cells, measurement, group, nreps=1000):
    """Calculate bootstrap sample means for a measurement column in 'cells'.
    
    Breaks up cells into groups based on the value in the 'group' column of the
    'cells' dataframe and individually performs bootstrap resampling for each
    group. Groups typically comprise a large number of cells (e.g. all cells in
    one image or all cells treated with the same condition), and measurement is
    something that was individually measured for each cell (e.g. intensity,
    number of clusters, etc.). Both 'group' and 'measurement' must be valid
    names of columns in 'cells'.
    
    Args:
        cells (pandas dataframe): Collection of cell-by-cell properties.
        images (pandas dataframe): Collection of image properties. Must contain
            an 'ImageNumber' column as well as a column matching the argument
            'prop'
        measurement (str): Name of the cell-specific measurement to be 
            bootstrapped.
        group (str): Property that unites many cells into sensible categories.
        nreps (int): Number of bootstrap iterations. Defaults to 1000.
            
    Returns:
         df (pandas dataframe): bootstrap results per group (each group is a
            separate column in the dataframe). The individual values in each
            column are mean values of 'measurement' from individual resampling
            runs.
    """
    result = {}
    groups_unique = cells[group].unique()
    
    for grp in groups_unique:
        
        cells_in_group = cells.loc[cells[group] == grp]
        total_cells = len(cells_in_group)
        
        # Bootstrap the samples to estimate uncertainties
        metric = []
        for i in range(nreps):
            subsamp_cells = np.random.choice(list(cells_in_group[measurement]),
                                             total_cells)
            metric_mean_curr_rep = np.mean(subsamp_cells)
            metric.append(metric_mean_curr_rep)
        
        result.update({grp : metric})
    
    df = pd.DataFrame(data=result)
    
    return df