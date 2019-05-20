# -*- coding: utf-8 -*-
"""

This module contains a collection of functions for extracting and processing
Cell Profiler outputs.

Written by Vladislav Belyy (UCSF)

"""
import numpy as np
import pandas as pd
import h5py


def get_data_cp_hdf5 (file, obj, data_fields=[]):
    """Extract data from a Cell Profiler HDF5 file as a Pandas dataframe.
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
            
            """
            # This block performs faster but only works if all indexes 
            # share the same order
            
            # calculate total number of data points
            num_data = np.sum(index[:,2] - index[:,1])
            # Check for consistent order of images in index
            print(field)
            if not np.any(prev_imagenumber):
                prev_imagenumber = index[:,0]
            elif not np.array_equal(prev_imagenumber, index[:,0]):
                raise Exception("Images not listed in the same order!")
                
            # Use index to parse data in cases when some points are not used
            if num_data != len(data):
                data_trunc = []
                for row in index:
                    data_trunc.extend(data[row[1]:row[2]])
                data = data_trunc
            data_dict.update({field : data})
            """       
        
        # Convert data to dataframe and return
        df = pd.DataFrame(data_dict)
    return df

def get_data_cp_csv (file, data_fields=None):
    """Extract data from a Cell Profiler CSV file as a Pandas dataframe.
    Returns data from all available fields unless specific ones are specified
    
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