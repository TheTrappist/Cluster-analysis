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
    """Extract data from a Cell Profiler HDF5 output file as a Pandas dataframe.
    Returns data from all available fields unless specific ones are specified
    
    Args:
        file (str): Full path to HDF5 file written by Cell Profiler.
        obj (str): Name of Cell Profiler object for which the HDF5 file 
            contains valid measurements.
        data_fields (list of str): Names of fields for which to get data. If left
            empty, function returns all available fields. Defaults to an empty tuple.
    Returns:
        df (pandas dataframe): all requested data from obj_name, as a pandas df.       
    """

    with h5py.File(file, 'r') as f:
        measurement_keys = list(f['Measurements'].keys())

        if len(measurement_keys) > 1:
            raise Exception("More than one dataset in file! Not sure which one to use.")
        else:
            key = measurement_keys[0]

        if not data_fields: # Specific fields not requested; grab all of them
            data_fields = list(f['Measurements'][key][obj])
            
        # Extract data and save as a dictionary
        data_dict= {} # dictionary to store parsed data
        #prev_imagenumber = [] # To check whether all images are listed in the same order
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
                print("Note: ", field, "has been replaced with NaNs due to a length mismatch")
            
            data_dict.update({field : data_trunc})
            
            """
            # This block performs faster but only works if all indexes share the same order
            
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
