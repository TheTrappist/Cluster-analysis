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

def add_image_prop_to_objects (objects, images, prop, fast=True):
    """Add a column for an image-specific property to the 'objects' dataframe.
    
    This function is useful for plotting the distribution of individual cell 
    measurements against properties that are common to the entire image, such
    as timepoint, type of treatment, etc. The 'objects' dataframe is appended
    with a new column that contains the value of the specified property 'prop'
    for each cell. 'prop' must be a valid field in the 'images' dataframe.
    
    Args:
        objects (pandas dataframe): Collection of cell-by-cell properties. Must 
            contain the 'ImageNumber' column.
        images (pandas dataframe): Collection of image properties. Must contain
            an 'ImageNumber' column as well as a column matching the argument
            'prop'
        prop (str): Name of the image-wide property to be added to the 
            'objects' dataframe. Must be a valid column name in 'images'.
        fast (bool): If True, uses a much faster vectorized algorithm. Defaults
            to True.

            
    Returns:
        None (the dataframe 'objects' gets modified in place).      
    """
    
    if prop in objects.columns:
        print('Warning: column named '+prop+' already exists in objects!')
   
    if fast: # Optimized method
        # Get image numbers for rows of objects and images
        obj_ids = objects['ImageNumber'].values
        img_ids = images['ImageNumber'].values
        
         # Find corresponding row in images for each row in objects
        index = np.argsort(img_ids)
        sorted_img = img_ids[index]
        sorted_index = np.searchsorted(sorted_img, obj_ids)
        img_ind = index[sorted_index] # array of image index for each object
        
        # Pull desired image property for each object
        objects[prop] = images.loc[img_ind, prop].values
    
    else: # old version of the method; slow
        property_values = []    
        for index, obj in objects.iterrows():
            curr_value = images.loc[images['ImageNumber'] == obj['ImageNumber'], 
                                    prop].item()
            property_values.append(curr_value)
        objects[prop] = property_values
    
    return None

def add_parent_prop (children, parents, prop, rel_col, result_name):
    """Add a column for a parent-derived property in the "children" dataframe.
    
    This function is much faster than add_parent_prop_thorough by vectorizing
    the sort and search operations. Each child must have exactly one parent,
    and both 'ImageNumber' and 'ObjectNumber' fields must all be nonnegative
    integers.
    
    Args:
        children (pandas dataframe): Collection of child object properties.
        parents (pandas dataframe): Collection of parent object properties.
        prop (str): Name of the property to be added to the 'children'
            dataframe. Must be a valid column name in 'parents'.
        rel_col (str): Name of the column in 'children' that contains the
            object ID of the parent.
        result_name (str): Name of the new column in 'children' holding the 
            parent-derived property
            
    Returns:
        None (the dataframe 'children' gets modified in place).      
    """
      
    # Get image and object numbers for parents and children
    p_ids = parents[['ImageNumber', 'ObjectNumber']].values
    c_ids = children[['ImageNumber', rel_col]].values
        
    # Determine minimum required multiplier for combining image and object
    # IDs into a single identifier
    all_ids = np.concatenate((p_ids, c_ids))
    objnum_mult = np.max(all_ids[:,1]) + 1
    
    # Assign new single-number IDs to parents and children. These are unique
    # values for parents, but multiple children all share their parent's ID
    p_uids = p_ids[:,0] * objnum_mult + p_ids[:,1] # parent unique IDs
    c_uids = c_ids[:,0] * objnum_mult + c_ids[:,1] # matching child IDs
    
    # Use the new single IDs indices to locate parent for each child
    index = np.argsort(p_uids)
    sorted_p = p_uids[index]
    sorted_index = np.searchsorted(sorted_p, c_uids)
    p_ind = index[sorted_index] # array containing parent index for each child
    
    # Pull desired parent property for each child
    children[result_name] = parents.loc[p_ind, prop].values
    
    return None

def add_parent_prop_thorough (children, parents, prop, rel_col, result_name):
    """Add a column for a parent-derived property in the "children" dataframe.
    
    This is a much slower but more careful version of this function that relies
    on the inefficient .iterrows() method. Consider using add_parent_prop 
    instead.
    
    Args:
        children (pandas dataframe): Collection of child object properties.
        parents (pandas dataframe): Collection of parent object properties.
        prop (str): Name of the property to be added to the 'children'
            dataframe. Must be a valid column name in 'parents'.
        rel_col (str): Name of the column in 'children' that contains the
            object ID of the parent.
        result_name (str): Name of the new column in 'children' holding the 
            parent-derived property
            
    Returns:
        None (the dataframe 'children' gets modified in place).      
    """
    
    if result_name in children.columns:
        print('Warning: column named '+result_name+' already exists!')
    
    property_values = []
    nan_flag = False
    
    for index, child in children.iterrows():
        parents_in_image = parents['ImageNumber'] == child['ImageNumber']
        parents_with_obj_num = child[rel_col] == parents['ObjectNumber']
        parent = parents_in_image & parents_with_obj_num
        try:
            curr_value = parents.loc[parent, prop].item()
        except:
            curr_value = np.nan
            nan_flag = True

        property_values.append(curr_value)
    
    children[result_name] = property_values
    
    if nan_flag:
        print("WARNING: Some children did not have valid parents!")
    
    return None

def add_child_prop_to_parents (parents, children, prop, rel_col, result_name,
                             statistic='mean'):
    """Add a column for a statistic of child objects to the parents dataframe.
    
    Finds all child objects for each object in parents, calculates the desired
    statistic of the chosen property ('prop') for that parent, and adds the 
    statistic to the parents dataframe as a new column. A use case would be to
    calculate the mean intensity of all clusters in each cell.
    
    Args:
        parents (pandas dataframe): Collection of parent properties.
        children (pandas dataframe): Collection of properties of the child
            objects. Must contain a column with a name that matches 'prop'.
        prop (str): Name of the property to be added to the 'parents'
            dataframe.
        rel_col (str): Name of the column in 'children' that contains the
            object ID of the parent.
        result_name (str): Name of the new 'parents' column holding the 
            calculated statistic.
        statistic (str): type of statistic to calculate for the chilren. Valid
            values are 'mean', 'median', and 'sum'. Defaults to 'mean'.

    Returns:
        None (the dataframe 'parents' gets modified in place).      
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
    for index, parent in parents.iterrows():
        
        # Find indices of children that match both ImageNumber and ObjectNumber
        # of the current parent
        children_in_image = children['ImageNumber'] == parent['ImageNumber']
        children_with_parent_id = children[rel_col] == parent['ObjectNumber']
        children_in_parent = children_in_image & children_with_parent_id
        
        # Extract the desired property from children of the current parent
        prop_values = children.loc[children_in_parent, prop]
        
        if prop_values.empty:
            curr_value = np.nan     
        else:
            curr_value = stat_fun[statistic](prop_values)
        
        property_values.append(curr_value)
        
    parents[result_name] = property_values
    
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

def norm_clust_time_by_track (cells, col_num_clust, time_step, min_clust=5):
    """Find valid single-cell cluster timecourses and normalize them by time.
    
    Goes through each unique trajectory in 'cells' (specified by the required 
    column 'Track_and_group', which serves as a unique identifier of a tracking
    trajectory), finds the first frame where clusters appear in the cell and
    the last frame before they disappear, then normalizes the trajectory time
    such that appearance of clusters corresponds to t=0 and disappearance of
    clusters corresponds to t=1.0. Returns the valid cells as a new dataframe
    with the additional columns 'Time_Norm', which is the normalized time as
    defined above, and 'Time_Aligned_hrs', which is a timestamp in hours where
    all trajectories are aligned by the starting point only without re-scaling.
    
    Args:
        cells (pandas dataframe): Collection of cell-by-cell properties.
        col_num_clust (str): Name of the column in 'cells' that stores the
            number of clusters er cell.
        time_step (float): Duration of one frame, in minutes
        min_clust (int): Minimum number of clusters that need to be formed in
            the cell at some point in the timecourse for it to be considered
            a valid trajectory.
            
    Returns:
         cells_filt (pandas dataframe): A truncated copy of 'cells' with all
             cells that do not fit in a valid trajectory removed. Contains the 
             two additional columns 'Time_Norm' and 'Time_Aligned_hrs.
    """
    
    trajectories = []
    cells_filt = cells.copy()
    cells_filt['Time_Norm'] = np.nan # will be used to store normalized values
    
    for traj_id in cells_filt['Track_and_group'].unique():
        traj = cells_filt.loc[cells_filt['Track_and_group'] == traj_id, 
                              col_num_clust]
        
        # Only keep trajectories that begin and end with zero clusters
        # but have min_clust clusters in the midle
        if traj.iloc[0] > 0 or traj.iloc[-1] > 0 or max(traj) < min_clust:
            continue
        
        #Find frame-to-frame differences in number of clusters
        deltas = (traj.iloc[1:].values - traj.iloc[:-1].values).astype(bool)
        
        # Normalize time scaling
        first_clust_frame = np.argmax(deltas)
        last_clust_frame = len(deltas) - np.argmax(np.flip(deltas))
        frame_interval = 1.0 / (last_clust_frame - first_clust_frame) 
        
        # Add rescaled time to cells_filt
        cells_in_traj = cells_filt['Track_and_group'] == traj_id
        n_time = cells_filt.loc[cells_in_traj, 'Metadata_Frame'].values
        offset_time = n_time - first_clust_frame
        n_time = offset_time * frame_interval
        cells_filt.loc[cells_in_traj, 
                       'Time_Aligned_hrs'] = offset_time * time_step / 60
        cells_filt.loc[cells_in_traj, 'Time_Norm'] = n_time
        trajectories.append(traj_id)
    
    # remove columns that don't match valid trajectories from cells_filt.
    rejected_cells = ~cells_filt['Track_and_group'].isin(set(trajectories))
    cells_filt.drop(cells_filt[rejected_cells].index, inplace=True)
    
    return cells_filt