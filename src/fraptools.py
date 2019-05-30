# -*- coding: utf-8 -*-
"""

This module contains a collection of functions for analyzing FRAP data.

Written by Vladislav Belyy (UCSF)

"""
import numpy as np
import matplotlib.pylab as plt
from scipy.optimize import curve_fit
import json

def frap_fun(x, a, b):
    # FRAP recovery formula
    return a*(1-np.exp(-b*x))

def read_ij_intensity (file):
    """Read a file that was written by the save_background.py ImageJ script.
    
    Args:
        file (str): Full path to .json file written by save_background.py
        
    Returns:
        intensity_data (dict): All data stored in the data file     
    """
    with open(file) as data_file:    
        intensity_data = json.load(data_file)
    
    return intensity_data

def get_traces_from_df_list(df_list, file_names, exclude=()):
    """Extract singe traces from Extract_two_radii_TrackMate.ijm data.
    
    Args:
        df_list (list of pandas dataframes): List of parsed data from
            the output file of the "Extract_two_radii_TrackMate.ijm"
            ImageJ macro.
        file_names (list of str): List of file names that the traces
            are coming from
        exclude (list of tuples): Traces to exclude. Each item in the list
            is a tuple containing file name and track ID of the trace to
            be excluded. Defaults to an emptry tuple.
    
    Returns:
        result_df (list of pandas dataframes): One dataframe per trace.
            Contains all the same columns as df_list; the only difference
            is that this one is split into smaller single-trace chunks.
        corr_int (list of numpy arrays): list of one array per trace,
            containing just the background-corrected intensity values.
        trace_IDs (list of tuples): list of filenames corresponding 
            to each trace and trace ID. Same length as result_df.
        
    """
    
    result_df = []
    result_np = []
    trace_IDs = []
    
    for df, file_name in zip(df_list, file_names):
        df.sort_values(by=['track_IDs', 'frames'])
        
        # Find all the transitions between tracks by changes in track_ID
        diff = df['track_IDs'][1:].values - df['track_IDs'][:-1].values
        track_start_points = np.append([0],(np.flatnonzero(diff)+1))
        track_slices = np.append(track_start_points, [len(diff)+1])            
        
        # break up data into smaller data frames, one per trace
        for i in range(len(track_slices)-1):
            curr_df = df[track_slices[i]:track_slices[i+1]]
            bkgnd_corr_int = curr_df['mean_int_inner']-curr_df['bknd_int']
            curr_np = bkgnd_corr_int.values

            # save the trace if it's not excluded
            curr_trace_number = curr_df['track_IDs'].iloc[0]
            trace_ID = (file_name, curr_trace_number)
            
            if trace_ID not in exclude:
            
                result_df.append(curr_df)
                result_np.append(curr_np)   
                trace_IDs.append(trace_ID)
        
    return result_df, result_np, trace_IDs

def read_nonclust_frap_data(df_list, file_names, exclude=()):
    """Extract singe traces from Manual_FRAP_ROI.ijm data.
    
    Args:
        df_list (list of pandas dataframes): List of parsed data from
            the output file of the "Manual_FRAP_ROI.ijm"
            ImageJ macro.
        file_names (list of str): List of file names that the traces
            are coming from
        exclude (list of tuples): Traces to exclude. Each item in the list
            is a tuple containing file name and track ID of the trace to
            be excluded. Defaults to an empty tuple.
    
    Returns:
        result_df (list of pandas dataframes): One dataframe per trace.
            Contains two columns, "ROI_intensity" and "Bkgnd_intensity".
        corr_int (list of numpy arrays): list of one array per trace,
            containing just the background-corrected intensity values.
        trace_IDs (list of tuples): list of filenames corresponding 
            to each trace and trace ID. Same length as result_df.
        
    """
    
    result_df = []
    result_np = []
    trace_IDs = []
    
    for df, file_name in zip(df_list, file_names):
        
        # Find the number of FRAP/background trace pairs in file
        num_roi_pairs = int((len(df.columns) - 1) / 4)
        
        # break up data into smaller data frames, one per trace
        for i in range(num_roi_pairs):
            
            # there are four columns per trace
            curr_df = df.iloc[:, i*4+1:i*4+5]
            bkgnd_corr_int = curr_df.iloc[:,1] - curr_df.iloc[:,3]
            curr_np = bkgnd_corr_int.values - min(bkgnd_corr_int.values)

            # save the trace if it's not excluded
            trace_ID = (file_name, i)
            
            if trace_ID not in exclude:
            
                result_df.append(curr_df)
                result_np.append(curr_np)   
                trace_IDs.append(trace_ID)
        
    return result_df, result_np, trace_IDs

def correct_background (bkgnd_data, track_coords, track_intensities):
    """Subtract background intensity from track_coords.
    
    Args:
        bkgnd_data (dict): Background intensities, e.g. from read_ij_intensity
        track_coords (list): numpy arrays of x,y,t spot positions per track
        track_intensities (list): List of float numpy arrays, one per track, 
            containing the fluorescence intensity values of each tracked spot
    Returns:
        bkgnd_int (list): List of numpy arrays (one for each track in 
            track_intensities) storing background intensity values
        bkgnd_corr_int (list): List of numpy arrays (one for each track in 
            track_intensities) storing background-subtracted track intensities     
    """
    # Subtract background intensities from track intensities
    bkgnd_raw = bkgnd_data['intensities']
    bkgnd_int = [] # background intensities per track
    bkgnd_corr_int = [] # spot intensities corrected for background
    
    for i, track in enumerate(track_coords):
        start = int(track[0,2])
        n_frames = track.shape[0]
        bkgnd_temp = bkgnd_raw[start:(start+n_frames)]
        corrected_int = track_intensities[i]-bkgnd_temp
        bkgnd_int.append(bkgnd_temp)
        bkgnd_corr_int.append(corrected_int)
    
    return bkgnd_int, bkgnd_corr_int

def fit_frap(bkgnd_corr_int, frame_interval=1.0, bleach_n_frames=1):
    """Normalize and fit each fluorescence trace to a standard FRAP curve.
    
    Args:
        bkgnd_corr_int (list): List of numpy arrays (one for each track in 
            track_intensities) storing background-subtracted track intensities 
        frame_interval (float): time, in seconds, between subsequent frames.
        bleach_n_frames (int): number of frames where bleaching beam is 
            visible. Defaults to 1.

    Returns:
        fit_results (list): List of fit parameters and y-values for the fit 
            curve. One element for each input track.
        norm_data (list): Normalized x and y data used for fitting. One element
            for each input track.            
    """
    traces = [list(x) for x in bkgnd_corr_int] # convert to lists
    fit_results = []    
    norm_data = []
    for trace in traces:
        # First, find the bleach frame by looking for the minimum intensity
        lowest = min(trace)
        bleach_frame = trace.index(lowest)
        bleach_time = bleach_frame * frame_interval
        trace_norm = np.array(trace-lowest)
        
        # Then find mean intensity prior to bleaching
        last_prebleach_frame = bleach_frame - bleach_n_frames
        last_prebleach_frame = max([last_prebleach_frame, 1])
        y_prebleach = trace_norm[:last_prebleach_frame]
        int_prebleach = np.mean(y_prebleach)
        
        trace_norm = trace_norm / int_prebleach
        
        # get subset of data for fitting
        y_data = trace_norm[bleach_frame:]
        x_data = np.array(range(len(y_data))) * frame_interval
        x_whole_trace = np.array(range(len(trace_norm))) * frame_interval
        
        # set up fit bounds
        a_upper = max(y_data)*2
        if a_upper <= 0:
            a_upper = 0.001 # catches traces that are fully lower than zero
        b_upper = 4 / frame_interval
        fit_bounds = (0, [a_upper, b_upper])
        
        # fit recovery curve
        popt, pcov = curve_fit(frap_fun, x_data, y_data, bounds=fit_bounds)
        fit_result_y = frap_fun(x_data, *popt)
        thalf = np.log(2) / popt[1]
        mobile_fraction = popt[0]
        
        # package fit results
        fit_result_curr = {'popt':popt, 'pcov':pcov, 'thalf':thalf,
                           'mobile_fraction':mobile_fraction}
        
        # package normalized data:
        t_trunc = x_data + bleach_time
        norm_d = {'t_trunc':t_trunc, 't_trunc_0':x_data, 'y_trunc':y_data,
                  't_full':x_whole_trace, 'y_full':trace_norm, 
                  'fit_curve':fit_result_y}
        
        fit_results.append(fit_result_curr)
        norm_data.append(norm_d)
        
    return fit_results, norm_data

def fit_frap_smart(traces, frame_interval=1.0, max_bleach_frames=6, 
                   prebleach_frames=10, skip_frames=2, min_bleach_jump=0.2):
    """Fit each trace to a FRAP curve with smart normalization.
    
    This function finds the bleach point of the trace using a maximum-slope
    method and is more robust at identifying the bleach and pre-bleach points
    than fit_frap in complex or noisy traces.
    
    Args:
        traces (list): List of numpy arrays (one for each track in 
            track_intensities) storing background-subtracted track intensities
        frame_interval (float): time, in seconds, between subsequent frames.
            Defaults to 1.0.
        max_bleach_frames (int): how many frames can the bleaching pulse last?
            Defaults to 6.
        prebleach_frames (int): Desired number of frames to average to find
            the pre-bleach intensity level. Defaults to 10.
        skip_frames (int): number of frames to skip at the beginning of the 
            trace (the first few frames of each trace can be bleachy and 
            unreliable). Defaults to 2.
        min_bleach_jump (float): how large of a jump in a single frame, 
            expressed as a fraction of initial intensity, counts as the start
            of a bleaching frame? Defaults to 0.2.

    Returns:
        fit_results (list): List of fit parameters and y-values for the fit 
            curve. One element for each input trace.
        norm_data (list): Normalized x and y data used for fitting. One element
            for each input track.            
    """
    fit_results = []    
    norm_data = []

    for trace in traces:
        # Remove the skipped frames
        trace = trace[skip_frames:]
        
        # Find the frame where recovery begins.
        # First, locate the point with the largest negative slope.
        trace_pre_norm = trace / trace[0]
        slopes = trace_pre_norm[1:] - trace_pre_norm[:-1]
        min_slope = min(slopes)
        min_slope_index = np.where(slopes == min_slope)[0][0] + 1
        
        # Find the local minimum around the minimum slope point:
        try:
            local_start_index = min_slope_index - max_bleach_frames
            local_end_index = min_slope_index + max_bleach_frames
            trace_local = trace[local_start_index:local_end_index]
            min_index_local = np.argmin(trace_local)
            min_index = min_index_local + local_start_index
        except:
            print("Good local minimum not found; using minimum slope.")
            min_index = min_slope_index
        
        recovery_start_index = min_index
        
        # Find pre-bleach intensity level to normalize the trace.
        # Use the mean of the first prebleach_frames data points or however 
        # many frames occur before a large (>min_bleach_jump) jump in intensity
        abs_slope = np.absolute(slopes)
        prebleach_frames = min(prebleach_frames, len(trace))
        try:
            big_jump_index = np.where(abs_slope > min_bleach_jump)[0][0] + 1
            last_prebleach_frame = min(prebleach_frames, big_jump_index)
        except:
            last_prebleach_frame = prebleach_frames 
               
        # Normalize the trace such that the pre-bleach level is 1 and lowest
        # intensity after bleaching is 0.
        lowest = trace[recovery_start_index]
        bleach_frame = recovery_start_index
        bleach_time = bleach_frame * frame_interval
        trace_norm = np.array(trace-lowest)
        y_prebleach = trace_norm[:last_prebleach_frame]
        int_prebleach = np.mean(y_prebleach)
        
        trace_norm = trace_norm / int_prebleach
        
        # get subset of data for fitting
        y_data = trace_norm[bleach_frame:]
        x_data = np.array(range(len(y_data))) * frame_interval
        x_whole_trace = np.array(range(len(trace_norm))) * frame_interval
        
        # set up fit bounds
        a_upper = max(y_data)*2
        if a_upper <= 0:
            a_upper = 0.001 # catches traces that are fully lower than zero
        b_upper = 4 / frame_interval
        fit_bounds = (0, [a_upper, b_upper])
        
        # fit recovery curve
        popt, pcov = curve_fit(frap_fun, x_data, y_data, bounds=fit_bounds)
        fit_result_y = frap_fun(x_data, *popt)
        thalf = np.log(2) / popt[1]
        mobile_fraction = popt[0]
        
        # package fit results
        fit_result_curr = {'popt':popt, 'pcov':pcov, 'thalf':thalf,
                           'mobile_fraction':mobile_fraction}
        
        # package normalized data:
        t_trunc = x_data + bleach_time
        norm_d = {'t_trunc':t_trunc, 't_trunc_0':x_data, 'y_trunc':y_data,
                  't_full':x_whole_trace, 'y_full':trace_norm, 
                  'fit_curve':fit_result_y}
        
        fit_results.append(fit_result_curr)
        norm_data.append(norm_d)
        
    return fit_results, norm_data

def plot_fit_results(fit_result, data, rc_params=[4,2]):
    """Plot the results of a FRAP fit.
    
    Args:
        fit_result (dict): List of fit parameters and y-values for the fit 
            curve.
        data (dict): Normalized x and y data used in fitting.
        rc_params (list): matplotlib figsize setting. Defaults to [18,8].

    Returns:
        fig (matplotib figure) : Output figure
        axarr (matplotlib axes array): Axes of the output figure (2 elements)          
    """
    plt.rcParams["figure.figsize"] = rc_params
    fig, axarr = plt.subplots(1,2)
    fig.tight_layout(pad=2)

    axarr[0].set_title('Full trace')
    axarr[0].set_xlabel('Time (s)')
    axarr[0].set_ylabel('Normalized intensity')
    axarr[0].plot(data['t_full'], data['y_full'], 'b-', label='data')
    axarr[0].plot(data['t_trunc'], data['fit_curve'], 'r-', label='fit')
    axarr[0].legend(loc='upper right', shadow=True)

    axarr[1].set_title('Recovery and fit')
    axarr[1].set_xlabel('Time post-bleaching (s)')
    axarr[1].set_ylabel('Normalized intensity')
    axarr[1].plot(data['t_trunc_0'], data['y_trunc'], 'b-', label='data')
    axarr[1].plot(data['t_trunc_0'], data['fit_curve'], 'r-', 
                         label='fit')
    axarr[1].legend(loc='lower right', shadow=True)
    
    # Adjust the y-limits of axis 1
    ylims_fit = axarr[1].get_ylim()
    y_max = max(ylims_fit[1], 1.5)
    ylims = [ylims_fit[0], y_max]
    
    axarr[0].set_ylim(ylims)
    
    return fig, axarr

