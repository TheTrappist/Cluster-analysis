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

def plot_fit_results(fit_result, data, rc_params=[18,8]):
    """Plot the results of a FRAP fit.
    
    Args:
        fit_result (dict): List of fit parameters and y-values for the fit 
            curve.
        data (dict): Normalized x and y data used in fitting.
        rc_params (list): matplotlib figsize setting. Defaults to [18,8].

    Returns:
        fig (matplotib figue) : Output figure
        axarr (matplotlib axes array): Axes of the output figure (2 elements)          
    """
    plt.rcParams["figure.figsize"] = rc_params
    fig, axarr = plt.subplots(1,2)

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
    
    return fig, axarr

