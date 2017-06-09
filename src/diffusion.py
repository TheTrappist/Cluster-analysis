# -*- coding: utf-8 -*-
"""

This module contains a collection of functions for analyzing diffusion data.

Written by Vladislav Belyy (UCSF)

"""
import os
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib.pylab as plt
#from scipy.optimize import curve_fit


def ax_fun(x, a):
    # Line with zero y-intercept
    return a*x

def calc_msd (track_coords, time_step):
    """Calculate and return mean squared displacement from provided tracks.
    
    Args:
        track_coords (list): numpy arrays of x,y,t spot positions, one array for
           each track.
        time_step (float): duration, in seconds, of each frame of the movie.
       
    Returns:
        t_dsq (list): numpy arrays of time vs. squared displacement, by track.
        msd_data (tuple): Contains the following data.
            time (numpy array): time vector for MSD plot
            msd (numpy array): mean squared displacements
            std (numpy array): standard deviations of MSD
            sterr (numpy array): standard errors of MSD 
    """
    t_dsq = list() # time vs sq. displacement for each particle
    for track in track_coords:
        # 'track' is an array of x,y,t

        track_norm = track - track[0,:] # set x,y,t=0
        xy_sq = np.square(track_norm[:,0:2])
        d_sq = np.sum(xy_sq, axis=1)
        t = track_norm[:,2] * time_step
        tds = np.vstack((t, d_sq))
        tds = np.transpose(tds)

        t_dsq.append(tds)

    # rebuild results as a padded numpy array for calculating means
    l_max = max((len(el) for el in t_dsq))

    d_sq = np.empty((l_max, len(t_dsq)))
    t_total = np.empty((l_max, len(t_dsq)))
    for i, x in enumerate(t_dsq):
        len1 = x.shape[0]
        pad_len = l_max - len1
        nanpad = np.zeros([pad_len,2])*np.nan
        x_padded = np.concatenate((x,nanpad), axis=0)

        d_sq[:,i] = x_padded[:,1]
        t_total[:,i] = x_padded[:,0]

    time = np.nanmean(t_total, axis=1)
    msd = np.nanmean(d_sq, axis=1)
    std_dsq = np.nanstd(d_sq, axis=1)
    sterr_dsq = stats.sem(d_sq, axis=1, nan_policy='omit',ddof=0)
    
    msd_data = (time, msd, std_dsq, sterr_dsq)
    
    return t_dsq, msd_data

def fit_diffusion_const(msd_data, dim=2, nframes=10):
    """Fit MSD data to a simple linear regression.
    
    Args:
        msd_data (tuple): Contains the following data.
            time (numpy array): time vector for MSD plot
            msd (numpy array): mean squared displacements
            std (numpy array): standard deviations of MSD
            sterr (numpy array): standard errors of MSD
        dim (int): dimensionality of the diffusion. Defaults to 2.
        nframes (int): number of frames from the start of the movie to be used
            for fitting. Defaults to 10.
    
    Returns:
        fit_params (dict): parameters of the fit.        
    """
    time, mean_dsq, __, __ = msd_data
    t_temp = time[0:nframes]
    msd_temp = mean_dsq[0:nframes]
    # get the slope:
    popt, pcov = curve_fit(ax_fun, t_temp, msd_temp)
    slope = popt[0]
    stdev = np.sqrt(np.diag(pcov))
    #slope = np.dot(t_temp, msd_temp) / np.dot(t_temp, t_temp)
    dc = slope / (2 * dim) # diffusion constant
    fit_params = {}
    fit_params['slope'] = slope
    fit_params['dc'] = dc
    fit_params['stdev'] = stdev
    
    return fit_params

def plot_msd(t_dsq, msd_data, fit_params=None, rc_params=[18,8]):
    """Plot the results of an MSD calculation and individual trajectories.
    
    Args:
        t_dsq (list): numpy arrays of time vs. squared displacement, by track.
        msd_data (tuple): Contains the following data.
            time (numpy array): time vector for MSD plot
            msd (numpy array): mean squared displacements
            std (numpy array): standard deviations of MSD
            sterr (numpy array): standard errors of MSD
        fit_params (dict): parameters of fit to be plotted, if desired
        rc_params (list): matplotlib figsize setting. Defaults to [18,8].

    Returns:
        fig (matplotib figue) : Output figure
        axarr (matplotlib axes array): Axes of the output figure (2 elements)          
    """
    plt.rcParams["figure.figsize"] = rc_params
    
    time, mean_dsq, std_dsq, sterr_dsq = msd_data
    
    fig, axarr = plt.subplots(1,2)
    axarr[0].set_title('All tracks')
    axarr[0].set_xlabel('Time (s)')
    axarr[0].set_ylabel(r'Squared displacement ($\mu m^2$)')
    axarr[1].set_title('Mean squared displacements')
    axarr[1].set_xlabel('Time (s)')
    axarr[1].set_ylabel(r'MSD ($\mu m^2$)')
    
    # Plot individual trajectories
    for x in t_dsq:
        axarr[0].plot(x[:,0], x[:,1])
        
    # Plot MSD with errors
    axarr[1].errorbar(time, mean_dsq, yerr=sterr_dsq, fmt='o')
    
    if fit_params is not None:
        label = r'D = {0:.5f} $\mu m^2/s$'.format(fit_params['dc'])     
        t_fit = np.arange(0,max(time), max(time)/5)
        y_fit = t_fit * fit_params['slope']
        axarr[1].plot(t_fit,y_fit,label=label)
        axarr[1].legend(loc='upper left', shadow=True)

    
    return fig, axarr

def make_fig_save_dirs(save_dir_root):
    """Create directories for saving figures in svg and png formats.
    
    Args:
        save_dir_root (string): root directory for saving figures

    Returns:
        dir_svg (string): valid directory
        dir_png (string): valid directory          
    """
    dir_png = os.path.join(save_dir_root, 'png')
    dir_svg = os.path.join(save_dir_root, 'svg')
    if not os.path.exists(dir_png):
        os.makedirs(dir_png)
    if not os.path.exists(dir_svg):
        os.makedirs(dir_svg)
        
    return dir_svg, dir_png

def save_fig(fig, file_name, dir_svg, dir_png):
    """Save figure in svg and png formats.
    
    Args:
        fig (matplotlib figure): figure to be saved
        file_name (string): desired name of figure file
        dir_svg (string): valid directory
        dir_png (string): valid directory  

    Returns:
        None         
    """
    fig_filename_svg = os.path.join(dir_svg, (file_name+'.svg'))
    fig_filename_png = os.path.join(dir_png, (file_name+'.png'))
    fig.savefig(fig_filename_svg)
    fig.savefig(fig_filename_png)
    
    return