""" 

Used for measuring the intensity of a selected ROI throughout an
entire stack and saving it (along with useful metadata) as a json
file so that it can be used, e.g., for background subtraction in FRAP 
experiments.

by Vladislav Belyy (UCSF)

"""

import java.awt.Color as Color
from ij import WindowManager as WindowManager
from ij.plugin.frame import RoiManager as RoiManager
from ij.process import ImageStatistics as ImageStatistics
from ij.measure import Measurements as Measurements
from ij import IJ as IJ
from ij.measure import CurveFitter as CurveFitter
from ij.gui import Plot as Plot
from ij.gui import PlotWindow as PlotWindow
from ij.io import SaveDialog as SaveDialog
import os
import math
import json
 
# Get ROI
roi = IJ.getImage().getRoi()

# Specify time units
frame_interval = 1 # data will be in frames
time_units = 'frames'
 
# Get current image plus and image processor
current_imp  = WindowManager.getCurrentImage()
stack        = current_imp.getImageStack()
calibration  = current_imp.getCalibration()


# Get current dir and filename
stack_title = current_imp.getTitle()
file_info = current_imp.getOriginalFileInfo()
file_dir = file_info.directory
n_slices = stack.getSize()

#############################################

# Get details about the ROI (for saving as metadata)
roi_bounds = roi.getBounds()
roi_type = roi.getTypeAsString()
roi_area = current_imp.getStatistics().area

roi_bounds_list = {'x':roi_bounds.x, 'y':roi_bounds.y, 'w':roi_bounds.width, 'h':roi_bounds.height}
 
# Collect intensity values
 
# Create empty lists of numbers
intensities = []  # Intensity values
 
# Loop over each slice of the stack
for i in range(0, n_slices):
  
    # Get the current slice 
    ip = stack.getProcessor(i+1)
  
    # Put the ROI on it
    ip.setRoi(roi)
  
    # Make a measurement in it
    stats = ImageStatistics.getStatistics(ip, Measurements.MEAN, calibration);
    mean  = stats.mean
  
    # Store the measurement in the list
    intensities.append( mean  )
 
IJ.log('For image ' + current_imp.getTitle() )
IJ.log('Time interval is ' + str(frame_interval) + ' ' + time_units)

# Build plot   
x = [i * frame_interval for i in range( n_slices ) ] 
y = intensities

  
plot = Plot("Backgrouncurve " + current_imp.getTitle(), "Time ("+time_units+')', "NU", [], [])
plot.setLimits(0, max(x), 0, max(y))
plot.setLineWidth(2)
 
 
plot.setColor(Color.BLACK)
plot.addPoints(x, y, Plot.LINE)
plot.addPoints(x,y,PlotWindow.X)
 
plot.setColor(Color.black)
plot_window =  plot.show()

###############################

# Save data as a json file

###############################

# Ask for filename
savename_temp = os.path.splitext(stack_title)[0] + '_cell_XX'
save_file = SaveDialog('Please choose a location to save results', file_dir, savename_temp, '.json')

# Prepare data for export
file_info = {'source_file':stack_title, 'source_file_dir':file_dir}
roi_info = {'roi_type':roi_type, 'roi_bounds_pixel':roi_bounds_list, 'roi_area':roi_area}
calb = {'pixel_size':calibration.pixelHeight, 'unit':calibration.getXUnit()}

data_out = {'file_info':file_info, 'roi_info':roi_info, 
	'intensities':intensities, 'calibration':calb}

# Save file
save_file_name = save_file.getFileName()
save_file_dir = save_file.getDirectory()
save_file_full = save_file_dir + save_file_name
with open(save_file_full, 'w') as outfile:
	json.dump(data_out, outfile)
