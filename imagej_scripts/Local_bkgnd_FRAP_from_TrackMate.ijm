// Takes the output of a TrackMate track and extracts mean intensities of all pixels 
// contained within an inner and outer radius of each spot. Additionally, measures
// the mean intensity of two larger circles centered around the spot, with inner radius
// r_ring_in and outer radius r_ring_out. This information can be used to create a 
// local background trace in FRAP experiments and compare fluorescence recovery within
// a small cluster located in the center of the spot and a larger diffuse area around
// the small cluster.

// To run this macro, a TrackMate file must be loaded with the accompanying image,
// the image stack file must be selected as ImageJ's currently active image,
// and TrackMate's analyis windows must be open (click "Analysis" in the main
// TrackMate window). This macro loads spot data directly from TrackMate's
// output table called 'Spots in tracks statistics'.

// Written by Vladislav Belyy on 2019-02-21.

//// PARAMETERS CONFIGURED BY USER ///////////////////////////////////////////////
// The variables in this section can be modified as needed

// Set the inner and outer radii of the spots
inner_radius = 1.2;
outer_radius = 1.9;

// Set the inner and outer radii of the larger ring used for background correction

r_ring_in = 2;
r_ring_out = 3;

// Saving params
save_output = true;
save_dir = 'spot_radii';
file_name_suffix = '_FRAP_with_local_bkgnd.txt'

debug_mode = true; // If true, runs slower but with more visual output

//// CODE STARTS HERE - DON'T MODIFY UNLESS YOU KNOW WHAT YOU'RE DOING //////////

if(!debug_mode) {
	setBatchMode(true); 
}
if(save_output) {
	parent_dir = getDirectory("image");
	save_dir_full = parent_dir + save_dir;
	File.makeDirectory(save_dir_full);
}

getPixelSize(unit, pw, ph, pd);
inner_px = inner_radius / pw;
outer_px = outer_radius / pw;
inner_ring_px = r_ring_in / pw;
outer_ring_px = r_ring_out / pw;

img_name = getTitle(); // Title of current image
dot_index = indexOf(img_name, "."); 
title = substring(img_name, 0, dot_index); 

spot_table = 'Spots in tracks statistics'; // TrackMate's name for the spot table

// This would be nice but table sorting is currently glitchy:
//Table.sort('FRAME', spot_table);
//Table.update(spot_table);

frames = Table.getColumn('FRAME', spot_table);
spot_IDs = Table.getColumn('ID', spot_table);
track_IDs = Table.getColumn('TRACK_ID', spot_table);
pos_x = Table.getColumn('POSITION_X', spot_table);
pos_y = Table.getColumn('POSITION_Y', spot_table);

// Remove any existing ROIs and clear results window
roiManager("reset");
run("Clear Results");

// Draw inner and outer circle ROIs for all relevant slices
run("Set Measurements...", "area mean redirect=None decimal=3");
for (i=0; i<spot_IDs.length; i++) {
	setSlice(frames[i]+1);
	x_inner = (pos_x[i] + pw/2 - inner_radius) / pw;
	y_inner = (pos_y[i] + pw/2  - inner_radius) / pw;
	makeOval(x_inner, y_inner, inner_px*2, inner_px*2);
	roiManager("add")

	// measure outer diameter intensity
	x_outer = (pos_x[i] + pw/2 - outer_radius) / pw;
	y_outer = (pos_y[i] + pw/2  - outer_radius) / pw;
	makeOval(x_outer, y_outer, outer_px*2, outer_px*2);
	roiManager("add")

	// maeasure inner and outer ring intensities
	x_inner_ring = (pos_x[i] + pw/2 - r_ring_in) / pw;
	y_inner_ring = (pos_y[i] + pw/2  - r_ring_in) / pw;
	makeOval(x_inner_ring, y_inner_ring, inner_ring_px*2, inner_ring_px*2);
	roiManager("add")

	x_outer_ring = (pos_x[i] + pw/2 - r_ring_out) / pw;
	y_outer_ring = (pos_y[i] + pw/2  - r_ring_out) / pw;
	makeOval(x_outer_ring, y_outer_ring, outer_ring_px*2, outer_ring_px*2);
	roiManager("add")
	
	
}

// Measure intensities contained within the inner and outer radii
inner_indexes = Array.getSequence(spot_IDs.length);
outer_indexes = newArray(spot_IDs.length);
inner_ring_indexes = newArray(spot_IDs.length);
outer_ring_indexes = newArray(spot_IDs.length);

for (i = 0; i < inner_indexes.length; i++) {
	inner_indexes[i] = inner_indexes[i] * 4;
	outer_indexes[i] = inner_indexes[i] + 1;
	inner_ring_indexes[i] = inner_indexes[i] + 2;
	outer_ring_indexes[i] = inner_indexes[i] + 3;
}

roiManager("select", inner_indexes);
roiManager("measure");
mean_int_inner = Table.getColumn("Mean", "Results");
area_inner = Table.get("Area", 0, "Results");
run("Clear Results");
roiManager("select", outer_indexes);
roiManager("measure");
mean_int_outer = Table.getColumn("Mean", "Results");
area_outer = Table.get("Area", 0, "Results");
run("Clear Results");
roiManager("select", inner_ring_indexes);
roiManager("measure");
mean_int_ring_inner = Table.getColumn("Mean", "Results");
area_ring_inner = Table.get("Area", 0, "Results");
run("Clear Results");
roiManager("select", outer_ring_indexes);
roiManager("measure");
mean_int_ring_outer = Table.getColumn("Mean", "Results");
area_ring_outer = Table.get("Area", 0, "Results");
run("Clear Results");

roiManager("deselect");

//roiManager("measure");

// Build array of background intensities for each spot:
small_ring_int = newArray(spot_IDs.length);
local_bkgnd_int = newArray(spot_IDs.length);
for (i = 0; i < mean_int_inner.length; i++) {
	small_ring_int[i] = (area_outer * mean_int_outer[i] - area_inner * mean_int_inner[i]) /
						( area_outer - area_inner);
	local_bkgnd_int[i] = (area_ring_outer * mean_int_ring_outer[i] -
						area_ring_inner * mean_int_ring_inner[i]) /
						( area_ring_outer - area_ring_inner);
}

setBatchMode(false);

Table.showArrays('Extraction', spot_IDs, track_IDs, frames,pos_x,pos_y, mean_int_inner, 
	mean_int_outer, small_ring_int, local_bkgnd_int);


if(save_output) {
	save_file = save_dir_full + File.separator + title + file_name_suffix;
	Table.save(save_file, 'Extraction');
}

// Clean up all unnecessary windows if not in debug mode
if(!debug_mode) {
	 close(spot_table);
	 close("Links in tracks statistics");
	 close("Track statistics");
	 close("Results");
	 close("Extraction");
	 close("ROI Manager");
}
