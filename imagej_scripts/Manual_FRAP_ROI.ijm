// This macro can be used to quickly select bleached ROIs and corresponding
// background ROIs in a typical FRAP experiment. Once the radius is set,
// a left click draws a FRAP ROI, and a subsequent click draws a corresponding
// background ROI. The macro ends when the user presses the right mouse button.


// Written by Vladislav Belyy at UCSF on 2019-02-22.

// User-adjustable parameters ////////////////////////
save_output = true;
save_dir = 'manual_ROI_FRAP';
default_radius = 2;

// Main routine //////////////////////////////////////

radius = get_radius(default_radius);

// populate ROI manager with user-selected circles
roiManager("Reset");
roiManager("Show All");
draw_circles_on_mouse_clicks(radius);

// measure intensities in all ROIs
if(save_output) {
	setBatchMode(true); 
}
roiManager("deselect");
run("Clear Results");
run("Set Measurements...", "area mean redirect=None decimal=3");
roiManager("Multi Measure");

// save results and clean up
if(save_output) {
	save_results(save_dir);
}

// Sub-routines //////////////////////////////////////

function save_results(save_dir) {
	parent_dir = getDirectory("image");
	img_name = getTitle(); // Title of current image
	dot_index = indexOf(img_name, "."); 
	title = substring(img_name, 0, dot_index); 
	
	save_dir_base = parent_dir + File.separator + save_dir + File.separator;
	save_dir_ROIs = save_dir_base + File.separator + "ROI_files";
	save_dir_results = save_dir_base + File.separator + "intensities";
	File.makeDirectory(save_dir_base);
	File.makeDirectory(save_dir_ROIs);
	File.makeDirectory(save_dir_results);

	roiManager("Show None");
	roiManager("save", save_dir_ROIs + File.separator + title + "_ROIs.zip");
	Table.save(save_dir_results + File.separator + title + 
		"_manual_FRAP_results.csv", "Results");
	close("Results");
	close("ROI Manager");

	setBatchMode(false);
}

function get_radius (default_r) {
    Dialog.create("Settings"); 
    Dialog.addNumber("Set radius of circle", default_r); 
    Dialog.show(); 
    radius = Dialog.getNumber();
    return radius;
}

function draw_circles_on_mouse_clicks (radius) {

	setOption("DisablePopupMenu", true); 
	getPixelSize(unit, pixel_width, pixel_height); 
	setTool("rectangle"); 
	left_button=16; 
	right_button=4; 
        
	height = 2*radius/pixel_height; 
	width = 2*radius/pixel_width; 
	x2=-1; y2=-1; z2=-1; flags2=-1; 
	getCursorLoc(x, y, z, flags); 
	was_left_pressed = false;
	even_label = false;
	while (flags&right_button==0){ 
		getCursorLoc(x, y, z, flags); 
                
		if (flags&left_button!=0) { 
			// Wait for it to be released 
			was_left_pressed = true; 
		} else if (was_left_pressed) { 
			was_left_pressed = false; 
			if (x!=x2 || y!=y2 || z!=z2 || flags!=flags2) { 
				x = x - width/2; 
				y = y - height/2; 
				makeOval(x, y, width, height); 
				roiManager("Add");
				if (even_label) {color="red";} else {color="green";}
				even_label = !even_label;
				Overlay.addSelection(color, 5);                      
			}
		} 
	} 
	setOption("DisablePopupMenu", false);
	Overlay.clear;
} 