// Open a set of Micro-Manager z-stacks from an acquisition folder and save them as average or max intensity
// projection tiffs, split by channel. The main data folder should contain sub-folders (each sub-folder in considered
// a separate condition). The names of all sub-folders should start with the user-defined parameter
// "template_name", which can be modified below. Each sub-folder can contain any number of multi-color
// TIFF files, one per position. Works both for single-timepoint images and movies.
//
// The output is placed in a new directory, which can be set below by changing save_dir_name, and
// consists of one sub-folder per condition, with one sub-folder for each individual file within
// each condition sub-folder containing one z-stack .tif file for each channel.


// Originally written by Vlad Belyy on 2019-04-30.

// These need to be set manually
template_name = "vVB_"; // what do all the data folders start with?
save_results = true;
save_dir_name = "Max_intensity_stacks";

z_stack = true; // make a z-stack?
// options for projection: "Average Intensity", "Max Intensity"
projection_type = "Max Intensity" // what type of projection to use?

// Ask user to point to data file directory
dir_working = getDirectory("Choose the parent directory of uManager data files");


// Create a list of directories containing valid acquisition data
valid_dir_list = newArray();
base_dir_list = getFileList(dir_working);
for (i = 0; i < base_dir_list.length; i++) {
    if (startsWith(base_dir_list[i], template_name)==1) {
    	valid_dir_list = Array.concat(valid_dir_list, base_dir_list[i]);
    }
}

Array.sort(valid_dir_list);
num_conditions = valid_dir_list.length;

// Iterate over conditions (assuming one directory per condition) and process data
setBatchMode(true);
for (i = 0; i < num_conditions; i++) {
	close("*"); // Clear all open windows

	// Prepare save dir for current condition
	if (save_results) {
		curr_condition_name = valid_dir_list[i];
		save_dir_base = dir_working + File.separator + save_dir_name + File.separator;
		File.makeDirectory(save_dir_base);
		save_dir_condition = save_dir_base + curr_condition_name + File.separator;
		File.makeDirectory(save_dir_condition);
	}

	// Build list of files for current condition
	curr_dir = dir_working + valid_dir_list[i];
	// only use .ome.tif files
	files_in_curr_dir_tmp = getFileList(curr_dir);
	files_in_curr_dir = newArray();
	for (i_file = 0; i_file < files_in_curr_dir_tmp.length; i_file++) {
		if (endsWith(files_in_curr_dir_tmp[i_file], ".ome.tif") == 1) {
			files_in_curr_dir = Array.concat(files_in_curr_dir, files_in_curr_dir_tmp[i_file]);
		}
	}

	// Open all the stacks(files) from the current condition one by one and process them
	for (ii = 0; ii < files_in_curr_dir.length; ii++) {		
		curr_file = files_in_curr_dir[ii]; 
		curr_file_path = curr_dir + curr_file;
		open(curr_file_path);

		// Prepare save dir for current file
		if (save_results) {
			save_dir_file = save_dir_condition + curr_file + File.separator;
			File.makeDirectory(save_dir_file);
		}

		// Make z-stack, if needed
		title = getTitle();
		if (z_stack) {
			run("Z Project...", "projection=["+projection_type+"] all");
			close(title);
		}

		// separate file into channels and close original file
		title = getTitle();
		getDimensions(width, height, channels, slices, frames);
		run("Split Channels");
		close(title);

		print(title); // for debugging

		// Save each channel individually
		for (iii = 0; iii < channels; iii++) {
			curr_title = "C"+parseInt(iii+1)+"-"+title;
			print(curr_title);
			selectImage(curr_title);
			run("Enhance Contrast", "saturated=0.35");
			
			if (save_results){
				saveAs("tiff", save_dir_file + curr_title);
			}
			close(curr_title);
		}
	}
	
	
}
setBatchMode(false);
close("Log");