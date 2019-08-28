// Open a set of Micro-Manager z-stacks from an acquisition folder and save them as average or max intensity
// projection tiffs, split by channel. This version of the macro is capable of processing files in arbitrarily
// deep nested folders. Works both for single-timepoint images and movies.
//
// The output is placed in a single new directory, which can be set below by changing save_dir_name.
// One .tif file for each channel is saved. The directory structure is not preserved, so files
// _MUST HAVE UNIQUE NAMES_ for the macro to work correctly.


// Originally written by Vlad Belyy on 2019-08-27.

// Options that should be set manually
save_results = true;
save_dir_name = "Avg_intensity_stacks";

z_stack = true; // make a z-stack?
// options for projection: "Average Intensity", "Max Intensity"
projection_type = "Average Intensity" // what type of projection to use?

// Ask user to point to data file directory
dir_working = getDirectory("Choose the parent directory of all uManager data files");


files = newArray();
files = listFiles(dir_working, files); 

function listFiles(dir, files) {
	list = getFileList(dir);
	for (i=0; i<list.length; i++) {
		if (endsWith(list[i], "/")) {
			files = listFiles(""+dir+list[i], files);
		}
		else {
			file = dir + list[i];
			files = Array.concat(files, file);
		}
	}
	return files;
}


// Prepare save dir
if (save_results) {
	save_dir_base = dir_working + File.separator + save_dir_name + File.separator;
	File.makeDirectory(save_dir_base);
}

print(files.length);
Array.sort(files);

// Iterate over files and process data
setBatchMode(true);
print("starting loop");
print(files.length);
for (i = 0; i < files.length; i++) {
	close("*"); // Clear all open windows
	file = files[i];
	print(file);

	// only use .ome.tif files
	if(endsWith(file, ".ome.tif") == 1) {
		open(file);
	
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
			selectImage(curr_title);
			run("Enhance Contrast", "saturated=0.35");
				
			if (save_results){
				saveAs("tiff", save_dir_base + curr_title);
			}
			close(curr_title);
		}
	}
}
print("finished loop");
	
setBatchMode(false);
close("Log");