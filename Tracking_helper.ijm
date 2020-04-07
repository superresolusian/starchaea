// file path stuff
dir = getDirectory("Choose base directory containing image and stardist folders");

images_dir = dir+"registered data";
rois_dir = dir+"stardist results";
tracks_dir = dir+"tracking results";

im_list = getFileList(images_dir);
rois_list = getFileList(rois_dir);

File.makeDirectory(tracks_dir); // prepare target folder for saving trackmate results

// variables to help loops and string manipulation
n_images = lengthOf(im_list);

prefix = "DRIFTCORRECTED_";
prefix_length = lengthOf(prefix);

// set correct measurements for trackmate script
run("Set Measurements...", "area center stack redirect=None decimal=3");

for(i=0; i<n_images; i++){

	// open first drift corrected image and get core filename
	im_name = im_list[i];
	im_stem = substring(im_name, prefix_length, lengthOf(im_name)-4);

	open(images_dir+"/"+im_name);
	run("RGB Color", "frames");

	// find corresponding rois file
	found = false;
	j = 0;
	
	while(!found){
		roi_string = rois_list[j];
		if(startsWith(roi_string, im_stem) && endsWith(roi_string, ".zip")) found = true;
		j++;
	}

	open(rois_dir+"/"+roi_string);

	// roi handling
	roiManager("Show All without labels");
	roiManager("Associate", "true");

	n_rois = roiManager("count");
	roi_index_array = newArray(n_rois);
	for(r=0; r<n_rois; r++) roi_index_array[r] = r;

	roiManager("Select", roi_index_array);
	roiManager("Measure");

	print("gonna save to "+tracks_dir);
	//tracks_dir_ = replace(tracks_dir, " ", "\\^");
	run("my tracking", "save_dir=["+tracks_dir+"]");

	// close everything for next iteration
	selectWindow("Results");
	run("Close");
	run("Close All");
	roiManager("reset");
}
