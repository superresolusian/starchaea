// make dialog
makeDialog();

// grab values from dialog
is_drift_corrected = Dialog.getCheckbox();
is_two_colour = Dialog.getCheckbox();
tracking_channel = Dialog.getNumber();

// get directory
dir = getDirectory("Choose base directory containing image and stardist folders");

if(is_drift_corrected){
	images_dir = dir+"registered data";
}
else{
	images_dir = dir+"raw data";
}

rois_dir = dir+"stardist results";

if(is_two_colour){
	rois_dir = rois_dir+"/channel "+tracking_channel;
}

tracks_dir = dir+"tracking results";

IJ.log(rois_dir)

im_list = getFileList(images_dir);
rois_list = getFileList(rois_dir);

File.makeDirectory(tracks_dir); // prepare target folder for saving trackmate results

// variables to help loops and string manipulation
n_images = lengthOf(im_list);

prefix = "";
if(is_drift_corrected){
	prefix = "DRIFTCORRECTED_";
}

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

	print("Saving Trackmate output to "+tracks_dir);
	//tracks_dir_ = replace(tracks_dir, " ", "\\^");
	run("my tracking", "save_dir=["+tracks_dir+"]");

	// close everything for next iteration
	selectWindow("Results");
	run("Close");
	run("Close All");
	roiManager("reset");
}

function makeDialog(){
	Dialog.create("Run Trackmate on Stardist detections");

	Dialog.addCheckbox("Drift-corrected data?", 0);
	Dialog.addCheckbox("Two-colour data?", 0);
	Dialog.addNumber("Channel for running tracking (only relevant for two-colour data", 1);

	Dialog.show();
}