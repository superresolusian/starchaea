// file path stuff
dir = getDirectory("Choose results>crop directory for this dataset");

// set correct measurements for trackmate script
run("Set Measurements...", "stack redirect=None decimal=3");

n_rois = roiManager("count");
roi_index_array = newArray(n_rois);
for(r=0; r<n_rois; r++) roi_index_array[r] = r;

if(isOpen("Results")) run("Close");

roiManager("Select", roi_index_array);
roiManager("Measure");

saveAs("Results", dir+"curated_results.csv");