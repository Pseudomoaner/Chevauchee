/*
 * Pre processes images for PIV analysis
 * Expects a folder containing individual tifs of just the brightfield for each timepoint
 * Blurs each frame, then takes a median projection, then divides the original stack by this projection
 * Eliminates uneven illumination and shadows from condensation
 * Based on Macro template to process multiple images in a folder
 */

#@ File (label = "Input directory", style = "directory") input
output = input;
#@ String (label = "File suffix", value = ".tif") suffix

processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
}

function processFile(input, output, file) {
	print(input+File.separator+file);
	open(input+File.separator+file);
	original_title=getTitle;
	original_title=substring(original_title,0,indexOf(original_title,"_BF-block.tif"));
	

	rename("WIP");
	run("Select All");
	run("Duplicate...", "duplicate");
	run("Gaussian Blur...", "sigma=5 stack");
	rename("blurred_stack");
	run("Z Project...", "projection=Median");
	rename("MED");

	selectWindow("blurred_stack");
	close();

	imageCalculator("Divide create 32-bit stack", "WIP","MED");
	selectWindow("WIP");
	close();
	selectWindow("MED");
	close();
	run("Invert", "stack");
	run("16-bit");

	newfolder = input+File.separator+original_title;
	File.makeDirectory(newfolder);
	rename(original_title+"_BG-subtracted_");
	run("Image Sequence... ", "select="+newfolder+File.separator+" dir="+newfolder+File.separator+" format=TIFF");
	close();
}


