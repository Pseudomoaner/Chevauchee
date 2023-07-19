/*
 * Takes 3 channel timelapses and makes a separate folder for each file, and each channel then saves
 * each timepoint as Frame%0d.tif
 * Prepares them for input into matlab analysis
 */

#@ File (label = "Input directory", style = "directory") input
output = input
#@ String (label = "File suffix", value = ".tif") suffix

processFolder(input);

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
	open(input+File.separator+file);
	original_title=getTitle;
	original_title=substring(original_title,0,indexOf(original_title,"_F-block.tif"));
	newfolder = input+File.separator+original_title;
	File.makeDirectory(newfolder);
	rename("WIP");
	run("Split Channels");
	selectWindow("C1-WIP");
	File.makeDirectory(newfolder+File.separator+"Channel_1");
	File.makeDirectory(newfolder+File.separator+"Channel_2");
	File.makeDirectory(newfolder+File.separator+"Channel_3");
	run("Image Sequence... ", "select="+newfolder+File.separator+"Channel_1"+File.separator+" dir="+newfolder+File.separator+"Channel_1"+File.separator+" format=TIFF name=Frame");
	//run("Image Sequence... ", "select=D:/Sean/SurfaceColonyPIV/Fluorescence_Blocks/Feb2_LB-f20_WT-R-C1TA-Y_OD01_0.5-12h_F-block/Channel_1/ dir=D:/Sean/SurfaceColonyPIV/Fluorescence_Blocks/Feb2_LB-f20_WT-R-C1TA-Y_OD01_0.5-12h_F-block/Channel_1/ format=TIFF name=Frame");
	close();
	selectWindow("C2-WIP");
	run("Image Sequence... ", "select="+newfolder+File.separator+"Channel_2"+File.separator+" dir="+newfolder+File.separator+"Channel_2"+File.separator+" format=TIFF name=Frame");
	close();
	selectWindow("C3-WIP");
	run("Image Sequence... ", "select="+newfolder+File.separator+"Channel_3"+File.separator+" dir="+newfolder+File.separator+"Channel_3"+File.separator+" format=TIFF name=Frame");
	close();
}
