// Written by Hiroyuki Hakozaki hhakozaki@ucsd.edu 4/5/2020
// This macro split channel and save each channel as tiff format
// Use Name #2, Name #3, Name #4 tag to identify staining
// Save DAPI and STAT1 channel - _DAPI, _STAT1
// It also perform backgrouns subtraciton

macro "Batch Channel Splitter" {
	// Check if any image is opened.
	// By selecting OK button, macro will close all opened images
	if (nImages != 0) {
		Dialog.create("Batch Channel Splitter Macro");
		Dialog.addMessage("There are images open.\nThis macro require to close all images\nDo you want to close all Images?");
		Dialog.show();
		close("*");
	}

	// Enable Bio-Format Importer
	run("Bio-Formats Macro Extensions");

	// Get source and destination directory
	SrcDir = getDirectory("Choose Source Image Directory");
	DstDir = getDirectory("Choose Distination Directory");
	
	// Create Dialog for recursive mode selection.
	Dialog.create("Batch Tiff Converter");
	Dialog.addCheckbox("Process Sub Folders:", true);
	// If not selected, it won't create subfolder and all image will be saved into one folder
	Dialog.addCheckbox("Create Sub Folders in Destination:", true);
	Dialog.show();
	// Get recursive option selection
	Recursive = Dialog.getCheckbox();
	SubFolder = Dialog.getCheckbox();

	// Enable batch mode
	setBatchMode(true);
	// If recursive mode is selected, call recursive function
	if (Recursive == true) BatchSplitter_Recursive(SrcDir, DstDir, SubFolder);
	// If not recursive, directly call BatchSplitter function
	else BatchSplitter(SrcDir, DstDir);
	// Disable batch mode
	setBatchMode(false);
}

//Recursively process into subfolders
function BatchSplitter_Recursive(SrcDir, DstDir, SubFolder) {

	// Process primary directory first
	BatchSplitter(SrcDir, DstDir);

	// Look for the directory and do recursive processing
	FileList = getFileList(SrcDir);
	for (i=0; i<FileList.length; i++) {
		if (File.isDirectory(SrcDir + FileList[i])) {
			// If SubFolder option is true, create subfolder in destination folder
 			if (SubFolder == true) {
				print("Found Directory : "+FileList[i]);
				NewDstDir = DstDir+FileList[i];
				if (!File.exists(NewDstDir)) {
					print("Create Distination Directory : "+NewDstDir);
					File.makeDirectory(NewDstDir);
				}
				// Recursive call to process children folder
				BatchSplitter_Recursive(SrcDir+FileList[i], NewDstDir, SubFolder);
 			}
 			else {
 				// SubFolder option is false. No subfolder creation
 				// Just recursive call to process children holder
 				BatchSplitter_Recursive(SrcDir+FileList[i], DstDir, SubFolder);				
 			}
		}
	}
}

// Main function to perform Z projection
function BatchSplitter(SrcDir, DstDir) {

	// Get file list in the folder
	FileList = getFileList(SrcDir);
	// Process each file
	for(i=0; i<FileList.length; i++) {
		// Check if Bio-Format Importer can open the image or not.
		FormatType = false;
		// If file isn't directory, get format type
		if (!File.isDirectory(SrcDir + FileList[i])) {
			Ext.isThisType(SrcDir+FileList[i], FormatType);
		}
		// If Fileformat can be opened by Bio-Format Importer
		// Process the file
		if (FormatType!=0) {
			// Initialize Bio-Format Importer to read number of series
			Ext.setId(SrcDir+FileList[i]);
			// Read number of series contained in the file
			Ext.getSeriesCount(seriesCount);

			// Process each series
			for(s=1; s<=seriesCount; s++) {
				filename = FileList[i];
				// Output log for open image
				print("Open Image :" + filename + ", Series :" + s);

				// Set command string to open the image
				cmd = "open=[" + SrcDir + FileList[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+d2s(s, 0);
				// Open image with Bio-Formats Importer
				run("Bio-Formats Importer", cmd);

				// Rename image for easy handling
				rename("Org");

				// Array to store channel name in metadata
				chname = newArray(4);
				// Array to store channel number - 0:DAPI, 1:STAT1, 2:
				chnum = newArray(3);

				// Get meta information
				metatxt = getImageInfo();

				// Find "Name #2 = " location and read parameter
				// This is for channel 1
				start = indexOf(metatxt, "Name #2 = ")+10;
				end = indexOf(metatxt, "\n", start);
				chname[0] = substring(metatxt, start, end);

				// Find "Name #3 = " location and read parameter
				// This is for channel 2
				start = indexOf(metatxt, "Name #3 = ")+10;
				end = indexOf(metatxt, "\n", start);
				chname[1] = substring(metatxt, start, end);
				
				// Find "Name #4 = " location and read parameter
				// This is for channel 3
				start = indexOf(metatxt, "Name #4 = ")+10;
				end = indexOf(metatxt, "\n", start);
				chname[2] = substring(metatxt, start, end);

				// Set DAPI and STAT1 channel number
				for(j=0; j<3; j++) {
					if (chname[j].contains("405")) dapiCh=j+1;
					if (chname[j].contains("561")) stat1Ch=j+1;				
				}

				// Split Channel
				run("Split Channels");

				// Set new file name base
				// If there is more than 1 series, put series number
				if (seriesCount != 1) NewFileName = File.nameWithoutExtension() + "_S" + d2s(s, 0);
				// Only one series, no series number
				else  NewFileName = File.nameWithoutExtension();

				// Save DAPI Channel
				dapiTitle = "C"+dapiCh+"-Org";
				selectWindow(dapiTitle);
				// Subtract background
				run("Subtract Background...", "rolling=200 stack");
				saveAs("Tiff", DstDir+NewFileName+"_DAPI");
				
				// Save STAT1 Channel
				stat1Title = "C"+stat1Ch+"-Org";
				selectWindow(stat1Title);
				// Subtract background
				run("Subtract Background...", "rolling=200 stack");
				saveAs("Tiff", DstDir+NewFileName+"_STAT1");

				// Close all open image
				close("*");
			}
		}
	}
}




