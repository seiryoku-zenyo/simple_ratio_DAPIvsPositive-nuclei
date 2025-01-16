//This script aims at batch processing dual channel data, where nuclei are tested positive/nbegative for a second marker.
//There are 3 main functions, one for segmenting the structures, a second one for analysing and writing a .csv file and a third to create representative images.

//In image analysis it's always important to keep in mind that sofware updates and the nature of different datasets, may result in suboptimal performance
//and the need to code maintenance. For any issues contact vasco.fachada@gmail.com

//Code: Vasco Fachada
//Principal investigator:Riikka Kivelä 
//Experimentists: Laura Ylä-Outinen; Kalle Kolari
//Other researchers involved: Sakari Mäntyselkä, Erik Niemi

while (nImages>0) { 
		selectImage(nImages); 
        close();    				  
	} 
	
list = getList("window.titles");
	for (i=0; i<list.length; i++){
    	winame = list[i];
      	selectWindow(winame);
      	print(winame);
     	run("Close");
     }


print("Opening data file...");
//Finding the files
studydir = getDirectory("Select your studiy´s main directory");
dir = studydir+File.separator+"data";
//dir = studydir+File.separator+"testing";
files = getFileList(dir);
start = getTime();
//Creating folder where ready images and raw data file will be be saved
savingDir = studydir+File.separator+"results";
File.makeDirectory(savingDir);
datafile = savingDir+File.separator+"datafile.csv";

if (File.exists(datafile)){
	Dialog.create("There is already a datafile.csv in your resultas folder!");
	items=newArray(2);
	items[0]="Sure, those were crappy results anyways, let's make new ones!";
	items[1]="NO! I don't want to lose my previous results.";
	Dialog.addRadioButtonGroup("Continuing will replace your previous datafile with a new one. \nAre you sure you want to continue?\n", items, 2, 1, items[1]);
	Dialog.show();
	if (Dialog.getRadioButton == items[1]){
		exit("yup, better backup your precious results somewhere else ;)");
	}else{
		File.delete(datafile);	
	}
}

//setBatchMode(true);
run("Bio-Formats Macro Extensions");
for(f=0; f<files.length; f++) {
	
	//get information from file, treatment, well and sample names into variables:
	if(endsWith(files[f], ".nd2")) { 		
		id = dir+File.separator+files[f];
		Ext.setId(id);
		Ext.getSeriesCount(seriesCount);
		filename=substring(files[f], 0, lengthOf(files[f])-4);
		treatment=substring(files[f], 0, 1);
		well=substring(files[f], 1, 4);
		rep=substring(files[f], 4, 5);

		print('opening file: '+filename);
		print('Treatment: '+treatment);
		print('Well: '+well);
		print('Replicate (n): '+rep);
		

		//OPENING FILE AND READING METADATA
		
		for (i=0; i<seriesCount; i++) {
			run("Bio-Formats Importer", "open=["+id+"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+(i+1));
			rename("opened");
			getVoxelSize(pWidth, pHeight, depth, unit); //detecting the voxel size of the current raw file, so the final image also has the same scale.
			Stack.getDimensions(width, height, channels, slices, frames);
			run("Properties...", "channels="+channels+" slices="+slices+" frames="+frames+" unit="+unit+" pixel_width="+pWidth+" pixel_height="+pHeight+" voxel_depth="+depth+" global");
			run("Split Channels");

			classes = newArray("all_nuc", "pos_nuc");
			bins = newArray(2);
			for (c=0; c<classes.length; c++){
				selectWindow("C"+c+1+"-opened");
				//run("Images to Stack", "name="+classes[c]+" title=C"+c+1);
				rename(classes[c]);
				run("HiLo");
				bins[c] = MLseg(classes[c]);
			}
			Analysis(bins);
			saveImgs();
			
		}
	
		while (nImages>0) { 
			selectImage(nImages); 
        	close();   
		} 
		call("java.lang.System.gc");
	}
}

end = getTime();
//saveAs("Text", savingDir+File.separator+"analysis_Log.txt");
Dialog.create("DONE!");
Dialog.addMessage("The analysis was completed within "+(end-start)/1000+" seconds ("+(end-start)/1000/60+" minutes).\n \nI can tell you, these numbers were hard to crunch, hope you get a Nobel for this ;)");
Dialog.show();

//This function uses two trained deep learning models to segment and classify all the detected nuclei + the nuclei positive for the second marker.
function MLseg(classifier){
	//setBatchMode(false);
	selectWindow(classifier);
	run("Duplicate...", "title="+classifier+"_bin duplicate");
	run("Trainable Weka Segmentation");
	wait(1000);
	//selectWindow("Trainable Weka Segmentation v3.2.34");
	//print(dir+File.separator+"code"+File.separator+classifier+".model");
	call("trainableSegmentation.Weka_Segmentation.loadClassifier", studydir+"code"+File.separator+classifier+".model");
	//call("trainableSegmentation.Weka_Segmentation.loadClassifier", "C:\\Users\\Vasco\\Desktop\\Laura\\code\\all_nuc.model");
	wait(2000);
	call("trainableSegmentation.Weka_Segmentation.getProbability");
	close("Trainable Weka Segmentation v3.2.34");

	selectWindow("Probability maps");
	run("Delete Slice", "delete=channel");
	run("Convert to Mask", "method=Otsu background=Light calculate");
	run("Despeckle", "stack");
	run("Remove Outliers...", "radius=2 threshold=50 which=Dark");
	run("Watershed", "stack");
	rename(classifier+"_bin");
	bins = getImageID();
	close("Probability maps");
	return bins;
}


//This function counts the number and measures the area of both total and postive-only nuclei. The ration between these variables is also calculated.
//Finally the function creates and writes a .csv file with these measurements.
function Analysis(bins){	
	//Create datafile in case there is none
	if (!File.exists(datafile)){
    	// Create a header
    	File.append("Image,Treatment,Well,Parallels,Series,All nuclei number,All nuclei total area,All nuclei AvgSize,Positive nuclei number,Positive nuclei total area,Positive nuclei AvgSize,%Positive nuclei number,%Positive nuclei area", datafile);
	}

	totCount = newArray(2);
	avgSize = newArray(2);
	totSize = newArray(2);
	for (b=0; b<classes.length; b++){
		//print(bins[b]);
		//selectWindow(classes[ci]+"_bin");
		selectImage(bins[b]);
		run("Analyze Particles...", "size=60-Infinity show=Outlines exclude summarize");
		run("8-bit");
		run("Invert");
		selectWindow("Summary");
		IJ.renameResults("Summary","Results");
		totCount[b] = getResult("Count", 0);
		avgSize[b] = getResult("Average Size", 0);
		totSize[b] = getResult("Total Area", 0);
		close("Results");
	}
	ratioNum = (totCount[1]*100)/totCount[0];
	ratioArea = (totSize[1]*100)/totSize[0];

    // insert data into table
	File.append(filename+"_series"+i+1+","+treatment+","+well+","+rep+","+i+1+","+totCount[0]+","+totSize[0]+","+avgSize[0]+","+totCount[1]+","+totSize[1]+","+avgSize[1]+","+ratioNum+","+ratioArea+",", datafile);

	//run("Clear Results");
}

//This function creates montage figures representative of the performed analysis, both for inspection and use in publications.
function saveImgs(){	
	setForegroundColor(255, 255, 255);
	run("Images to Stack", "name=bin title=Drawing use");
	run("Make Montage...", "columns=2 rows=1 scale=1 border=3 label use");
	run("RGB Color");
	
	run("Images to Stack", "name=original title=_nuc use");
	run("8-bit");
	run("Fire");
	run("Make Montage...", "columns=2 rows=1 scale=1 border=3 label use");
	run("RGB Color");

	run("Images to Stack", "method=[Copy (center)] name=final_montage title=Montage");
	run("Reverse");
	run("Make Montage...", "columns=1 rows=2 scale=1 border=6 use");
	saveAs("tiff", savingDir+File.separator+filename+"_series"+i+1+".tif");
}

