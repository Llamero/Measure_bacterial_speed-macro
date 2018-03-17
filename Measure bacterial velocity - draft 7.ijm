close("*");
run("Set Measurements...", "min redirect=None decimal=3");
dir = getDirectory("Please select the input directory.");
fileList = getFileList(dir);
map = getBoolean("Do you want to generate a velocity map for each sample? (NOTE: This greatly slows the analysis)");
nBins = 100; //Variable to set number of bins in histogram
reduce = 1; //Interval at which to sample stack (every x frames);
setBatchMode(true);

//Create an ouput directory to save the results
outDir = dir + "Output\\";
File.makeDirectory(outDir);

//Find all nd2 files in directory
nd2Count = 0;
nd2Array = newArray(0);
for(b=0; b<2; b++){	
	for(a=0; a<fileList.length; a++){
		if(endsWith(fileList[a], ".nd2")){
			if(b) nd2Array[nd2Count] = fileList[a];
			nd2Count++;
		}
	}
	if(!b){
		nd2Array = newArray(nd2Count);
		nd2Count = 0;
	}
}

//Process each nd2 file in the list
run("Bio-Formats Macro Extensions");
for(a=0; a<nd2Array.length; a++){
	//Get metadata from the file
	Ext.setId(dir + nd2Array[a]);
	Ext.getSizeT(nFrames); //Number of frames
	Ext.getSizeZ(sizeZ); //Number of slices
	Ext.getSizeC(sizeC); //Number of channels
	Ext.getMetadataValue("dAvgPeriodDiff", dt); //Time step per frame in ms
	Ext.getMetadataValue("dCalibration", dx); //Pixel size in microns
	dt /= 1000; //convert dt unit to seconds - from some reason the line dt = dt / 1000; does not work with the bioformats ext code

	//If the nd2 file is a single channel time-lapse, then analyze it
	if(nFrames > 1 && sizeZ == 1 && sizeC == 1){
		Ext.openImagePlus(dir + nd2Array[a]);

		//Down sample time if desired
		run("Reduce...", "reduction=" + reduce);
		nFrames = nSlices;
		dt = dt * reduce;

		//Crop the analysis to only the cell
		cropToCell(nd2Array[a]);

		//Process the image for analysis
		processStack(nd2Array[a], 3, 0.5); //image, high-pass value, minimum area in sq. microns

		//Measure the velocity of each object in the stack per frame
		measureVelocity(nd2Array[a], nFrames, dx, dt);

		//Save the results
		resultName = replace(nd2Array[a], ".nd2$", " - Velocity Histogram.csv");
		mapName = replace(nd2Array[a], ".nd2$", " - Velocity Map.tif");
		rawVel = replace(nd2Array[a], ".nd2$", " - Raw Velocities.tif");
		saveAs("Results", outDir + resultName);
		selectWindow("1");
		saveAs("Tiff", outDir + rawVel);
		close(rawVel);
		if(map){
			selectWindow("Stack");
			saveAs("Tiff", outDir + mapName);
			close(mapName);
		}
	}
	close("*");
}
setBatchMode("exit and display");

function processStack(i, hp, minArea){
	//De-noise stack
	run("Median...", "radius=1 stack");
	
	//Run high-pass filter to remove background and out-of-focus bacteria
	selectWindow(i);
	run("Duplicate...", "title=1 duplicate");
	selectWindow("1");
	run("Gaussian Blur...", "sigma=" + hp + " stack");
	run("Gaussian Blur...", "sigma=1 stack");
	imageCalculator("Subtract stack", i,"1");
	close("1");

	//Create mask
	setBatchMode("show");
	setAutoThreshold("Huang dark stack");
	run("Threshold...");
	waitForUser("Please set threshold");
	setBatchMode("hide");

	//Remove objects smaller than cutoff
	run("Analyze Particles...", "size=" + minArea + "-Infinity show=Masks clear stack");
	close(i);
	selectWindow("Mask of " + i);
	rename(i);

	//Create distance map of mask
	run("Duplicate...", "title=1 duplicate");
	selectWindow("1");
	run("Invert", "stack");
	run("Distance Map", "stack");

	
	//Create dx/dt stack
	selectWindow("1");
	run("Reverse");
	setSlice(nSlices);
	run("Add Slice");
	run("Reverse");
	run("Add...", "value=1 stack"); //Add 1 to distance so velocities of 0 can also be recorded

	selectWindow(i);
	setSlice(nSlices);
	run("Add Slice");
	run("Divide...", "value=255 stack");

	imageCalculator("Multiply stack", i,"1");
	close("1");
}

function cropToCell(i){
	selectWindow(i);
	getDimensions(width, height, channels, slices, frames);
	newImage("1", "8-bit black", width, height, 1);
	selectWindow(i);
	setBatchMode("show");
	select = true;
	while(select){
		setTool("polygon");
		waitForUser("Draw selection around cell");
		if(selectionType() == 2) select = false;
		else showMessage("Error: Selection needs to be a polygon.");
	}
	setBatchMode("hide");
	selectWindow("1");
	run("Restore Selection");
	setForegroundColor(255, 255, 255);
	run("Fill", "slice");
	saveAs("Tiff", outDir + replace(i, ".nd2$", " - selection.tif"));
	selectWindow(i);
	run("Crop");
	
}

function measureVelocity(i, slices, dx, dt){
	selectWindow(i);
	setThreshold(1, 255);
	width = 1;
	newImage("1", "32-bit black", width, slices, 1); //Create matrix to store velocities
	countArray = newArray(slices); //Create an array to store the number of objects in each slices
	
	//Measure the displacement of each object
	for(a=1; a<=slices; a++){
		showProgress((a-1),slices);
		selectWindow(i);
		setSlice(a);
		run("Analyze Particles...", "  show=[Count Masks] display  clear slice");

		//Store the results in the corresponding matrix
		selectWindow("1");
		countArray[a-1] = nResults;
		if(width < nResults){  //Adjust output matrix size if needed
			run("Canvas Size...", "width=" + nResults + " height=" + slices + " position=Center-Left zero");
			width = nResults;
		}
		for(b=0; b<nResults; b++){
			vel = (getResult("Max",b) - 1)*dx/dt;
			setPixel(b,a-1,vel);	
		}

		//Create mask of pixel displacements
		selectWindow("Count Masks of " + i);
		if(map){
			for(b=0; b<nResults; b++){
				max = getResult("Max",b);
				changeValues(b+1, b+1, max);
			}
			run("Select None");
			if(a==1) rename("Stack");
			else run("Concatenate...", "  title=Stack image1=Stack image2=[Count Masks of " + i + "] image3=[-- None --]");
		}
		else close("Count Masks of " + i);
	}

	//Fill all empty points on matrix with -1 so that they are not part of histogram
	selectWindow("1");
	countSum = 0;
	for(a=0; a<slices; a++){
		for(b=countArray[a]; b<width; b++) setPixel(b,a,-1);
		countSum += countArray[a];
	}

	//Generate histogram of velocities
	maxBins = (nBins + 1)*dx/dt;
	halfBin = 0.5*dx/dt;
	getHistogram(values, counts, nBins, 0-halfBin, maxBins-halfBin); //Offset the histogram by a half-bin so pixels are counted in the correct bin
	run("Clear Results");
	for(a=0; a<nBins; a++){
		setResult("Velocity (um/s)", a, values[a] + halfBin);
		setResult("Percent", a, counts[a]/countSum);
	}
	
}
