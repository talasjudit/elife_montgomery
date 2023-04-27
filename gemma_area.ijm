input = "/Volumes/groups/berger/user/sean.montgomery/Documents/Fertility/Gemmae/input/";

function action(input, filename) {
	open(input + filename);
	run("Split Channels");
	run("Gaussian Blur...", "sigma=4");
	setAutoThreshold("Default");
	//run("Threshold...");
	//setThreshold(0, 110);
	setOption("BlackBackground", false);
	run("Convert to Mask");
	run("Fill Holes");
	roiManager("Deselect");
	run("Analyze Particles...", "size=0.50-Infinity display exclude add");
	run("Close");
	run("Close");
	run("Close");
}

//setBatchMode(true); 
list = getFileList(input);
for (i = 0; i < list.length; i++){
	if (matches(list[i],"(.*)png") == 1) {
        action(input, list[i]);
	}
}
//setBatchMode(false);

//saveAs("Results", "/Volumes/groups/berger/user/sean.montgomery/Documents/Fertility/Gemmae/input/gemmae_areas_tmp.xls");