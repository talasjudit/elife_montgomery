input1 = "/Volumes/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/Archegonia_TH/H3K9me1-H3-DAPI/2daf/input/";
output1 = "/Volumes/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/Archegonia_TH/H3K9me1-H3-DAPI/2daf/output/";

input2 = "/Volumes/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/Archegonia_TH/H3K9me1-H3-DAPI/3daf/input/";
output2 = "/Volumes/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/Archegonia_TH/H3K9me1-H3-DAPI/3daf/output/";

input3 = "/Volumes/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/Archegonia_TH/H3K27me3-H3-DAPI/2daf/input/";
output3 = "/Volumes/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/Archegonia_TH/H3K27me3-H3-DAPI/2daf/output/";

input4 = "/Volumes/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/Archegonia_TH/H3K27me3-H3-DAPI/3daf/input/nostack/";
output4 = "/Volumes/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/Archegonia_TH/H3K27me3-H3-DAPI/3daf/output/";

input5 = "/Volumes/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/Archegonia_TH/H3K27me3-H3-DAPI/3daf/input/stack/";
output5 = "/Volumes/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/Archegonia_TH/H3K27me3-H3-DAPI/3daf/output/";

function actionK9(input, output, filename) {
		//run("Bio-Formats Importer", "open=["+ input + filename +"] autoscale color_mode=Custom view=Hyperstack stack_order=XYCZT series_0_channel_0_red=0 series_0_channel_0_green=255 series_0_channel_0_blue=0 series_0_channel_1_red=0 series_0_channel_1_green=255 series_0_channel_1_blue=255");
		run("Bio-Formats Importer", "open=["+ input + filename +"] autoscale color_mode=Custom view=Hyperstack stack_order=XYCZT series_0_channel_0_red=255 series_0_channel_0_green=0 series_0_channel_0_blue=0 series_0_channel_1_red=0 series_0_channel_1_green=255 series_0_channel_1_blue=0 series_0_channel_2_red=0 series_0_channel_2_green=0 series_0_channel_2_blue=255");
		//Maximum intensity projection
		run("Z Project...", "projection=[Max Intensity]");
		//Make and save merged image
		run("Make Composite");
		run("Scale Bar...", "width=10 height=12 font=72 color=White background=None location=[Lower Right] bold overlay");
		//saveAs("Jpeg", output + filename + "merged.jpg");
		saveAs("PNG", output + filename + "merged.png");
		//Make and save single channel images
		run("Split Channels");
		run("Scale Bar...", "width=10 height=0 font=0 color=Black background=None location=[Lower Right] bold overlay");
		//saveAs("Jpeg", output + filename + "dapi.jpg");
		saveAs("PNG", output + filename + "dapi.png");
		close();
		run("Scale Bar...", "width=10 height=0 font=0 color=Black background=None location=[Lower Right] bold overlay");
		//saveAs("Jpeg", output + filename + "immuno1.jpg");
		saveAs("PNG", output + filename + "immuno1.png");
		close();
		run("Scale Bar...", "width=10 height=0 font=0 color=Black background=None location=[Lower Right] bold overlay");
		//saveAs("Jpeg", output + filename + "immuno2.jpg");
		saveAs("PNG", output + filename + "immuno2.png");
		close();
		//Close remaining windows
		close(); 
}

function actionK27(input, output, filename) {
		//run("Bio-Formats Importer", "open=["+ input + filename +"] autoscale color_mode=Custom view=Hyperstack stack_order=XYCZT series_0_channel_0_red=0 series_0_channel_0_green=255 series_0_channel_0_blue=0 series_0_channel_1_red=0 series_0_channel_1_green=255 series_0_channel_1_blue=255");
		run("Bio-Formats Importer", "open=["+ input + filename +"] autoscale color_mode=Custom view=Hyperstack stack_order=XYCZT series_0_channel_0_red=0 series_0_channel_0_green=255 series_0_channel_0_blue=0 series_0_channel_1_red=255 series_0_channel_1_green=0 series_0_channel_1_blue=0 series_0_channel_2_red=0 series_0_channel_2_green=0 series_0_channel_2_blue=255");
		//Maximum intensity projection
		run("Z Project...", "projection=[Max Intensity]");
		//Make and save merged image
		run("Make Composite");
		run("Scale Bar...", "width=10 height=12 font=72 color=White background=None location=[Lower Right] bold overlay");
		//saveAs("Jpeg", output + filename + "merged.jpg");
		saveAs("PNG", output + filename + "merged.png");
		//Make and save single channel images
		run("Split Channels");
		run("Scale Bar...", "width=10 height=0 font=0 color=Black background=None location=[Lower Right] bold overlay");
		//saveAs("Jpeg", output + filename + "dapi.jpg");
		saveAs("PNG", output + filename + "dapi.png");
		close();
		run("Scale Bar...", "width=10 height=0 font=0 color=Black background=None location=[Lower Right] bold overlay");
		//saveAs("Jpeg", output + filename + "immuno1.jpg");
		saveAs("PNG", output + filename + "immuno1.png");
		close();
		run("Scale Bar...", "width=10 height=0 font=0 color=Black background=None location=[Lower Right] bold overlay");
		//saveAs("Jpeg", output + filename + "immuno2.jpg");
		saveAs("PNG", output + filename + "immuno2.png");
		close();
		//Close remaining windows
		close(); 
}

function actionnostack(input, output, filename) {
		//run("Bio-Formats Importer", "open=["+ input + filename +"] autoscale color_mode=Custom view=Hyperstack stack_order=XYCZT series_0_channel_0_red=0 series_0_channel_0_green=255 series_0_channel_0_blue=0 series_0_channel_1_red=0 series_0_channel_1_green=255 series_0_channel_1_blue=255");
		run("Bio-Formats Importer", "open=["+ input + filename +"] autoscale color_mode=Custom view=Hyperstack stack_order=XYCZT series_0_channel_0_red=255 series_0_channel_0_green=0 series_0_channel_0_blue=0 series_0_channel_1_red=0 series_0_channel_1_green=255 series_0_channel_1_blue=0 series_0_channel_2_red=0 series_0_channel_2_green=0 series_0_channel_2_blue=255");
		//Make and save merged image
		run("Make Composite");
		run("Scale Bar...", "width=10 height=12 font=72 color=White background=None location=[Lower Right] bold overlay");
		//saveAs("Jpeg", output + filename + "merged.jpg");
		saveAs("PNG", output + filename + "merged.png");
		//Make and save single channel images
		run("Split Channels");
		run("Scale Bar...", "width=10 height=0 font=0 color=Black background=None location=[Lower Right] bold overlay");
		//saveAs("Jpeg", output + filename + "dapi.jpg");
		saveAs("PNG", output + filename + "dapi.png");
		close();
		run("Scale Bar...", "width=10 height=0 font=0 color=Black background=None location=[Lower Right] bold overlay");
		//saveAs("Jpeg", output + filename + "immuno1.jpg");
		saveAs("PNG", output + filename + "immuno1.png");
		close();
		run("Scale Bar...", "width=10 height=0 font=0 color=Black background=None location=[Lower Right] bold overlay");
		//saveAs("Jpeg", output + filename + "immuno2.jpg");
		saveAs("PNG", output + filename + "immuno2.png");
		close();
}

setBatchMode(true); 
list = getFileList(input1);
for (i = 0; i < list.length; i++){
	if (matches(list[i],"(.*)czi") == 1) {
        actionK9(input1, output1, list[i]);
	}
}
setBatchMode(false);

setBatchMode(true); 
list = getFileList(input2);
for (i = 0; i < list.length; i++){
	if (matches(list[i],"(.*)czi") == 1) {
        actionK9(input2, output2, list[i]);
	}
}
setBatchMode(false);

setBatchMode(true); 
list = getFileList(input3);
for (i = 0; i < list.length; i++){
	if (matches(list[i],"(.*)czi") == 1) {
        actionK27(input3, output3, list[i]);
	}
}
setBatchMode(false);

setBatchMode(true); 
list = getFileList(input4);
for (i = 0; i < list.length; i++){
	if (matches(list[i],"(.*)czi") == 1) {
        actionnostack(input4, output4, list[i]);
	}
}
setBatchMode(false);

setBatchMode(true); 
list = getFileList(input5);
for (i = 0; i < list.length; i++){
	if (matches(list[i],"(.*)czi") == 1) {
        actionK9(input5, output5, list[i]);
	}
}
setBatchMode(false);