input1 = "/Volumes/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/output/Mutant/H3K27me3-488/slide19/";
output1 = "/Volumes/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/output/Mutant/H3K27me3-488/pngs/";

input2 = "/Volumes/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/output/Mutant/H3K27me3-488/slide20/";
output2 = "/Volumes/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/output/Mutant/H3K27me3-488/pngs/";

input3 = "/Volumes/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/output/Tak2xTak1/H3K27me3_488/";
output3 = "/Volumes/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/output/Tak2xTak1/pngs/";

input4 = "/Volumes/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/output/Cam2ez23xCam1/H3K27me3-488/";
output4 = "/Volumes/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/output/Cam2ez23xCam1/pngs/";

input5 = "/Volumes/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/output/Cam2xCak1/H3K27me3-488/";
output5 = "/Volumes/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/output/Cam2xCam1/pngs/";

input6 = "/Volumes/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/output/Tak2xCam1/H3K27me3-488/";
output6 = "/Volumes/groups/berger/user/sean.montgomery/Documents/Results/NI immuno/output/Tak2xCam1/pngs/";

function action(input, output, filename) {
		//run("Bio-Formats Importer", "open=["+ input + filename +"] autoscale color_mode=Custom view=Hyperstack stack_order=XYCZT series_0_channel_0_red=0 series_0_channel_0_green=255 series_0_channel_0_blue=0 series_0_channel_1_red=0 series_0_channel_1_green=255 series_0_channel_1_blue=255");
		run("Bio-Formats Importer", "open=["+ input + filename +"] autoscale color_mode=Custom view=Hyperstack stack_order=XYCZT series_0_channel_0_red=255 series_0_channel_0_green=0 series_0_channel_0_blue=0 series_0_channel_1_red=0 series_0_channel_1_green=0 series_0_channel_1_blue=255");
		//Maximum intensity projection
		run("Z Project...", "projection=[Max Intensity]");
		//Make and save merged image
		run("Make Composite");
		run("Scale Bar...", "width=1 height=12 font=72 color=White background=None location=[Lower Right] bold overlay");
		//saveAs("Jpeg", output + filename + "merged.jpg");
		saveAs("PNG", output + filename + "merged.png");
		//Make and save single channel images
		run("Split Channels");
		run("Scale Bar...", "width=1 height=0 font=0 color=Black background=None location=[Lower Right] bold overlay");
		//saveAs("Jpeg", output + filename + "dapi.jpg");
		saveAs("PNG", output + filename + "dapi.png");
		close();
		run("Scale Bar...", "width=1 height=0 font=0 color=Black background=None location=[Lower Right] bold overlay");
		//saveAs("Jpeg", output + filename + "immuno.jpg");
		saveAs("PNG", output + filename + "immuno.png");
		close();
		//Close remaining windows
		close(); 
}


setBatchMode(true); 
list = getFileList(input1);
for (i = 0; i < list.length; i++){
	if (matches(list[i],"(.*)ics") == 1) {
        action(input1, output1, list[i]);
	}
}
setBatchMode(false);

setBatchMode(true); 
list = getFileList(input2);
for (i = 0; i < list.length; i++){
	if (matches(list[i],"(.*)ics") == 1) {
        action(input2, output2, list[i]);
	}
}
setBatchMode(false);

setBatchMode(true); 
list = getFileList(input3);
for (i = 0; i < list.length; i++){
	if (matches(list[i],"(.*)ics") == 1) {
        action(input3, output3, list[i]);
	}
}
setBatchMode(false);

setBatchMode(true); 
list = getFileList(input4);
for (i = 0; i < list.length; i++){
	if (matches(list[i],"(.*)ics") == 1) {
        action(input4, output4, list[i]);
	}
}
setBatchMode(false);

setBatchMode(true); 
list = getFileList(input5);
for (i = 0; i < list.length; i++){
	if (matches(list[i],"(.*)ics") == 1) {
        action(input5, output5, list[i]);
	}
}
setBatchMode(false);

setBatchMode(true); 
list = getFileList(input6);
for (i = 0; i < list.length; i++){
	if (matches(list[i],"(.*)ics") == 1) {
        action(input6, output6, list[i]);
	}
}
setBatchMode(false);