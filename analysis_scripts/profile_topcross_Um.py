# Script for running velocity profiles in time. 
#
# Note: you must create destination folders before beginning! 
#

for num in range(1,682):

	OpenDatabase("localhost:/Users/Bosque/IBAMR/peri-gPC/newcode_runs/viz_IB2d"+str(num)+"/lag_data.visit", 0)
	#OpenDatabase("localhost:/Users/Spectre/Dropbox/peri-gPC/tess"+str(num).zfill(2)+"/viz_IB2d/lag_data.visit",0)
	AddPlot("Mesh", "heart_race_512_vertices", 1, 0)
	AddPlot("Mesh", "heart_tube_512_mesh", 1, 0)
	AddPlot("Mesh", "heart_tube_512_vertices", 1, 0)
	DrawPlots()
	OpenDatabase("localhost:/Users/Bosque/IBAMR/peri-gPC/newcode_runs/viz_IB2d"+str(num)+"/dumps.visit", 0)
	#OpenDatabase("localhost:/Users/Spectre/Dropbox/peri-gPC/tess01/viz_IB2d/dumps.visit", 0)
	HideActivePlots()
	AddPlot("Pseudocolor", "U_magnitude", 1, 0)
	SetActivePlots(4)
	SetActivePlots(4)
	DrawPlots()
	globalLineout = GetGlobalLineoutAttributes()
	globalLineout.Dynamic=1;
	globalLineout.curveOption=0;
	SetGlobalLineoutAttributes(globalLineout)
	print globalLineout
	SetTimeSliderState(0)
	SetQueryFloatFormat("%g")
	n = TimeSliderGetNStates()
	for i in range(TimeSliderGetNStates()):
		SetActiveWindow(1)
		TimeSliderNextState()
		Query("Lineout", end_point=(0, 0.201, 0), num_samples=50, start_point=(0, 0.099, 0), use_sampling=0)
		SetActiveWindow(2)
		SetActivePlots(0)
		SaveWindowAtts = SaveWindowAttributes()
		SaveWindowAtts.outputToCurrentDirectory = 0
		SaveWindowAtts.outputDirectory = "/Users/Bosque/IBAMR/peri-gPC/newcode_runs/newcode"+str(num)+"_profiles"
		SaveWindowAtts.fileName = "newcode"+str(num).zfill(3)+"_Um_profile"
#		SaveWindowAtts.fileName = "tess01_topcross_Um_profile"
		SaveWindowAtts.family = 1
		SaveWindowAtts.format = SaveWindowAtts.CURVE  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
		SaveWindowAtts.width = 1024
		SaveWindowAtts.height = 1024
		SaveWindowAtts.screenCapture = 0
		SaveWindowAtts.saveTiled = 0
		SaveWindowAtts.quality = 80
		SaveWindowAtts.progressive = 0
		SaveWindowAtts.binary = 0
		SaveWindowAtts.stereo = 0
		SaveWindowAtts.compression = SaveWindowAtts.PackBits  # None, PackBits, Jpeg, Deflate
		SaveWindowAtts.forceMerge = 0
		SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions  # NoConstraint, EqualWidthHeight, ScreenProportions
		SaveWindowAtts.advancedMultiWindowSave = 0
		SetSaveWindowAttributes(SaveWindowAtts)
		SaveWindow()
		DeleteActivePlots()
	SetActiveWindow(1)
	DeleteActivePlots()
	DeleteActivePlots()
	DeleteActivePlots()
	DeleteActivePlots()
	DeleteActivePlots()
	#OpenDatabase("localhost:/Users/Spectre/Dropbox/peri-gPC/tess"+str(num).zfill(2)+"/viz_IB2d/dumps.visit")
	#OpenDatabase("localhost:/Users/Spectre/Dropbox/peri-gPC/tess"+str(num).zfill(2)+"/viz_IB2d/lag_data.visit")


exit()
