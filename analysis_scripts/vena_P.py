##### Top Cross, Um ####
#
# Instructions: Make the following edits:
# 
# Line 12: OpenDatabase: add file path to lag_data.visit
# Line 17: OpenDatabase: add file path to dumps.visit
# Line 87: File name: change file name for saving .curve file
#
# Then copy and paste text into VisIt > Controls > Command... window and hit "Execute"
#

for num in range(1,682):

	OpenDatabase("localhost:/Users/Bosque/IBAMR/peri-gPC/newcode_runs/viz_IB2d"+str(num)+"/lag_data.visit", 0)
	AddPlot("Mesh", "heart_race_512_vertices", 1, 0)
	AddPlot("Mesh", "heart_tube_512_mesh", 1, 0)
	AddPlot("Mesh", "heart_tube_512_vertices", 1, 0)
	DrawPlots()
	OpenDatabase("localhost:/Users/Bosque/IBAMR/peri-gPC/newcode_runs/viz_IB2d"+str(num)+"/dumps.visit", 0)
	HideActivePlots()
	AddPlot("Pseudocolor", "P", 1, 0)
	AddOperator("Box", 0)
	SetActivePlots(4)
	SetActivePlots(4)
	BoxAtts = BoxAttributes()
	BoxAtts.amount = BoxAtts.Some  # Some, All
	BoxAtts.minx = -0.23
	BoxAtts.maxx = -0.22
	BoxAtts.miny = -0.20
	BoxAtts.maxy = -0.10
	BoxAtts.minz = 0
	BoxAtts.maxz = 1
	BoxAtts.inverse = 0
	SetOperatorOptions(BoxAtts, 0)
	DrawPlots()
	SetActivePlots(4)
	SetActivePlots(4)
	QueryOverTimeAtts = GetQueryOverTimeAttributes()
	QueryOverTimeAtts.timeType = QueryOverTimeAtts.DTime  # Cycle, DTime, Timestep
	QueryOverTimeAtts.startTimeFlag = 0
	QueryOverTimeAtts.startTime = 0
	QueryOverTimeAtts.endTimeFlag = 0
	QueryOverTimeAtts.endTime = 1
	QueryOverTimeAtts.strideFlag = 0
	QueryOverTimeAtts.stride = 1
	QueryOverTimeAtts.createWindow = 1
	QueryOverTimeAtts.windowId = 2
	SetDefaultQueryOverTimeAttributes(QueryOverTimeAtts)
	QueryOverTimeAtts = GetQueryOverTimeAttributes()
	QueryOverTimeAtts.timeType = QueryOverTimeAtts.DTime  # Cycle, DTime, Timestep
	QueryOverTimeAtts.startTimeFlag = 0
	QueryOverTimeAtts.startTime = 0
	QueryOverTimeAtts.endTimeFlag = 0
	QueryOverTimeAtts.endTime = 1
	QueryOverTimeAtts.strideFlag = 0
	QueryOverTimeAtts.stride = 1
	QueryOverTimeAtts.createWindow = 1
	QueryOverTimeAtts.windowId = 2
	SetQueryOverTimeAttributes(QueryOverTimeAtts)
	SetQueryFloatFormat("%g")
	QueryOverTime("Max", end_time=50, start_time=0, stride=1, use_actual_data=1)
	QueryOverTime("Average Value", end_time=50, start_time=0, stride=1)
	SetActiveWindow(2)
	ViewCurveAtts = ViewCurveAttributes()
	ViewCurveAtts.domainCoords = (-0.2625, 2.7625)
	ViewCurveAtts.rangeCoords = (-0.0384302, 0.404432)
	ViewCurveAtts.viewportCoords = (0.2, 0.95, 0.15, 0.95)
	ViewCurveAtts.domainScale = ViewCurveAtts.LINEAR  # LINEAR, LOG
	ViewCurveAtts.rangeScale = ViewCurveAtts.LINEAR  # LINEAR, LOG
	SetViewCurve(ViewCurveAtts)
	# End spontaneous state
	
	# Begin spontaneous state
	ViewCurveAtts = ViewCurveAttributes()
	ViewCurveAtts.domainCoords = (-0.580125, 3.08013)
	ViewCurveAtts.rangeCoords = (-0.0849308, 0.450933)
	ViewCurveAtts.viewportCoords = (0.2, 0.95, 0.15, 0.95)
	ViewCurveAtts.domainScale = ViewCurveAtts.LINEAR  # LINEAR, LOG
	ViewCurveAtts.rangeScale = ViewCurveAtts.LINEAR  # LINEAR, LOG
	SetViewCurve(ViewCurveAtts)
	# End spontaneous state
	
	HideActivePlots()
	SetActivePlots((0, 1))
	SetActivePlots(2)
	#HideActivePlots()
	SaveWindowAtts = SaveWindowAttributes()
	SaveWindowAtts.outputToCurrentDirectory = 1
	SaveWindowAtts.outputDirectory = "/Users/Bosque/IBAMR/peri-gPC/newcode_runs/curvefiles"
	SaveWindowAtts.fileName = "newcode"+str(num).zfill(3)+"_vena_P_max"
	SaveWindowAtts.family = 0
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
	
	HideActivePlots()
	SetActivePlots((0, 0))
	SetActivePlots(1)
	#HideActivePlots()
	SaveWindowAtts = SaveWindowAttributes()
	SaveWindowAtts.outputToCurrentDirectory = 1
	SaveWindowAtts.outputDirectory = "/Users/Bosque/IBAMR/peri-gPC/newcode_runs/curvefiles"
	SaveWindowAtts.fileName = "newcode"+str(num).zfill(3)+"_vena_P_avg"
	SaveWindowAtts.family = 0
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

	
	DeleteAllPlots()
	SetActiveWindow(1)
	DeleteAllPlots()
	
exit()



