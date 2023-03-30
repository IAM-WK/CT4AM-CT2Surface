/*************************************************************************************************
**************************************************************************************************
 *  CT2Surface
 *  Analyse CT stack for surface parameters
 *	An STL is needed to define the surfaces to be used for analysis
 *	Additionally subsurface parameters can be found and output of height maps can be performed for
 *	external analysis
 *
 *
 *  Copyright 2023 Lukas Englert
 *
 *
 *  Institute of Materials Science and Engineering, <http://www.iam.kit.edu/wk/english/index.php>
 *  Karlsruhe Institute of Technology
 *
 *  If you intend to use this work for your scientific publication please cite the
 *  appropriate publications listed on https://sourceforge.net/projects/ctfam/
 *
 **************************************************************************************************
 **************************************************************************************************/

#include "Misc.h"
#include "Geom.h"
#include "Surface.h"
#include "Subsurface.h"
#include "vtkWriter.h"
#include "progress.h"
 //#include "ShapeDetection.h"

 // boost libs
#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc, char* argv[])
{
	//**********************************************************************************************************************************************
	// INPUT
	//**********************************************************************************************************************************************
	std::string translationstring;
	std::string ImageFilename, MeshFilename, HeightMapFilename;

	std::vector<double> translationvec;
	itk::Vector<double, 3> BuildDirection;
	double MeasurementThickness, SubSurfaceDepth, SubDivisionSize, WindowDiameter;// , SamplingDiameter;
	unsigned int MinMeasurementAreaPixels;
	int MeasurementMargin, EdgeMargin;
	float MaxMergeAngle; double CreaseAngle;
	bool verboseOutput = false, descriptorOutput = true, Masking = true, ShapeRemoval = true, vtkOutput = true, plausibilityMasking = true;
	bool maskAir = false, writeStacks = false, writeMasks = false, writeShapes = false, writeASCII = false, writeHMs = false, writeSubSurfs = false, writeAirMask = false;
	bool slidingMode = false, HandleEdgeMarginsSeperately = false;
	//bool doClassification = false, doCleanup = true;
	bool sourceImageisBinary = true, detectSourceImageType = false, useSubDivision = false, RestrictToImageArea = false;
	std::string BuildDirectionString;
	std::size_t imageMode;

	try {
		// Declare the supported options
		po::options_description general("general");
		general.add_options()
			("help,h", "produce help message")
			("ImageFilename,i", po::value<std::string>(&ImageFilename)->required(), "filename of Image - .mhd! (.mhd should be in mm)")
			("MeshFilename,m", po::value<std::string>(&MeshFilename)->required(), "filename of Mesh - .stl! (should be in mm)")
			("HMFilename,o", po::value<std::string>(&HeightMapFilename)->required(), "filename of heightmaps - will be appended an index for each plane - .mhd")
			("translation,u", po::value< std::string >(&translationstring)->required(), "Vector to translate CT image coords to fit to stl file (as string): \"x y z\" (in mm)")
			;
		po::options_description AnalysisModes("LocalAnalysisModes");
		AnalysisModes.add_options()
			("SubDivisionSize,q", po::value<double>(&SubDivisionSize)->default_value(-1.0), "Size of analysis patches (mm^2); -1 = use max size; usage enables analysis subdivision")
			("WindowDiameter,w", po::value<double>(&WindowDiameter)->default_value(-1.0), "Diameter of analysis window (mm); usage enables sliding analysis window around each triangle")
			;
		po::options_description MeasurementOptions("MeasurementOptions");
		MeasurementOptions.add_options()
			("MeasurementThickness,t", po::value<double>(&MeasurementThickness)->default_value(0.5), "Thickness of projection stack for height measurement (in mm)")
			("SubSurfaceDepth,s", po::value<double>(&SubSurfaceDepth)->default_value(0.45), "Search depth for near surface porosity calculation (in mm)")
			("MaxMergeAngle,a", po::value<float>(&MaxMergeAngle)->default_value(5.f), "relative angle to pivot face of faces to be measured at once")
			("MeasurementMargin,z", po::value<int>(&MeasurementMargin)->default_value(0), "Margin for mask. Positive values = increase mask size, negative values = decrease size")
			("CreaseAngle,c", po::value<double>(&CreaseAngle), "Angle above which model edges are handled differently in masking")
			("EdgeMargin,e", po::value<int>(&EdgeMargin), "Margin for mask at model edges > crease angle. Positive values = increase mask size, negative values = decrease size")
			("MinMeasurementArea,p", po::value<unsigned int>(&MinMeasurementAreaPixels)->default_value(1000), "min area in voxels to be analysed")
			("RestrictToImageArea,r", po::value<bool>(&RestrictToImageArea)->default_value(true), "Only analyse triangles with midpoint inside CT image region.")
			;
		po::options_description StackProcessing("StackProcessing");
		StackProcessing.add_options()
			("MaskSurroundingAir,k", "Mask the air surrounding the object to not mistake it for porosity.")
			("noMasking", "do not mask out of plane values")
			("noPlausibilityMasking", "do not extend mask with non-plausible values")
			("noShapeRemoval", "do not remove shape")
			;
		po::options_description ImageOutput("ImageOutput");
		ImageOutput.add_options()
			("writeStacks", "output the projected image stacks")
			("writeASCII", "output ASCII file of height maps, usable in surface analysis software")
			("writeHMs", "output height maps as .mhd images")
			("writeSubSurfs", "output subsurface maps as .mhd images")
			("writeShapes", "output shape maps as .mhd images")
			("writeMasks", "output mask maps as .mhd images")
			;
		po::options_description FileOutput("FileOutput");
		FileOutput.add_options()
			("noDescriptorOutput,n", "do not create descriptor csv output")
			("noVTKOutput", "do not create vtk output")
			;
		po::options_description Misc("Misc");
		Misc.add_options()
			("SourceImageMode,b", po::value<std::size_t>(&imageMode)->default_value(0), "Type of source image, 0 = binarized, 1 = gray value, 2 = auto detection")
			("BuildDirection,j", po::value<std::string>(&BuildDirectionString), "Build direction wrt to STL as comma seperated values, default: \"0,0,1\"")
			("verboseOutput,v", "create verbose console output")
			; //end add options
		po::options_description all("Allowed options");
		all.add(general).add(AnalysisModes).add(MeasurementOptions).add(StackProcessing).add(ImageOutput).add(FileOutput).add(Misc);
		po::variables_map vm; // declare vm object
		po::store(po::parse_command_line(argc, argv, all), vm); // let vm contain options from command line

		if (vm.count("help")) {
			std::cout << all << "\n";
			std::cout << "sample invocation: \n";
			std::cout << ".\\CT2Surface.exe -i Image.mhd -m mesh.stl -u \"-12.4 -12.45 -0.5\"  -o \"HeightMap.mhd\" --MaxMergeAngle 10 --MinMeasurementArea 3600 --verboseOutput\n";
			std::cout << "Prerequisites: \n";
			std::cout << " - rotate your CT Image to align to part orientation in stl file\n";
			std::cout << " - Image has to be .mhd format (and binarised, if -b 1 or 2 is not set)! \n";
			std::cout << " - set image to correct origin in .mhd or provide correct translation vector with -u option! \n";
			return 1;
		}
		// place to insert warnings aside from help
		po::notify(vm); // run notify after dealing with help, otherwise required options will throw exception

		if (vm.count("MaskSurroundingAir")) { maskAir = true; }
		if (vm.count("writeStacks")) { writeStacks = true; }
		if (vm.count("writeHMs")) { writeHMs = true; }
		if (vm.count("writeASCII")) { writeASCII = true; }
		if (vm.count("writeSubSurfs")) { writeSubSurfs = true; }
		if (vm.count("writeShapes")) { writeShapes = true; }
		if (vm.count("writeMasks")) { writeMasks = true; }
		if (vm.count("verboseOutput")) { verboseOutput = true; }
		if (vm.count("noDescriptorOutput")) { descriptorOutput = false; }
		if (vm.count("noVTKOutput")) { vtkOutput = false; }
		if (vm.count("noMasking")) { Masking = false; }
		if (vm.count("noPlausibilityMasking")) { plausibilityMasking = false; }
		if (vm.count("noShapeRemoval")) { ShapeRemoval = false; }
		if (vm.count("SubDivisionSize") && SubDivisionSize != -1.0) { useSubDivision = true; }
		if (vm.count("WindowDiameter") && WindowDiameter != -1.0) { slidingMode = true; }

		if (((vm.count("EdgeMargin")) || (vm.count("CreaseAngle"))) && EdgeMargin != MeasurementMargin) {
			HandleEdgeMarginsSeperately = true; descriptiveIndices::s_EdgeCalculationRequested = true;
		}
		else { // sane defaults...
			CreaseAngle = 45;
			EdgeMargin = MeasurementMargin;
		}
		descriptiveIndices::s_creaseAngle = CreaseAngle;

		if (useSubDivision && slidingMode) {
			std::cerr << "SubDivision and sliding analysis window are not meant to be used in combination! \n";
			std::cerr << "...exiting!" << std::endl;
		}

		if (MeasurementMargin > 0 && ShapeRemoval) {
			std::cerr << "Positive image margins are currently not supported with shape removal!\n";
			std::cerr << "The usage will lead to missing shape removal at the margins!\n";
			std::cerr << "Continuing the calculation...\n";
		}

		if (!Masking && ShapeRemoval) {
			std::cerr << "Shape removal without masking currently unsupported!\n";
		}

		if (imageMode == 0) {
			sourceImageisBinary = true;
		}
		else if (imageMode == 1) {
			sourceImageisBinary = false;
		}
		else if (imageMode == 2) {
			detectSourceImageType = true;
		}
		else {
			std::cerr << "SourceImageMode has to be set to 0 = binarized, 1 = gray value or 2 = auto detection!" << std::endl;
		}

		size_t start = 0, end = 0, tokens = 0; // token counting is necessary to catch last value
		while ((end = translationstring.find(' ', start)) != std::string::npos || (tokens < 3)) {
			translationvec.push_back(stod(translationstring.substr(start, end - start)));
			start = end + 1;
			++tokens;
		}
		if (!(ImageFilename.substr(ImageFilename.length() - 4, ImageFilename.length()) == ".mhd")) {
			ImageFilename = ImageFilename + ".mhd";
		}

		if (!vm["BuildDirection"].empty()) {
			size_t startBdir = 0, endBdir = 0;
			for (size_t i = 0; i < BuildDirection.size(); ++i) {
				if ((endBdir = BuildDirectionString.find(',', startBdir)) != std::string::npos) {
					BuildDirection[i] = stod(BuildDirectionString.substr(startBdir, endBdir - startBdir));
					startBdir = endBdir + 1;
				}
			}
		}
		else {
			BuildDirection[0] = 0.0;
			BuildDirection[1] = 0.0;
			BuildDirection[2] = 1.0;
		}
	}
	catch (std::exception& e)
	{
		std::cerr << "Error: " << e.what() << "\n";
		return EXIT_FAILURE;
	}
	// store command line arguments for documentation
	std::vector<std::string> argList(argv + 1, argv + argc);
	//**********************************************************************************************************************************************
	// END OF INPUT
	//**********************************************************************************************************************************************
	// Todos future
	//
	// - no overload for leveling (masked), just multiply with mask (as in subsurf)
	// 
	// - implement some volume parameters
	//		- abbott curve derived
	//
	// - improve the global analysis code path
	//	-> fix cleanupsmallgroups to avoid leaving unassigned triangles
	//
	// - use level set segmentation (sub pixel precision)
	//
	// - offer optional height measurement normal to shape, not in projection direction
	// 
	// - detect and develop simple shapes
	//		- zero gaussian: cones or cylinders could be detected as example
	//		- algorithm develops height maps in concise arithmetic to contiguous HM
	//		- similar developed elevation angle map can be correlated to the HM
	//
	// - detect shapes and save for each triangle to correlate surface and subsurface params 
	// 
	// - find projection height automatically -> area proportions, area fraction curve slope...
	//
	//**********************************************************************************************************************************************
	//**********************************************************************************************************************************************
	auto timeatstart = std::chrono::steady_clock::now();
	//**********************************************************************************************************************************************
	// read stack image and align it to defined position
	//**********************************************************************************************************************************************
	StackReaderType::Pointer StackReader = StackReaderType::New();
	StackReader->SetFileName(ImageFilename);
	try
	{
		StackReader->Update();
	}
	catch (itk::ExceptionObject& err)
	{
		std::cerr << "error reading image! " << err << std::endl;
		return 1;
	}
	std::cout << "Image reading finished!" << "\n";
	// in this "pre" version the registration will happen virtually - we can also include an automatic registration here
	//*********************************** align image (ITK) ****************************************************************************************
	// image alignment through ITK (changeinformationfilter)
	// translate origin according to values from user
	// this makes up for translation registration of image to STL
	// rotation alignment should be done externally!
	typedef itk::ChangeInformationImageFilter<StackImageType> ChangeInfoFilterType;
	ChangeInfoFilterType::Pointer ChangeInfofilter = ChangeInfoFilterType::New();
	ChangeInfofilter->SetInput(StackReader->GetOutput());
	StackImageType::PointType origin = StackReader->GetOutput()->GetOrigin();
	origin[0] += translationvec[0];
	origin[1] += translationvec[1];
	origin[2] += translationvec[2];
	ChangeInfofilter->SetOutputOrigin(origin);
	ChangeInfofilter->ChangeOriginOn();
	ChangeInfofilter->CenterImageOff();
	ChangeInfofilter->Update();
	std::cout << "Image translation finished!" << "\n";
	//**********************************************************************************************************************************************
	// end stack reading
	//**********************************************************************************************************************************************

	//**********************************************************************************************************************************************
	// read stl mesh and calculate normals to obtain surface orientations
	//**********************************************************************************************************************************************
	// register STLMeshIOFactory
	// itk 5.3 may provide simpler reader functs
	itk::STLMeshIOFactory::RegisterOneFactory();
	typedef itk::MeshFileReader < MeshType > MeshReaderType;
	MeshReaderType::Pointer MeshReader = MeshReaderType::New();
	MeshReader->SetFileName(MeshFilename);
	try
	{
		MeshReader->Update();
	}
	catch (itk::ExceptionObject& err)
	{
		std::cerr << "error reading stl! " << err << std::endl;
		return 1;
	}
	MeshType::Pointer STLMesh = MeshReader->GetOutput();
	std::cout << "Mesh reading finished!" << "\n";
	auto timeafterinput = std::chrono::steady_clock::now();
	if (verboseOutput) {
		std::cout << "Input reading finished " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(timeafterinput - timeatstart).count()) << "ms\n";
	}
	// calculate normals as itk doesn't seem to read the normals from the stl file
	auto normalcalculator = NormalFilterType::New();
	normalcalculator->SetInput(STLMesh);
	normalcalculator->SetWeight(weight_type);
	normalcalculator->Update();
	NormalQEMeshType::Pointer STLnormalmesh = normalcalculator->GetOutput();
	NormalQEMeshType::CellDataContainerPointer FaceNormals = STLnormalmesh->GetCellData();
	NormalQEMeshType::PointDataContainerPointer VertexNormals = STLnormalmesh->GetPointData();

	auto timeafternormals = std::chrono::steady_clock::now();

	//**********************************************************************************************************************************************
	// create FaceGroups by comparing normals
	//**********************************************************************************************************************************************
	// create datastructure to aggregate faces with equal (or similar) normals
	// FaceGroupsCellNums[GroupNum] = vector of cell nums with equal normals
	std::vector<std::tuple<std::vector<itk::IdentifierType>, std::vector<itk::IdentifierType>, itk::CovariantVector< double, StackMeshDimensionality>>> FaceGroupsCellNums;

	if (!slidingMode) {
		// aggregate Faces with Normals more similar than MaxMergeAngle
		aggregateFaces(STLMesh, FaceGroupsCellNums, FaceNormals, MaxMergeAngle);
	}
	else {
		// aggregate Faces in a WindowDiameter/2 radius
		aggregateFacesAroundCell(STLMesh, FaceGroupsCellNums, FaceNormals, MaxMergeAngle, WindowDiameter / 2.0);
	}
	std::cout << "Face identification finished!" << "\n";
	if (verboseOutput) {
		std::cout << "STL contains " << STLMesh->GetNumberOfCells() << " faces!\n";
		std::cout << "Number of grouped faces found: " << FaceGroupsCellNums.size() << "\n";
		auto timeafteraggregation = std::chrono::steady_clock::now();
		if (verboseOutput) {
			std::cout << "Time elapsed for face aggregation " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(timeafteraggregation - timeafternormals).count()) << "ms\n";
		}
	}
	if (useSubDivision) {
		std::cout << "Sub dividing elements... \n";
		subDivideLargeGroups(STLMesh, FaceGroupsCellNums, SubDivisionSize, ((MinMeasurementAreaPixels + 1) * ChangeInfofilter->GetOutput()->GetSpacing()[0] * ChangeInfofilter->GetOutput()->GetSpacing()[1]));
		if (verboseOutput) {
			std::cout << "Number of grouped faces after subdivision: " << FaceGroupsCellNums.size() << "\n";
			auto timeaftersubdiv = std::chrono::steady_clock::now();
			if (verboseOutput) {
				std::cout << "Time elapsed for face aggregation including subdiv" << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(timeaftersubdiv - timeafternormals).count()) << "ms\n";
			}
		}
	}

	// DEBUG
	//if (doCleanup) {
	//	std::cout << "cleaning up face groups ... ";
	//	// use empirical guess for clean up
	//	// 2 times min area
	//	const double MaxCleanupAreaPhysicalUnits = static_cast<double>(StackReader->GetOutput()->GetSpacing()[0] * StackReader->GetOutput()->GetSpacing()[1] * 2.0 * MinMeasurementArea);
	//	cleanUpSmallGroups(STLMesh, FaceGroupsCellNums, FaceNormals, MaxMergeAngle, MaxCleanupAreaPhysicalUnits);
	//	std::cout << "Finished!\n ";
	//	if (verboseOutput) {
	//		std::cout << "Number of grouped faces after clean up: " << FaceGroupsCellNums.size() << "\n";
	//	}
	//}
	// DEBUG

	//**********************************************************************************************************************************************
	// detect source image type if requested
	//**********************************************************************************************************************************************
	StackImageType::Pointer ImageStack = ChangeInfofilter->GetOutput();
	ImageStack->DisconnectPipeline();

	if (detectSourceImageType) {
		// check if image is a binary image or greyvalue image
		typedef itk::Statistics::ImageToHistogramFilter<StackImageType> HistoFilterType;
		HistoFilterType::Pointer HistoFilter = HistoFilterType::New();
		HistoFilter->SetInput(ImageStack);
		HistoFilter->SetAutoMinimumMaximum(true);
		HistoFilterType::HistogramSizeType HistoSize(1);
		HistoSize[0] = 255;
		HistoFilter->SetHistogramSize(HistoSize);
		HistoFilter->Update();
		if ((HistoFilter->GetOutput()->GetFrequency(0) + HistoFilter->GetOutput()->GetFrequency(254)) == HistoFilter->GetOutput()->GetTotalFrequency()) {
			std::cout << "Source Image is binary!\n";
			std::cout << "Continuing with NearestNeighbourInterpolation! \n";
			sourceImageisBinary = true;
		}
		else {
			std::cout << "Source image is not binarized! \n";
			std::cout << "Continuing with LinearInterpolation and Huang thresholding!  \n";
			sourceImageisBinary = false;
		}
	}

	//**********************************************************************************************************************************************
	// if not already, create a binarized copy
	// label the image stack
	// find the label which is outside air
	// use as a reference image later, to check if projected dark pixels are outside air or porosity
	//**********************************************************************************************************************************************

	auto timebeforeAirMask = std::chrono::steady_clock::now();
	StackImageType::Pointer LabelCopyTemp = StackImageType::New();
	LabelCopyTemp->Graft(ImageStack);
	LabelCopyTemp->DisconnectPipeline();
	LabelBinStackImageType::Pointer LabelCopy = LabelBinStackImageType::New();
	typedef itk::BinaryShapeKeepNObjectsImageFilter<LabelBinStackImageType> BinaryOpeningFilterType;
	BinaryOpeningFilterType::Pointer BinaryAirFilter = BinaryOpeningFilterType::New();
	LabelBinStackImageType::Pointer OutSideAirMask = LabelBinStackImageType::New();
	//**********************************************************************************************************************************************
	// do segmentation if image is not binarized
	//**********************************************************************************************************************************************
	if (maskAir) {
		// we need to threshold the present image
		typedef itk::HuangThresholdImageFilter<StackImageType, LabelBinStackImageType> HuangThresholdFilterType;
		HuangThresholdFilterType::Pointer HuangFilterStack = HuangThresholdFilterType::New();
		HuangFilterStack->SetAutoMinimumMaximum(true);
		HuangFilterStack->SetInsideValue(255);
		HuangFilterStack->SetOutsideValue(0);
		HuangFilterStack->SetInput(LabelCopyTemp);
		HuangFilterStack->SetNumberOfHistogramBins(256);
		HuangFilterStack->UpdateLargestPossibleRegion();
		LabelCopy->Graft(HuangFilterStack->GetOutput());

		BinaryAirFilter->SetInput(LabelCopy);
		BinaryAirFilter->SetBackgroundValue(0);
		BinaryAirFilter->SetNumberOfObjects(1);
		BinaryAirFilter->SetAttribute(BinaryOpeningFilterType::LabelObjectType::NUMBER_OF_PIXELS_ON_BORDER);
		BinaryAirFilter->Update();
		OutSideAirMask->Graft(BinaryAirFilter->GetOutput());

		if (writeAirMask) {
			BinStackWriterType::Pointer AirMaskWriter = BinStackWriterType::New();
			AirMaskWriter->SetInput(OutSideAirMask);
			std::string outFileName = ImageFilename + "AirMask.mhd";
			AirMaskWriter->SetFileName(outFileName);
			try
			{
				AirMaskWriter->Update();
			}
			catch (itk::ExceptionObject& err)
			{
				std::cerr << err << std::endl;
				return 1;
			}

			if (verboseOutput) {
				std::cout << "Finished writing " << outFileName << "\n";
			}
		}
		auto timeafterAirMask = std::chrono::steady_clock::now();
		if (verboseOutput) {
			std::cout << "Time elapsed for creating surronding air mask " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(timeafterAirMask - timebeforeAirMask).count()) << "ms\n";
		}
	}

	//**********************************************************************************************************************************************
	// preparation for actual measurements
	//**********************************************************************************************************************************************
	typedef itk::ResampleImageFilter<StackImageType, StackImageType> ResampleFilterType;
	typedef itk::ResampleImageFilter<LabelBinStackImageType, LabelBinStackImageType> ResampleLabelFilterType;
	typedef itk::AffineTransform<double, StackMeshDimensionality> TransformType; // double means coordinates are calculated as doubles
	typedef itk::LinearInterpolateImageFunction<StackImageType, double> GVInterpolatorType;
	typedef itk::NearestNeighborInterpolateImageFunction<StackImageType, double> BinInterpolatorType;
	typedef itk::NearestNeighborInterpolateImageFunction<LabelBinStackImageType, double> LabelInterpolatorType;

	// each entry contains parameters of facegroup ID in array val [Sa, Sq, Sz,...]
	std::vector<std::array<double, 7> > SurfaceGroupsParameters;
	// each entry contains parameters of facegroup ID in array val [max porosity, clustersPerArea, max cluster size, avgClusterSize, ProjPorosity, VolPorosity]
	std::vector<std::array<double, 6> > SubSurfaceGroupsParameters;
	// each entry contains parameters of facegroup ID in array val [max porosity, clustersPerArea, max cluster size, avgClusterSize, ProjPorosity, VolPorosity]
	std::vector<std::array<double, 6> > SubSurfaceClosedGroupsParameters;
	// each entry contains parameters of facegroup ID in array val [max porosity, clustersPerArea, max cluster size, avgClusterSize, ProjPorosity, VolPorosity]
	std::vector<std::array<double, 6> > SubSurfaceOpenGroupsParameters;
	// each entry contains a state variable to indicate if face was characterised successfully and level sucess
	// [0] = validity (0 is valid, 1 is invalid, 2 is invalid bc out of image), [1] = levelSuccess  (0 is sucess, 1 is sucess after stretch and 10 is fail, -1 for invalid surfaces)
	std::vector<std::array<int, 2>> SurfaceValidity;
	// contains size of each facegroup
	std::vector<std::size_t> SizeofFaceGroups, HMlevelingSuccess;
	SizeofFaceGroups.resize(FaceGroupsCellNums.size(), 0);
	HMlevelingSuccess.resize(FaceGroupsCellNums.size(), 0);
	SurfaceGroupsParameters.resize(FaceGroupsCellNums.size(), { 0,0,0,0,0,0,0 });
	SubSurfaceGroupsParameters.resize(FaceGroupsCellNums.size(), { 0,0,0,0,0,0 });
	SubSurfaceClosedGroupsParameters.resize(FaceGroupsCellNums.size(), { 0,0,0,0,0,0 });
	SubSurfaceOpenGroupsParameters.resize(FaceGroupsCellNums.size(), { 0,0,0,0,0,0 });
	SurfaceValidity.resize(FaceGroupsCellNums.size(), { 0,0 });
	std::vector<std::shared_future<std::array<double, 7>>> pending_SurfParamComputations;
	pending_SurfParamComputations.resize(FaceGroupsCellNums.size(), std::shared_future<std::array<double, 7>>());
	std::vector<std::shared_future<std::array<double, 6>>> pending_SubSurfParamComputations;
	pending_SubSurfParamComputations.resize(FaceGroupsCellNums.size(), std::shared_future<std::array<double, 6>>());
	std::vector<std::shared_future<std::array<double, 6>>> pending_SubSurfClosedParamComputations;
	pending_SubSurfClosedParamComputations.resize(FaceGroupsCellNums.size(), std::shared_future<std::array<double, 6>>());
	std::vector<std::shared_future<std::array<double, 6>>> pending_SubSurfOpenParamComputations;
	pending_SubSurfOpenParamComputations.resize(FaceGroupsCellNums.size(), std::shared_future<std::array<double, 6>>());

	// precompute indices if not already done
	precomputeIndices(STLMesh, FaceNormals);
	precomputeAggregates(FaceGroupsCellNums);

	std::cout << "Starting Surface analysis!" << "\n";
	//**********************************************************************************************************************************************
	// loop through grouped faces and generate projected images, height maps, subsurface maps, roughness parameters...
	//**********************************************************************************************************************************************
	for (int curProjGroup = 0; curProjGroup < FaceGroupsCellNums.size(); ++curProjGroup) {
		//**********************************************************************************************************************************************
		// parametrize ResampleFilter
		//**********************************************************************************************************************************************
		// Instantiation of input variables seems to be necessary in each loop iteration, else weird behaviour of pipeline in writer
		ResampleFilterType::Pointer ResampleFilter = ResampleFilterType::New();
		TransformType::Pointer Transformation = TransformType::New();
		ResampleFilter->SetTransform(Transformation);
		BinInterpolatorType::Pointer BinInterpolator = BinInterpolatorType::New();
		GVInterpolatorType::Pointer GVInterpolator = GVInterpolatorType::New();
		if (sourceImageisBinary) {
			ResampleFilter->SetInterpolator(BinInterpolator);
		}
		else {
			ResampleFilter->SetInterpolator(GVInterpolator);
		}
		ResampleFilter->SetDefaultPixelValue(0); // explicit. If pixels are outside
		ResampleFilter->SetInput(ImageStack);
		// pixel spacing in millimeters along X and Y is needed - use spacing of input
		const itk::Vector<double, StackMeshDimensionality> ResampleSpacing = { ImageStack->GetSpacing() };
		ResampleFilter->SetOutputSpacing(ResampleSpacing);
		ResampleFilter->UpdateLargestPossibleRegion();

		// Instantiation of input variables seems to be necessary in each loop iteration, else weird behaviour of pipeline in writer
		//==========================================================================================================
		// Notes on Resampling
		//==========================================================================================================
		// transformations are defined from output to input
		// so we need to use inverse cosine directions
		StackImageType::DirectionType direction = calculateOrientedDirection(FaceGroupsCellNums, curProjGroup);

		ResampleFilter->SetOutputDirection(direction.GetInverse().as_matrix());
		itk::Vector<double, StackMeshDimensionality> facenormal;
		facenormal[0] = direction(2, 0);
		facenormal[1] = direction(2, 1);
		facenormal[2] = direction(2, 2);

		//**********************************************************************************************************************************************
		// define size and origin
		//**********************************************************************************************************************************************
		// first select size based on oriented bounding box about vertices of local face group
		// we first multiply each point with the DCM to get them in alignment with the future axis
		// then find the bounding box of the rotated point set - z thickness should be small, x and y extent will define the size
		// origin should be selected base on local face group, e.g. point with coordinates that are the most in the direction against the direction cosines 1 and 2

		float xmin = std::numeric_limits<float>::max(), xmax = -std::numeric_limits<float>::max(), ymin = std::numeric_limits<float>::max(), ymax = -std::numeric_limits<float>::max(), zmin = std::numeric_limits<float>::max(), zmax = -std::numeric_limits<float>::max();

		for (std::size_t CellIteratorinGroup = 0; CellIteratorinGroup < std::get<1>(FaceGroupsCellNums[curProjGroup]).size(); ++CellIteratorinGroup) {
			itk::IdentifierType cellindexinGroup = std::get<1>(FaceGroupsCellNums[curProjGroup])[CellIteratorinGroup];
			for (size_t i = 0; i < descriptiveIndices::s_cellVerticesPoints[cellindexinGroup].size(); ++i)
			{
				itk::Point<double, 3> ProjPoint = (direction * descriptiveIndices::s_cellVerticesPoints[cellindexinGroup][i]);
				if (ProjPoint[0] < xmin) {
					xmin = ProjPoint[0];
				}
				if (ProjPoint[0] > xmax) {
					xmax = ProjPoint[0];
				}

				if (ProjPoint[1] < ymin) {
					ymin = ProjPoint[1];
				}
				if (ProjPoint[1] > ymax) {
					ymax = ProjPoint[1];
				}

				if (ProjPoint[2] < zmin) {
					zmin = ProjPoint[2];
				}
				if (ProjPoint[2] > zmax) {
					zmax = ProjPoint[2];
				}
			}
		}
		// code was to find the bounding box -> use std::minmax_element later on
		unsigned int sizeX = static_cast<unsigned int>((static_cast<double>(xmax) - static_cast<double>(xmin)) / ImageStack->GetSpacing()[0]);
		unsigned int sizeY = static_cast<unsigned int>((static_cast<double>(ymax) - static_cast<double>(ymin)) / ImageStack->GetSpacing()[1]);
		// Size regularization
		// skip small areas
		SizeofFaceGroups[curProjGroup] = static_cast<itk::ImageIORegion::SizeValueType>(sizeX * sizeY);
		if ((SizeofFaceGroups[curProjGroup] < MinMeasurementAreaPixels)) {
			SurfaceValidity[curProjGroup] = { 1,-1 };
			continue;
		}
		// define size
		StackImageType::SizeType ResampleSize;
		ResampleSize[0] = sizeX;
		ResampleSize[1] = sizeY;
		ResampleSize[2] = (MeasurementThickness + SubSurfaceDepth) / StackReader->GetOutput()->GetSpacing()[2]; //CHANGED
		ResampleFilter->SetSize(ResampleSize);

		// define origin
		StackImageType::DirectionType InvDirection = direction.GetInverse().as_matrix();
		// Physical space coordinate of origin for X and Y
		// origin is translated about MeasurementThickness in projection direction
		itk::Vector<double, StackMeshDimensionality> ResampleOrigin;
		ResampleOrigin[0] = xmin;
		ResampleOrigin[1] = ymin;
		ResampleOrigin[2] = zmin;
		ResampleOrigin = InvDirection * ResampleOrigin;
		// shift origin after inversion so that one measurementthickness can be projected in both directions from STL face
		ResampleOrigin[0] = ResampleOrigin[0] - ((MeasurementThickness * 0.5) + SubSurfaceDepth) * facenormal[0]; //CHANGED
		ResampleOrigin[1] = ResampleOrigin[1] - ((MeasurementThickness * 0.5) + SubSurfaceDepth) * facenormal[1]; // 0.5 because we shift half the thickness
		ResampleOrigin[2] = ResampleOrigin[2] - ((MeasurementThickness * 0.5) + SubSurfaceDepth) * facenormal[2]; //
		ResampleFilter->SetOutputOrigin(ResampleOrigin); // origin = pixel whose indexes are all 0

		//**********************************************************************************************************************************************
		// check if curProjGroup should be excluded from analysis by bounding box
		//**********************************************************************************************************************************************
		if (RestrictToImageArea) {
			// the cell could still touch the CT Image, just not the midpoint!
			itk::Point<double, 3> pivotMid = descriptiveIndices::s_cellMidPoints[std::get<0>(FaceGroupsCellNums[curProjGroup])[0]]; // get<0> = pivot cell
			StackImageType::IndexType pivotMidAsIndex; //
			bool CellIsInside = ImageStack->TransformPhysicalPointToIndex(pivotMid, pivotMidAsIndex);
			if (!CellIsInside) {
				SurfaceValidity[curProjGroup] = { 2,-1 };
				continue;
			}
		}

		//**********************************************************************************************************************************************
		// actual resampling of image
		//**********************************************************************************************************************************************
		ResampleFilter->Update();

		StackImageType::Pointer projectedImageTemp = ResampleFilter->GetOutput();
		StackImageType::Pointer projectedImage = StackImageType::New();
		projectedImageTemp->DisconnectPipeline();

		//==========================================================================================================
		// Instantiate Label Resampling for Air identification
		//==========================================================================================================
		LabelBinStackImageType::Pointer projectedLabelStack = LabelBinStackImageType::New();
		ResampleLabelFilterType::Pointer ResampleLabelFilter = ResampleLabelFilterType::New();
		TransformType::Pointer LabelTransformation = TransformType::New();
		LabelInterpolatorType::Pointer LabelInterpolator = LabelInterpolatorType::New();
		if (maskAir) {
			ResampleLabelFilter->SetTransform(Transformation);
			ResampleLabelFilter->SetInterpolator(LabelInterpolator);
			ResampleLabelFilter->SetDefaultPixelValue(251); // 1 should be the label for surrounding air , 251 -> out of image domain
			ResampleLabelFilter->SetInput(OutSideAirMask);
			// pixel spacing in millimeters along X and Y is needed - use spacing of input
			const itk::Vector<double, StackMeshDimensionality> LabelResampleSpacing = { ImageStack->GetSpacing() };
			ResampleLabelFilter->SetOutputSpacing(LabelResampleSpacing);
			ResampleLabelFilter->UpdateLargestPossibleRegion();
			// resampling of labelImage
			ResampleLabelFilter->SetOutputDirection(direction.GetInverse().as_matrix());
			ResampleLabelFilter->SetSize(ResampleSize);
			ResampleLabelFilter->SetOutputOrigin(ResampleOrigin); // origin = pixel whose indexes are all 0
			ResampleLabelFilter->Update();
			projectedLabelStack->Graft(ResampleLabelFilter->GetOutput());
		}
		//**********************************************************************************************************************************************
		// do segmentation if image is not binarized
		//**********************************************************************************************************************************************
		if (!sourceImageisBinary) {
			// we need to threshold the present image
			// TODO: make this better with level set segmentation or similar
			// TODO: Make it possible to get threshold value from user input
			typedef itk::HuangThresholdImageFilter<StackImageType, StackImageType> HuangThresholdFilterType;
			HuangThresholdFilterType::Pointer HuangFilter = HuangThresholdFilterType::New();
			HuangFilter->SetInsideValue(0);
			HuangFilter->SetOutsideValue(255);
			HuangFilter->SetInput(projectedImageTemp);
			HuangFilter->SetNumberOfHistogramBins(256);
			HuangFilter->UpdateLargestPossibleRegion();
			projectedImage->Graft(HuangFilter->GetOutput());
		}
		else {
			projectedImage->Graft(projectedImageTemp);
		}
		//**********************************************************************************************************************************************
		// write projected stacks if requested
		//**********************************************************************************************************************************************
		if (writeStacks) {
			StackWriterType::Pointer writer = StackWriterType::New();
			writer->SetInput(projectedImage);
			std::string outFileName = "Stackfile_no" + std::to_string(curProjGroup) + "_Face.mhd";
			writer->SetFileName(outFileName);
			try
			{
				writer->Update();
			}
			catch (itk::ExceptionObject& err)
			{
				std::cerr << err << std::endl;
				return 1;
			}
		}

		//**********************************************************************************************************************************************
		// project stacks to height mappings
		//**********************************************************************************************************************************************
		HeightMapImageType::IndexType HMstartIndex = { 0,0 };
		HeightMapImageType::RegionType HMRegion;
		HeightMapImageType::SizeType HMsize;
		HMsize[0] = projectedImage->GetLargestPossibleRegion().GetSize()[0];
		HMsize[1] = projectedImage->GetLargestPossibleRegion().GetSize()[1];
		HMRegion.SetIndex(HMstartIndex);
		HMRegion.SetSize(HMsize);
		const double xSpacing = projectedImage->GetSpacing()[0];
		const double ySpacing = projectedImage->GetSpacing()[1];
		const double zSpacing = projectedImage->GetSpacing()[2];
		// set HeightMap to region with params
		HeightMapImageType::Pointer HMImageMap = HeightMapImageType::New();
		HMImageMap->SetRegions(HMRegion);
		HMImageMap->Allocate();
		HMImageMap->FillBuffer(itk::NumericTraits<HeightMapPixelType>::Zero); // catch as signal for non-valid height values later. Min realistic value == subsurf height
		// set spacing for mhd output
		itk::Vector<itk::ImageBase<2>::SpacingValueType, 2> spacing;
		spacing[0] = xSpacing;
		spacing[1] = ySpacing;
		HMImageMap->SetSpacing(spacing);

		// doing the trick
		auto Stack2HMThread = std::thread(Stack2HeightMap, std::ref(HMImageMap), std::ref(projectedImage), SubSurfaceDepth);

		//**********************************************************************************************************************************************
		// generate Mask
		//**********************************************************************************************************************************************
		MaskImageType::Pointer HMMask = MaskImageType::New();
		MaskImageType::RegionType MaskRegion;
		MaskImageType::SizeType MaskSize = HMsize;
		itk::Vector<double, 3> ProjectedOrigin = direction * ResampleOrigin;
		MaskImageType::IndexType MaskstartIndex = { static_cast<itk::IndexValueType>(ProjectedOrigin[0] / spacing[0]), static_cast<itk::IndexValueType>(ProjectedOrigin[1] / spacing[1]) };
		MaskRegion.SetIndex(MaskstartIndex);
		MaskRegion.SetSize(MaskSize);
		HMMask->SetRegions(MaskRegion);
		HMMask->SetSpacing(spacing);
		HMMask->Allocate();
		HMMask->FillBuffer(itk::NumericTraits<MaskPixelType>::Zero);
		// Initialize valid region with size of unmasked HeightMap. This way value can be used either with or without masking -> for example in subsurfcomputation
		unsigned int MaskValidAreaPixels = HMsize[0] * HMsize[1];// HMImageMap->GetLargestPossibleRegion().GetNumberOfPixels();
		if (Masking) {
			// generate Mask from geometric information
			generateMask(HMMask, STLMesh, FaceGroupsCellNums, curProjGroup, direction);

			//**********************************************************************************************************************************************
			// remove non plausible points from mask
			//**********************************************************************************************************************************************
			if (plausibilityMasking && !maskAir) {
				if (Stack2HMThread.joinable()) {
					Stack2HMThread.join();
				}
				// measured height == maxheight/minheight to mask means non plausible
				extendMaskbyPlausibilityHeight(HMMask, HMImageMap, MeasurementThickness, SubSurfaceDepth, zSpacing);
			}
			if (plausibilityMasking && maskAir) {
				// no transition from material to air or highest pt in material means non plausible
				extendMaskbyPlausibilityAir(HMMask, projectedLabelStack, SubSurfaceDepth, zSpacing);
			}

			int currentMargin = MeasurementMargin;
			// check if seperate margin for edges is to be used
			if (HandleEdgeMarginsSeperately && aggregatePredicates::s_FacegrouptouchesEdge[curProjGroup]) {
				currentMargin = EdgeMargin;
			}
			if (currentMargin > 0) {
				// use dilate on the mask
				itk::FlatStructuringElement<2>::RadiusType Radius;
				Radius.Fill(static_cast<itk::Size<2U>::SizeValueType>(currentMargin));
				itk::FlatStructuringElement<2> StructElement = itk::FlatStructuringElement<2>::Ball(Radius);
				typedef itk::BinaryDilateImageFilter<MaskImageType, MaskImageType, itk::FlatStructuringElement<2>> dilateFilterType;
				auto dilateFilter = dilateFilterType::New();
				dilateFilter->SetInput(HMMask);
				dilateFilter->SetKernel(StructElement);
				dilateFilter->SetForegroundValue(1); // Value to dilate
				dilateFilter->Update();
				HMMask = dilateFilter->GetOutput();
			}
			if (currentMargin < 0) {
				// use erode on the mask
				itk::FlatStructuringElement<2>::RadiusType Radius;
				Radius.Fill((-1) * static_cast<itk::Size<2U>::SizeValueType>(currentMargin));
				itk::FlatStructuringElement<2> StructElement = itk::FlatStructuringElement<2>::Ball(Radius);
				typedef itk::BinaryErodeImageFilter<MaskImageType, MaskImageType, itk::FlatStructuringElement<2>> erodeFilterType;
				auto erodeFilter = erodeFilterType::New();
				erodeFilter->SetInput(HMMask);
				erodeFilter->SetKernel(StructElement);
				erodeFilter->SetForegroundValue(1); // Value to erode
				erodeFilter->SetBoundaryToForeground(false);
				erodeFilter->Update();
				HMMask = erodeFilter->GetOutput();
			}

			if (Stack2HMThread.joinable()) {
				Stack2HMThread.join();
			}
			// mask out of image values for maskAir plausibility masking and no (plausibility) masking -> always wanted
			if ((plausibilityMasking && maskAir) || !plausibilityMasking || !Masking) {
				maskOutofImageValues(HMMask, HMImageMap, zSpacing);
			}

			// check if masking makes the image too small for MinMeasurementArea
			using StatisticsImageFilterType = itk::StatisticsImageFilter<MaskImageType>;
			auto statisticsImageFilter = StatisticsImageFilterType::New();
			statisticsImageFilter->SetInput(HMMask);
			statisticsImageFilter->Update();
			MaskValidAreaPixels = statisticsImageFilter->GetSum();
			if (MaskValidAreaPixels < MinMeasurementAreaPixels) {
				SizeofFaceGroups[curProjGroup] = MaskValidAreaPixels;
				SurfaceValidity[curProjGroup] = { 1,-1 };
				if (plausibilityMasking && maskAir) {
					// join the thread to avoid segfaults..
					//std::cout << "continue bc too small after erode " << std::endl;
					//if (Stack2HMThread.joinable()) {
					//	Stack2HMThread.join();
					//}
				}
				continue;
			}
		}
		else {
			if (Stack2HMThread.joinable()) {
				Stack2HMThread.join();
			}
		}

		if (writeMasks && Masking) {
			MaskwriterType::Pointer Maskwriter = MaskwriterType::New();
			Maskwriter->SetInput(HMMask);
			Maskwriter->SetFileName("Mask" + std::to_string(curProjGroup) + ".mhd");
			try
			{
				Maskwriter->Update();
			}
			catch (itk::ExceptionObject& err)
			{
				std::cerr << "Error while writing Maskwriter images!" << std::endl;
				std::cerr << err << std::endl;
				return 1;
			}
		}

		//**********************************************************************************************************************************************
		// generate subsurface map
		//**********************************************************************************************************************************************
		// typing, size etc is identical to HeightMap
		HeightMapImageType::Pointer SubSurfaceMapping = HeightMapImageType::New();
		SubSurfaceMapping->SetRegions(HMRegion);
		SubSurfaceMapping->SetSpacing(spacing);
		SubSurfaceMapping->Allocate();
		SubSurfaceMapping->FillBuffer(itk::NumericTraits<HeightMapPixelType>::Zero);
		// with AirMask: generate overall, closed, open porosity
		// essentially just needed with maskAir true, but simpler this way
		// will just waste memory
		HeightMapImageType::Pointer SubSurfaceMappingOpen = HeightMapImageType::New();
		SubSurfaceMappingOpen->SetRegions(HMRegion);
		SubSurfaceMappingOpen->SetSpacing(spacing);
		SubSurfaceMappingOpen->Allocate();
		SubSurfaceMappingOpen->FillBuffer(itk::NumericTraits<HeightMapPixelType>::Zero);
		HeightMapImageType::Pointer SubSurfaceMappingClosed = HeightMapImageType::New();
		SubSurfaceMappingClosed->SetRegions(HMRegion);
		SubSurfaceMappingClosed->SetSpacing(spacing);
		SubSurfaceMappingClosed->Allocate();
		SubSurfaceMappingClosed->FillBuffer(itk::NumericTraits<HeightMapPixelType>::Zero);
		if (maskAir) {
			Stack2SubSurface(HMImageMap, projectedImage, SubSurfaceMapping, SubSurfaceMappingClosed, SubSurfaceMappingOpen, SubSurfaceDepth, projectedLabelStack);
		}
		else {
			Stack2SubSurface(HMImageMap, projectedImage, SubSurfaceMapping, SubSurfaceDepth);
		}

		// apply masking if requested
		if (Masking) {
			MaskSubSurface(SubSurfaceMapping, HMMask);
			if (maskAir) {
				MaskSubSurface(SubSurfaceMappingClosed, HMMask);
				MaskSubSurface(SubSurfaceMappingOpen, HMMask);
			}
		}

		if (writeSubSurfs) {
			HMWriterType::Pointer SubSurfwriter = HMWriterType::New();
			SubSurfwriter->SetInput(SubSurfaceMapping);
			std::string SubSurfFileName = HeightMapFilename.substr(0, HeightMapFilename.length() - 4) + "_subsurf_" + std::to_string(curProjGroup) + ".mhd";
			SubSurfwriter->SetFileName(SubSurfFileName);
			try
			{
				SubSurfwriter->Update();
			}
			catch (itk::ExceptionObject& err)
			{
				std::cerr << "Error while writing SubSurfacePorosity images!" << std::endl;
				std::cerr << err << std::endl;
				return 1;
			}
			if (maskAir) {
				// also write closed pores
				HMWriterType::Pointer SubSurfClosedwriter = HMWriterType::New();
				SubSurfClosedwriter->SetInput(SubSurfaceMappingClosed);
				std::string SubSurfClosedFileName = HeightMapFilename.substr(0, HeightMapFilename.length() - 4) + "_subsurfclosed_" + std::to_string(curProjGroup) + ".mhd";
				SubSurfClosedwriter->SetFileName(SubSurfClosedFileName);
				try
				{
					SubSurfClosedwriter->Update();
				}
				catch (itk::ExceptionObject& err)
				{
					std::cerr << "Error while writing closed SubSurfacePorosity images!" << std::endl;
					std::cerr << err << std::endl;
					return 1;
				}
				// also write open pores
				HMWriterType::Pointer SubSurfOpenwriter = HMWriterType::New();
				SubSurfOpenwriter->SetInput(SubSurfaceMappingOpen);
				std::string SubSurfOpenFileName = HeightMapFilename.substr(0, HeightMapFilename.length() - 4) + "_subsurfopen_" + std::to_string(curProjGroup) + ".mhd";
				SubSurfOpenwriter->SetFileName(SubSurfOpenFileName);
				try
				{
					SubSurfOpenwriter->Update();
				}
				catch (itk::ExceptionObject& err)
				{
					std::cerr << "Error while writing open SubSurfacePorosity images!" << std::endl;
					std::cerr << err << std::endl;
					return 1;
				}
			}
		}
		//**********************************************************************************************************************************************
		// remove Shape
		//**********************************************************************************************************************************************
		// do after subsurf analysis, because otherwise these functions would be much more complicated
		// idea is similar to masking: triangles of FG are projected
		// then, interpolate height of the triangles
		// subtract height of each triangle pt
		HeightMapImageType::Pointer ShapeMap = HeightMapImageType::New();
		//MaskRegion.SetIndex(MaskstartIndex);
		//MaskRegion.SetSize(MaskSize);
		ShapeMap->SetRegions(MaskRegion);
		ShapeMap->SetSpacing(spacing);
		ShapeMap->Allocate();
		ShapeMap->FillBuffer(itk::NumericTraits<MaskPixelType>::Zero);

		if (Masking && ShapeRemoval) {
			removeShape(HMImageMap, ShapeMap, HMMask, STLMesh, FaceGroupsCellNums, curProjGroup, direction);
			// if unsupported combination -> reported above, too late down here
		}

		if (writeShapes && ShapeRemoval) {
			HMWriterType::Pointer Shapewriter = HMWriterType::New();
			Shapewriter->SetInput(ShapeMap);
			Shapewriter->SetFileName("shape" + std::to_string(curProjGroup) + ".mhd");
			try
			{
				Shapewriter->Update();
			}
			catch (itk::ExceptionObject& err)
			{
				std::cerr << "Error while writing Shapewriter images!" << std::endl;
				std::cerr << err << std::endl;
				return 1;
			}
		}

		//**********************************************************************************************************************************************
		// level height maps
		//**********************************************************************************************************************************************
		if (Masking) {
			HMlevelingSuccess[curProjGroup] = levelHeightMap(HMImageMap, HMMask);
		}
		else {
			HMlevelingSuccess[curProjGroup] = levelHeightMap(HMImageMap);
		}

		SurfaceValidity[curProjGroup] = { 0,static_cast<int>(HMlevelingSuccess[curProjGroup]) };

		if (writeASCII) {
			std::ofstream ASCIIfile;
			std::string HeightMapFilenameTXT = HeightMapFilename.substr(0, HeightMapFilename.length() - 4) + "_face_" + std::to_string(curProjGroup) + ".txt";
			ASCIIfile.open(HeightMapFilenameTXT);
			HeightMapConstIteratorTypeIndexed HMIterator(HMImageMap, HMImageMap->GetLargestPossibleRegion());
			HMIterator.GoToBegin();
			HeightMapImageType::IndexType oldIndex;
			oldIndex = HMIterator.GetIndex();
			while (!HMIterator.IsAtEnd()) {
				if (HMIterator.GetIndex()[1] != oldIndex[1]) {
					// reached new image line
					ASCIIfile << "\n";
				}
				ASCIIfile << HMIterator.Get() << "\t";
				oldIndex = HMIterator.GetIndex();
				++HMIterator;
			}
			ASCIIfile.close();
		}

		if (writeHMs) {
			HMWriterType::Pointer HMwriter = HMWriterType::New();
			HMwriter->SetInput(HMImageMap);
			std::string HeightMapFilenameMHD = HeightMapFilename.substr(0, HeightMapFilename.length() - 4) + "_face_" + std::to_string(curProjGroup) + ".mhd";
			HMwriter->SetFileName(HeightMapFilenameMHD);
			try
			{
				HMwriter->Update();
			}
			catch (itk::ExceptionObject& err)
			{
				std::cerr << "Error while writing HeightMap images!" << std::endl;
				std::cerr << err << std::endl;
				return 1;
			}
		}
		// Sa, Sq, Sz, Sp, Sv, Ssk, Sku
		if (Masking) {
			pending_SurfParamComputations[curProjGroup] = std::async(std::launch::async, computeHeightParamsMasked, HMImageMap, HMMask);
		}
		else {
			pending_SurfParamComputations[curProjGroup] = std::async(std::launch::async, computeHeightParams, HMImageMap);
		}
		// compute Subsurf parameters
		// max porosity, clustersPerArea, max cluster size, avgClusterSize, ProjPorosity, VolPorosity
		pending_SubSurfParamComputations[curProjGroup] = std::async(std::launch::async, computeSubSurfParams, SubSurfaceMapping, MaskValidAreaPixels, SubSurfaceDepth, zSpacing);
		pending_SubSurfClosedParamComputations[curProjGroup] = std::async(std::launch::async, computeSubSurfParams, SubSurfaceMappingClosed, MaskValidAreaPixels, SubSurfaceDepth, zSpacing);
		pending_SubSurfOpenParamComputations[curProjGroup] = std::async(std::launch::async, computeSubSurfParams, SubSurfaceMappingOpen, MaskValidAreaPixels, SubSurfaceDepth, zSpacing);

		progress::progressbar(curProjGroup, FaceGroupsCellNums.size());
	}
	progress::progressbar(FaceGroupsCellNums.size(), FaceGroupsCellNums.size()); // if last facegroups are all too small, progress could be <100%
	std::cout << "\nFinished Surface analysis!" << "\n";
	auto timeafterProjections = std::chrono::steady_clock::now();
	if (verboseOutput) {
		std::cout << "Time elapsed after projections " << std::to_string(static_cast<size_t>(std::round(std::chrono::duration_cast<std::chrono::milliseconds>(timeafterProjections - timeafternormals).count() / 1000.0))) << " sec\n";
	}
	std::cout << "\nFinishing Parameter computation ..." << "\n";
	for (int curProjGroup = 0; curProjGroup < FaceGroupsCellNums.size(); ++curProjGroup) {
		progress::progressbar(curProjGroup, FaceGroupsCellNums.size());
		// check if a value has been calculated (not being done if e.g. size was too small)
		if (pending_SurfParamComputations[curProjGroup].valid()) {
			SurfaceGroupsParameters[curProjGroup] = pending_SurfParamComputations[curProjGroup].get();
		}
		if (pending_SubSurfParamComputations[curProjGroup].valid()) {
			SubSurfaceGroupsParameters[curProjGroup] = pending_SubSurfParamComputations[curProjGroup].get();
		}
		if (pending_SubSurfClosedParamComputations[curProjGroup].valid()) {
			SubSurfaceClosedGroupsParameters[curProjGroup] = pending_SubSurfClosedParamComputations[curProjGroup].get();
		}
		if (pending_SubSurfOpenParamComputations[curProjGroup].valid()) {
			SubSurfaceOpenGroupsParameters[curProjGroup] = pending_SubSurfOpenParamComputations[curProjGroup].get();
		}
	}
	progress::progressbar(FaceGroupsCellNums.size(), FaceGroupsCellNums.size()); // set progress to 100%
	auto timeafterParameters = std::chrono::steady_clock::now();
	if (verboseOutput) {
		std::cout << "\nTime elapsed for finishing parameter computation " << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(timeafterParameters - timeafterProjections).count()) << " ms\n";
	}
	//**********************************************************************************************************************************************
	// write CSV descriptor file
	//**********************************************************************************************************************************************
	if (descriptorOutput) {
		std::ofstream CSVFile;
		std::string HeightMapDescriptorFilename = HeightMapFilename.substr(0, HeightMapFilename.length() - 4) + "_description_" + ".csv";
		CSVFile.open(HeightMapDescriptorFilename);
		// write header
		CSVFile << "Facegroup; Filename ; Inclination Angle ; levelingSuccess ; type of surface ; Sa ; Sq ; Sz; Sp ; Sv; Ssk ; Sku ; TotalavgPoreArea ; TotalmaxPoreArea; TotalPoresPerArea ; TotalVolSubSurfPorosity ; TotalProjSubSurfPorosity ; TotalmaxSubSurfPorosity ; ClosedavgPoreArea ; ClosedmaxPoreArea; ClosedPoresPerArea ; ClosedVolSubSurfPorosity ; ClosedProjSubSurfPorosity ; ClosedmaxSubSurfPorosity ; OpenavgPoreArea ; OpenmaxPoreArea; OpenPoresPerArea ; OpenVolSubSurfPorosity ; OpenProjSubSurfPorosity ; OpenmaxSubSurfPorosity ; PivotMidpointX; PivotMidpointY; PivotMidpointZ ; FaceIDsPivot ; FaceIDs " << "\n";

		for (size_t IDVar = 0; IDVar < FaceGroupsCellNums.size(); ++IDVar) {
			if ((SizeofFaceGroups[IDVar] < MinMeasurementAreaPixels) || (SurfaceValidity[IDVar][0] != 0)) {
				continue;
			}
			// calculate inclination angle
			double inclinationAngle = computeBuildingAngle(std::get<2>(FaceGroupsCellNums[IDVar]), BuildDirection);
			std::string SurfType;
			if (inclinationAngle <= 90.5 && inclinationAngle >= 89.5) {
				SurfType = "vertical";
			}
			else if (inclinationAngle <= 180.5 && inclinationAngle >= 179.5) {
				SurfType = "horizontal down";
			}
			else if (inclinationAngle <= 0.5 && inclinationAngle >= 0) {
				SurfType = "horizontal up";
			}
			else if (inclinationAngle < 90) {
				SurfType = "upskin";
			}
			else if (inclinationAngle > 90) {
				SurfType = "downskin";
			}
			else {
				SurfType = "ERROR";
			}
			// generate comma seperated string of faces in group
			std::string involvedFaces;
			for (const auto& facenum : std::get<1>(FaceGroupsCellNums[IDVar])) {
				involvedFaces += std::to_string(facenum) + ",";
			}
			// generate comma seperated string of faces in pivotgroup
			std::string involvedFacesPivot;
			itk::Point<double, 3> pivotMidPoint({ 0,0,0 });
			for (const auto& facenum : std::get<0>(FaceGroupsCellNums[IDVar])) {
				involvedFacesPivot += std::to_string(facenum) + ",";
				pivotMidPoint[0] += descriptiveIndices::s_cellMidPoints[facenum][0];
				pivotMidPoint[1] += descriptiveIndices::s_cellMidPoints[facenum][1];
				pivotMidPoint[2] += descriptiveIndices::s_cellMidPoints[facenum][2];
			}
			if (std::get<0>(FaceGroupsCellNums[IDVar]).size() != 0) {
				pivotMidPoint[0] = pivotMidPoint[0] / std::get<0>(FaceGroupsCellNums[IDVar]).size();
				pivotMidPoint[1] = pivotMidPoint[1] / std::get<0>(FaceGroupsCellNums[IDVar]).size();
				pivotMidPoint[2] = pivotMidPoint[2] / std::get<0>(FaceGroupsCellNums[IDVar]).size();
			}
			// filename regeneration
			std::string HeightMapFilenameperID = HeightMapFilename.substr(0, HeightMapFilename.length() - 4) + "_face_" + std::to_string(IDVar) + ".mhd";
			// output string
			CSVFile << std::to_string(IDVar) << ";" << HeightMapFilenameperID << ";" << inclinationAngle << ";" << HMlevelingSuccess[IDVar] << ";" << SurfType << ";" << SurfaceGroupsParameters[IDVar][0] << ";" << SurfaceGroupsParameters[IDVar][1] << ";" << SurfaceGroupsParameters[IDVar][2] << ";" << SurfaceGroupsParameters[IDVar][3] << ";" << SurfaceGroupsParameters[IDVar][4] << ";" << SurfaceGroupsParameters[IDVar][5] << ";" << SurfaceGroupsParameters[IDVar][6] << ";"
				<< SubSurfaceGroupsParameters[IDVar][3] << ";" << SubSurfaceGroupsParameters[IDVar][2] << ";" << SubSurfaceGroupsParameters[IDVar][1] << ";" << SubSurfaceGroupsParameters[IDVar][5] << ";" << SubSurfaceGroupsParameters[IDVar][4] << ";" << SubSurfaceGroupsParameters[IDVar][0] << ";"
				<< SubSurfaceClosedGroupsParameters[IDVar][3] << ";" << SubSurfaceClosedGroupsParameters[IDVar][2] << ";" << SubSurfaceClosedGroupsParameters[IDVar][1] << ";" << SubSurfaceClosedGroupsParameters[IDVar][5] << ";" << SubSurfaceClosedGroupsParameters[IDVar][4] << ";" << SubSurfaceClosedGroupsParameters[IDVar][0] << ";"
				<< SubSurfaceOpenGroupsParameters[IDVar][3] << ";" << SubSurfaceOpenGroupsParameters[IDVar][2] << ";" << SubSurfaceOpenGroupsParameters[IDVar][1] << ";" << SubSurfaceOpenGroupsParameters[IDVar][5] << ";" << SubSurfaceOpenGroupsParameters[IDVar][4] << ";" << SubSurfaceOpenGroupsParameters[IDVar][0] << ";"
				<< pivotMidPoint[0] << ";" << pivotMidPoint[1] << ";" << pivotMidPoint[2] << ";" << involvedFacesPivot << ";" << involvedFaces << "\n";
		}
		CSVFile.close();
		std::cout << "Finished writing .csv file!" << "\n";
	}

	if (vtkOutput) {
		if (slidingMode) {
			vtkWriter::s_slidingModeenabled = true;
		}
		vtkWriter::s_BuildDirection = BuildDirection;
		vtkWriter::prepareData(STLMesh, FaceGroupsCellNums, FaceNormals, SurfaceValidity, argList);
		const std::string vtkFilename = HeightMapFilename.substr(0, HeightMapFilename.length() - 4) + ".vtk";
		vtkWriter::writerWorker(vtkFilename, SurfaceGroupsParameters, SubSurfaceGroupsParameters, SubSurfaceClosedGroupsParameters, SubSurfaceOpenGroupsParameters);
	}
	return EXIT_SUCCESS;
}
