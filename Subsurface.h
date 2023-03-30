#pragma once
#include "Misc.h"
#include "itkMultiplyImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkConnectedComponentImageFilter.h"

// generate subsurface porosity - simple
void Stack2SubSurface(const HeightMapImageType::Pointer HeightMap, const StackImageType::Pointer projectedImage, HeightMapImageType::Pointer SubsurfaceMap, const double searchDepth) {
	// get Height of Surface from Heightmap
	// count dark pixels below height and calulate rel. porosity below surface
	// generate a double map for every point in the image with the porosity

	HeightMapConstIteratorTypeIndexed HMIterator(HeightMap, HeightMap->GetLargestPossibleRegion());
	const double xSpacing = projectedImage->GetSpacing()[0];
	const double ySpacing = projectedImage->GetSpacing()[1];
	const double zSpacing = projectedImage->GetSpacing()[2];

	const InputPixelType MaterialThres = 127; // Threshold to seperate material from porosity, image should be binarized at this point, with 255 being the material and 0 the air

	// if air is not masked, just count dark pixels
	HMIterator.GoToBegin();
	// iterate through image
	while (!HMIterator.IsAtEnd()) {
		size_t HeightValueinVoxels = static_cast<size_t>((HMIterator.Get()) / zSpacing); // add additional height in this deeper projection
		if (HeightValueinVoxels == 0) {
			SubsurfaceMap->SetPixel(HMIterator.GetIndex(), 0.0);
			++HMIterator;
			continue;
		}
		else if ((static_cast<int>(HeightValueinVoxels) - static_cast<int>((searchDepth / zSpacing))) < 0) {
			SubsurfaceMap->SetPixel(HMIterator.GetIndex(), -1); // should not be possible to happen
			++HMIterator;
			continue;
		}
		size_t SubsurfacePorosityCount = 0;
		size_t searchStartIndex = HeightValueinVoxels - static_cast<size_t>((searchDepth / zSpacing));
		// check subsurface porosity below this height value
		for (size_t height = searchStartIndex; height < HeightValueinVoxels; ++height) {
			itk::Index<3> IndexInProjectedImage;
			IndexInProjectedImage[0] = HMIterator.GetIndex()[0];
			IndexInProjectedImage[1] = HMIterator.GetIndex()[1];
			IndexInProjectedImage[2] = height;
			if (projectedImage->GetPixel(IndexInProjectedImage) < MaterialThres) {
				++SubsurfacePorosityCount;
			}
		}
		double SubsurfacePorosity = static_cast<double>(SubsurfacePorosityCount) / (searchDepth / zSpacing);
		SubsurfaceMap->SetPixel(HMIterator.GetIndex(), SubsurfacePorosity);
		++HMIterator;
	}
}

// generate subsurface porosity
// air masked case - we can calculate open, closed and overall porosity
// threading interferes with const pointers....
void Stack2SubSurface(const HeightMapImageType::Pointer HeightMap, const StackImageType::Pointer projectedImage, HeightMapImageType::Pointer SubsurfaceMap, HeightMapImageType::Pointer SubsurfaceMapClosed, HeightMapImageType::Pointer SubsurfaceMapOpen, const double searchDepth, const LabelBinStackImageType::Pointer projectedAirMask) {
	// get Height of Surface from Heightmap
	// count dark pixels below height and calulate rel. porosity below surface
	// generate a double map for every point in the image with the porosity

	HeightMapConstIteratorTypeIndexed HMIterator(HeightMap, HeightMap->GetLargestPossibleRegion());
	const double xSpacing = projectedImage->GetSpacing()[0];
	const double ySpacing = projectedImage->GetSpacing()[1];
	const double zSpacing = projectedImage->GetSpacing()[2];

	const InputPixelType MaterialThres = 127; // Threshold to seperate material from porosity, image should be binarized at this point, with 255 being the material and 0 the air

	HMIterator.GoToBegin();
	// iterate through image
	while (!HMIterator.IsAtEnd()) {
		size_t HeightValueinVoxels = static_cast<size_t>((HMIterator.Get()) / zSpacing); // add additional height in this deeper projection
		if (HeightValueinVoxels == 0) {
			SubsurfaceMap->SetPixel(HMIterator.GetIndex(), 0.0);
			++HMIterator;
			continue;
		}
		else if ((static_cast<int>(HeightValueinVoxels) - static_cast<int>((searchDepth / zSpacing))) < 0) {
			SubsurfaceMap->SetPixel(HMIterator.GetIndex(), -1); // should not be possible to happen
			++HMIterator;
			continue;
		}
		size_t SubsurfacePorosityCount = 0;
		size_t SubsurfaceClosedPorosityCount = 0;
		size_t SubsurfaceOpenPorosityCount = 0;
		size_t outOfImageValuesCount = 0;
		size_t searchStartIndex = HeightValueinVoxels - static_cast<size_t>((searchDepth / zSpacing));
		// check subsurface porosity below this height value
		for (size_t height = searchStartIndex; height < HeightValueinVoxels; ++height) {
			itk::Index<3> IndexInProjectedImage;
			IndexInProjectedImage[0] = HMIterator.GetIndex()[0];
			IndexInProjectedImage[1] = HMIterator.GetIndex()[1];
			IndexInProjectedImage[2] = height;
			if (projectedImage->GetPixel(IndexInProjectedImage) < MaterialThres && projectedAirMask->GetPixel(IndexInProjectedImage) != 251) {		// 251 -> out of image domain) {
				++SubsurfacePorosityCount;
				if (projectedAirMask->GetPixel(IndexInProjectedImage) == 0) {	// count as closed if not surrounding air
					++SubsurfaceClosedPorosityCount;
				}
				else {															// count as open if not surrounding air
					++SubsurfaceOpenPorosityCount;
				}
			}
			else if (projectedAirMask->GetPixel(IndexInProjectedImage) == 251) {
				++outOfImageValuesCount;
			}
		}
		double AvailPixelCount = (searchDepth / zSpacing) - outOfImageValuesCount;
		if (AvailPixelCount <= 0) {
			SubsurfaceMap->SetPixel(HMIterator.GetIndex(), 0.0);
			SubsurfaceMapClosed->SetPixel(HMIterator.GetIndex(), 0.0);
			SubsurfaceMapOpen->SetPixel(HMIterator.GetIndex(), 0.0);
		}
		else {
			double SubsurfacePorosity = static_cast<double>(SubsurfacePorosityCount) / (AvailPixelCount);
			double SubsurfacePorosityClosed = static_cast<double>(SubsurfaceClosedPorosityCount) / (AvailPixelCount);
			double SubsurfacePorosityOpen = static_cast<double>(SubsurfaceOpenPorosityCount) / (AvailPixelCount);
			SubsurfaceMap->SetPixel(HMIterator.GetIndex(), SubsurfacePorosity);
			SubsurfaceMapClosed->SetPixel(HMIterator.GetIndex(), SubsurfacePorosityClosed);
			SubsurfaceMapOpen->SetPixel(HMIterator.GetIndex(), SubsurfacePorosityOpen);
		}
		++HMIterator;
	}
}

// apply mask by multiplying subsurface with the mask
// we can either change the pointer of SSM by passing the reference or Graft the output
// graft the output for now
void MaskSubSurface(HeightMapImageType::Pointer SubsurfaceMap, const MaskImageType::Pointer Mask) {
	typedef itk::MultiplyImageFilter<HeightMapImageType, MaskImageType, HeightMapImageType> MaskMultiplyerType;
	auto Maskingfilter = MaskMultiplyerType::New();
	itk::ImageRegion<2> Region;
	Region = Mask->GetLargestPossibleRegion();
	SubsurfaceMap->SetRegions(Region);
	Maskingfilter->SetInput1(SubsurfaceMap);
	Maskingfilter->SetInput2(Mask);
	//Maskingfilter->SetReleaseDataFlag(true);
	//Maskingfilter->InPlaceOn();
	Maskingfilter->Update();
	//SubsurfaceMap = Maskingfilter->GetOutput();
	SubsurfaceMap->Graft(Maskingfilter->GetOutput());
}

std::array<double, 6> computeSubSurfParams(const HeightMapImageType::Pointer SubSurfMapping, const unsigned int MeasurementAreaPixels, const double SubSurfDepth, const double zSpacing) {
	std::array<double, 6> SubSurfParams;
	double maxPoro = 0.0;			// maximum subsurface porosity value found
	double clustersPerArea = 0.0;	// count of clusters found in subsurface
	double maxClusterSize = 0.0;	// max size of cluster found in subsurface
	double avgClusterSize = 0.0;	// avg size of cluster found in subsurface
	double poroRatio = 0.0;			// percentage of area in subsurface with porosity

	if (AreFloatsEqual(MeasurementAreaPixels, 0.0)) {
		std::cerr << "MeasurementAreainPixels was zero in SubSurfaceComputation! " << std::endl;
	}
	if (AreFloatsEqual(zSpacing, 0.0)) {
		std::cerr << "zSpacing was zero in SubSurfaceComputation! " << std::endl;
	}

	using LabelImageType = itk::Image<unsigned int, 2>;
	using ConnectedComponentImageFilterType = itk::ConnectedComponentImageFilter<HeightMapImageType, LabelImageType>;
	auto connectedComponentImageFilter = ConnectedComponentImageFilterType::New();
	connectedComponentImageFilter->SetInput(SubSurfMapping);
	connectedComponentImageFilter->SetBackgroundValue(0);
	connectedComponentImageFilter->Update();

	using LabelStatisticsImageFilterType = itk::LabelStatisticsImageFilter<HeightMapImageType, LabelImageType>;
	using LabelPixelType = LabelStatisticsImageFilterType::LabelPixelType;
	auto labelStatisticsImageFilter = LabelStatisticsImageFilterType::New();
	labelStatisticsImageFilter->SetLabelInput(connectedComponentImageFilter->GetOutput());
	labelStatisticsImageFilter->SetInput(SubSurfMapping);
	labelStatisticsImageFilter->Update();

	const double sizeofAreaPhysical = MeasurementAreaPixels * SubSurfMapping->GetSpacing()[0] * SubSurfMapping->GetSpacing()[1];

	clustersPerArea = static_cast<double>(labelStatisticsImageFilter->GetNumberOfLabels() - 1) / sizeofAreaPhysical; // -1 because of background
	for (auto vIt = labelStatisticsImageFilter->GetValidLabelValues().begin(); vIt != labelStatisticsImageFilter->GetValidLabelValues().end(); ++vIt)
	{
		if (labelStatisticsImageFilter->HasLabel(*vIt) && (*vIt) != 0) // label 0 is background
		{
			LabelPixelType labelValue = *vIt;
			if (maxPoro < labelStatisticsImageFilter->GetMaximum(labelValue)) {
				maxPoro = labelStatisticsImageFilter->GetMaximum(labelValue);
			}
			if (maxClusterSize < labelStatisticsImageFilter->GetCount(labelValue)) {
				maxClusterSize = labelStatisticsImageFilter->GetCount(labelValue);
			}
			avgClusterSize += labelStatisticsImageFilter->GetCount(labelValue);
		}
	}

	poroRatio = avgClusterSize / MeasurementAreaPixels;
	avgClusterSize = avgClusterSize / static_cast<double>(labelStatisticsImageFilter->GetNumberOfLabels());

	// calculation of volumetric porosity
	// count of defect voxels in depth  = percentage * depth
	// sum of percentage*depth over image / MeasurementAreaPixels*depth = volumetric porosity
	itk::ImageRegionConstIterator<HeightMapImageType> imageIterator(SubSurfMapping, SubSurfMapping->GetLargestPossibleRegion());
	const size_t SubSurfDepthVoxels = static_cast<size_t>(SubSurfDepth / zSpacing);
	size_t defectiveVoxels = 0;
	while (!imageIterator.IsAtEnd())
	{
		defectiveVoxels += static_cast<size_t>(imageIterator.Get() * SubSurfDepthVoxels);
		++imageIterator;
	}
	const size_t validVoxelsVolumetric = MeasurementAreaPixels * SubSurfDepthVoxels;
	double subSurfPorosityVolumetric = static_cast<double>(defectiveVoxels) / static_cast<double>(validVoxelsVolumetric);

	SubSurfParams[0] = maxPoro;
	SubSurfParams[1] = clustersPerArea;
	SubSurfParams[2] = maxClusterSize * SubSurfMapping->GetSpacing()[0] * SubSurfMapping->GetSpacing()[1];
	SubSurfParams[3] = avgClusterSize * SubSurfMapping->GetSpacing()[0] * SubSurfMapping->GetSpacing()[1];
	SubSurfParams[4] = poroRatio;
	SubSurfParams[5] = subSurfPorosityVolumetric;

	return SubSurfParams;
}
