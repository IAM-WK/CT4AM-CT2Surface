#pragma once
#include "Misc.h"
#include <CGAL/MP_Float.h>

// "magic value" - value deemed plausible for angular registration misfit
constexpr double plausibleFitAngle = 5.0;

// level height map based on arithmetic mean
// return value 0: successful normal leveling, 1: stretched leveling used, 2: only height correction
int levelHeightMap(HeightMapImageType::Pointer HeightMap, const MaskImageType::Pointer Mask) {
	int levelingMode = 0;
	// subtract the plane from the actual height map
	using ConstIteratorTypeIndexed = itk::ImageRegionConstIteratorWithIndex<itk::Image<double, 2 >>;
	using IteratorTypeIndexed = itk::ImageRegionIteratorWithIndex<itk::Image<double, 2 >>;
	using ConstIteratorType = itk::ImageRegionConstIterator<itk::Image<double, 2 >>;
	using IteratorType = itk::ImageRegionIterator<itk::Image<double, 2 >>;
	using MaskConstIteratorType = itk::ImageRegionConstIterator<MaskImageType>;

	ConstIteratorTypeIndexed HMCIndexediterator(HeightMap, HeightMap->GetLargestPossibleRegion());
	IteratorTypeIndexed HMIndexediterator(HeightMap, HeightMap->GetLargestPossibleRegion());
	ConstIteratorType HMCiterator(HeightMap, HeightMap->GetLargestPossibleRegion());
	IteratorType HMiterator(HeightMap, HeightMap->GetLargestPossibleRegion());
	MaskConstIteratorType MaskCiterator(Mask, Mask->GetLargestPossibleRegion());

	// subtract least squares fitted plane
	// we assume the plane is well conditioned for later subtraction, i.e. value defined anywhere on the HM set
	std::vector<Point3D> HMpts;
	size_t CountPixelsinMask = 0;
	// fit the plane
	Plane FittedPlane;
	HMCIndexediterator.GoToBegin();
	MaskCiterator.GoToBegin();
	while (!HMCIndexediterator.IsAtEnd()) {
		if (MaskCiterator.Get() == 1) {
			HMpts.push_back(Point3D(HMCIndexediterator.GetIndex()[0] * HeightMap->GetSpacing()[0], HMCIndexediterator.GetIndex()[1] * HeightMap->GetSpacing()[1], HMCIndexediterator.Get()));
			++CountPixelsinMask;
		}
		++HMCIndexediterator;
		++MaskCiterator;
	}
	// returned double is quality of fit, 0 = equal variance, 1 = null variance orthogonal to fit
	// idea: could be stored in a later version for distortion or Registration purposes

	CGAL::linear_least_squares_fitting_3(HMpts.begin(), HMpts.end(), FittedPlane, CGAL::Dimension_tag<0>());
	double fitangle = CGAL::approximate_angle(FittedPlane.orthogonal_vector(), { 0.0,0.0,1.0 });

	if (fitangle > 90.0) {
		fitangle = 180 - fitangle; // select smaller angle if (by chance) the large angle was selected
	}

	if (fitangle < plausibleFitAngle) {
		// do the actual plane subtraction
		// idea: do a reprojection with the lsq plane instead of subtraction
		HMIndexediterator.GoToBegin();
		MaskCiterator.GoToBegin();
		while (!HMIndexediterator.IsAtEnd()) {
			// subtract z coord of projection of index coords to fitted plane from heightmap
			if (MaskCiterator.Get() == 1) {
				HMIndexediterator.Set(HMIndexediterator.Get() - (FittedPlane.a() * HMIndexediterator.GetIndex()[0] * HeightMap->GetSpacing()[0] + FittedPlane.b() * HMIndexediterator.GetIndex()[1] * HeightMap->GetSpacing()[1] + FittedPlane.d()) / (-1.0 * FittedPlane.c()));
			}
			else {
				// we could set height to arbitrary value here for other measurement indices...
				HMIndexediterator.Set(0.0);
			}
			++HMIndexediterator;
			++MaskCiterator;
		}
	}
	else {
		// usually higher fit angles occur because of very small measurement patches (if registration is correct)
		// try fitting a double stretched Version
		// if this version manages to find a small fit angle, apply this angle, to get the points at least a bit better
		std::vector<Point3D> stretchedHMpts;
		for (auto& HMpoint : HMpts) {
			stretchedHMpts.push_back(Point3D(HMpoint[0] * 2, HMpoint[1] * 2, HMpoint[2]));
		}
		linear_least_squares_fitting_3(stretchedHMpts.begin(), stretchedHMpts.end(), FittedPlane, CGAL::Dimension_tag<0>());
		double stretchfitangle = CGAL::approximate_angle(FittedPlane.orthogonal_vector(), { 0.0,0.0,1.0 });
		if (stretchfitangle > 90.0) {
			stretchfitangle = 180 - stretchfitangle; // select smaller angle if (by chance) the large angle was selected
		}
		if (stretchfitangle < plausibleFitAngle) {
			levelingMode = 1;
			// subtract stretched plane
			// in theory, this should be too little subtracted, but at least correct the points a bit
			HMIndexediterator.GoToBegin();
			MaskCiterator.GoToBegin();
			while (!HMIndexediterator.IsAtEnd()) {
				// subtract z coord of projection of index coords to fitted plane from heightmap
				// resubstitute doubled coords
				if (MaskCiterator.Get() == 1) {
					HMIndexediterator.Set(HMIndexediterator.Get() - (FittedPlane.a() * HMIndexediterator.GetIndex()[0] * 2 * HeightMap->GetSpacing()[0] + FittedPlane.b() * HMIndexediterator.GetIndex()[1] * 2 * HeightMap->GetSpacing()[1] + FittedPlane.d()) / (-1.0 * FittedPlane.c()));
				}
				else {
					// we could set height to arbitrary value here for other measurement indices...
					// use min val?
					HMIndexediterator.Set(0.0);
				}
				++HMIndexediterator;
				++MaskCiterator;
			}
		}
	}

	// make sure height map is wrt to mid plane by subtracting avg value
	// normally this should be the case after lsqf
	// fallback if lsqf fails
	// find avg pixel val
	HMCiterator.GoToBegin();
	MaskCiterator.GoToBegin();
	double avgPixelVal = 0;
	while (!HMCiterator.IsAtEnd()) {
		if (MaskCiterator.Get() == 1) {
			avgPixelVal += HMCiterator.Get();
		}
		++HMCiterator;
		++MaskCiterator;
	}
	if (CountPixelsinMask != 0) {
		avgPixelVal = avgPixelVal / CountPixelsinMask;
	}
	// if avgPixel height error after lsqf is larger than 10% of pixelspacing, issue warning and subtract the value
	// if avgPixelVal > 10% pixelspacing -> this should only happen if leveling was unsuccessful
	if (avgPixelVal > (HeightMap->GetSpacing()[0] * 0.1)) {
		//std::cerr << "Height map leveling encountered problems with leveling! \n";
		HMiterator.GoToBegin();
		MaskCiterator.GoToBegin();
		while (!HMiterator.IsAtEnd()) {
			if (MaskCiterator.Get() == 1) {
				HMiterator.Set(HMiterator.Get() - avgPixelVal);
			}
			else {
				// we could set height to arbitrary value here for other measurement indices...
				HMiterator.Set(0.0);
			}
			++HMiterator;
			++MaskCiterator;
		}
		//std::cerr << "An attempt to correct the height map was undertaken! \n";
		levelingMode += 10;
	}
	return levelingMode;
}

// level height map based on arithmetic mean
// non mask overload
// return value 0: successful normal leveling, 1: stretched leveling used, 10: only height correction
int levelHeightMap(HeightMapImageType::Pointer HeightMap) {
	int levelingMode = 0;
	// subtract the plane from the actual height map
	using ConstIteratorTypeIndexed = itk::ImageRegionConstIteratorWithIndex<itk::Image<double, 2 >>;
	using ConstIteratorType = itk::ImageRegionConstIterator<itk::Image<double, 2 >>;
	using IteratorType = itk::ImageRegionIterator<itk::Image<double, 2 >>;

	ConstIteratorTypeIndexed HMCIndexediterator(HeightMap, HeightMap->GetLargestPossibleRegion());
	ConstIteratorType HMCiterator(HeightMap, HeightMap->GetLargestPossibleRegion());
	IteratorType HMiterator(HeightMap, HeightMap->GetLargestPossibleRegion());

	// subtract least squares fitted plane
	std::vector<Point3D> HMpts;
	HMCIndexediterator.GoToBegin();

	// fit the plane
	Plane FittedPlane;
	while (!HMCIndexediterator.IsAtEnd()) {
		HMpts.push_back(Point3D(HMCIndexediterator.GetIndex()[0] * HeightMap->GetSpacing()[0], HMCIndexediterator.GetIndex()[1] * HeightMap->GetSpacing()[1], HMCIndexediterator.Get()));
		++HMCIndexediterator;
	}
	linear_least_squares_fitting_3(HMpts.begin(), HMpts.end(), FittedPlane, CGAL::Dimension_tag<0>());

	double fitangle = CGAL::approximate_angle(FittedPlane.orthogonal_vector(), { 0,0,1 });
	if (fitangle > 90.0) {
		fitangle = 180 - fitangle; // select smaller angle if (by chance) the large angle was selected
	}

	if (fitangle < plausibleFitAngle) {
		HMCIndexediterator.GoToBegin();
		HMiterator.GoToBegin();
		while (!HMCIndexediterator.IsAtEnd()) {
			// subtract z coord of projection of index coords to fitted plane from heightmap
			HMiterator.Set(HMiterator.Get() - (FittedPlane.a() * HMCIndexediterator.GetIndex()[0] * HeightMap->GetSpacing()[0] + FittedPlane.b() * HMCIndexediterator.GetIndex()[1] * HeightMap->GetSpacing()[1] + FittedPlane.d()) / (-1.0 * FittedPlane.c()));
			++HMiterator;
			++HMCIndexediterator;
		}
	}
	else {
		// usually higher fit angles occur because of very small measurement patches (if registration is correct)
		// try fitting a double stretched Version
		// if this version manages to find a small fit angle, apply this angle, to get the points at least a bit better
		// if not, only height is corrected
		std::vector<Point3D> stretchedHMpts;
		for (auto& HMpoint : HMpts) {
			stretchedHMpts.push_back(Point3D(HMpoint[0] * 2, HMpoint[1] * 2, HMpoint[2]));
		}
		linear_least_squares_fitting_3(stretchedHMpts.begin(), stretchedHMpts.end(), FittedPlane, CGAL::Dimension_tag<0>());
		double stretchfitangle = CGAL::approximate_angle(FittedPlane.orthogonal_vector(), { 0,0,1 });
		if (stretchfitangle > 90.0) {
			stretchfitangle = 180 - stretchfitangle; // select smaller angle if (by chance) the large angle was selected
		}
		if (stretchfitangle < plausibleFitAngle) {
			levelingMode = 1;
			// subtract stretched plane
			HMCIndexediterator.GoToBegin();
			HMiterator.GoToBegin();
			while (!HMCIndexediterator.IsAtEnd()) {
				// subtract z coord of projection of index coords to fitted plane from heightmap
				HMiterator.Set(HMiterator.Get() - (FittedPlane.a() * HMCIndexediterator.GetIndex()[0] * 2 * HeightMap->GetSpacing()[0] + FittedPlane.b() * HMCIndexediterator.GetIndex()[1] * 2 * HeightMap->GetSpacing()[1] + FittedPlane.d()) / (-1.0 * FittedPlane.c()));
				++HMiterator;
				++HMCIndexediterator;
			}
		}
	}

	// make sure height map is wrt to mid plane by subtracting avg value
	HMCiterator.GoToBegin();
	double avgPixelVal = 0;
	while (!HMCiterator.IsAtEnd()) {
		avgPixelVal += HMCiterator.Get();

		++HMCiterator;
	}
	avgPixelVal = avgPixelVal / HeightMap->GetLargestPossibleRegion().GetNumberOfPixels();
	// if avgPixel height error after lsqf is larger than 10% of pixelspacing, issue warning and subtract the value
	// if avgPixelVal > 10% pixelspacing -> this should only happen if leveling was unsuccessful
	if (avgPixelVal > (HeightMap->GetSpacing()[0] * 0.1)) {
		//std::cerr << "Height map leveling encountered problems with leveling! \n";
		HMiterator.GoToBegin();

		while (!HMiterator.IsAtEnd()) {
			HMiterator.Set(HMiterator.Get() - avgPixelVal);

			++HMiterator;
		}
		//std::cerr << "An attempt to correct the height map was undertaken! \n";
		levelingMode += 10;
	}
	return levelingMode;
}

// compute Sa, Sq and Sz value of given Height Map
// ISO 25178-2 Height parameters
// input is assumed to be in mm, output will be given in µm
// only pixel positions inside the mask are used for parameters calculation
std::array<double, 7> computeHeightParamsMasked(const HeightMapImageType::Pointer HeightMap, const MaskImageType::Pointer Mask) {
	std::array<double, 7> HeightParams;
	using ConstIteratorType = itk::ImageRegionConstIterator<itk::Image<double, 2 >>;
	ConstIteratorType HMiterator(HeightMap, HeightMap->GetLargestPossibleRegion());
	using MaskConstIteratorType = itk::ImageRegionConstIterator<MaskImageType>;
	MaskConstIteratorType MaskCiterator(Mask, Mask->GetLargestPossibleRegion());

	CGAL::MP_Float Sa = CGAL::MP_Float(0.0);
	CGAL::MP_Float Sq = CGAL::MP_Float(0.0);
	double Sz = 0.0;
	double lowestPt = std::numeric_limits<double>::max();
	double highestPt = -std::numeric_limits<double>::max();
	CGAL::MP_Float SskPrecise = CGAL::MP_Float(0.0);
	CGAL::MP_Float SkuPrecise = CGAL::MP_Float(0.0);

	double HeightValue = 0.0;
	const double xSpacing = HeightMap->GetSpacing()[0];
	const double ySpacing = HeightMap->GetSpacing()[1];
	size_t CountPixelsinMask = 0;

	// accumulate abs value of every pixel
	HMiterator.GoToBegin();
	MaskCiterator.GoToBegin();
	while (!HMiterator.IsAtEnd()) {
		if (MaskCiterator.Get() == 1) {
			HeightValue = HMiterator.Get();
			Sa += CGAL::MP_Float(abs(HeightValue) * xSpacing * ySpacing);
			Sq += CGAL::MP_Float(HeightValue * HeightValue * xSpacing * ySpacing);
			if (HeightValue > highestPt) {
				highestPt = HeightValue;
			}
			if (HeightValue < lowestPt) {
				lowestPt = HeightValue;
			}
			SskPrecise += CGAL::MP_Float(HeightValue * HeightValue * HeightValue * xSpacing * ySpacing);
			SkuPrecise += CGAL::MP_Float(HeightValue * HeightValue * HeightValue * HeightValue * xSpacing * ySpacing);
			++CountPixelsinMask;
		}

		++MaskCiterator;
		++HMiterator;
	}
	// divide value by area
	// area = number of pixels * area of pixel
	double imageArea = 0;
	imageArea = CountPixelsinMask * xSpacing * ySpacing;

	// make sure imageArea != 0
	//if (imageArea != 0) {
	// avoid FP comparison with zero. imageArea should be larger or equal than one pixel area
	if (imageArea > (xSpacing * ySpacing)) {
		Sa = CGAL::approximate_division(Sa, imageArea);
		Sq = approximate_sqrt(CGAL::approximate_division(Sq, imageArea));
		Sz = highestPt - lowestPt;
		// might happen for cut positions in images
		//if (Sq != 0) {// avoid FP comparison with zero. imageArea should be larger or equal than one pixel area
		//if (Sq > (xSpacing * ySpacing)) {
		if (!AreFloatsEqual(CGAL::to_double(Sq), 0.0)) {
			SskPrecise = CGAL::approximate_division(SskPrecise, CGAL::MP_Float((Sq * Sq * Sq * imageArea)));
			SkuPrecise = CGAL::approximate_division(SkuPrecise, CGAL::MP_Float((Sq * Sq * Sq * Sq * imageArea)));
		}
		else {
			// sane default, if Sq is zero evaluated, Ssk and Sku should be too
			SskPrecise = CGAL::MP_Float(0.0);
			SkuPrecise = CGAL::MP_Float(0.0);
		}
	}

	HeightParams[0] = CGAL::to_double(Sa) * 1000;
	HeightParams[1] = CGAL::to_double(Sq) * 1000;
	HeightParams[2] = Sz * 1000;
	HeightParams[3] = highestPt * 1000;
	HeightParams[4] = abs(lowestPt) * 1000;
	HeightParams[5] = CGAL::to_double(SskPrecise);
	HeightParams[6] = CGAL::to_double(SkuPrecise);

	return HeightParams;
}

// no mask
std::array<double, 7> computeHeightParams(const HeightMapImageType::Pointer HeightMap) {
	std::array<double, 7> HeightParams;
	using ConstIteratorType = itk::ImageRegionConstIterator<itk::Image<double, 2 >>;
	ConstIteratorType HMiterator(HeightMap, HeightMap->GetLargestPossibleRegion());

	CGAL::MP_Float Sa = CGAL::MP_Float(0.0);
	CGAL::MP_Float Sq = CGAL::MP_Float(0.0);
	double Sz = 0.0;
	double lowestPt = std::numeric_limits<double>::max();
	double highestPt = -std::numeric_limits<double>::max();
	CGAL::MP_Float Ssk = CGAL::MP_Float(0.0);
	CGAL::MP_Float Sku = CGAL::MP_Float(0.0);

	double HeightValue = 0.0;
	const double xSpacing = HeightMap->GetSpacing()[0];
	const double ySpacing = HeightMap->GetSpacing()[1];
	// accumulate abs value of every pixel
	HMiterator.GoToBegin();
	while (!HMiterator.IsAtEnd()) {
		HeightValue = HMiterator.Get();
		Sa += CGAL::MP_Float(abs(HeightValue) * xSpacing * ySpacing);
		Sq += CGAL::MP_Float(HeightValue * HeightValue * xSpacing * ySpacing);
		if (HeightValue > highestPt) {
			highestPt = HeightValue;
		}
		if (HeightValue < lowestPt) {
			lowestPt = HeightValue;
		}
		Ssk += CGAL::MP_Float(HeightValue * HeightValue * HeightValue * xSpacing * ySpacing);
		Sku += CGAL::MP_Float(HeightValue * HeightValue * HeightValue * HeightValue * xSpacing * ySpacing);

		++HMiterator;
	}
	// divide value by area
	// area = number of pixels * area of pixel
	double imageArea = HeightMap->GetLargestPossibleRegion().GetNumberOfPixels() * xSpacing * ySpacing;

	// make sure imageArea != 0
	if (imageArea != 0) {
		Sa = CGAL::approximate_division(Sa, imageArea);
		Sq = approximate_sqrt(CGAL::approximate_division(Sq, imageArea));
		Sz = highestPt - lowestPt;
		// might happen for cut positions in images
		if (Sq != 0) {
			Ssk = CGAL::approximate_division(Ssk, CGAL::MP_Float((Sq * Sq * Sq * imageArea)));
			Sku = CGAL::approximate_division(Sku, CGAL::MP_Float((Sq * Sq * Sq * Sq * imageArea)));
		}
		else {
			// sane default, if Sq is zero evaluated, Ssk and Sku should be too
			Ssk = CGAL::MP_Float(0.0);
			Sku = CGAL::MP_Float(0.0);
		}
	}

	HeightParams[0] = CGAL::to_double(Sa) * 1000;
	HeightParams[1] = CGAL::to_double(Sq) * 1000;
	HeightParams[2] = Sz * 1000;
	HeightParams[3] = highestPt * 1000;
	HeightParams[4] = abs(lowestPt) * 1000;
	HeightParams[5] = CGAL::to_double(Ssk);
	HeightParams[6] = CGAL::to_double(Sku);
	return HeightParams;
}

// generate height maps
void Stack2HeightMap(HeightMapImageType::Pointer HeightMap, StackImageType::Pointer projectedImage, const double SubSurfDepth) {
	// create iterator to iterate through stack
	StackConstIteratorTypeIndexed StackIterator(projectedImage, projectedImage->GetLargestPossibleRegion());
	const double xSpacing = projectedImage->GetSpacing()[0];
	const double ySpacing = projectedImage->GetSpacing()[1];
	const double zSpacing = projectedImage->GetSpacing()[2];

	StackIterator.GoToBegin();
	itk::Index<3> StartIndex;
	StartIndex.Fill(0);
	StartIndex[2] = static_cast<size_t>(SubSurfDepth / zSpacing);
	StackIterator.SetIndex(StartIndex);
	// iterate through image
	while (!StackIterator.IsAtEnd()) {
		if (StackIterator.Get() != 0) {
			// pixel value is non-zero --> is object
			// check if higher value at this x,y has already been found
			if ((StackIterator.GetIndex()[2] * zSpacing) > HeightMap->GetPixel({ StackIterator.GetIndex()[0], StackIterator.GetIndex()[1] })) {
				// write z value at pos x,y to HeightMap
				HeightMap->SetPixel({ StackIterator.GetIndex()[0], StackIterator.GetIndex()[1] }, (StackIterator.GetIndex()[2] * zSpacing));
			}
		}
		++StackIterator;
	}
}

// generate mask for determination of valid surface pts
// mask to be filled, Mesh, FaceGroupVector, InvDirection matrix to calculate projection of triangles
void generateMask(MaskImageType::Pointer ImageMask, const MeshType::Pointer STLMesh, std::vector<std::tuple<std::vector<itk::IdentifierType>, std::vector<itk::IdentifierType>, itk::CovariantVector< double, 3>>>& FaceGroupsCellNums, const std::size_t curProjGroup, const StackImageType::DirectionType Direction) {
	// create mask
	// project each triangle from facegroup to HeightMap domain
	//
	std::vector<Triangle2D> TrianglesToCheck;
	TrianglesToCheck.reserve(std::get<1>(FaceGroupsCellNums[curProjGroup]).size());
	for (std::size_t CellIteratorinGroup = 0; CellIteratorinGroup < std::get<1>(FaceGroupsCellNums[curProjGroup]).size(); ++CellIteratorinGroup) {
		// get vertices of cell
		itk::IdentifierType cellindexinGroup = std::get<1>(FaceGroupsCellNums[curProjGroup])[CellIteratorinGroup];
		// project to image region
		itk::QuadEdgeMeshPoint<MeshCoordType, StackMeshDimensionality> a = Direction * descriptiveIndices::s_cellVerticesPoints[cellindexinGroup][0];
		itk::QuadEdgeMeshPoint<MeshCoordType, StackMeshDimensionality> b = Direction * descriptiveIndices::s_cellVerticesPoints[cellindexinGroup][1];
		itk::QuadEdgeMeshPoint<MeshCoordType, StackMeshDimensionality> c = Direction * descriptiveIndices::s_cellVerticesPoints[cellindexinGroup][2];

		// convert to CGAL vector of triangles
		Point2D PointA(static_cast<double>(a[0]), static_cast<double>(a[1]));
		Point2D PointB(static_cast<double>(b[0]), static_cast<double>(b[1]));
		Point2D PointC(static_cast<double>(c[0]), static_cast<double>(c[1]));
		TrianglesToCheck.push_back(Triangle2D(PointA, PointB, PointC));
	}

	size_t lastSuccessfullIndex = 0;
	// setup mask image iteration
	MaskImageIteratorTypeIndexed MaskIterator(ImageMask, ImageMask->GetLargestPossibleRegion());
	MaskIterator.GoToBegin();
	while (!MaskIterator.IsAtEnd()) {
		itk::QuadEdgeMeshPoint<MeshCoordType, HeightMapDimensionality> IndexPointinPhysicalCoords;
		ImageMask->TransformIndexToPhysicalPoint(MaskIterator.GetIndex(), IndexPointinPhysicalCoords);
		Point2D PointinMask(IndexPointinPhysicalCoords[0], IndexPointinPhysicalCoords[1]);
		bool PointIsInsideMesh = false;

		// starting at lastSuccessfullIndex should save a few CPU cycles
		// doesn't really happen in practice/release mode?
		// how can cycles be saved if point is not in any triangle?
		size_t numTriangles = TrianglesToCheck.size();
		for (size_t indexofTriangle = lastSuccessfullIndex; indexofTriangle < (numTriangles + lastSuccessfullIndex); ++indexofTriangle) {//for (const auto& Tri : TrianglesToCheck) {
			if (TrianglesToCheck[indexofTriangle % numTriangles].has_on_bounded_side(PointinMask)) {//if (Tri.has_on_bounded_side(PointinMask)) {
				PointIsInsideMesh = true;
				lastSuccessfullIndex = indexofTriangle % numTriangles;
				break; // if it is in one triangle, we can end our search
			}
		}

		if (PointIsInsideMesh) {
			MaskIterator.Set(1);
		}
		else {
			MaskIterator.Set(0);
		}
		++MaskIterator;
	}
}

// since maximum/minimum possible measured values are MeasurementThickness/2, this values can be found easily
// with this function, they are considered non plausible
void extendMaskbyPlausibilityHeight(MaskImageType::Pointer ImageMask, const HeightMapImageType::Pointer HeightMap, const double MeasurementThickness, const double SubSurfaceDepth, const double heightSpacing) {
	// highest and lowest possible values can deviate from "analytical" values because of the discretization of the images
	// calculate nearest possible value with zSpacing
	const double nonPlausibleValueHigh = MeasurementThickness + SubSurfaceDepth;
	const double nonPlausibleValueLow = SubSurfaceDepth;
	HeightMapConstIteratorTypeIndexed HMIterator(HeightMap, HeightMap->GetLargestPossibleRegion());
	MaskImageIteratorTypeIndexed MaskIterator(ImageMask, ImageMask->GetLargestPossibleRegion());
	HMIterator.GoToBegin();
	MaskIterator.GoToBegin();
	while (!HMIterator.IsAtEnd()) {
		if (MaskIterator.Get() == 0) {
			// no need to check already out of mask values...
			++HMIterator;
			++MaskIterator;
			continue;
		}
		// zero is also non plausible and can happen for out of image values
		if ((std::abs(HMIterator.Get() - nonPlausibleValueLow) < 2.0 * heightSpacing) || (std::abs(HMIterator.Get() - nonPlausibleValueHigh) < 2.0 * heightSpacing) || (std::abs(HMIterator.Get()) < heightSpacing)) {
			// set Mask value to 0
			// height map is just naively placed in physical domain
			// mask is placed in physical coordinates
			//MaskIterator.SetIndex(HMIterator.GetIndex(),0);
			MaskIterator.Set(0);
		}
		++HMIterator;
		++MaskIterator;
	}
}
// improved plausibility function -> did and air<->material transition happen? if not -> not plausible value!
void extendMaskbyPlausibilityAir(MaskImageType::Pointer ImageMask, const LabelBinStackImageType::Pointer projectedAirMask, const double SubSurfaceDepth, const double zSpacing) {
	// material to air transition has to happen in measurement thickness
	size_t SubSurfHeightinVoxels = static_cast<size_t>(SubSurfaceDepth / zSpacing); // height of SubSurf <-> Border to Surf
	size_t MeasurementEndIndex = static_cast<size_t>(projectedAirMask->GetLargestPossibleRegion().GetSize()[2]); // largest possible height

	MaskImageIteratorTypeIndexed MaskIterator(ImageMask, ImageMask->GetLargestPossibleRegion());
	MaskIterator.GoToBegin();
	// since mask is in physical domain, we have to calculate offsets for "naive" domain airMask
	itk::Index<3> MaskOffsetIndex;
	MaskOffsetIndex[0] = MaskIterator.GetIndex()[0]; // valid at first iterator point
	MaskOffsetIndex[1] = MaskIterator.GetIndex()[1];
	MaskOffsetIndex[2] = 0;
	while (!MaskIterator.IsAtEnd()) {
		if (MaskIterator.Get() == 0) {
			// no need to check already out of mask values...
			++MaskIterator;
			continue;
		}
		bool changeHasHappend = false;
		bool lastStatewasMaterial = false; //
		bool lastStatewasAir = false; // both are checked for changes. This way the first change doesn't count
		for (size_t height = SubSurfHeightinVoxels; height < MeasurementEndIndex; ++height) {
			itk::Index<3> IndexInProjectedImage;
			IndexInProjectedImage[0] = MaskIterator.GetIndex()[0] - MaskOffsetIndex[0];
			IndexInProjectedImage[1] = MaskIterator.GetIndex()[1] - MaskOffsetIndex[1];
			IndexInProjectedImage[2] = height;
			unsigned char AirMaskPixel = projectedAirMask->GetPixel(IndexInProjectedImage);
			if (AirMaskPixel != 251) {		// 251 -> out of image domain) {
				if (AirMaskPixel == 0) {
					lastStatewasMaterial = true;
					lastStatewasAir = false;
				}
				else if (AirMaskPixel == 255) {
					if (lastStatewasMaterial && !lastStatewasAir) {
						changeHasHappend = true;
						break;
					}
					lastStatewasMaterial = false;
					lastStatewasAir = true;
				}
			}
		}
		// height map is just naively placed in physical domain
		// mask is placed in physical coordinates
		if (!changeHasHappend) {
			// when no transition happens -> invalid!
			MaskIterator.Set(0);
			++MaskIterator;
			continue; // in this case the following special case does not have to be checked
		}
		// when a transition happens but uppermost pixel is material, we also expect non-plausible height values!
		itk::Index<3> HighestIndexInProjectedImage;
		HighestIndexInProjectedImage[0] = MaskIterator.GetIndex()[0] - MaskOffsetIndex[0];
		HighestIndexInProjectedImage[1] = MaskIterator.GetIndex()[1] - MaskOffsetIndex[1];
		HighestIndexInProjectedImage[2] = MeasurementEndIndex - 1; // end index = size, so -1 for max index
		if (changeHasHappend && (projectedAirMask->GetPixel(HighestIndexInProjectedImage) == 0)) {
			MaskIterator.Set(0);
		}
		++MaskIterator;
	}
}
// remove shape by projection of triangles to HeightMap domain and subtraction of interpolated heights
void removeShape(HeightMapImageType::Pointer HeightMap, HeightMapImageType::Pointer ShapeMap, const MaskImageType::Pointer ImageMask, const MeshType::Pointer STLMesh, std::vector<std::tuple<std::vector<itk::IdentifierType>, std::vector<itk::IdentifierType>, itk::CovariantVector< double, 3>>>& FaceGroupsCellNums, const std::size_t curProjGroup, StackImageType::DirectionType Direction) {
	std::vector<Triangle2D> TrianglesToCheck;
	std::vector<Triangle3D> ProjectedTriangles;
	for (std::size_t CellIteratorinGroup = 0; CellIteratorinGroup < std::get<1>(FaceGroupsCellNums[curProjGroup]).size(); ++CellIteratorinGroup) {
		// get vertices of cell
		itk::IdentifierType cellindexinGroup = std::get<1>(FaceGroupsCellNums[curProjGroup])[CellIteratorinGroup];
		// project to image region
		itk::QuadEdgeMeshPoint<MeshCoordType, StackMeshDimensionality> a = Direction * descriptiveIndices::s_cellVerticesPoints[cellindexinGroup][0]; // STLMesh->GetPoint(PointIDsofCurCell[0]);
		itk::QuadEdgeMeshPoint<MeshCoordType, StackMeshDimensionality> b = Direction * descriptiveIndices::s_cellVerticesPoints[cellindexinGroup][1]; // STLMesh->GetPoint(PointIDsofCurCell[1]);
		itk::QuadEdgeMeshPoint<MeshCoordType, StackMeshDimensionality> c = Direction * descriptiveIndices::s_cellVerticesPoints[cellindexinGroup][2]; // STLMesh->GetPoint(PointIDsofCurCell[2]);

		// convert to CGAL vector of triangles
		Point2D PointA(static_cast<double>(a[0]), static_cast<double>(a[1]));
		Point2D PointB(static_cast<double>(b[0]), static_cast<double>(b[1]));
		Point2D PointC(static_cast<double>(c[0]), static_cast<double>(c[1]));
		TrianglesToCheck.push_back(Triangle2D(PointA, PointB, PointC));

		Point3D PointA3D(static_cast<double>(a[0]), static_cast<double>(a[1]), static_cast<double>(a[2]));
		Point3D PointB3D(static_cast<double>(b[0]), static_cast<double>(b[1]), static_cast<double>(b[2]));
		Point3D PointC3D(static_cast<double>(c[0]), static_cast<double>(c[1]), static_cast<double>(c[2]));
		ProjectedTriangles.push_back(Triangle3D(PointA3D, PointB3D, PointC3D));
	}
	// triangle height should be normalized
	// this of course only works for uniform triangle sizes
	// should work "good enough" in practice to let the leveling do the rest later on
	size_t numofTriangles = TrianglesToCheck.size();
	double zAvg = 0.0;
	for (size_t indexofTriangle = 0; indexofTriangle < numofTriangles; ++indexofTriangle) {
		zAvg += ProjectedTriangles[indexofTriangle].vertex(0).z() + ProjectedTriangles[indexofTriangle].vertex(1).z() + ProjectedTriangles[indexofTriangle].vertex(2).z();
	}
	zAvg = zAvg / (static_cast<double>(numofTriangles) * 3.0);

	// setup mask image iteration
	MaskImageIteratorTypeIndexed MaskIterator(ImageMask, ImageMask->GetLargestPossibleRegion());
	HeightMapIteratorTypeIndexed ShapeIterator(ShapeMap, ShapeMap->GetLargestPossibleRegion());
	HeightMapIteratorTypeIndexed HMIterator(HeightMap, HeightMap->GetLargestPossibleRegion());
	MaskIterator.GoToBegin();
	ShapeIterator.GoToBegin();
	HMIterator.GoToBegin();
	size_t lastSuccessfullIndex = 0;
	size_t indexTri = 0;

	// additional bool check, since erode/dilate possibly makes unmasked areas where no triangle is available
	// -> maybe positive margins should issue a warning with remove shape?
	bool TriFound = false;

	while (!MaskIterator.IsAtEnd()) {
		itk::QuadEdgeMeshPoint<MeshCoordType, HeightMapDimensionality> IndexPointinPhysicalCoords;
		ImageMask->TransformIndexToPhysicalPoint(MaskIterator.GetIndex(), IndexPointinPhysicalCoords);
		Point2D PointinMask(IndexPointinPhysicalCoords[0], IndexPointinPhysicalCoords[1]);
		// bool PointIsInsideMesh = false;

		// ask Maskiterator if inside Mask
		if (MaskIterator.Get() == 1) {
			//size_t numTriangles = TrianglesToCheck.size();

			// wrap around loop starting at last successfull triangle
			// saves a few CPU cycles
			for (size_t indexofTriangle = lastSuccessfullIndex; indexofTriangle < (numofTriangles + lastSuccessfullIndex); ++indexofTriangle) {// const auto & Tri : TrianglesToCheck) {
				if (TrianglesToCheck[indexofTriangle % numofTriangles].has_on_bounded_side(PointinMask)) {
					//PointIsInsideMesh = true;
					lastSuccessfullIndex = indexofTriangle % numofTriangles;
					indexTri = indexofTriangle % numofTriangles;
					TriFound = true;
					break; // if it is in one triangle, we can end our search
				}
			}
			// TODO: What if normal of triangle ~= 0? --> should not be possible for reasonable mergeangles and projections
			if (TriFound) {
				// calculate height over this pt
				Vector3D NormalofTriangle = ProjectedTriangles[indexTri].supporting_plane().orthogonal_vector();
				NormalofTriangle = NormalofTriangle / sqrt(NormalofTriangle.squared_length());
				// catching the case where the triangle is not hit is not necessary because this is checked through the projection above
				// catching the case where the triangle intersection is a line should also not be necessary because of above check
				// it would then be exactly on boundary and not on_bounded_side
				// numerically possible to happen ?
				double ZatCurrentPoint = ProjectedTriangles[indexTri][0].z() - ((IndexPointinPhysicalCoords[0] - ProjectedTriangles[indexTri][0].x()) * NormalofTriangle.x() + (IndexPointinPhysicalCoords[1] - ProjectedTriangles[indexTri][0].y()) * NormalofTriangle.y()) / NormalofTriangle.z();
				ZatCurrentPoint = ZatCurrentPoint - zAvg; // to keep values near 0
				ShapeIterator.Set(ZatCurrentPoint);
				HMIterator.Set(HMIterator.Get() - ZatCurrentPoint);
			}
			else {
				ShapeIterator.Set(0.0);
			}
		}
		else {
			ShapeIterator.Set(0.0);
		}
		++MaskIterator;
		++ShapeIterator;
		++HMIterator;
	}
}

// since maximum/minimum possible measured values are MeasurementThickness/2, this values can be found easily
// with this function, they are considered non plausible
void maskOutofImageValues(MaskImageType::Pointer ImageMask, const HeightMapImageType::Pointer HeightMap, const double heightSpacing) {
	// mask height values which are 0 or -1
	HeightMapConstIteratorTypeIndexed HMIterator(HeightMap, HeightMap->GetLargestPossibleRegion());
	MaskImageIteratorTypeIndexed MaskIterator(ImageMask, ImageMask->GetLargestPossibleRegion());
	HMIterator.GoToBegin();
	MaskIterator.GoToBegin();
	while (!HMIterator.IsAtEnd()) {
		if (MaskIterator.Get() == 0) {
			// no need to check already out of mask values...
			++HMIterator;
			++MaskIterator;
			continue;
		}
		// zero is also non plausible and can happen for out of image values
		if (HMIterator.Get() < heightSpacing) {
			MaskIterator.Set(0);
		}
		++HMIterator;
		++MaskIterator;
	}
}
