#pragma once
#include "Misc.h"

namespace descriptiveIndices
{
	// set bool to true after indices have been calulated and to false if breaking changes were made
	bool s_indicesValid = false;
	double s_creaseAngle = 0;
	bool s_EdgeCalculationRequested = false;
	size_t s_numCells = 0;

	std::vector<double> s_cellAreaIndex;
	std::vector< itk::Point<double, 3>> s_cellMidPoints;
	std::vector< itk::Vector<double, 3>> s_cellNormals;
	std::vector< std::array<itk::Point<double, 3>, 3> > s_cellVerticesPoints;
	std::vector< std::array<MeshType::PointIdentifier, 3>> s_cellVerticesIDs;

	std::vector<bool> s_isanEdgeCell; // cell is located at an edge > XX deg
}

namespace aggregatePredicates {
	bool s_predicatesValid = false;
	size_t s_numFacegroups = 0;

	std::vector<bool> s_FacegrouptouchesEdge; // FG contains a cell which isanEdgeCell
}
//
//namespace aggregates
//{
//	// set bool to true after indices have been calculated and to false if breaking changes were made
//	bool s_aggregatesValid = false;
//	size_t s_numFaceGroups = 0;
//
//	std::vector<std::array<double, 6>> s_BoundingBoxes;
//	std::vector< StackImageType::DirectionType> s_Directions;
//	std::vector< itk::Vector<double, StackMeshDimensionality>> s_ResampleOrigins;
//	std::vector< itk::Vector<double, StackMeshDimensionality>> s_ProjectedOrigins;
//	std::vector< StackImageType::SizeType> s_ResampleSizes;
//}

// takes normals, deviation < maxRelDiff
bool AreFacesParallel(const itk::Vector<double, 3> a, const itk::Vector<double, 3> b, const float maxRelDiff = FLT_EPSILON) noexcept
{
	bool Xequal = false, Yequal = false, Zequal = false;
	// Calculate the difference.
	const double diffX = fabs(a[0] - b[0]);
	// Find the largest of A/B abs. and assign to largest
	const double largestX = (fabs(b[0]) > fabs(a[0])) ? fabs(b[0]) : fabs(a[0]);
	// Calculate the difference.
	const double diffY = fabs(a[1] - b[1]);
	// Find the largest of A/B abs. and assign to largest
	const double largestY = (fabs(b[1]) > fabs(a[1])) ? fabs(b[1]) : fabs(a[1]);
	// Calculate the difference.
	const double diffZ = fabs(a[2] - b[2]);
	// Find the largest of A/B abs. and assign to largest
	const double largestZ = (fabs(b[2]) > fabs(a[2])) ? fabs(b[2]) : fabs(a[2]);
	if (diffX <= largestX * maxRelDiff) {
		Xequal = true;
	}
	if (diffY <= largestY * maxRelDiff) {
		Yequal = true;
	}
	if (diffZ <= largestZ * maxRelDiff) {
		Zequal = true;
	}
	if (Xequal && Yequal && Zequal) {
		return true;
	}
	return false;
}

bool AreFloatsEqual(const float A, const float B, const float maxRelDiff = FLT_EPSILON) noexcept
{
	// Calculate the difference.
	const float diff = fabs(A - B);
	// Find the largest of A/B abs. and assign to largest
	const float largest = (fabs(B) > fabs(A)) ? fabs(B) : fabs(A);
	if (diff <= largest * maxRelDiff) {
		return true;
	}
	return false;
}

// compute if two normals are pointing in similar direction, deviation <= maxRelAngle1
bool AreNormalsSimilar(itk::Vector<double> a, itk::Vector<double> b, const float maxRelAngle = 10.0f)
{
	// make sure vectors are normed
	a.Normalize();
	b.Normalize();

	// replace acos with atan2(det crossprod, scalarprod) for numeric reasons
	double dotprod = std::clamp((a[0] * b[0] + a[1] * b[1] + a[2] * b[2]), -1., 1.);
	itk::Vector c = itk::CrossProduct(a, b);
	double relAngle = atan2(c.GetNorm(), dotprod) * 180 / itk::Math::pi;

	if (relAngle <= maxRelAngle) {
		return true;
	}
	return false;
}

// center point based on three vertices
itk::Point<double, 3> getCenterPointOfCell(const itk::Point<double, 3>& Pone, const itk::Point<double, 3>& Ptwo, const itk::Point<double, 3>& Pthree) noexcept
{
	itk::Point<double, 3> center;
	center[0] = (Pone[0] + Ptwo[0] + Pthree[0]) / 3.0;
	center[1] = (Pone[1] + Ptwo[1] + Pthree[1]) / 3.0;
	center[2] = (Pone[2] + Ptwo[2] + Pthree[2]) / 3.0;
	return center;
}

double sqdistBetweenITKPoints(const itk::Point<double, 3>& PointA, const itk::Point<double, 3>& PointB) noexcept
{
	double sqDist = (PointA[0] - PointB[0]) * (PointA[0] - PointB[0]) + (PointA[1] - PointB[1]) * (PointA[1] - PointB[1]) + (PointA[2] - PointB[2]) * (PointA[2] - PointB[2]);
	return sqDist;
}

bool isCellinThresRadius(const itk::Point<double, 3>& Pmid, const double DistThres, const itk::Point<double, 3>& Pone, const itk::Point<double, 3>& Ptwo, const itk::Point<double, 3>& Pthree) noexcept
{
	double sqDistThres = DistThres * DistThres;
	if ((sqdistBetweenITKPoints(Pmid, Pone) < sqDistThres) && (sqdistBetweenITKPoints(Pmid, Ptwo) < sqDistThres) && (sqdistBetweenITKPoints(Pmid, Pthree) < sqDistThres)) {
		return true;
	}
	return false;
}

double sqdistBetweenITKPoints2D(const itk::Point<double, 2>& PointA, const itk::Point<double, 2>& PointB) noexcept
{
	double sqDist = (PointA[0] - PointB[0]) * (PointA[0] - PointB[0]) + (PointA[1] - PointB[1]) * (PointA[1] - PointB[1]);
	return sqDist;
}

// calculate physical area of triangles in FG
double calculateSizeofFaceGroup(const std::tuple<std::vector<itk::IdentifierType>, std::vector<itk::IdentifierType>, itk::CovariantVector< double, 3>>& Facegroup) noexcept
{
	double GroupArea = 0;
	for (const itk::IdentifierType CellElem : std::get<1>(Facegroup)) {
		GroupArea += descriptiveIndices::s_cellAreaIndex[CellElem];
	}
	return GroupArea;
}

StackImageType::DirectionType calculateSimpleDirection(const std::vector<std::tuple<std::vector<itk::IdentifierType>, std::vector<itk::IdentifierType>, itk::CovariantVector< double, 3>>>& FaceGroupsCellNums, const size_t requestedFaceGroupIndex) {
	//
	StackImageType::DirectionType direction;
	std::size_t cellID = std::get<1>(FaceGroupsCellNums[requestedFaceGroupIndex])[0];
	if (cellID > descriptiveIndices::s_numCells) {
		std::cerr << "calculateDirection: RequestedCellID doesn't exist!" << std::endl;
		throw;
	}

	// yep 3 dimensions hardcoded here..
	// calculate cross product of dir 0 and dir 2 -> at least well defined orthogonal system
	itk::CovariantVector<double, StackMeshDimensionality> thirdVec;
	itk::Vector<double, StackMeshDimensionality> firstedge;
	itk::Vector<double, StackMeshDimensionality> facenormal;
	facenormal[0] = std::get<2>(FaceGroupsCellNums[requestedFaceGroupIndex])[0];
	facenormal[1] = std::get<2>(FaceGroupsCellNums[requestedFaceGroupIndex])[1];
	facenormal[2] = std::get<2>(FaceGroupsCellNums[requestedFaceGroupIndex])[2];
	// set face normal as 3rd columns
	direction(2, 0) = facenormal[0];
	direction(2, 1) = facenormal[1];
	direction(2, 2) = facenormal[2];

	// use first edge of cell as dir1
	//itk::Array<itk::IdentifierType> PointIDsofFirstCell = CAPtoFirstCell->GetPointIdsContainer();
	firstedge[0] = static_cast<double>((descriptiveIndices::s_cellVerticesPoints[cellID][2])[0]) - static_cast<double>((descriptiveIndices::s_cellVerticesPoints[cellID][1])[0]);
	firstedge[1] = static_cast<double>((descriptiveIndices::s_cellVerticesPoints[cellID][2])[1]) - static_cast<double>((descriptiveIndices::s_cellVerticesPoints[cellID][1])[1]);
	firstedge[2] = static_cast<double>((descriptiveIndices::s_cellVerticesPoints[cellID][2])[2]) - static_cast<double>((descriptiveIndices::s_cellVerticesPoints[cellID][1])[2]);

	firstedge.Normalize(); // returns the norm if needed
	direction(0, 0) = firstedge[0];
	direction(0, 1) = firstedge[1];
	direction(0, 2) = firstedge[2];

	// dir2 is result of crossproduct
	itk::CrossProduct(thirdVec, facenormal, firstedge);
	thirdVec.Normalize(); // should not be neccessary
	direction(1, 0) = thirdVec[0];
	direction(1, 1) = thirdVec[1];
	direction(1, 2) = thirdVec[2];
	//
	return direction;
}

// call calculateSimpleDirection to get a preliminary direction
// afterwards the oriented bounding box of preliminaryDirection * points is calculated
// longest axis of OBB sets new direction
StackImageType::DirectionType calculateOrientedDirection(const std::vector<std::tuple<std::vector<itk::IdentifierType>, std::vector<itk::IdentifierType>, itk::CovariantVector< double, 3>>>& FaceGroupsCellNums, const size_t requestedFaceGroupIndex) {
	StackImageType::DirectionType preliminaryDirection = calculateSimpleDirection(FaceGroupsCellNums, requestedFaceGroupIndex);
	StackImageType::DirectionType preliminaryInverseDirection = preliminaryDirection.GetInverse().as_matrix();
	std::vector<Point2D> ProjectedPointsofFG;
	std::vector<Point3D> ProjectedPointsofFG3D;

	for (std::size_t CellIteratorinGroup = 0; CellIteratorinGroup < std::get<1>(FaceGroupsCellNums[requestedFaceGroupIndex]).size(); ++CellIteratorinGroup) {
		itk::IdentifierType cellindexinGroup = std::get<1>(FaceGroupsCellNums[requestedFaceGroupIndex])[CellIteratorinGroup];

		descriptiveIndices::s_cellVerticesPoints[cellindexinGroup];
		for (size_t i = 0; i < descriptiveIndices::s_cellVerticesPoints[cellindexinGroup].size(); ++i)
		{
			ProjectedPointsofFG.push_back(Point2D((preliminaryDirection * descriptiveIndices::s_cellVerticesPoints[cellindexinGroup][i])[0], (preliminaryDirection * descriptiveIndices::s_cellVerticesPoints[cellindexinGroup][i])[1]));
			ProjectedPointsofFG3D.push_back(Point3D((preliminaryDirection * descriptiveIndices::s_cellVerticesPoints[cellindexinGroup][i])[0], (preliminaryDirection * descriptiveIndices::s_cellVerticesPoints[cellindexinGroup][i])[1], (preliminaryDirection * descriptiveIndices::s_cellVerticesPoints[cellindexinGroup][i])[2]));
		}
	}
	std::vector<Point2D> ConvexHull;
	// Projected points have to be brought into a convex polygon
	try {
		// CGALs AKL/Touissant is tainted for this input sometimes -> somehow unable to catch?
		//CGAL::convex_hull_2(ProjectedPointsofFG.begin(), ProjectedPointsofFG.end(), std::back_inserter(ConvexHull));
		CGAL::ch_bykat(ProjectedPointsofFG.begin(), ProjectedPointsofFG.end(), std::back_inserter(ConvexHull));
	}
	// general catch and try some other convex hull?
	catch (const CGAL::Postcondition_exception e) {
		// TODO: catch something else?
		std::cerr << "Postcondition failed convex hull at requestedFaceGroupIndex: " << requestedFaceGroupIndex << std::endl;
		// skip silently
	}
	catch (const CGAL::Failure_exception& e) {
		std::cerr << "Generic CGAL failure of convex hull at requestedFaceGroupIndex: " << requestedFaceGroupIndex << std::endl;
		std::cerr << e.what() << std::endl;
	}
	catch (const std::exception& e) {
		std::cerr << "std exception thrown by convex hull at requestedFaceGroupIndex: " << requestedFaceGroupIndex << std::endl;
		std::cerr << e.what() << std::endl;
	}
	catch (...) {
		std::cerr << "Unknown failure of convex hull at requestedFaceGroupIndex: " << requestedFaceGroupIndex << std::endl;
	}
	// uncomment and switch with the catch statement to see postcondition errors in debug mode
	//	catch (const CGAL::Postcondition_exception e) {
	//		// output failing postchecks in debug mode
	//		std::cerr << "Postcondition error in CGAL::convex_hull_2 at FaceGroup " << requestedFaceGroupIndex << "! \n";
	//		std::cerr << e.what() << std::endl;
	//}
	// The points denoted by the range [points_begin, points_end) form the boundary of a simple convex polygon P in counterclockwise orientation.
	Polygon_2 EnclosingRectangle;
	CGAL::min_rectangle_2(
		ConvexHull.begin(), ConvexHull.end(), std::back_inserter(EnclosingRectangle));

	// redefine direction to respect facegroup orientation
	StackImageType::DirectionType OrientedDirection;

	// keep old normal
	itk::CovariantVector<double, StackMeshDimensionality> thirdVec;
	itk::Vector<double, StackMeshDimensionality> tempfirstedge; // before projection
	itk::Vector<double, StackMeshDimensionality> firstedge;
	itk::Vector<double, StackMeshDimensionality> facenormal;
	facenormal[0] = std::get<2>(FaceGroupsCellNums[requestedFaceGroupIndex])[0];
	facenormal[1] = std::get<2>(FaceGroupsCellNums[requestedFaceGroupIndex])[1];
	facenormal[2] = std::get<2>(FaceGroupsCellNums[requestedFaceGroupIndex])[2];

	// set face normal as 3rd columns
	OrientedDirection(2, 0) = facenormal[0];
	OrientedDirection(2, 1) = facenormal[1];
	OrientedDirection(2, 2) = facenormal[2];

	// use main edge of oriented rectangle as firstedge
	// multiply direction "in-plane" with inverted preliminary direction
	// define longest axis as "main-axis"
	if (EnclosingRectangle.edge(0).squared_length() > EnclosingRectangle.edge(1).squared_length()) {
		tempfirstedge[0] = EnclosingRectangle.edge(0).direction().dx();
		tempfirstedge[1] = EnclosingRectangle.edge(0).direction().dy();
	}
	else {
		tempfirstedge[0] = EnclosingRectangle.edge(1).direction().dx();
		tempfirstedge[1] = EnclosingRectangle.edge(1).direction().dy();
	}
	tempfirstedge[2] = 0;
	// now rotate the found axis back
	tempfirstedge.Normalize(); // returns the norm if needed

	firstedge[0] = (preliminaryInverseDirection * tempfirstedge)[0];
	firstedge[1] = (preliminaryInverseDirection * tempfirstedge)[1];
	firstedge[2] = (preliminaryInverseDirection * tempfirstedge)[2];

	firstedge.Normalize(); // returns the norm if needed
	OrientedDirection(0, 0) = firstedge[0];
	OrientedDirection(0, 1) = firstedge[1];
	OrientedDirection(0, 2) = firstedge[2];

	// dir2 is result of crossproduct
	itk::CrossProduct(thirdVec, facenormal, firstedge);
	thirdVec.Normalize(); // should not be neccessary
	OrientedDirection(1, 0) = thirdVec[0];
	OrientedDirection(1, 1) = thirdVec[1];
	OrientedDirection(1, 2) = thirdVec[2];

	return OrientedDirection;
}

std::array<double, 4> getFaceGroupBoundingBox(const MeshType::Pointer STLMesh, const std::vector<std::tuple<std::vector<itk::IdentifierType>, std::vector<itk::IdentifierType>, itk::CovariantVector< double, 3>>>& FaceGroupsCellNums, const size_t requestedFaceGroupIndex, const StackImageType::DirectionType& direction) {
	std::vector<Point2D> VerticesinFaceGroup;
	VerticesinFaceGroup.reserve(std::get<1>(FaceGroupsCellNums[requestedFaceGroupIndex]).size() * 3);
	for (const itk::IdentifierType CellElem : std::get<1>(FaceGroupsCellNums[requestedFaceGroupIndex])) {
		itk::Point<double> Pone = direction * descriptiveIndices::s_cellVerticesPoints[CellElem][0];
		itk::Point<double> Ptwo = direction * descriptiveIndices::s_cellVerticesPoints[CellElem][1];
		itk::Point<double> Pthree = direction * descriptiveIndices::s_cellVerticesPoints[CellElem][2];
		VerticesinFaceGroup.push_back(Point2D(Pone[0], Pone[1]));
		VerticesinFaceGroup.push_back(Point2D(Ptwo[0], Ptwo[1]));
		VerticesinFaceGroup.push_back(Point2D(Pthree[0], Pthree[1]));
	}

	Rectangle_2D BRect = CGAL::bounding_box(VerticesinFaceGroup.begin(), VerticesinFaceGroup.end());

	// lets just use the AABB for now
	// later: aligned parallelogram?
	std::array<double, 4> BBox{};
	BBox[0] = BRect.xmin();
	BBox[1] = BRect.xmax();
	BBox[2] = BRect.ymin();
	BBox[3] = BRect.ymax();

	return BBox;
}

// compute build angle of face
double computeBuildingAngle(itk::CovariantVector<double> FaceNormal, itk::Vector<double, 3> BuildDirection)
{
	// make sure that vectors are normalized
	FaceNormal.Normalize();
	BuildDirection.Normalize();

	itk::Vector<double, 3> NormalasVector;
	NormalasVector[0] = FaceNormal[0];
	NormalasVector[1] = FaceNormal[1];
	NormalasVector[2] = FaceNormal[2];

	double dotprod = std::clamp((FaceNormal[0] * BuildDirection[0] + FaceNormal[1] * BuildDirection[1] + FaceNormal[2] * BuildDirection[2]), -1., 1.);
	itk::Vector c = itk::CrossProduct(NormalasVector, BuildDirection);
	double relAngle = atan2(c.GetNorm(), dotprod) * 180 / itk::Math::pi;
	return relAngle;
}

void precomputeIndices(const MeshType::Pointer STLMesh, const NormalQEMeshType::CellDataContainerPointer FaceNormals) {
	if (descriptiveIndices::s_indicesValid) {
		return;
	}
	descriptiveIndices::s_numCells = STLMesh->GetNumberOfCells();
	descriptiveIndices::s_cellAreaIndex.resize(descriptiveIndices::s_numCells, 0.0);
	descriptiveIndices::s_cellMidPoints.resize(descriptiveIndices::s_numCells, 0.0);;
	descriptiveIndices::s_cellNormals.resize(descriptiveIndices::s_numCells, 0.0);;
	descriptiveIndices::s_cellVerticesPoints.resize(descriptiveIndices::s_numCells, std::array<itk::Point<double, 3>, 3>{});
	descriptiveIndices::s_cellVerticesIDs.resize(descriptiveIndices::s_numCells, std::array<MeshType::PointIdentifier, 3>{});

	descriptiveIndices::s_isanEdgeCell.resize(descriptiveIndices::s_numCells, false);

	for (size_t cellIndex = 0; cellIndex < descriptiveIndices::s_numCells; ++cellIndex) {
		// get vertices
		MeshType::CellType::CellAutoPointer CAPtoFace;
		STLMesh->GetCell(cellIndex, CAPtoFace);
		itk::Array<itk::IdentifierType> PointIDsofCell = CAPtoFace->GetPointIdsContainer();
		itk::Point<double> Pone = STLMesh->GetPoint(PointIDsofCell[0]);
		itk::Point<double> Ptwo = STLMesh->GetPoint(PointIDsofCell[1]);
		itk::Point<double> Pthree = STLMesh->GetPoint(PointIDsofCell[2]);

		// area computation
		itk::Vector<double> SideOne = Ptwo - Pone;
		itk::Vector<double> SideTwo = Pthree - Pone;
		itk::Vector<double> CrossProd = itk::CrossProduct(SideOne, SideTwo);
		descriptiveIndices::s_cellAreaIndex[cellIndex] = 0.5 * sqrt(CrossProd.GetSquaredNorm());
		descriptiveIndices::s_cellMidPoints[cellIndex] = getCenterPointOfCell(Pone, Ptwo, Pthree);
		descriptiveIndices::s_cellNormals[cellIndex] = FaceNormals->GetElement(cellIndex);
		descriptiveIndices::s_cellVerticesPoints[cellIndex] = { Pone, Ptwo, Pthree };
		descriptiveIndices::s_cellVerticesIDs[cellIndex] = { PointIDsofCell[0],PointIDsofCell[1],PointIDsofCell[2] };

		// set isAnEdgeCell to false. it will either be set to true in the following or be changed
		descriptiveIndices::s_isanEdgeCell[cellIndex] = false;
		if (descriptiveIndices::s_EdgeCalculationRequested) {
			for (const auto& pointID : PointIDsofCell) {
				// define all vertex neighbours as neighbours
				MeshType::QEType* QEElement = STLMesh->FindEdge(pointID);
				MeshType::QEType* TempQEElement = QEElement;
				do
				{
					if (!AreNormalsSimilar(descriptiveIndices::s_cellNormals[cellIndex], FaceNormals->GetElement(TempQEElement->GetLeft()), descriptiveIndices::s_creaseAngle)) {
						// found a neighbour which has a normal which differs more than creaseAngle from normal of cell
						// set the edge predicate to true and break the search
						descriptiveIndices::s_isanEdgeCell[cellIndex] = true;
						break;
					}
					TempQEElement = TempQEElement->GetOnext();
				} while (QEElement != TempQEElement); // ring?
			}
		}
	}
	//
	// angle and distance conditions could be implemented as a deferred async operation in a symmetric data storage
	// if requested, conditions will be evaluated, check for dist or angle in reverse index configuration will be cheaper the second time
	//
	//
	//// distances and angles
	//for (size_t cellIndex = 0; cellIndex < numCells; ++cellIndex) {
	//	for (size_t othercellIndex = 0; othercellIndex < numCells; ++othercellIndex) {
	//		cellsqDistanceIndex[cellIndex][othercellIndex] = (sqdistBetweenITKPoints(CellMidPoints[cellIndex], CellMidPoints[othercellIndex]) < sqMergeDist);
	//		cellAngleIndex[cellIndex][othercellIndex] = AreNormalsSimilar(CellNormals[cellIndex], CellNormals[othercellIndex], MaxMergeAngle);
	//	}
	//}
	descriptiveIndices::s_indicesValid = true;
}

void precomputeAggregates(const std::vector<std::tuple<std::vector<itk::IdentifierType>, std::vector<itk::IdentifierType>, itk::CovariantVector< double, StackMeshDimensionality>>>& FaceGroups) {
	if (aggregatePredicates::s_predicatesValid) {
		return; // nothing to do here
	}
	if (!descriptiveIndices::s_indicesValid) {
		std::cout << "precompute mesh indices before aggregates!" << std::endl;
		return;
	}
	aggregatePredicates::s_numFacegroups = FaceGroups.size();
	aggregatePredicates::s_FacegrouptouchesEdge.resize(aggregatePredicates::s_numFacegroups, false);
	//
	for (int curProjGroup = 0; curProjGroup < aggregatePredicates::s_numFacegroups; ++curProjGroup) {
		size_t numCellsinGroup = std::get<1>(FaceGroups[curProjGroup]).size();
		for (size_t cellID = 0; cellID < numCellsinGroup; ++cellID) {
			if (descriptiveIndices::s_isanEdgeCell[std::get<1>(FaceGroups[curProjGroup])[cellID]]) {
				aggregatePredicates::s_FacegrouptouchesEdge[curProjGroup] = true;
				break;
			}
		}
	}

	aggregatePredicates::s_predicatesValid = true;
}

//void precomputeAggregates(std::vector<std::tuple<std::vector<itk::IdentifierType>, std::vector<itk::IdentifierType>, itk::CovariantVector< double, 3>>>& FaceGroupsCellNums, const double VoxelSpacing, const double MeasurementThickness, const double SubSurfaceDepth, const bool Masking) {
//	if (aggregates::s_aggregatesValid) {
//		return;
//	}
//	aggregates::s_numFaceGroups = FaceGroupsCellNums.size();
//	aggregates::s_BoundingBoxes.resize(aggregates::s_numFaceGroups);
//	aggregates::s_Directions.resize(aggregates::s_numFaceGroups);
//	aggregates::s_ResampleOrigins.resize(aggregates::s_numFaceGroups);
//	aggregates::s_ProjectedOrigins.resize(aggregates::s_numFaceGroups);
//	aggregates::s_ResampleSizes.resize(aggregates::s_numFaceGroups);
//
//	for (size_t groupIndex = 0; groupIndex < aggregates::s_numFaceGroups; ++groupIndex) {
//		aggregates::s_Directions[groupIndex] = calculateOrientedDirection(FaceGroupsCellNums, groupIndex);
//		itk::Vector<double, StackMeshDimensionality> facenormal;
//		facenormal[0] = aggregates::s_Directions[groupIndex](2, 0);
//		facenormal[1] = aggregates::s_Directions[groupIndex](2, 1);
//		facenormal[2] = aggregates::s_Directions[groupIndex](2, 2);
//
//		float xmin = std::numeric_limits<float>::max(), xmax = -std::numeric_limits<float>::max(), ymin = std::numeric_limits<float>::max(), ymax = -std::numeric_limits<float>::max(), zmin = std::numeric_limits<float>::max(), zmax = -std::numeric_limits<float>::max();
//
//		for (std::size_t CellIteratorinGroup = 0; CellIteratorinGroup < std::get<1>(FaceGroupsCellNums[groupIndex]).size(); ++CellIteratorinGroup) {
//			itk::IdentifierType currentCellIndex = std::get<1>(FaceGroupsCellNums[groupIndex])[CellIteratorinGroup];
//			//CellAutoPointer CAPofCurrentCell;
//			//STLMesh->GetCell(cellindexinGroup, CAPofCurrentCell);
//			//itk::Array<MeshType::PointIdentifier> PointIDsofCurCell = CAPofCurrentCell->GetPointIdsContainer();CAPofCurrentCell->GetNumberOfPoints()
//			for (size_t i = 0; i < descriptiveIndices::s_cellVerticesPoints[currentCellIndex].size(); ++i)
//			{
//				if ((aggregates::s_Directions[groupIndex] * descriptiveIndices::s_cellVerticesPoints[currentCellIndex][i])[0] < xmin) {
//					xmin = (aggregates::s_Directions[groupIndex] * descriptiveIndices::s_cellVerticesPoints[currentCellIndex][i])[0];
//				}
//
//				if ((aggregates::s_Directions[groupIndex] * descriptiveIndices::s_cellVerticesPoints[currentCellIndex][i])[0] > xmax) {
//					xmax = (aggregates::s_Directions[groupIndex] * descriptiveIndices::s_cellVerticesPoints[currentCellIndex][i])[0];
//				}
//
//				if ((aggregates::s_Directions[groupIndex] * descriptiveIndices::s_cellVerticesPoints[currentCellIndex][i])[1] < ymin) {
//					ymin = (aggregates::s_Directions[groupIndex] * descriptiveIndices::s_cellVerticesPoints[currentCellIndex][i])[1];
//				}
//
//				if ((aggregates::s_Directions[groupIndex] * descriptiveIndices::s_cellVerticesPoints[currentCellIndex][i])[1] > ymax) {
//					ymax = (aggregates::s_Directions[groupIndex] * descriptiveIndices::s_cellVerticesPoints[currentCellIndex][i])[1];
//				}
//
//				if ((aggregates::s_Directions[groupIndex] * descriptiveIndices::s_cellVerticesPoints[currentCellIndex][i])[2] < zmin) {
//					zmin = (aggregates::s_Directions[groupIndex] * descriptiveIndices::s_cellVerticesPoints[currentCellIndex][i])[2];
//				}
//
//				if ((aggregates::s_Directions[groupIndex] * descriptiveIndices::s_cellVerticesPoints[currentCellIndex][i])[2] > zmax) {
//					zmax = (aggregates::s_Directions[groupIndex] * descriptiveIndices::s_cellVerticesPoints[currentCellIndex][i])[2];
//				}
//			}
//		}
//		// code was to find the bounding box -> use std::minmax_element later on
//		unsigned int sizeX = static_cast<unsigned int>((static_cast<double>(xmax) - static_cast<double>(xmin)) / VoxelSpacing);
//		unsigned int sizeY = static_cast<unsigned int>((static_cast<double>(ymax) - static_cast<double>(ymin)) / VoxelSpacing);
//		aggregates::s_ResampleSizes[groupIndex][0] = sizeX;
//		aggregates::s_ResampleSizes[groupIndex][1] = sizeY;
//		aggregates::s_ResampleSizes[groupIndex][2] = (MeasurementThickness + SubSurfaceDepth) / VoxelSpacing; //CHANGED
//
//		itk::Vector<double, 3> OriginTemp;
//		OriginTemp[0] = xmin;
//		OriginTemp[1] = ymin;
//		OriginTemp[2] = zmin;
//		StackImageType::DirectionType InvDirection = aggregates::s_Directions[groupIndex].GetInverse().as_matrix();
//		OriginTemp = InvDirection * OriginTemp;
//		// shift origin after inversion so that one measurementthickness can be projected in both directions from STL face
//		aggregates::s_ResampleOrigins[groupIndex] = OriginTemp[0] - ((MeasurementThickness * 0.5) + SubSurfaceDepth) * facenormal[0]; //CHANGED
//		aggregates::s_ResampleOrigins[groupIndex] = OriginTemp[1] - ((MeasurementThickness * 0.5) + SubSurfaceDepth) * facenormal[1]; // 0.5 because we shift half the thickness
//		aggregates::s_ResampleOrigins[groupIndex] = OriginTemp[2] - ((MeasurementThickness * 0.5) + SubSurfaceDepth) * facenormal[2]; //
//		aggregates::s_ProjectedOrigins[groupIndex] = aggregates::s_Directions[groupIndex] * aggregates::s_ResampleOrigins[groupIndex];
//		//**********************************************************************************************************************************************
//		// generate Mask
//		//**********************************************************************************************************************************************
//		MaskImageType::Pointer HMMask = MaskImageType::New();
//		MaskImageType::RegionType MaskRegion;
//		MaskImageType::SizeType MaskSize = { aggregates::s_ResampleSizes[groupIndex][0],aggregates::s_ResampleSizes[groupIndex][1] };
//		MaskImageType::IndexType MaskstartIndex = { static_cast<itk::IndexValueType>(aggregates::s_ProjectedOrigins[groupIndex][0] / VoxelSpacing), static_cast<itk::IndexValueType>(aggregates::s_ProjectedOrigins[groupIndex][1] / VoxelSpacing) };
//		MaskRegion.SetIndex(MaskstartIndex);
//		MaskRegion.SetSize(MaskSize);
//		HMMask->SetRegions(MaskRegion);
//		itk::Vector<itk::ImageBase<2>::SpacingValueType, 2> spacing;
//		spacing[0] = VoxelSpacing;
//		spacing[1] = VoxelSpacing;
//		HMMask->SetSpacing(spacing);
//		HMMask->Allocate();
//		HMMask->FillBuffer(itk::NumericTraits<MaskPixelType>::Zero);
//		// Initialize valid region with size of unmasked HeightMap. This way value can be used either with or without masking
//		unsigned int MaskValidAreaPixels = aggregates::s_ResampleSizes[groupIndex][0] * aggregates::s_ResampleSizes[groupIndex][1];// HMImageMap->GetLargestPossibleRegion().GetNumberOfPixels();
//		if (Masking) {
//			// generate Mask from geometric information
//			generateMask(HMMask, STLMesh, FaceGroupsCellNums, curProjGroup, direction);
//			if (MeasurementMargin > 0) {
//				// use dilate on the mask
//				itk::FlatStructuringElement<2>::RadiusType Radius;
//				Radius.Fill(static_cast<itk::Size<2U>::SizeValueType>(MeasurementMargin));
//				itk::FlatStructuringElement<2> StructElement = itk::FlatStructuringElement<2>::Ball(Radius);
//				typedef itk::BinaryDilateImageFilter<MaskImageType, MaskImageType, itk::FlatStructuringElement<2>> dilateFilterType;
//				auto dilateFilter = dilateFilterType::New();
//				dilateFilter->SetInput(HMMask);
//				dilateFilter->SetKernel(StructElement);
//				dilateFilter->SetForegroundValue(1); // Value to dilate
//				dilateFilter->Update();
//				HMMask = dilateFilter->GetOutput();
//			}
//			if (MeasurementMargin < 0) {
//				// use erode on the mask
//				itk::FlatStructuringElement<2>::RadiusType Radius;
//				Radius.Fill((-1) * static_cast<itk::Size<2U>::SizeValueType>(MeasurementMargin));
//				itk::FlatStructuringElement<2> StructElement = itk::FlatStructuringElement<2>::Ball(Radius);
//				typedef itk::BinaryErodeImageFilter<MaskImageType, MaskImageType, itk::FlatStructuringElement<2>> erodeFilterType;
//				auto erodeFilter = erodeFilterType::New();
//				erodeFilter->SetInput(HMMask);
//				erodeFilter->SetKernel(StructElement);
//				erodeFilter->SetForegroundValue(1); // Value to erode
//				erodeFilter->SetBoundaryToForeground(false);
//				erodeFilter->Update();
//				HMMask = erodeFilter->GetOutput();
//			}
//		}
//		aggregates::s_aggregatesValid = true;
//	}
//}
// aggregateFaces starting from a pivot cell
// pivot cell is selected always as the next cell which isn't already assigned
// aggregation is limited by MaxMergeAngle
void aggregateFaces(const MeshType::Pointer STLMesh, std::vector<std::tuple<std::vector<itk::IdentifierType>, std::vector<itk::IdentifierType>, itk::CovariantVector< double, 3>>>& FaceGroupsCellNums, const NormalQEMeshType::CellDataContainerPointer FaceNormals, const float MaxMergeAngle) {
	// datastructure of FaceIDs not assigned yet
	std::vector<itk::IdentifierType> FaceIDs(STLMesh->GetNumberOfCells()); // vector with numCells entries
	std::iota(std::begin(FaceIDs), std::end(FaceIDs), 0); // Fill with 0, 1, ..., numCells-1
	std::size_t currentGroup = 0;
	std::set<itk::IdentifierType> IDsalreadyUsed;
	std::vector<bool> NeighboursMarkedForDeletion;
	NeighboursMarkedForDeletion.resize(STLMesh->GetNumberOfCells(), false);
	// iterate through FaceIDs and aggregate them into groups of same direction
	precomputeIndices(STLMesh, FaceNormals);
	auto cellIt = FaceIDs.begin();
	while (cellIt != FaceIDs.end())
	{
		if (NeighboursMarkedForDeletion[*cellIt]) {
			// this
			++cellIt;
			continue;
		}

		itk::IdentifierType currentCellNum = *cellIt;
		itk::IdentifierType FirstCellNumofGroup = currentCellNum;
		itk::CovariantVector< double, 3> initNormal(0.0);
		std::vector<itk::IdentifierType> initVector;
		initVector.push_back(currentCellNum);
		FaceGroupsCellNums.push_back(std::make_tuple(initVector, initVector, initNormal));
		// remove currentCellNum from FaceIDs, it's not available anymore.
		FaceIDs.erase(std::remove(FaceIDs.begin(), FaceIDs.end(), currentCellNum), FaceIDs.end());
		// additional set to keep track which elements are added to facegroup, vector is used in facegroup to keep order
		IDsalreadyUsed.insert(currentCellNum);

		// helper vars for finding all faces neighbouring currentCellNum
		std::set<itk::IdentifierType> currentNeighbours;// = new std::set<itk::IdentifierType>;
		itk::IdentifierType currentIDinGroup = 0;

		do {
			CellAutoPointer CAPofcurrentCell;
			STLMesh->GetCell(currentCellNum, CAPofcurrentCell);
			itk::Array<MeshType::PointIdentifier> pointIDofcurrentCell = CAPofcurrentCell->GetPointIdsContainer();

			for (const auto& pointID : pointIDofcurrentCell) {
				// since GetCellNeighbors doesn't work, we use the find edges approach from itk tutorials, and define all vertex neighbours as neighbours
				MeshType::QEType* QEElement = STLMesh->FindEdge(pointID);
				MeshType::QEType* TempQEElement = QEElement;
				do
				{
					currentNeighbours.insert(TempQEElement->GetLeft());
					TempQEElement = TempQEElement->GetOnext();
				} while (QEElement != TempQEElement); // ring?
			}

			// check normals of neighbouring faces
			for (const auto& neighbour : (currentNeighbours)) {
				if (AreNormalsSimilar(FaceNormals->GetElement(FirstCellNumofGroup), FaceNormals->GetElement(neighbour), MaxMergeAngle)) {
					// add only neighbouring faces with "similar" normal
					// check if face is not already inside
					if (IDsalreadyUsed.count(neighbour) != 1) {
						std::get<1>(FaceGroupsCellNums[currentGroup]).push_back(neighbour);
						IDsalreadyUsed.insert(neighbour);
					}
					// delete from FaceIDs since this ID has already been assigned
					NeighboursMarkedForDeletion[neighbour] = true;
				}
			}

			if (++currentIDinGroup < std::get<1>(FaceGroupsCellNums[currentGroup]).size()) {
				currentCellNum = std::get<1>(FaceGroupsCellNums[currentGroup]).at(currentIDinGroup);
			}
		} while (currentIDinGroup < std::get<1>(FaceGroupsCellNums[currentGroup]).size());

		itk::CovariantVector normalofFaceGroup(0.0);
		for (const auto& curCellinGroup : std::get<1>(FaceGroupsCellNums[currentGroup])) {
			// add weighted normal to normaloffacegroup
			normalofFaceGroup[0] = normalofFaceGroup[0] + FaceNormals->GetElement(curCellinGroup)[0] * descriptiveIndices::s_cellAreaIndex[curCellinGroup];
			normalofFaceGroup[1] = normalofFaceGroup[1] + FaceNormals->GetElement(curCellinGroup)[1] * descriptiveIndices::s_cellAreaIndex[curCellinGroup];
			normalofFaceGroup[2] = normalofFaceGroup[2] + FaceNormals->GetElement(curCellinGroup)[2] * descriptiveIndices::s_cellAreaIndex[curCellinGroup];
		}
		// normalize normal :D
		normalofFaceGroup.Normalize();
		std::get<2>(FaceGroupsCellNums[currentGroup]) = normalofFaceGroup;

		++currentGroup;
		// reset iterator after element deletion
		cellIt = FaceIDs.begin();
	}
}

// subdivide groups with Area larger than or equal MinDivisibleArea if resulting Groups are larger than MinAreaAfterDivision
// Groups are divided in the middle between the largest extending dimension of their bounding box
void subDivideLargeGroups(const MeshType::Pointer STLMesh, std::vector<std::tuple<std::vector<itk::IdentifierType>, std::vector<itk::IdentifierType>, itk::CovariantVector< double, 3>>>& FaceGroupsCellNums, const double MinDivisibleArea, const double MinAreaAfterDivision) {
	// first calculate areal size of each FaceGroup
	std::vector<double> sizeofFaceGroups;
	for (const std::tuple<std::vector<itk::IdentifierType>, std::vector<itk::IdentifierType>, itk::CovariantVector< double, 3>>&GroupElem : FaceGroupsCellNums) {
		sizeofFaceGroups.push_back(calculateSizeofFaceGroup(GroupElem));
	}

	bool subdividedAGroup = false;

	for (size_t groupID = 0; groupID < sizeofFaceGroups.size(); ++groupID) {
		if (subdividedAGroup) {
			--groupID; // last group was deleted, so continue there
		}
		// check if FaceGroup division is applicable
		if (sizeofFaceGroups[groupID] >= MinDivisibleArea && sizeofFaceGroups[groupID] >= 2 * MinAreaAfterDivision) { // even in the best case, the area has to be at least double the requested area to make it possible for the divided groups to be large enough
			StackImageType::DirectionType direction = calculateOrientedDirection(FaceGroupsCellNums, groupID);
			//StackImageType::DirectionType invDir = direction.GetInverse().as_matrix();
			std::array<double, 4> BBoxFG = getFaceGroupBoundingBox(STLMesh, FaceGroupsCellNums, groupID, direction);
			// divide FG in middle of larger BBox Span
			double xspan = std::abs(BBoxFG[1] - BBoxFG[0]);
			double yspan = std::abs(BBoxFG[3] - BBoxFG[2]);
			double xmiddle = 0.5 * (BBoxFG[1] + BBoxFG[0]);
			double ymiddle = 0.5 * (BBoxFG[3] + BBoxFG[2]);

			bool xDivision = xspan > yspan;

			// candidates for future FaceGroups

			std::tuple<std::vector<itk::IdentifierType>, std::vector<itk::IdentifierType>, itk::CovariantVector< double, 3>> FG_A = std::make_tuple(std::get<0>(FaceGroupsCellNums[groupID]), std::vector<itk::IdentifierType>(), std::get<2>(FaceGroupsCellNums[groupID]));
			std::tuple<std::vector<itk::IdentifierType>, std::vector<itk::IdentifierType>, itk::CovariantVector< double, 3>> FG_B = std::make_tuple(std::get<0>(FaceGroupsCellNums[groupID]), std::vector<itk::IdentifierType>(), std::get<2>(FaceGroupsCellNums[groupID]));

			for (const auto& cellElem : std::get<1>(FaceGroupsCellNums[groupID])) {
				MeshType::CellType::CellAutoPointer CAPtoCell;
				STLMesh->GetCell(cellElem, CAPtoCell);
				itk::Array<itk::IdentifierType> PointIDsofCell = CAPtoCell->GetPointIdsContainer();
				itk::Point<double> Pone = (direction * STLMesh->GetPoint(PointIDsofCell[0]));
				itk::Point<double> Ptwo = (direction * STLMesh->GetPoint(PointIDsofCell[1]));
				itk::Point<double> Pthree = (direction * STLMesh->GetPoint(PointIDsofCell[2]));
				// check for every face if each vertex  is below x middle. If yes, add to new facegroup A else to Facegroup B
				if (xDivision) {
					// divide along x-axis
					if (Pone[0] < xmiddle && Ptwo[0] < xmiddle && Pthree[0] < xmiddle) {
						// add to FG_A
						std::get<1>(FG_A).push_back(cellElem);
					}
					else if ((Pone[0] < xmiddle || AreFloatsEqual(Pone[0], xmiddle, 3 * FLT_EPSILON)) && (Ptwo[0] < xmiddle || AreFloatsEqual(Ptwo[0], xmiddle, 3 * FLT_EPSILON)) && (Pthree[0] < xmiddle || AreFloatsEqual(Pthree[0], xmiddle, 3 * FLT_EPSILON))) {
						// add to FG_A
						std::get<1>(FG_A).push_back(cellElem);
					}
					else {
						std::get<1>(FG_B).push_back(cellElem);
					}
				}
				else {
					// divide along y-axis
					if (Pone[1] < ymiddle && Ptwo[1] < ymiddle && Pthree[1] < ymiddle) {
						// add to FG_A
						std::get<1>(FG_A).push_back(cellElem);
					}
					else if ((Pone[1] < ymiddle || AreFloatsEqual(Pone[1], ymiddle, 3 * FLT_EPSILON)) && (Ptwo[1] < ymiddle || AreFloatsEqual(Ptwo[1], ymiddle, 3 * FLT_EPSILON)) && (Pthree[1] < ymiddle || AreFloatsEqual(Pthree[1], ymiddle, 3 * FLT_EPSILON))) {
						// add to FG_A
						std::get<1>(FG_A).push_back(cellElem);
					}
					else {
						std::get<1>(FG_B).push_back(cellElem);
					}
				}
			}

			// afterwards calculcate size of Facegroups A and B. If above MinAreaAfterDivision -> add to FaceGroupCellNums and delete GroupID from FaceGroupCellNums
			// then update sizeofFaceGroups and restart the loop
			double sizeofFG_A = calculateSizeofFaceGroup(FG_A);
			double sizeofFG_B = calculateSizeofFaceGroup(FG_B);
			if (sizeofFG_A > MinAreaAfterDivision && sizeofFG_B > MinAreaAfterDivision) {
				// was a sucessfull division
				// delete groupID and add FG_A and FG_B
				FaceGroupsCellNums.erase(FaceGroupsCellNums.begin() + groupID);
				FaceGroupsCellNums.push_back(FG_A);
				FaceGroupsCellNums.push_back(FG_B);

				sizeofFaceGroups.erase(sizeofFaceGroups.begin() + groupID);
				sizeofFaceGroups.push_back(sizeofFG_A);
				sizeofFaceGroups.push_back(sizeofFG_B);
				subdividedAGroup = true;
			}
			else {
				subdividedAGroup = false;
				// groupID is non divisible
			}
		}
		else {
			subdividedAGroup = false;
			// groupID is non divisible
		}
	}
}
//
//// TODO: function not finished
//// small groups of cells might be remaining between two larger groups
//// add this groups to the existing ones if relaxed condition is met
//// example: merge angle is measured to average normal instead of pivot normal
//// caveat: if a cell is added to another facegroup, the added size might make this facegroup large enough to survive
//void cleanUpSmallGroups(const MeshType::Pointer STLMesh, std::vector<std::pair<std::vector<itk::IdentifierType>, itk::CovariantVector< double, 3>>>& FaceGroupsCellNums, const NormalQEMeshType::CellDataContainerPointer FaceNormals, const float MaxMergeAngle, const double MaxCleanupArea) {
//	// clean possible too small elements which remained
//	// first calculate areal size of each FaceGroup
//	std::vector<double> sizeofFaceGroups;
//	for (const std::pair<std::vector<itk::IdentifierType>, itk::CovariantVector< double, 3>>&GroupElem : FaceGroupsCellNums) {
//		double GroupArea = 0;
//		for (const itk::IdentifierType CellElem : GroupElem.first) {
//			GroupArea += descriptiveIndices::s_cellAreaIndex[CellElem];
//		}
//		sizeofFaceGroups.push_back(GroupArea);
//	}
//
//	for (size_t groupID = 0; groupID < sizeofFaceGroups.size(); ++groupID) {
//		// check area vs thres
//		if (sizeofFaceGroups[groupID] < MaxCleanupArea) {
//			// start cleanup attempt
//			// we search neighbouring faces and if there are, check their normals
//			// we add to the neighbour with the smallest normal difference if this difference angle is smaller than MaxMergeAngle/2
//
//			// 01: Find neighbouring FaceGroups
//			// iterate through the cells in the too small group
//			for (std::size_t cellNuminGroup = 0; cellNuminGroup < FaceGroupsCellNums[groupID].first.size(); ++cellNuminGroup)// FaceGroupsCellNums[groupID].first.end() && (cellIt + 1 != FaceGroupsCellNums[groupID].first.end()))
//			{
//				// helper vars for finding all faces neighbouring currentCellNum
//				std::set<itk::IdentifierType> currentNeighbours;// = new std::set<itk::IdentifierType>;
//
//				CellAutoPointer CAPofcurrentCell;
//				STLMesh->GetCell(FaceGroupsCellNums[groupID].first[cellNuminGroup], CAPofcurrentCell);
//				itk::Array<MeshType::PointIdentifier> pointIDsofcurrentCell = CAPofcurrentCell->GetPointIdsContainer();
//
//				for (const auto& pointID : pointIDsofcurrentCell) {
//					// since GetCellNeighbors doesn't work, we use the find edges approach from itk tutorials, and define all vertex neighbours as neighbours
//					MeshType::QEType* QEElement = STLMesh->FindEdge(pointID);
//					MeshType::QEType* TempQEElement = QEElement;
//					do
//					{
//						currentNeighbours.insert(TempQEElement->GetLeft());
//						TempQEElement = TempQEElement->GetOnext();
//					} while (QEElement != TempQEElement); // ring?
//				}
//
//				for (const auto& neighbour : (currentNeighbours)) {
//					// check which facegroup the neighbours are
//					std::size_t FaceGroupIDcontainingNeighbour = 0;
//					std::vector<itk::IdentifierType>::iterator searchIter;
//					for (size_t facegroupNum = 0; facegroupNum < FaceGroupsCellNums.size(); ++facegroupNum)
//					{
//						searchIter = std::find(FaceGroupsCellNums[facegroupNum].first.begin(), FaceGroupsCellNums[facegroupNum].first.end(), neighbour);
//						if (searchIter != FaceGroupsCellNums[facegroupNum].first.end() && (sizeofFaceGroups[facegroupNum] > MaxCleanupArea)) {
//							FaceGroupIDcontainingNeighbour = facegroupNum;
//							break;
//						}
//						else {
//							continue;
//						}
//					}
//					// now we have a facegroupID of this neighbour
//					itk::Vector<double, 3> NormalofNeighbourGroup;
//					NormalofNeighbourGroup[0] = FaceGroupsCellNums[FaceGroupIDcontainingNeighbour].second[0];
//					NormalofNeighbourGroup[1] = FaceGroupsCellNums[FaceGroupIDcontainingNeighbour].second[1];
//					NormalofNeighbourGroup[2] = FaceGroupsCellNums[FaceGroupIDcontainingNeighbour].second[2];
//					// compare if normal of the neighbouring group and the lonely cell are similar
//					if (AreNormalsSimilar(NormalofNeighbourGroup, FaceNormals->GetElement(FaceGroupsCellNums[groupID].first[cellNuminGroup]), 2 * MaxMergeAngle)) {
//						// add to this facegroup
//						FaceGroupsCellNums[FaceGroupIDcontainingNeighbour].first.push_back(FaceGroupsCellNums[groupID].first[cellNuminGroup]);
//						// remove currentCellNum from FaceIDs, it's not available anymore.
//						FaceGroupsCellNums[groupID].first.erase(std::remove(FaceGroupsCellNums[groupID].first.begin(), FaceGroupsCellNums[groupID].first.end(), FaceGroupsCellNums[groupID].first[cellNuminGroup]), FaceGroupsCellNums[groupID].first.end());
//
//						// correct normals
//						itk::CovariantVector<double, 3> normalofFaceGroupAddedTo(0.0);
//						for (const auto& curCellinGroup : FaceGroupsCellNums[FaceGroupIDcontainingNeighbour].first) {
//							// add weighted normal to normaloffacegroup
//							normalofFaceGroupAddedTo[0] = normalofFaceGroupAddedTo[0] + FaceNormals->GetElement(curCellinGroup)[0] * descriptiveIndices::s_cellAreaIndex[curCellinGroup];
//							normalofFaceGroupAddedTo[1] = normalofFaceGroupAddedTo[1] + FaceNormals->GetElement(curCellinGroup)[1] * descriptiveIndices::s_cellAreaIndex[curCellinGroup];
//							normalofFaceGroupAddedTo[2] = normalofFaceGroupAddedTo[2] + FaceNormals->GetElement(curCellinGroup)[2] * descriptiveIndices::s_cellAreaIndex[curCellinGroup];
//							// normalize normal :D
//							normalofFaceGroupAddedTo.Normalize();
//							FaceGroupsCellNums[FaceGroupIDcontainingNeighbour].second = normalofFaceGroupAddedTo;
//						}
//
//						// correct size of facegroup
//						double GroupAreaNew = 0;
//						for (const itk::IdentifierType CellElem : FaceGroupsCellNums[FaceGroupIDcontainingNeighbour].first) {
//							GroupAreaNew += descriptiveIndices::s_cellAreaIndex[CellElem];
//						}
//						sizeofFaceGroups[FaceGroupIDcontainingNeighbour] = GroupAreaNew;
//						// we achieved the goal, try to fix next cell
//						break;
//					}
//				}
//			}
//		}
//	}
//	// remove facegroup if empty
//	for (std::size_t checkedIndex = 0; checkedIndex < FaceGroupsCellNums.size(); ++checkedIndex) {
//		if (FaceGroupsCellNums[checkedIndex].first.size() == 0) {
//			FaceGroupsCellNums.erase(FaceGroupsCellNums.begin() + checkedIndex);
//		}
//	}
//}

// aggregates Faces around a pivot cell keeping only cells within a certain distance to pivot cell
// changes to aggregateFaces
// Faces are not (virtually) deleted and can be used for the group of the next cell
// so faces have to be deleted for each run and reactivated in the next run
// the pivot cell is not just the "next" cell but every cell
//void aggregateFacesAroundCell(const MeshType::Pointer STLMesh, std::vector<std::tuple<std::vector<itk::IdentifierType>, std::vector<itk::IdentifierType>, itk::CovariantVector< double, 3>>>& FaceGroupsCellNums, const NormalQEMeshType::CellDataContainerPointer FaceNormals, const float MaxMergeAngle, const double MaxMergeDistance) {
//	// datastructure of FaceIDs
//	const size_t numCells = STLMesh->GetNumberOfCells();
//	std::vector<itk::IdentifierType> FaceIDs(numCells); // vector with numCells entries
//	std::iota(std::begin(FaceIDs), std::end(FaceIDs), 0); // Fill with 0, 1, ..., numCells-1
//	std::size_t currentGroup = 0;
//	std::vector<bool> NeighboursMarkedForDeletion;
//	NeighboursMarkedForDeletion.resize(numCells, false);
//	// preparation work
//	// generate a cellAreaIndex with areas for each cell, so they don't have to be recomputed
//
//	std::vector<double> cellAreaIndex(numCells, 0.0);
//	for (size_t cellIndex = 0; cellIndex < numCells; ++cellIndex) {
//		cellAreaIndex[cellIndex] = getFaceArea(cellIndex, STLMesh);
//	}
//
//	// iterate through FaceIDs and aggregate them into groups of same direction
//	auto cellIt = FaceIDs.begin();
//	while (cellIt != FaceIDs.end())
//	{
//		// reset NeighboursMarkedForDeletion
//		std::fill(NeighboursMarkedForDeletion.begin(), NeighboursMarkedForDeletion.end(), false);
//		// set is placed inside loop to be reset for each group
//		std::set<itk::IdentifierType> IDsalreadyUsed;
//
//		itk::IdentifierType pivotCellNum = *cellIt;
//		itk::IdentifierType currentCellNum = pivotCellNum;
//		itk::CovariantVector< double, 3> initNormal(0.0);
//		std::vector<itk::IdentifierType> initVector; // empty vector
//		initVector.push_back(pivotCellNum);
//		FaceGroupsCellNums.push_back(std::make_tuple(initVector, initVector, initNormal));
//
//		// helper vars for finding all faces neighbouring currentCellNum
//		std::set<itk::IdentifierType> currentNeighbours;// = new std::set<itk::IdentifierType>;
//		itk::IdentifierType currentIDinGroup = 0;
//
//		// calculate midpoint of pivotCell
//		CellAutoPointer CAPofpivotCell;
//		STLMesh->GetCell(pivotCellNum, CAPofpivotCell);
//		itk::Array<MeshType::PointIdentifier> pointIDofpivotCell = CAPofpivotCell->GetPointIdsContainer();
//		itk::Point<double, 3> pivotMidPoint = getCenterPointOfCell(STLMesh->GetPoint(pointIDofpivotCell[0]), STLMesh->GetPoint(pointIDofpivotCell[1]), STLMesh->GetPoint(pointIDofpivotCell[2]));
//
//		do {
//			// find (vertex) neighbours of current cell inside this growing region
//			CellAutoPointer CAPofcurrentCell;
//			STLMesh->GetCell(currentCellNum, CAPofcurrentCell);
//			itk::Array<MeshType::PointIdentifier> pointIDofcurrentCell = CAPofcurrentCell->GetPointIdsContainer();
//			for (const auto& pointID : pointIDofcurrentCell) {
//				MeshType::QEType* QEElement = STLMesh->FindEdge(pointID);
//				MeshType::QEType* TempQEElement = QEElement;
//				do
//				{
//					currentNeighbours.insert(TempQEElement->GetLeft());
//					TempQEElement = TempQEElement->GetOnext();
//				} while (QEElement != TempQEElement); // ring?
//			}
//
//			// two checks for neighbouring faces (neighbouring to current cell in growing region)
//			// 1. distance to pivot face
//			// 2. normal deviation to pivot face
//			for (const auto& neighbour : (currentNeighbours)) {
//				if (NeighboursMarkedForDeletion[neighbour]) {
//					// this neighbour has already been assigned
//					// skip to save processing time
//					continue;
//				}
//				// calculate midpoint of neighbourCell
//				CellAutoPointer CAPofneighbourCell;
//				STLMesh->GetCell(neighbour, CAPofneighbourCell);
//				itk::Array<MeshType::PointIdentifier> pointIDsofneighbourCell = CAPofneighbourCell->GetPointIdsContainer();
//				//itk::Point<double, 3> neighbourMidPoint = getCenterPointOfCell(STLMesh->GetPoint(pointIDsofneighbourCell[0]), STLMesh->GetPoint(pointIDsofneighbourCell[1]), STLMesh->GetPoint(pointIDsofneighbourCell[2]));
//				// distance
//
//				//bool distanceCondition = sqrt(sqdistBetweenITKPoints(pivotMidPoint, neighbourMidPoint)) < MaxMergeDistance;
//				bool distanceCondition = isCellinThresRadius(pivotMidPoint, MaxMergeDistance, STLMesh->GetPoint(pointIDsofneighbourCell[0]), STLMesh->GetPoint(pointIDsofneighbourCell[1]), STLMesh->GetPoint(pointIDsofneighbourCell[2])); //sqrt(sqdistBetweenITKPoints(pivotMidPoint, neighbourMidPoint)) < MaxMergeDistance;
//				bool angleCondition = AreNormalsSimilar(FaceNormals->GetElement(pivotCellNum), FaceNormals->GetElement(neighbour), MaxMergeAngle);
//
//				if (distanceCondition && angleCondition) {
//					// add only neighbouring faces with "similar" normal and inside the radius
//					// check if face is not already inside
//					if (IDsalreadyUsed.count(neighbour) != 1) {
//						std::get<1>(FaceGroupsCellNums[currentGroup]).push_back(neighbour);
//						IDsalreadyUsed.insert(neighbour);
//					}
//				}
//				// delete from FaceIDs since this ID either already been assigned or didn't match the conditions
//				NeighboursMarkedForDeletion[neighbour] = true;
//			}
//
//			if (++currentIDinGroup < std::get<1>(FaceGroupsCellNums[currentGroup]).size()) {
//				currentCellNum = std::get<1>(FaceGroupsCellNums[currentGroup]).at(currentIDinGroup);
//			}
//		} while (currentIDinGroup < std::get<1>(FaceGroupsCellNums[currentGroup]).size());
//
//		itk::CovariantVector normalofFaceGroup(0.0);
//		for (const auto& curCellinGroup : std::get<1>(FaceGroupsCellNums[currentGroup])) {
//			// add weighted normal to normaloffacegroup
//			normalofFaceGroup[0] = normalofFaceGroup[0] + FaceNormals->GetElement(curCellinGroup)[0] * cellAreaIndex[curCellinGroup];
//			normalofFaceGroup[1] = normalofFaceGroup[1] + FaceNormals->GetElement(curCellinGroup)[1] * cellAreaIndex[curCellinGroup];
//			normalofFaceGroup[2] = normalofFaceGroup[2] + FaceNormals->GetElement(curCellinGroup)[2] * cellAreaIndex[curCellinGroup];
//		}
//		// normalize normal :D
//		normalofFaceGroup.Normalize();
//		std::get<2>(FaceGroupsCellNums[currentGroup]) = normalofFaceGroup;
//
//		// it can be assumed that all neighbours have been found
//		// thus, increase currentGroup and search next
//		++currentGroup;
//		++cellIt;
//	}
//}

// aggregates Faces around a pivot cell keeping only cells within a certain distance to pivot cell
// changes to aggregateFaces
// Faces are not (virtually) deleted and can be used for the group of the next cell
// so faces have to be deleted for each run and reactivated in the next run
// the pivot cell is not just the "next" cell but every cell
void aggregateFacesAroundCell(const MeshType::Pointer STLMesh, std::vector<std::tuple<std::vector<itk::IdentifierType>, std::vector<itk::IdentifierType>, itk::CovariantVector< double, 3>>>& FaceGroupsCellNums, const NormalQEMeshType::CellDataContainerPointer FaceNormals, const float MaxMergeAngle, const double MaxMergeDistance) {
	// datastructure of FaceIDs
	const size_t numCells = STLMesh->GetNumberOfCells();
	std::vector<itk::IdentifierType> FaceIDs(numCells); // vector with numCells entries
	std::iota(std::begin(FaceIDs), std::end(FaceIDs), 0); // Fill with 0, 1, ..., numCells-1

	// iterate through FaceIDs and aggregate them into groups of same direction
	FaceGroupsCellNums.resize(numCells, std::tuple<std::vector<itk::IdentifierType>, std::vector<itk::IdentifierType>, itk::CovariantVector< double, 3 >>{});
	precomputeIndices(STLMesh, FaceNormals);

	for (size_t cellID = 0; cellID < FaceIDs.size(); ++cellID)
	{
		itk::IdentifierType pivotCellNum = cellID;
		itk::IdentifierType currentCellNum = pivotCellNum;
		itk::CovariantVector<double, 3> initNormal(0.0);
		std::vector<itk::IdentifierType> initVector; // empty vector
		initVector.push_back(pivotCellNum);
		FaceGroupsCellNums[cellID] = std::make_tuple(initVector, initVector, initNormal);

		// helper vars for finding all faces neighbouring currentCellNum
		// set is placed inside loop to be reset for each group
		std::set<itk::IdentifierType> currentNeighbours;
		std::vector<itk::IdentifierType> curNeighboursOrdered;
		itk::IdentifierType currentIDinGroup = 0;

		// calculate midpoint of pivotCell
		itk::Point<double, 3> pivotMidPoint = descriptiveIndices::s_cellMidPoints[pivotCellNum];

		size_t lastindexfinished = 0;
		do {
			// find neighbours
			for (const auto& pointID : descriptiveIndices::s_cellVerticesIDs[currentCellNum]) {
				MeshType::QEType* QEElement = STLMesh->FindEdge(pointID); // last remaining call on STLMesh
				MeshType::QEType* TempQEElement = QEElement;
				do
				{
					if (currentNeighbours.insert(TempQEElement->GetLeft()).second) {
						curNeighboursOrdered.push_back(TempQEElement->GetLeft());
					}
					TempQEElement = TempQEElement->GetOnext();
				} while (QEElement != TempQEElement); // ring?
			}

			// two checks for neighbouring faces (neighbouring to current cell in growing region)
			// 1. distance to pivot face
			// 2. normal deviation to pivot face
			// start at last index, since all others have been checked
			for (size_t indexofNeighbour = lastindexfinished; indexofNeighbour < curNeighboursOrdered.size(); ++indexofNeighbour) {
				bool distanceCondition = isCellinThresRadius(pivotMidPoint, MaxMergeDistance, descriptiveIndices::s_cellVerticesPoints[curNeighboursOrdered[indexofNeighbour]][0], descriptiveIndices::s_cellVerticesPoints[curNeighboursOrdered[indexofNeighbour]][1], descriptiveIndices::s_cellVerticesPoints[curNeighboursOrdered[indexofNeighbour]][2]);// STLMesh->GetPoint(pointIDsofneighbourCell[0]), STLMesh->GetPoint(pointIDsofneighbourCell[1]), STLMesh->GetPoint(pointIDsofneighbourCell[2])); //sqrt(sqdistBetweenITKPoints(pivotMidPoint, neighbourMidPoint)) < MaxMergeDistance;
				bool angleCondition = AreNormalsSimilar(descriptiveIndices::s_cellNormals[pivotCellNum], descriptiveIndices::s_cellNormals[curNeighboursOrdered[indexofNeighbour]], MaxMergeAngle);//FaceNormals->GetElement(pivotCellNum), FaceNormals->GetElement(curNeighboursOrdered[indexofNeighbour]), MaxMergeAngle);

				if (distanceCondition && angleCondition) {
					// add only neighbouring faces with "similar" normal and inside the radius
					std::get<1>(FaceGroupsCellNums[cellID]).push_back(curNeighboursOrdered[indexofNeighbour]);
				}
			}
			lastindexfinished = curNeighboursOrdered.size();

			if (++currentIDinGroup < std::get<1>(FaceGroupsCellNums[cellID]).size()) {
				currentCellNum = std::get<1>(FaceGroupsCellNums[cellID]).at(currentIDinGroup);
			}
		} while (currentIDinGroup < std::get<1>(FaceGroupsCellNums[cellID]).size());

		itk::CovariantVector normalofFaceGroup(0.0);
		for (const auto& curCellinGroup : std::get<1>(FaceGroupsCellNums[cellID])) {
			// add weighted normal to normaloffacegroup
			normalofFaceGroup[0] = normalofFaceGroup[0] + descriptiveIndices::s_cellNormals[curCellinGroup][0] * descriptiveIndices::s_cellAreaIndex[curCellinGroup];
			normalofFaceGroup[1] = normalofFaceGroup[1] + descriptiveIndices::s_cellNormals[curCellinGroup][1] * descriptiveIndices::s_cellAreaIndex[curCellinGroup];
			normalofFaceGroup[2] = normalofFaceGroup[2] + descriptiveIndices::s_cellNormals[curCellinGroup][2] * descriptiveIndices::s_cellAreaIndex[curCellinGroup];
		}
		// normalize normal :D
		normalofFaceGroup.Normalize();
		std::get<2>(FaceGroupsCellNums[cellID]) = normalofFaceGroup;
		// it can be assumed that all neighbours have been found
	}
}
