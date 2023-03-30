//// CGAL includes.
//#include "Misc.h"
//
//#include <CGAL/property_map.h>
//#include <CGAL/IO/read_points.h>
//#include <CGAL/Point_with_normal_3.h>
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Shape_detection/Efficient_RANSAC.h>
//
//#include <CGAL/Timer.h>
//#include <map>
//
//// Type declarations.
//typedef CGAL::Exact_predicates_inexact_constructions_kernel  KernelEPIC;
//typedef KernelEPIC::FT                                           FType;
//typedef KernelEPIC::Point_3									Point3SD;
//typedef KernelEPIC::Vector_3								Vector3SD;
//
//typedef std::pair<KernelEPIC::Point_3, KernelEPIC::Vector_3>         Point_with_normal;
//typedef std::vector<Point_with_normal>                       Pwn_vector;
//typedef CGAL::First_of_pair_property_map<Point_with_normal>  Point_map;
//typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;
//typedef CGAL::Shape_detection::Efficient_RANSAC_traits
//<KernelEPIC, Pwn_vector, Point_map, Normal_map>             TraitsSD;
//typedef CGAL::Shape_detection::Efficient_RANSAC<TraitsSD> Efficient_ransac;
//typedef CGAL::Shape_detection::Cone<TraitsSD>             ConeSD;
//typedef CGAL::Shape_detection::Cylinder<TraitsSD>         CylinderSD;
//typedef CGAL::Shape_detection::Plane<TraitsSD>            PlaneSD;
//typedef CGAL::Shape_detection::Sphere<TraitsSD>           SphereSD;
//typedef CGAL::Shape_detection::Torus<TraitsSD>            TorusSD;
//
//Point3SD itkPt2CGALPt(const MeshType::PointType& a) {
//	Point3SD convertedPt = Point3SD(a[0], a[1], a[2]);
//	return convertedPt;
//}
//
//Vector3SD itkVec2CGALVec(const MeshType::VectorType& a) {
//	Vector3SD convertedVec = Vector3SD(a[0], a[1], a[2]);
//	return convertedVec;
//}
//
//int ShapeClassification(const MeshType::Pointer& STLMesh, const NormalQEMeshType::CellDataContainerPointer& FaceNormals, std::vector<int>& Classification) {
//	// read in Mesh File
//	// just reading the same Mesh seems easier than converting the itk datastructure
//
//	// set some values for search parameterisation
//	const size_t numofTrials = 10;
//
//	// itk mesh zu punkten umwandeln
//	// normals aus faces verwenden
//	// Points with normals.
//	Pwn_vector points;
//	// map to keep track which pointID belongs to which cellID
//	// 1. pointID, 2. cellID
//	std::map<std::size_t, std::size_t> repPointofCell;
//
//	for (size_t CellID = 0; CellID < STLMesh->GetNumberOfCells(); ++CellID) {
//		CellAutoPointer CAPofcurrentCell;
//		STLMesh->GetCell(CellID, CAPofcurrentCell);
//		itk::Array<MeshType::PointIdentifier> pointIDsofcurrentCell = CAPofcurrentCell->GetPointIdsContainer();
//
//		double curCentroidX = 0;
//		double curCentroidY = 0;
//		double curCentroidZ = 0;
//		// use vertices of face with cell normal
//		for (const auto& pointID : pointIDsofcurrentCell) {
//			points.push_back(std::make_pair(itkPt2CGALPt(STLMesh->GetPoint(pointID)), itkVec2CGALVec(FaceNormals->GetElement(CellID))));
//			//points.push_back(std::make_pair(itkPt2CGALPt(STLMesh->GetPoint(pointID)), itkVec2CGALVec(VertexNormals->GetElement(pointID))));
//			curCentroidX += itkPt2CGALPt(STLMesh->GetPoint(pointID)).x();
//			curCentroidY += itkPt2CGALPt(STLMesh->GetPoint(pointID)).y();
//			curCentroidZ += itkPt2CGALPt(STLMesh->GetPoint(pointID)).z();
//		}
//		curCentroidX = curCentroidX / 3.0; curCentroidY = curCentroidY / 3.0; curCentroidZ = curCentroidZ / 3.0;
//		// and the face barycenter with the face normal
//		// the classification of the barycenter will be used for the classification of the cell later on, since vertex points can be multiply assigned
//		points.push_back(std::make_pair(Point3SD(curCentroidX, curCentroidY, curCentroidZ), itkVec2CGALVec(FaceNormals->GetElement(CellID))));
//		// last entry of points is the rep Point for this cell
//		repPointofCell[points.size() - 1] = CellID;
//	}
//
//	//std::cout << points.size() << " points" << std::endl;
//	// Instantiate shape detection engine.
//	Efficient_ransac ransac;
//	// Provide input data.
//	ransac.set_input(points);
//	// Register shapes for detection.
//	ransac.add_shape_factory<PlaneSD>();
//	// ransac.add_shape_factory<SphereSD>();
//	ransac.add_shape_factory<CylinderSD>();
//	ransac.add_shape_factory<ConeSD>();
//
//	// Set parameters for shape detection.
//	Efficient_ransac::Parameters parameters;
//	// Set probability to miss the largest primitive at each iteration.
//	parameters.probability = 0.001;
//	// Detect shapes with at least 200 points.
//	parameters.min_points = 40;
//	// Set maximum Euclidean distance between a point and a shape.
//	parameters.epsilon = 0.1;
//	// Set maximum Euclidean distance between points to be clustered.
//	parameters.cluster_epsilon = 10.0; // TODO: sollte mind groesse der BBox haben
//	// Set maximum normal deviation.
//	// 0.9 < dot(surface_normal, point_normal);
//	parameters.normal_threshold = 0.999;
//	// Detect shapes.
//	CGAL::Timer time;
//	//ransac.detect(parameters);
//
//	// Perform detection several times and choose result with the highest coverage.
//	Efficient_ransac::Shape_range shapes = ransac.shapes();
//	FType best_coverage = 0;
//	FType best_numofShapes = FType(points.size());
//	std::cout << "Attempting to find a good shape detection solution... \n";
//	for (std::size_t Trial = 0; Trial < numofTrials; ++Trial) {
//		// Reset timer.
//		time.reset();
//		time.start();
//		// Detect shapes.
//		ransac.detect(parameters);
//		// Measure time after detection.
//		time.stop();
//		// Compute coverage, i.e. ratio of the points assigned to a shape.
//		FType coverage = FType(points.size() - ransac.number_of_unassigned_points()) / FType(points.size());
//		FType  numofShapes = ransac.shapes().end() - ransac.shapes().begin();
//		// Print number of assigned shapes and unassigned points.
//		std::cout << Trial << ": time: " << time.time() * 1000 << "ms -- ";
//		std::cout << numofShapes << " primitives, " << coverage << " coverage \n";
//
//		// Choose result with the highest coverage
//		if (coverage > best_coverage) {
//			best_coverage = coverage;
//			best_numofShapes = numofShapes;
//			// Efficient_ransac::shapes() provides
//			// an iterator range to the detected shapes.
//			shapes = ransac.shapes();
//		}// and lowest num of shapes (avoid oversegmentation) if indifferent
//		else if (coverage == best_coverage && numofShapes < best_numofShapes) {
//			best_coverage = coverage;
//			best_numofShapes = numofShapes;
//			// Efficient_ransac::shapes() provides
//			// an iterator range to the detected shapes.
//			shapes = ransac.shapes();
//		}
//	}
//	// Print number of detected shapes and unassigned points.
//	std::cout << best_numofShapes << " detected shapes, " << ransac.number_of_unassigned_points() << " unassigned points. \n";
//	// Efficient_ransac::shapes() provides
//	// an iterator range to the detected shapes.
//	//Efficient_ransac::Shape_range shapes = ransac.shapes();
//	Efficient_ransac::Shape_range::iterator it = shapes.begin();
//
//	//
//	//while (it != shapes.end()) {
//	//	// print found shapes
//	//	std::cout << (*it)->info();
//	//	// Sums distances of points to the detected shapes.
//	//	FType sum_distances = 0;
//	//	// Iterate through point indices assigned to each detected shape.
//	//	std::vector<std::size_t>::const_iterator index_it = (*it)->indices_of_assigned_points().begin();
//	//	while (index_it != (*it)->indices_of_assigned_points().end()) {
//	//		// Retrieve point.
//	//		const Point_with_normal& p = *(points.begin() + (*index_it));
//	//		// Adds Euclidean distance between point and shape.
//	//		sum_distances += CGAL::sqrt((*it)->squared_distance(p.first));
//	//		// Proceed with the next point.
//	//		index_it++;
//	//	}
//	//	// Compute and print the average distance.
//	//	FType average_distance = sum_distances / (*it)->indices_of_assigned_points().size();
//	//	std::cout << " average distance: " << average_distance << std::endl;
//	//	// Proceed with the next detected shape.
//	//	++it;
//	//}
//
//	// ich brauche eine liste der face center pts
//	// (*it)->indices_of_assigned_points().
//	// mit dem befehl krieg ich einen vector der point ids die zu dem shape gehören
//	// wenn ich oben den point vec aufstelle, konstruiere ich mir noch einen vector der mir auflöst, zu welcher cell id eine pt id gehört
//	// also vec[cellid] = pointcentroidID
//	// Dann durchsuche ich den cellid vector nach der pointcentroidID
//	// der richtige index ist die gesuchte cellid, das schreibe ich mir für einen vtk color field raus
//	// irgendwie: plane = 1, cyl = 2, cone = 3, unassigned = 0
//	// dann noch sowas wie eine durchgehende nummer, also zwei felder
//	Classification.resize(STLMesh->GetNumberOfCells(), 0);
//
//	it = shapes.begin();
//	while (it != shapes.end()) {
//		// determine shape
//		std::size_t shapeIdentifier = 0;
//		if (CylinderSD* cyl = dynamic_cast<CylinderSD*>(it->get())) {
//			shapeIdentifier = 2;
//		}
//		else if (PlaneSD* pl = dynamic_cast<PlaneSD*>(it->get())) {
//			shapeIdentifier = 1;
//		}
//		else if (ConeSD* csd = dynamic_cast<ConeSD*>(it->get())) {
//			shapeIdentifier = 3;
//		}
//		// search for point index in cellid Vec
//		// Iterate through point indices assigned to each detected shape.
//		//std::vector<std::size_t>::const_iterator indexC_it = (*it)->indices_of_assigned_points().begin();
//		//while (indexC_it != (*it)->indices_of_assigned_points().end()) {
//		//	// Retrieve point.
//		//	//const Point_with_normal& p = *(points.begin() + (*index_it));
//		//	if (repPointsperCell.count((*indexC_it))) {
//		//		Classification[repPointsperCell[*indexC_it]] = shapeIdentifier;
//		//	}
//		//	// Proceed with the next point.
//		//	++indexC_it;
//		//}
//		for (size_t i = 0; i < (*it)->indices_of_assigned_points().size(); ++i) {
//			if ((repPointofCell.count((*it)->indices_of_assigned_points()[i])) > 0) {
//				Classification[repPointofCell[(*it)->indices_of_assigned_points()[i]]] = shapeIdentifier;
//			}
//		}
//		++it;
//	}
//
//	//CGAL::IO::write_PLY(fullpath, polygon_mesh);
//
//	return EXIT_SUCCESS;
//}
