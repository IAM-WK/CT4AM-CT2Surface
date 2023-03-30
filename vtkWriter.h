#pragma once
#include "Misc.h"

// base class to provide basics for writing vtk files

namespace vtkWriter
{
	// anonymous namespace for encapsulated data storage
	namespace {
		bool s_datahasbeenprepared = false;
		bool s_contourFieldWritten = false;
		bool s_slidingModeenabled = false;
		size_t s_numberofVertices = 0;
		size_t s_numberofCells = 0;
		std::vector<std::string> s_argList = std::vector<std::string>();

		// CellID, FaceGroupID, PtIDs, Normal
		std::vector < std::tuple<std::size_t, std::size_t, std::vector<std::size_t>, NormalVectorType > > s_PolydataStruct;

		// index = ID, entry = x,y,z
		std::vector <itk::Vector<MeshCoordType, 3>> s_PointStruct;

		// build angles per cell
		std::vector <double> s_buildAngles;
		itk::Vector<double, 3> s_BuildDirection;

		// isAnEdgeCell per cell
		std::vector<int> s_isAnEdgeCell;

		// state vars
		std::vector<int> s_InvalidSurfSize; // size was large enough = 0
		std::vector<int> s_LevelingValidity; // leveling was successfull = 0
	}

	void prepareData(const MeshType::Pointer STLMesh, std::vector<std::tuple<std::vector<itk::IdentifierType>, std::vector<itk::IdentifierType>, itk::CovariantVector< double, 3>>>& FaceGroupsCellNums, const NormalQEMeshType::CellDataContainerPointer FaceNormals, const std::vector<std::array<int, 2>>& SurfaceValidity, const std::vector<std::string>& argList = std::vector<std::string>()) {
		s_numberofVertices = STLMesh->GetNumberOfPoints();
		s_numberofCells = STLMesh->GetNumberOfCells();
		s_argList = argList;

		for (size_t CellID = 0; CellID < s_numberofCells; ++CellID) {
			// find in which FaceGroup current cell number is
			std::size_t FaceGroupID = 0;
			bool CellIDfound = false; // for sanity check
			if (s_slidingModeenabled) {
				CellIDfound = true; // trivial
				FaceGroupID = CellID;
			}
			else {
				std::vector<itk::IdentifierType>::iterator searchIter;
				for (size_t facegroupNum = 0; facegroupNum < FaceGroupsCellNums.size(); ++facegroupNum)
				{
					searchIter = std::find(std::get<1>(FaceGroupsCellNums[facegroupNum]).begin(), std::get<1>(FaceGroupsCellNums[facegroupNum]).end(), CellID);
					if (searchIter != std::get<1>(FaceGroupsCellNums[facegroupNum]).end()) {
						FaceGroupID = facegroupNum;
						CellIDfound = true;
						break;
					}
				}
				if (!CellIDfound) {
					std::cerr << "Error in vtkWriter::prepareData: CellID " << CellID << " was not found in any facegroup! \n";
				}
			}

			// find PointIDs corresponding to this cell
			CellAutoPointer CAPtoCellID;
			STLMesh->GetCell(CellID, CAPtoCellID);
			itk::Array<MeshType::PointIdentifier> PointIDsofCellID = CAPtoCellID->GetPointIdsContainer();
			// we deal only with triangle meshs -- 3 pointIDs per cell
			std::vector<itk::IdentifierType> PointIDs = { PointIDsofCellID[0],PointIDsofCellID[1],PointIDsofCellID[2] };

			s_PolydataStruct.push_back(std::make_tuple(CellID, FaceGroupID, PointIDs, FaceNormals->GetElement(CellID)));
			itk::CovariantVector<double, 3> curNormal;
			curNormal[0] = FaceNormals->GetElement(CellID)[0];
			curNormal[1] = FaceNormals->GetElement(CellID)[1];
			curNormal[2] = FaceNormals->GetElement(CellID)[2];
			s_buildAngles.push_back(computeBuildingAngle(curNormal, s_BuildDirection));
			s_isAnEdgeCell.push_back(descriptiveIndices::s_isanEdgeCell[CellID]);

			s_InvalidSurfSize.push_back(SurfaceValidity[FaceGroupID][0]);
			s_LevelingValidity.push_back(SurfaceValidity[FaceGroupID][1]);
		}
		// fill Point struct
		for (size_t PointID = 0; PointID < s_numberofVertices; ++PointID) {
			itk::Vector<MeshCoordType, 3> pointCoords;
			pointCoords[0] = STLMesh->GetPoint(PointID)[0];
			pointCoords[1] = STLMesh->GetPoint(PointID)[1];
			pointCoords[2] = STLMesh->GetPoint(PointID)[2];
			s_PointStruct.push_back(pointCoords);
		}
		s_datahasbeenprepared = true;
	}

	// vtk output
	// this writes coordinates to a .vtk file
	void outputvtkfile(const std::string& filename) {
		// preparation:
		// count number of points
		//vtkWriter::s_numberofVertices = STLMesh->GetNumberOfPoints();
		if (!s_datahasbeenprepared) {
			std::cerr << "Error encountered: Please first prepare the data by invoking vtkWriter::prepareData(...) before invoking outputvtkfile!" << std::endl;
		}

		std::ofstream vtkfile(filename, std::ios::binary);

		// print header
		vtkfile << "# vtk DataFile Version 4.2 # "; // the command line comment has to be in the same line as the version string to preserve compatibility
		for (size_t j = 0; j < s_argList.size(); ++j) {
			vtkfile << s_argList[j] << " ";
		}
		vtkfile << "\n";
		vtkfile << "Mesh Coordinates\n";
		vtkfile << "ASCII\n";
		vtkfile << "DATASET POLYDATA\n";
		vtkfile << "POINTS " << s_numberofVertices << " float\n";
		for (size_t i = 0; i < s_PointStruct.size(); ++i) {
			vtkfile << s_PointStruct[i][0] << " " << s_PointStruct[i][1] << " " << s_PointStruct[i][2] << "\n";
		}
		vtkfile << "\n";
		vtkfile << "POLYGONS " << s_numberofCells << " " << (4 * s_numberofCells) << "\n"; // triangle = 3 + 3 pts = 4 numbers per cell

		for (size_t cellNum = 0; cellNum < s_numberofCells; ++cellNum) {
			vtkfile << "3";
			for (size_t pointNum = 0; pointNum < 3; ++pointNum) {
				vtkfile << " " << std::get<2>(s_PolydataStruct[cellNum])[pointNum];
			}
			vtkfile << "\n";
		}

		vtkfile.flush();
		vtkfile.close();
		std::cout << "successfully wrote vtk file! \n";
	}

	// this writes coordinates to a .vtk file
	template <typename T>
	void outputcontourvtkfile(const std::vector<T>& fieldvec, const std::string& filename, const std::string& fielddenominator) {
		std::ofstream vtkfile;
		vtkfile.open(filename, std::ios_base::app);
		// print color scalars
		vtkfile << "\n";
		//point data can only be specified once
		if (!s_contourFieldWritten) { // this is omitted if there were other fields written before
			s_contourFieldWritten = true;
			vtkfile << "CELL_DATA " << fieldvec.size() << " \n";
		}
		// just a consistency check
		if (s_numberofCells != fieldvec.size()) {
			std::cerr << "WARNING! contourfield size doesn't match cell size in field: " << fielddenominator << "!" << std::endl;
			std::cerr << "contourfield size = " << fieldvec.size() << " - while cell size is = " << s_numberofCells << "!" << std::endl;
		}
		// typeid().name may give mangled names on not MSVC compilers
		vtkfile << "SCALARS" << " " << fielddenominator << " " << typeid(T).name() << " 1 " << "\n";
		vtkfile << "LOOKUP_TABLE default \n";
		for (size_t i = 0; i < fieldvec.size(); ++i) {
			vtkfile << fieldvec[i] << " \n";
		}

		vtkfile.flush();
		vtkfile.close();
	}

	// worker process that outputs data
	void writerWorker(const std::string& vtkfilename, const std::vector<std::array<double, 7>>& RoughnessTable, const std::vector<std::array<double, 6>>& SubSurfTable, const std::vector<std::array<double, 6>>& SubSurfClosedTable, const std::vector<std::array<double, 6>>& SubSurfOpenTable) {
		std::cout << "writing vtk file!\n";
		outputvtkfile(vtkfilename);

		std::size_t tablelength = s_numberofCells;
		std::vector<double> Sa, Sq, Sz, Sp, Sv, Ssk, Sku, maxPoro, clustersPerArea, maxClusterSize, avgClusterSize, projectedPorosity, volumetricPorosity;
		std::vector<double>  ClosedmaxPoro, ClosedclustersPerArea, ClosedmaxClusterSize, ClosedavgClusterSize, ClosedprojectedPorosity, ClosedvolumetricPorosity;
		std::vector<double>  OpenmaxPoro, OpenclustersPerArea, OpenmaxClusterSize, OpenavgClusterSize, OpenprojectedPorosity, OpenvolumetricPorosity;
		Sa.reserve(tablelength);
		Sq.reserve(tablelength);
		Sz.reserve(tablelength);
		Sp.reserve(tablelength);
		Sv.reserve(tablelength);
		Ssk.reserve(tablelength);
		Sku.reserve(tablelength);
		maxPoro.reserve(tablelength);
		clustersPerArea.reserve(tablelength);
		maxClusterSize.reserve(tablelength);
		avgClusterSize.reserve(tablelength);
		projectedPorosity.reserve(tablelength);
		volumetricPorosity.reserve(tablelength);
		// closed
		ClosedmaxPoro.reserve(tablelength);
		ClosedclustersPerArea.reserve(tablelength);
		ClosedmaxClusterSize.reserve(tablelength);
		ClosedavgClusterSize.reserve(tablelength);
		ClosedprojectedPorosity.reserve(tablelength);
		ClosedvolumetricPorosity.reserve(tablelength);
		// open
		OpenmaxPoro.reserve(tablelength);
		OpenclustersPerArea.reserve(tablelength);
		OpenmaxClusterSize.reserve(tablelength);
		OpenavgClusterSize.reserve(tablelength);
		OpenprojectedPorosity.reserve(tablelength);
		OpenvolumetricPorosity.reserve(tablelength);

		std::vector<int> FGisAtEdge;
		FGisAtEdge.reserve(tablelength);

		std::size_t faceGroupNum = 0;
		for (size_t cellNum = 0; cellNum < s_numberofCells; ++cellNum) {
			faceGroupNum = std::get<1>(s_PolydataStruct[cellNum]);
			Sa.push_back(RoughnessTable[faceGroupNum][0]);
			Sq.push_back(RoughnessTable[faceGroupNum][1]);
			Sz.push_back(RoughnessTable[faceGroupNum][2]);
			Sp.push_back(RoughnessTable[faceGroupNum][3]);
			Sv.push_back(RoughnessTable[faceGroupNum][4]);
			Ssk.push_back(RoughnessTable[faceGroupNum][5]);
			Sku.push_back(RoughnessTable[faceGroupNum][6]);
			// and fill subsurf params
			maxPoro.push_back(SubSurfTable[faceGroupNum][0]);
			clustersPerArea.push_back(SubSurfTable[faceGroupNum][1]);
			maxClusterSize.push_back(SubSurfTable[faceGroupNum][2]);
			avgClusterSize.push_back(SubSurfTable[faceGroupNum][3]);
			projectedPorosity.push_back(SubSurfTable[faceGroupNum][4]);
			volumetricPorosity.push_back(SubSurfTable[faceGroupNum][5]);
			// and fill subsurf params closed
			ClosedmaxPoro.push_back(SubSurfClosedTable[faceGroupNum][0]);
			ClosedclustersPerArea.push_back(SubSurfClosedTable[faceGroupNum][1]);
			ClosedmaxClusterSize.push_back(SubSurfClosedTable[faceGroupNum][2]);
			ClosedavgClusterSize.push_back(SubSurfClosedTable[faceGroupNum][3]);
			ClosedprojectedPorosity.push_back(SubSurfClosedTable[faceGroupNum][4]);
			ClosedvolumetricPorosity.push_back(SubSurfClosedTable[faceGroupNum][5]);
			// and fill subsurf params open
			OpenmaxPoro.push_back(SubSurfOpenTable[faceGroupNum][0]);
			OpenclustersPerArea.push_back(SubSurfOpenTable[faceGroupNum][1]);
			OpenmaxClusterSize.push_back(SubSurfOpenTable[faceGroupNum][2]);
			OpenavgClusterSize.push_back(SubSurfOpenTable[faceGroupNum][3]);
			OpenprojectedPorosity.push_back(SubSurfOpenTable[faceGroupNum][4]);
			OpenvolumetricPorosity.push_back(SubSurfOpenTable[faceGroupNum][5]);
			// fill debug helpers
			FGisAtEdge.push_back(aggregatePredicates::s_FacegrouptouchesEdge[faceGroupNum]);
		}

		std::vector<int> FaceGroupsofCells;
		FaceGroupsofCells.reserve(s_PolydataStruct.size());
		for (const std::tuple<std::size_t, std::size_t, std::vector<std::size_t>, NormalVectorType >& cellInfos : s_PolydataStruct) {
			FaceGroupsofCells.push_back(std::get<1>(cellInfos));
		}

		std::cout << "writing mesh descriptor values to vtk file!\n";
		outputcontourvtkfile(FaceGroupsofCells, vtkfilename, "MFaceGroup");
		outputcontourvtkfile(s_InvalidSurfSize, vtkfilename, "MSizeIsInvalid");
		outputcontourvtkfile(s_LevelingValidity, vtkfilename, "MLevelingSuccess");
		outputcontourvtkfile(s_buildAngles, vtkfilename, "MBuildAngles");
		outputcontourvtkfile(s_isAnEdgeCell, vtkfilename, "MEdgeCellProperty");
		outputcontourvtkfile(FGisAtEdge, vtkfilename, "MEdgeFGProperty");
		std::cout << "writing Height parameter values to vtk file!\n";
		outputcontourvtkfile(Sa, vtkfilename, "Sa");
		outputcontourvtkfile(Sq, vtkfilename, "Sq");
		outputcontourvtkfile(Sz, vtkfilename, "Sz");
		outputcontourvtkfile(Sp, vtkfilename, "Sp");
		outputcontourvtkfile(Sv, vtkfilename, "Sv");
		outputcontourvtkfile(Ssk, vtkfilename, "Ssk");
		outputcontourvtkfile(Sku, vtkfilename, "Sku");
		std::cout << "writing SubSurf descriptors overall to vtk file!\n";
		outputcontourvtkfile(maxPoro, vtkfilename, "TotalmaxSubSurfPorosity");
		outputcontourvtkfile(clustersPerArea, vtkfilename, "TotalPoresPerArea");
		outputcontourvtkfile(maxClusterSize, vtkfilename, "TotalmaxPoreSize");
		outputcontourvtkfile(avgClusterSize, vtkfilename, "TotalavgPoreSize");
		outputcontourvtkfile(projectedPorosity, vtkfilename, "TotalAreaSubSurfPorosity");
		outputcontourvtkfile(volumetricPorosity, vtkfilename, "TotalVolSubSurfPorosity");
		std::cout << "writing SubSurf descriptors closed to vtk file!\n";
		outputcontourvtkfile(ClosedmaxPoro, vtkfilename, "ClosedmaxSubSurfPorosity");
		outputcontourvtkfile(ClosedclustersPerArea, vtkfilename, "ClosedPoresPerArea");
		outputcontourvtkfile(ClosedmaxClusterSize, vtkfilename, "ClosedmaxPoreSize");
		outputcontourvtkfile(ClosedavgClusterSize, vtkfilename, "ClosedavgPoreSize");
		outputcontourvtkfile(ClosedprojectedPorosity, vtkfilename, "ClosedAreaSubSurfPorosity");
		outputcontourvtkfile(ClosedvolumetricPorosity, vtkfilename, "ClosedVolSubSurfPorosity");
		std::cout << "writing SubSurf descriptors open to vtk file!\n";
		outputcontourvtkfile(OpenmaxPoro, vtkfilename, "OpenmaxSubSurfPorosity");
		outputcontourvtkfile(OpenclustersPerArea, vtkfilename, "OpenPoresPerArea");
		outputcontourvtkfile(OpenmaxClusterSize, vtkfilename, "OpenmaxPoreSize");
		outputcontourvtkfile(OpenavgClusterSize, vtkfilename, "OpenavgPoreSize");
		outputcontourvtkfile(OpenprojectedPorosity, vtkfilename, "OpenAreaSubSurfPorosity");
		outputcontourvtkfile(OpenvolumetricPorosity, vtkfilename, "OpenVolSubSurfPorosity");
	}
};
