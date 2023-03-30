#pragma once
#include <cmath>
#include <tuple>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <limits>
#include <algorithm>
#include <vector>
#include <string>

#include <functional>
#include <future>
#include <chrono>

#include <itkImage.h>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkChangeInformationImageFilter.h"

#include "itkMesh.h"
#include "itkQuadEdgeMesh.hxx"
#include "itkSTLMeshIOFactory.h"
#include "itkSTLMeshIO.h"
#include "itkMeshFileReader.h"
#include "itkNormalQuadEdgeMeshFilter.h"

#include "itkResampleImageFilter.h"
#include "itkAffineTransform.h"
#include "itkImageToHistogramFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h" // for binarized source images
#include "itkLinearInterpolateImageFunction.h" // for greylevel source images
#include "itkHuangThresholdImageFilter.h"

#include <itkVector.h>
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"

// to manipulate masking size
#include "itkFlatStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkStatisticsImageFilter.h"

// for outside air masking
#include "itkBinaryShapeKeepNObjectsImageFilter.h"

//***********************************************************************************************************************************************
// ITK
//***********************************************************************************************************************************************

typedef unsigned short InputPixelType;
typedef double HeightMapPixelType;
typedef unsigned char MaskPixelType;
typedef double MeshCoordType; // coords of points in STL

const unsigned int StackMeshDimensionality = 3; // Mesh (and image) dimension
const unsigned int HeightMapDimensionality = 2; // Height map dimension

typedef itk::QuadEdgeMesh <MeshCoordType, StackMeshDimensionality> MeshType;
typedef itk::Image<InputPixelType, StackMeshDimensionality > StackImageType;
typedef itk::Image<unsigned int, StackMeshDimensionality> LabelStackImageType;
typedef itk::Image<unsigned char, StackMeshDimensionality> LabelBinStackImageType;
typedef itk::Image<HeightMapPixelType, HeightMapDimensionality > HeightMapImageType;
typedef itk::Image<MaskPixelType, HeightMapDimensionality > MaskImageType;

typedef itk::ImageFileReader< StackImageType >  StackReaderType;
typedef itk::ImageFileWriter< StackImageType >  StackWriterType;
typedef itk::ImageFileWriter< LabelBinStackImageType >  BinStackWriterType;
typedef itk::ImageFileWriter< HeightMapImageType >  HMWriterType;
typedef itk::ImageFileWriter< MaskImageType >  MaskwriterType;
typedef itk::ImageRegionConstIteratorWithIndex < StackImageType > StackConstIteratorTypeIndexed;
typedef itk::ImageRegionConstIteratorWithIndex < HeightMapImageType > HeightMapConstIteratorTypeIndexed;
typedef itk::ImageRegionIteratorWithIndex < HeightMapImageType > HeightMapIteratorTypeIndexed;
typedef itk::ImageRegionIteratorWithIndex < MaskImageType > MaskImageIteratorTypeIndexed;
typedef itk::ImageRegionConstIteratorWithIndex < MaskImageType > MaskImageConstIteratorTypeIndexed;

typedef itk::Vector<MeshCoordType, StackMeshDimensionality> NormalVectorType; // type for normals
typedef itk::QuadEdgeMesh <NormalVectorType, StackMeshDimensionality> NormalQEMeshType;
typedef itk::NormalQuadEdgeMeshFilter<MeshType, NormalQEMeshType> NormalFilterType;
NormalFilterType::WeightEnum  weight_type = itk::NormalQuadEdgeMeshFilterEnums::Weight::GOURAUD; // we don't care for the vertex normals ...

typedef MeshType::CellType::CellAutoPointer CellAutoPointer;

//***********************************************************************************************************************************************
// CGAL
//***********************************************************************************************************************************************

// CGAL libs to create plane fit for leveling
#include <CGAL/Simple_cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Default_diagonalize_traits.h>
typedef double                      FT;
typedef CGAL::Simple_cartesian<FT>  K;
typedef K::Plane_3                  Plane;
typedef K::Point_3                  Point3D;
typedef K::Point_2                  Point2D;
// CGAL for masking
typedef K::Triangle_2				Triangle2D;
// CGAL for shape removal
#include <CGAL/intersections.h>
typedef K::Triangle_3			Triangle3D;
typedef K::Vector_3				Vector3D;
typedef K::Line_3				Line3D;

// getFaceGroupBoundingBox
#include <CGAL/bounding_box.h>
typedef K::Iso_rectangle_2 Rectangle_2D;

// getFaceGroupOptimalBoundingBox
#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/min_quadrilateral_2.h>
typedef CGAL::Polygon_2<K> Polygon_2;

//***********************************************************************************************************************************************
