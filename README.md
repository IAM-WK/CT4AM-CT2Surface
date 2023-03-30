# CT4AM Surface

CT4AM Surface is a C++ tool used to analyze the complete surface of an additively manufactured part from XCT scans. It provides a detailed analysis of the surface quality and subsurface porosity parameters of the part, which are mapped onto the input mesh to show how the surface quality is distributed over the part.


## Cite

To be done

## Prerequisites

- STL of the part: 
  - It is recommended to use a finely meshed, isotropic meshed version of the geometry, which can be obtained by remeshing the geometry in MeshLab using isotropic explicit remeshing. 
  - STL file should be in ASCII format!

- XCT scan of the part:
  - Should be in .mhd format.
  - Should be aligned in the rotational axis to the orientation of the STL.
  - The translation vector should be found using Paraview. It is needed to parametrize the program.

## Usage

To use CT4AM Surface, follow these steps:

./CT2Surface /path/to/stl_file /path/to/xct_file


## Troubleshooting

If you encounter any issues while using CT2Surface, please check the following:

- Ensure that the input files are in the correct format and meet the prerequisites listed above.

If you are still experiencing issues, please feel free to open an issue on this GitHub repository.

