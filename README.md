# CT4AM CT2Surface

![graphical abstract](./doc/GrabsV1.png)

CT4AM CT2Surface is a C++ tool used to analyze the complete surface of an additively manufactured part from XCT scans. It provides a detailed analysis of the surface quality and subsurface porosity parameters of the part, which are mapped onto the input mesh to show how the surface quality is distributed over the part.


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

## Download

[Download](https://github.com/IAM-WK/CT4AM_Surface/releases/download/v0.9/CT2Surface.exe)


## Usage
Exemplary call of CT4AM CT2Surface:

`.\CT2Surface.exe -i .\XCT_file.mhd -m .\STL_file.stl -u "-12.4 -3.75 -3.725" -o "output.mhd" -a 10 -p 400 -v -w 2 -k -e -10`

To get a list of options, type 

`.\CT2Surface.exe -h`

## Troubleshooting

If you encounter any issues while using CT2Surface, please check the following:

- Ensure that the input files are in the correct format and meet the prerequisites listed above.

If you are still experiencing issues, please feel free to open an issue on this GitHub repository.

