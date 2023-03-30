# CT4AM Surface 

Surface analysis tool to analyze surfaces in XCT images 

The tool allows to analyze the complete surface of an additively manufactured part containing undercuts and freeform surfaces.
The results are mapped to the triangles of the input mesh, providing information how the surface quality is distributed over the part.

Also, subsurface porosity parameters are calculated and mapped to the input mesh. 



Cite

To be done


Prerequisites

STL of the part: 
- it is recommended to use a finely meshed, isotropic meshed version of the geometry, which can be obtained by remeshing the geometry in MeshLab using isotropic explicit remeshing. 
- STL file should be in ASCII format!

XCT scan of the part
- Should be in .mhd format.
- should be aligned in the rotational axis to the orientation of the STL
- the translation vector should be found using Paraview. It is needed to parametrize the program. 




Troubleshooting


