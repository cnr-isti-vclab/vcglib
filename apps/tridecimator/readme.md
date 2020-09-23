
   VCGLib  http://www.vcglib.net 
    Copyright(C) 2005-2006             
   Visual Computing Lab  http://vcg.isti.cnr.it          
   ISTI - Italian National Research Council                 
   

                                                                       
This program is free software; you can redistribute it and/or modify      
it under the terms of the GNU General Public License as published by      
the Free Software Foundation; either version 2 of the License, or         
(at your option) any later version.                                       
                                                                          
This program is distributed in the hope that it will be useful,           
but WITHOUT ANY WARRANTY; without even the implied warranty of            
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             
GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          
for more details.                                                 

--- Synopsis ---

`tridecimator fileIn fileOut face_num [opt]`

Tridecimator is a commandline mesh simplifier based on a variant of the quadric error edge collapse strategy. 
It supports `PLY, OFF, OBJ` format for input and output and require a target number of faces. 
The following options are supported:

-e# QuadricError threshold  (range [0,inf) default inf) 
-b# Boundary Weight (default .5) 
-p# Quality quadric Weight (default .5) 
-q# Quality threshold (range [0.0, 0.866],  default .3 ) 
-n# Normal threshold  (degree range [0,180] default 90) 
-w# Quality weight factor  (10) 
-E# Minimal admitted quadric value (default 1e-15, must be >0) 
-Q[y|n]  Use or not Face Quality Threshold (default yes) 
-H[y|n]  Use or not HardQualityCheck (default no) 
-N[y|n]  Use or not Face Normal Threshold (default no) 
-P[y|n]  Add or not QualityQuadric (default no) 
-A[y|n]  Use or not Area Checking (default no) 
-O[y|n]  Use or not vertex optimal placement (default yes) 
-S[y|n]  Use or not Scale Independent quadric measure(default yes)  
-B[y|n]  Preserve or not mesh boundary (default no) 
-T[y|n]  Preserve or not Topology (default no) 
-W[y|n]  Use or not per vertex Quality to weight the quadric error (default no) 
-C       Before simplification, remove duplicate & unreferenced vertices 
    

This simplification tool employ a quadric error based edge collapse iterative approach. 
Each possible collapse is scored taking into account
- its quadric error
- its quadric boundary error
- the normal variation caused by the collapse itself
- the future quality (aspect ratio) of the triangles involved in the collapse

The simplification can ends when the desired number of triangles has been reached.

The 'Safe Heap Update' change the set of edges re-inserted in the heap
after a collapse resulting in a vertex v; this option put back in the
heap not only incident in the v, but all the edges of the triangles incident
in v. It slows down a lot.

The 'Scale Independent quadric' is useful when you want simplify 
two different meshes at the same 'error'; otherwise the quadric 
error used in the heap is somewhat normalized and independent 
of the size of the mesh.

The Normal threshold avoid to make large changes in variation of the normal
of the surfaces, but on the other hand it prevent the removal of small 'folded'
triangles that can be already present. Therefore in most cases is not very useful.

Cleaning the mesh is mandatory for some input format like STL that always
duplicates all the vertices.



