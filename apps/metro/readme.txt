
   VCGLib  http://vcg.sf.net                                         o o     
   Visual and Computer Graphics Library                            o     o   
                                                                  _   O  _   
   Copyright(C) 2004-2005                                                \/)\/    
   Visual Computing Lab  http://vcg.isti.cnr.it                    /\/|      
   ISTI - Italian National Research Council                           |      
                                                                      \      
   Metro 4.01 21/09/2004
   All rights reserved.                                                      
   

                                                                       
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

Metro is a tool designed to evaluate the difference between two triangular meshes. 
Metro adopts an approximated approach based on surface sampling and point-to-surface distance computation. 
Please, when using this tool, cite the following reference:

P. Cignoni, C. Rocchini and R. Scopigno
"Metro: measuring error on simplified surfaces"
Computer Graphics Forum, Blackwell Publishers, vol. 17(2), June 1998, pp 167-174

You can find some sample mesh to test in the 'Metro Sample dataset' package downloadable from sourceforge.

For any question about this software please contact:
Paolo Cignoni ( p.cignoni@isti.cnr.it )

--- General Info ---

Metro is a tool designed to evaluate the difference between two triangular meshes. 
Metro adopts an approximated approach based on surface sampling and point-to-surface distance computation. 
Three different surface sampling methods are implemented:

    *   Montecarlo sampling (pick k random samples in the interior of each face)
    *   Subdivision sampling (recursively subdivide each face along the longest edge and choose the sample in the center of each cell)
    *   Similar Triangles sampling (subdivide each face F in k polygons similar to F and sample the face in correspondence with the vertices of these polygons, internal to F)

Note that the three methods described above are used to sample only the interior of each face. 
A different scheme is used to sample vertices and edges: vertices are sampled in the straightforward manner, 
while edges are sampled by uniformly interleaving samples along each edge.

--- Basic usage ---

Metro is a command-line tool which allows the user to select among different sampling schemes. 
A list of the command-line parameters accepted by the tool is shown in the following.

Usage: Metro file1 file2 [opts]

where "file1" and "file2" are the input meshes in PLY or STL format, and opts can be:

   -v       disable vertex sampling
   -e       disable edge sampling
   -f       disable face sampling
   -u       does not ignore unreferred vertices (sample also unreferenced vertices 
            (useful for sampling point clouds against meshes)
   -Sx      set the face sampling mode
            where x can be:
             -S0  montecarlo sampling
             -S1  subdivision sampling
             -S2  similar triangles sampling (Default)
   -n#      set the required number of samples (overrides -a)
   -a#      set the required number of samples per area unit (overrides -n)
   -c       save computed error as vertex colour and quality in two ply files
   -C # #   Set the min/max values used for color mapping (useful for taking snapshot with coherent color ramp)
   -L       Remove duplicated and unreferenced vertices before processing to avoid 
