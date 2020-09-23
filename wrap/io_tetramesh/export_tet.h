/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *   
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

/**
@name Load and Save in Tetgen File format
*/
//@{

#ifndef __VCGLIB_TETRAEXPORT_TET
#define __VCGLIB_TETRAEXPORT_TET

#include<wrap/io_trimesh/precision.h>

#include <stdio.h>

namespace vcg {
namespace tetra {
namespace io {


template <class SaveMeshType>
class ExporterTET
{
public:
    typedef ::vcg::ply::PropDescriptor PropDescriptor ;
    typedef typename SaveMeshType::VertexPointer VertexPointer;
    typedef typename SaveMeshType::ScalarType ScalarType;
    typedef typename SaveMeshType::VertexType VertexType;
    typedef typename SaveMeshType::TetraType TetraType;
    typedef typename SaveMeshType::TetraPointer TetraPointer;
    typedef typename SaveMeshType::VertexIterator VertexIterator;
    typedef typename SaveMeshType::TetraIterator TetraIterator;

    static bool Save(SaveMeshType &m,  const char * filename)	// V1.0
    {
        FILE *  eleFile;
        FILE * nodeFile;

        const int DGT   = vcg::tri::io::Precision<ScalarType>::digits();

        std::string eleName(filename);
        std::string nodeName(filename);

        eleName.append(".ele");
        nodeName.append(".node");

        eleFile = fopen(eleName.c_str() , "wb");
        if (eleFile == NULL)
            return 1;
        
        nodeFile = fopen(nodeName.c_str() , "wb");
        if (nodeFile == NULL)
            return 1;
        
        fprintf(eleFile,
                "#VCG lib .ele generated file\n"
                "#This file contains the tetra definition in accord to Tetgen file format\n");
        
        fprintf(nodeFile,
                "#VCG lib .node generated file\n"
                "#This file contains the vertex definition in accord to Tetgen file format\n");
       

        fprintf(nodeFile, "%d 3 0 0\n", m.VN());
        int vi = 0;
        SimpleTempData<typename SaveMeshType::VertContainer, int> indices(m.vert);

        ForEachVertex(m, [&] (VertexType & v) {
            indices[&v] = vi;
            fprintf(nodeFile, "%d %g %g %g\n", vi++, double(v.P().X()), double(v.P().Y()), double(v.P().Z()));
        });

        fprintf(eleFile, "%d 4 0\n", m.TN());
        int ti = 0;
        ForEachTetra(m, [&] (SaveMeshType::TetraType & t) {
            fprintf(eleFile, "%d %d %d %d %d\n", ti++, indices[t.V(0)], indices[t.V(1)], indices[t.V(2)], indices[t.V(3)]);
        });


        fclose(eleFile);
        fclose(nodeFile);

        return 0;
    }




}; // end class



} // end namespace tetra
} // end namespace io
} // end namespace vcg

#endif
