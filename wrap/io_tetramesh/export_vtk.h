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

#ifndef __VCGLIB_TETRAEXPORT_VTK
#define __VCGLIB_TETRAEXPORT_VTK

#include <wrap/io_trimesh/precision.h>

#include <stdio.h>

namespace vcg
{
namespace tetra
{
namespace io
{

template <class SaveMeshType>
class ExporterVTK
{
public:
    typedef ::vcg::ply::PropDescriptor PropDescriptor;
    typedef typename SaveMeshType::VertexPointer VertexPointer;
    typedef typename SaveMeshType::ScalarType ScalarType;
    typedef typename SaveMeshType::VertexType VertexType;
    typedef typename SaveMeshType::TetraType TetraType;
    typedef typename SaveMeshType::TetraPointer TetraPointer;
    typedef typename SaveMeshType::VertexIterator VertexIterator;
    typedef typename SaveMeshType::TetraIterator TetraIterator;

    enum
    {
        VTK_VERTEX = 1,
        VTK_POLY_VERTEX,
        VTK_LINE,
        VTK_POLY_LINE,
        VTK_TRIANGLE,
        VTK_TRIANGLE_STRIP,
        VTK_POLYGON,
        VTK_PIXEL,
        VTK_QUAD,
        VTK_TETRA,
        VTK_VOXEL,
        VTK_HEXAHEDRON,
        VTK_WEDGE,
        VTK_PYRAMID
    };

    static bool Save(SaveMeshType &m, const char *filename) // V1.0
    {
        FILE *vtkFile;

        const int DGT    = vcg::tri::io::Precision<ScalarType>::digits();
        const char* vttp = vcg::tri::io::Precision<ScalarType>::typeName();

        std::string vtkName(filename);

        vtkFile = fopen(vtkName.c_str(), "wb");
        if (vtkFile == NULL)
            return 1;

        fprintf(vtkFile, "# vtk DataFile Version 2.0\n");
        fprintf(vtkFile, "VCG_mesh\n");
        fprintf(vtkFile, "ASCII\n");
        fprintf(vtkFile, "DATASET UNSTRUCTURED_GRID\n");
        fprintf(vtkFile, "POINTS %d %s\n", m.VN(), vttp);

        int vi = 0;
        SimpleTempData<typename SaveMeshType::VertContainer, int> indices(m.vert);

        int maxQ = 0;
        ForEachVertex(m, [&](VertexType &v) {
            indices[&v] = vi++;
            fprintf(vtkFile, "%g %g %g\n", v.P().X(), v.P().Y(), v.P().Z());
        });
        
        fprintf(vtkFile, "CELLS %d %d\n", m.TN(), m.TN()*5);

        int ti = 0;
        ForEachTetra(m, [&](SaveMeshType::TetraType &t) {
            ++ti;
            fprintf(vtkFile, "%d %d %d %d %d\n", 4, indices[t.V(0)], indices[t.V(1)], indices[t.V(2)], indices[t.V(3)]);
        });

        fprintf(vtkFile, "CELL_TYPES %d\n", m.TN());
        ForEachTetra(m, [&](SaveMeshType::TetraType &t) {
            fprintf(vtkFile, "%d\n", VTK_TETRA);
        });

        if( HasPerVertexQuality(m))
        {
            fprintf(vtkFile, "POINT_DATA %d\n", m.VN());
            fprintf(vtkFile, "SCALARS vquality %s 1\n", vttp);
            fprintf(vtkFile, "LOOKUP_TABLE default\n");
            ForEachVertex(m, [&](VertexType &v) {
                fprintf(vtkFile, "%g\n", v.Q());
            });
        }

        if( HasPerTetraQuality(m))
        {
            fprintf(vtkFile, "CELL_DATA %d\n", m.TN());
            fprintf(vtkFile, "SCALARS tquality %s 1\n", vttp);
            fprintf(vtkFile, "LOOKUP_TABLE default\n");
            ForEachTetra(m, [&](TetraType &t) {
                fprintf(vtkFile, "%g\n", t.Q());
            });
        }

        fclose(vtkFile);

        return 0;
    }

}; // end class

} // namespace io
} // namespace tetra
} // end namespace vcg

#endif
